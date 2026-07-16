//! The primary `compare` command: an on-rails pipeline that takes two PAFs and
//! produces, in one invocation, the per-read comparison table plus (optionally)
//! the per-set alninfo + readinfo tables.
//!
//! It owns its preconditions rather than trusting the user:
//!   1. **sorts** both inputs by `Query_Name` (identical deterministic rule), so
//!      grouping and matching order are guaranteed,
//!   2. **verifies** the two PAFs carry the same `Query_Name` set (O(1) check),
//!      erroring by default if they differ, then
//!   3. in a **single in-memory pass**, collapses both sorted PAFs in lock-step
//!      and feeds the merge-join directly — no readinfo written-then-reread. The
//!      alninfo + readinfo tables are tee'd out as side outputs as it goes
//!      (suppressible with `--no-alninfo` / `--no-readinfo`).
//!
//! This is the porcelain over the plumbing subcommands (`paf2tables`,
//! `compare-readinfo`, …): the comparison table is byte-identical to running
//! `compare-readinfo` on the sorted readinfo files, and the side outputs are
//! byte-identical to `paf2tables` on the sorted PAFs.
//!
//! Precondition (documented, not enforced): a `Query_Name` uniquely identifies a
//! single read/sequence — so sorting by name alone (no `Read_Len` secondary key)
//! is sufficient for the downstream `(Read_Name, Read_Len)` merge-join.

use std::collections::HashMap;
use std::fs;
use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::compare_junctions::{emit_compare_junctions_row, write_compare_junctions_header};
use crate::compare_streaming::{emit_compare_row, write_compare_header, CompareMode};
use crate::compare_summary::{classify, CompareSummary};
use crate::external_sort::{parse_mem, read_id_set_check, sort_paf_to_file};
use crate::find_query_diff::{self, FindQueryDiffArgs};
use crate::io_utils::{open_input, open_output};
use crate::paf_groups::PafGroups;
use crate::readinfo::{collapse_group, ReadInfoRow, READINFO_HEADER};
use crate::record::AlnInfo;

/// Reject stdin (`-`) for `compare`'s two-file inputs. `compare` always needs two
/// independent files, so piping a single stream in for one side doesn't make
/// sense — require a real path for both `--paf-a` and `--paf-b`.
fn require_paf_path(s: &str) -> Result<String, String> {
    if s == "-" {
        Err("stdin ('-') is not supported here; provide a path to a PAF file".to_string())
    } else {
        Ok(s.to_string())
    }
}

#[derive(clap::Args, Debug)]
pub struct CompareArgs {
    /// PAF for dataset A ('.gz' auto-decompressed).
    #[arg(short = 'a', long = "paf-a", value_name = "a.paf", value_parser = require_paf_path)]
    paf_a: String,

    /// PAF for dataset B ('.gz' auto-decompressed).
    #[arg(short = 'b', long = "paf-b", value_name = "b.paf", value_parser = require_paf_path)]
    paf_b: String,

    /// Label for dataset A.
    #[arg(long = "label-a", value_name = "LABEL", default_value = "SetA")]
    label_a: String,

    /// Label for dataset B.
    #[arg(long = "label-b", value_name = "LABEL", default_value = "SetB")]
    label_b: String,

    /// Output directory.
    #[arg(short = 'o', long = "outdir", value_name = "DIR")]
    outdir: String,

    /// Filename prefix for all outputs.
    #[arg(short = 'p', long = "prefix", value_name = "NAME")]
    prefix: String,

    /// Comparison view: `full` (94 cols) or `junctions` (47-col splice view).
    #[arg(long = "mode", value_enum, default_value_t = CompareMode::Full)]
    mode: CompareMode,

    /// In-memory sort buffer per file (K/M/G suffix, or plain bytes).
    #[arg(long = "sort-mem", value_name = "SIZE", default_value = "1G")]
    sort_mem: String,

    /// Temp directory for temp out of memory sort files (default: --outdir).
    #[arg(long = "tmp-dir", value_name = "DIR")]
    tmp_dir: Option<String>,

    /// Number of sort threads (default: 1).
    #[arg(long = "sort-threads", value_name = "N", default_value_t = 1)]
    sort_threads: usize,

    /// Compare the shared intersection of aligned sequences instead of erroring when the two PAFs do
    /// not carry the exact same "Query_Name" set.
    #[arg(long = "allow-id-mismatch")]
    allow_id_mismatch: bool,

    /// Keep the intermediate sorted PAFs instead of deleting them at the end.
    #[arg(long = "keep-sorted-paf")]
    keep_sorted_paf: bool,

    /// Skip the internal sort: assume both PAFs already contain the same reads,
    /// grouped by "Query_Name" and in the same relative order.
    /// Not combinable with --allow-id-mismatch or --keep-sorted-paf.
    #[arg(long = "presorted", conflicts_with_all = ["allow_id_mismatch", "keep_sorted_paf"])]
    presorted: bool,

    /// Do not write the per-set alninfo (35-col) tables.
    #[arg(long = "no-alninfo")]
    no_alninfo: bool,

    /// Do not write the per-set readinfo (33-col) tables.
    #[arg(long = "no-readinfo")]
    no_readinfo: bool,

    /// Skip the "find-query-diff" step run at the end of the comparison.
    #[arg(long = "skip-find-query-diff")]
    skip_find_query_diff: bool,
}

pub fn run(args: &CompareArgs) -> Result<()> {
    // ── Step 0: setup ─────────────────────────────────────────────────────────
    let outdir = Path::new(&args.outdir);
    fs::create_dir_all(outdir)
        .with_context(|| format!("cannot create --outdir '{}'", args.outdir))?;
    let tmp_dir = args
        .tmp_dir
        .as_deref()
        .map(PathBuf::from)
        .unwrap_or_else(|| outdir.to_path_buf());
    let _ = fs::create_dir_all(&tmp_dir);
    let mem = parse_mem(&args.sort_mem)?;
    let junctions = matches!(args.mode, CompareMode::Junctions);

    let path = |name: String| outdir.join(name).to_string_lossy().into_owned();
    let a_sorted = path(format!("{}.{}.sorted.paf.gz", args.prefix, args.label_a));
    let b_sorted = path(format!("{}.{}.sorted.paf.gz", args.prefix, args.label_b));
    let a_alninfo = path(format!("{}.{}.alninfo.tsv.gz", args.prefix, args.label_a));
    let b_alninfo = path(format!("{}.{}.alninfo.tsv.gz", args.prefix, args.label_b));
    let a_readinfo = path(format!("{}.{}.readinfo.tsv.gz", args.prefix, args.label_a));
    let b_readinfo = path(format!("{}.{}.readinfo.tsv.gz", args.prefix, args.label_b));
    let compare_out = path(format!(
        "{}.compare{}.tsv.gz",
        args.prefix,
        if junctions { ".junctions" } else { "" }
    ));
    let summary_out = path(format!(
        "{}.compare{}.summary.tsv",
        args.prefix,
        if junctions { ".junctions" } else { "" }
    ));

    // Inputs fed to the compare pass: the freshly sorted temp files by default,
    // or the user's PAFs directly under --presorted (no sort, no set-check).
    let (a_in, b_in): (String, String) = if args.presorted {
        // ── --presorted: skip sort (Step 1) and set-check (Step 2) ────────────
        // The lex set-check assumes byte-lex order, which we don't require here;
        // instead the lock-step compare pass verifies the two PAFs carry the same
        // reads in the same order, erroring on the first divergence.
        eprintln!(
            "[INFO] --presorted: skipping sort and read-ID set-check; \
             same read order is verified during the compare pass."
        );
        (args.paf_a.clone(), args.paf_b.clone())
    } else {
        // ── Step 1: sort both PAFs by Query_Name (consistent rule) ────────────
        eprintln!(
            "[INFO] Step 1/3 — sorting both PAFs by Query_Name (mem={} bytes, tmp={})",
            mem,
            tmp_dir.display()
        );
        sort_paf_to_file(&args.paf_a, &a_sorted, mem, &tmp_dir, Some(args.sort_threads))
            .with_context(|| format!("sorting PAF A ({})", args.paf_a))?;
        sort_paf_to_file(&args.paf_b, &b_sorted, mem, &tmp_dir, Some(args.sort_threads))
            .with_context(|| format!("sorting PAF B ({})", args.paf_b))?;

        // ── Step 2: read-ID set-equality check (O(1) memory), before any output ─
        eprintln!("[INFO] Step 2/3 — verifying the two PAFs share the same read-ID set...");
        let chk = read_id_set_check(&a_sorted, &b_sorted, 5)?;
        eprintln!(
            "  shared: {}   only in {}: {}   only in {}: {}",
            chk.shared, args.label_a, chk.only_a, args.label_b, chk.only_b
        );
        if chk.only_a > 0 || chk.only_b > 0 {
            if !args.allow_id_mismatch {
                if !args.keep_sorted_paf {
                    let _ = fs::remove_file(&a_sorted);
                    let _ = fs::remove_file(&b_sorted);
                }
                bail!(
                    "read-ID sets differ between the two PAFs: {shared} shared, \
                     {oa} only in {la} (e.g. {exa}), {ob} only in {lb} (e.g. {exb}). \
                     Re-run with --allow-id-mismatch to compare the shared intersection.",
                    shared = chk.shared,
                    oa = chk.only_a,
                    la = args.label_a,
                    exa = chk.examples_a.join(", "),
                    ob = chk.only_b,
                    lb = args.label_b,
                    exb = chk.examples_b.join(", "),
                );
            }
            eprintln!(
                "  WARNING: read-ID sets differ; proceeding on the shared intersection \
                 (--allow-id-mismatch)."
            );
        }
        (a_sorted.clone(), b_sorted.clone())
    };

    // ── Step 3: single in-memory lock-step pass (collapse + compare + tee) ────
    if args.presorted {
        eprintln!("[INFO] comparing in one pass ({compare_out})...");
    } else {
        eprintln!("[INFO] Step 3/3 — comparing in one pass ({compare_out})...");
    }
    let mut summary = CompareSummary::default();
    let counts = compare_sorted_pafs(
        &a_in,
        &b_in,
        &args.label_a,
        &args.label_b,
        &compare_out,
        if args.no_readinfo { None } else { Some(&a_readinfo) },
        if args.no_readinfo { None } else { Some(&b_readinfo) },
        if args.no_alninfo { None } else { Some(&a_alninfo) },
        if args.no_alninfo { None } else { Some(&b_alninfo) },
        junctions,
        args.allow_id_mismatch,
        &mut summary,
    );
    if let Err(e) = counts {
        // The compare pass can fail partway (e.g. --presorted inputs that are
        // not actually in the same order), having already written part of the
        // output. Remove the partial artifacts so the failure leaves nothing
        // half-written, then surface the error (with a hint under --presorted).
        let _ = fs::remove_file(&compare_out);
        if !args.no_alninfo {
            let _ = fs::remove_file(&a_alninfo);
            let _ = fs::remove_file(&b_alninfo);
        }
        if !args.no_readinfo {
            let _ = fs::remove_file(&a_readinfo);
            let _ = fs::remove_file(&b_readinfo);
        }
        if !args.presorted && !args.keep_sorted_paf {
            let _ = fs::remove_file(&a_sorted);
            let _ = fs::remove_file(&b_sorted);
        }
        return if args.presorted {
            Err(e).context(
                "--presorted requires both PAFs to contain the same reads in the \
                 same order (grouped by Query_Name); omit --presorted to sort them \
                 automatically",
            )
        } else {
            Err(e)
        };
    }

    // ── cleanup + summary ─────────────────────────────────────────────────────
    if !args.presorted && !args.keep_sorted_paf {
        let _ = fs::remove_file(&a_sorted);
        let _ = fs::remove_file(&b_sorted);
    }
    // Aggregate summary statistics → sidecar TSV + stderr block.
    summary.write_tsv(&summary_out, &args.label_a, &args.label_b)?;
    summary.render_stderr(&args.label_a, &args.label_b);
    eprintln!("Outputs in {}:", args.outdir);
    if !args.no_alninfo {
        eprintln!("  {a_alninfo}");
        eprintln!("  {b_alninfo}");
    }
    if !args.no_readinfo {
        eprintln!("  {a_readinfo}");
        eprintln!("  {b_readinfo}");
    }
    eprintln!("  {compare_out}");
    eprintln!("  {summary_out}");
    if args.keep_sorted_paf {
        eprintln!("  {a_sorted}");
        eprintln!("  {b_sorted}");
    }

    // ── find-query-diff (default-on): query-different reads + merged regions ──
    if !args.skip_find_query_diff {
        eprintln!(
            "[INFO] running find-query-diff on the comparison table \
             (skip with --skip-find-query-diff)..."
        );
        let fq_args = FindQueryDiffArgs::for_compare(
            compare_out.clone(),
            args.outdir.clone(),
            args.prefix.clone(),
        );
        find_query_diff::run(&fq_args)?;
    }

    Ok(())
}

/// Serialize a `ReadInfoRow` to its readinfo-TSV line (no trailing newline) so it
/// can be parsed into a by-name column map for the comparison emitters.
fn readinfo_line(ri: &ReadInfoRow) -> Result<String> {
    let mut buf: Vec<u8> = Vec::with_capacity(256);
    ri.write(&mut buf)?;
    let s = String::from_utf8(buf).expect("readinfo serialization is valid UTF-8");
    Ok(s.trim_end_matches(|c| c == '\n' || c == '\r').to_string())
}

/// Pull the next per-read group from `groups`, teeing its alignment rows to the
/// optional `alninfo` sink and writing the collapsed readinfo row to the optional
/// `readinfo` sink, returning the collapsed `ReadInfoRow` (or `None` at EOF).
fn pull<R: BufRead>(
    groups: &mut PafGroups<R>,
    alninfo: &mut Box<dyn Write>,
    readinfo: &mut Box<dyn Write>,
) -> Result<Option<ReadInfoRow>> {
    // Fresh per-call tee handle (borrow ends when this returns).
    let mut sink: Option<&mut dyn Write> = Some(alninfo.as_mut());
    match groups.next_group_tee(&mut sink)? {
        None => Ok(None),
        Some(mut group) => {
            let ri = collapse_group(&mut group);
            ri.write(readinfo.as_mut())?;
            Ok(Some(ri))
        }
    }
}

/// The fused pass: lock-step over the two sorted PAFs. Writes the comparison
/// table to `compare_out`, and (when the corresponding path is `Some`) the
/// per-set alninfo / readinfo tables. Returns (matched, a_only, b_only).
#[allow(clippy::too_many_arguments)]
fn compare_sorted_pafs(
    a_sorted: &str,
    b_sorted: &str,
    label_a: &str,
    label_b: &str,
    compare_out: &str,
    readinfo_a: Option<&str>,
    readinfo_b: Option<&str>,
    alninfo_a: Option<&str>,
    alninfo_b: Option<&str>,
    junctions: bool,
    allow_id_mismatch: bool,
    summary: &mut CompareSummary,
) -> Result<(u64, u64, u64)> {
    // Comparison output + header.
    let mut out = open_output(Some(compare_out))?;
    if junctions {
        write_compare_junctions_header(&mut out, label_a, label_b)?;
    } else {
        write_compare_header(&mut out, label_a, label_b)?;
    }

    // Per-set side outputs. A suppressed table writes to `io::sink()` (no file is
    // created and the bytes are discarded) — this keeps every writer a concrete
    // `Box<dyn Write>`, avoiding the `Option<&mut dyn Write>` lifetime pitfalls.
    let open_or_sink = |p: Option<&str>| -> Result<Box<dyn Write>> {
        Ok(match p {
            Some(path) => open_output(Some(path))?,
            None => Box::new(io::sink()),
        })
    };
    let mut ri_a = open_or_sink(readinfo_a)?;
    let mut ri_b = open_or_sink(readinfo_b)?;
    let mut al_a = open_or_sink(alninfo_a)?;
    let mut al_b = open_or_sink(alninfo_b)?;
    writeln!(ri_a, "{READINFO_HEADER}")?;
    writeln!(ri_b, "{READINFO_HEADER}")?;
    AlnInfo::write_header(al_a.as_mut())?;
    AlnInfo::write_header(al_b.as_mut())?;

    let mut groups_a = PafGroups::new(open_input(a_sorted)?, /* warn_unsorted = */ false);
    let mut groups_b = PafGroups::new(open_input(b_sorted)?, false);

    let header_cols: Vec<&str> = READINFO_HEADER.split('\t').collect();

    // The lock-step merge runs in a helper that borrows each `Box<dyn Write>` only
    // for the call (the v0.9.0 pattern), so the writers are free to flush after.
    let counts = run_merge(
        &mut groups_a,
        &mut groups_b,
        &mut out,
        &mut al_a,
        &mut ri_a,
        &mut al_b,
        &mut ri_b,
        &header_cols,
        junctions,
        allow_id_mismatch,
        summary,
    )?;

    out.flush()?;
    ri_a.flush()?;
    ri_b.flush()?;
    al_a.flush()?;
    al_b.flush()?;
    Ok(counts)
}

/// The lock-step merge of two sorted PAFs. Each `Box<dyn Write>` is borrowed only
/// for the call duration (the owning boxes live in the caller), so they're free
/// to flush afterward. Suppressed side outputs are `io::sink()` boxes.
#[allow(clippy::too_many_arguments)]
fn run_merge<R: BufRead>(
    groups_a: &mut PafGroups<R>,
    groups_b: &mut PafGroups<R>,
    out: &mut Box<dyn Write>,
    al_a: &mut Box<dyn Write>,
    ri_a: &mut Box<dyn Write>,
    al_b: &mut Box<dyn Write>,
    ri_b: &mut Box<dyn Write>,
    header_cols: &[&str],
    junctions: bool,
    allow_id_mismatch: bool,
    summary: &mut CompareSummary,
) -> Result<(u64, u64, u64)> {
    let mut pending_a = pull(groups_a, al_a, ri_a)?;
    let mut pending_b = pull(groups_b, al_b, ri_b)?;
    let mut n_matched: u64 = 0;
    let mut n_a_only: u64 = 0;
    let mut n_b_only: u64 = 0;

    loop {
        match (pending_a.take(), pending_b.take()) {
            (None, None) => break,

            (Some(ra), None) => {
                if !allow_id_mismatch {
                    bail!(
                        "PAF A has more reads than PAF B (B exhausted after {n_matched} \
                         matched; next unmatched A read is {:?}).",
                        ra.read_name
                    );
                }
                n_a_only += 1; // ra was already pulled (tables written)
                summary.note_a_only_id();
                while pull(groups_a, al_a, ri_a)?.is_some() {
                    n_a_only += 1;
                    summary.note_a_only_id();
                }
                break;
            }
            (None, Some(rb)) => {
                if !allow_id_mismatch {
                    bail!(
                        "PAF B has more reads than PAF A (A exhausted after {n_matched} \
                         matched; next unmatched B read is {:?}).",
                        rb.read_name
                    );
                }
                n_b_only += 1;
                summary.note_b_only_id();
                while pull(groups_b, al_b, ri_b)?.is_some() {
                    n_b_only += 1;
                    summary.note_b_only_id();
                }
                break;
            }

            (Some(ra), Some(rb)) => {
                if ra.read_name == rb.read_name {
                    let line_a = readinfo_line(&ra)?;
                    let line_b = readinfo_line(&rb)?;
                    let map_a: HashMap<&str, &str> =
                        header_cols.iter().copied().zip(line_a.split('\t')).collect();
                    let map_b: HashMap<&str, &str> =
                        header_cols.iter().copied().zip(line_b.split('\t')).collect();

                    // One pair of accessor closures, shared by the summary classifier
                    // and the row emitter (closures are Copy — they capture &map_*).
                    let get_a = |c: &str| *map_a.get(c).unwrap_or(&"");
                    let get_b = |c: &str| *map_b.get(c).unwrap_or(&"");
                    summary.observe(&classify(&get_a, &get_b));

                    if junctions {
                        emit_compare_junctions_row(out, &ra.read_name, ra.read_len, get_a, get_b)?;
                    } else {
                        emit_compare_row(out, &ra.read_name, ra.read_len, get_a, get_b)?;
                    }

                    n_matched += 1;
                    if n_matched % 100_000 == 0 {
                        eprintln!("[INFO]   compared {n_matched} reads...");
                    }
                    pending_a = pull(groups_a, al_a, ri_a)?;
                    pending_b = pull(groups_b, al_b, ri_b)?;
                } else if !allow_id_mismatch {
                    bail!(
                        "read-name mismatch at read #{}: A has {:?} but B has {:?}.",
                        n_matched + 1,
                        ra.read_name,
                        rb.read_name
                    );
                } else if ra.read_name < rb.read_name {
                    n_a_only += 1;
                    summary.note_a_only_id();
                    pending_b = Some(rb); // keep B; advance A
                    pending_a = pull(groups_a, al_a, ri_a)?;
                } else {
                    n_b_only += 1;
                    summary.note_b_only_id();
                    pending_a = Some(ra); // keep A; advance B
                    pending_b = pull(groups_b, al_b, ri_b)?;
                }
            }
        }
    }

    Ok((n_matched, n_a_only, n_b_only))
}
