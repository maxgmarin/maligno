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
//!      alninfo + readinfo tables are tee'd out as side outputs as it goes, when
//!      requested (`--emit-alninfo` / `--emit-readinfo`; off by default).
//!
//! This is the porcelain over the `compare-toolkit` plumbing subcommands
//! (`paf2tables`, `merge-readinfo`, …): the comparison table is byte-identical
//! to running `compare-toolkit merge-readinfo` on the sorted readinfo files,
//! and the side outputs are byte-identical to `compare-toolkit paf2tables` on
//! the sorted PAFs.
//!
//! By default the same pass also drives `find-aln-diff`'s core (via
//! `find_query_diff::AlnDiffAccumulator`) at its default settings (`--space
//! query --compare-by all`), so `compare` additionally emits the differing
//! reads + region tables without a second read of its own output table —
//! byte-identical to running standalone `compare-toolkit find-aln-diff`
//! against the emitted comparison table. `--skip-find-aln-diff` opts out;
//! other `--space`/`--compare-by` combinations still require the standalone
//! command.
//!
//! Precondition (documented, not enforced): a `Query_Name` uniquely identifies a
//! single read/sequence — so sorting by name alone (no `Read_Len` secondary key)
//! is sufficient for the downstream `(Read_Name, Read_Len)` merge-join.

use std::collections::HashMap;
use std::fs;
use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::compare_streaming::validate_set_label;
use crate::comparison_row::{write_compare_header, ComparisonRow};
use crate::parquet_out::{ComparisonParquetWriter, OutputFormat};
use crate::compare_summary::{classify, CompareSummary};
use crate::external_sort::{parse_mem, read_id_set_check, sort_paf_to_file};
use crate::find_query_diff::{AlnDiffAccumulator, AlnDiffStats, CompareBy, DiffSpace};
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

    /// Name for dataset A, recorded in the comparison table's `Label_A` column.
    #[arg(long = "label-a", value_name = "LABEL", default_value = "SetA",
          value_parser = validate_set_label)]
    label_a: String,

    /// Name for dataset B, recorded in the comparison table's `Label_B` column.
    #[arg(long = "label-b", value_name = "LABEL", default_value = "SetB",
          value_parser = validate_set_label)]
    label_b: String,

    /// Output directory.
    #[arg(short = 'o', long = "outdir", value_name = "DIR")]
    outdir: String,

    /// Filename prefix for all outputs.
    #[arg(short = 'p', long = "prefix", value_name = "NAME")]
    prefix: String,

    /// Output format for the comparison table.
    #[arg(long = "format", value_enum, default_value_t = OutputFormat::Both)]
    format: OutputFormat,

    /// Compare the shared intersection of aligned sequences instead of erroring when the two PAFs do
    /// not carry the exact same "Query_Name" set.
    #[arg(long = "allow-id-mismatch")]
    allow_id_mismatch: bool,

    /// Skip the internal sort: assume both PAFs already contain the same reads,
    /// grouped by "Query_Name" and in the same relative order.
    /// Not combinable with --allow-id-mismatch or --keep-sorted-paf.
    #[arg(long = "presorted", conflicts_with_all = ["allow_id_mismatch", "keep_sorted_paf"])]
    presorted: bool,

    /// Write the per-set alninfo (35-col) tables (off by default).
    #[arg(long = "emit-alninfo")]
    emit_alninfo: bool,

    /// Write the per-set readinfo (33-col) tables (off by default).
    #[arg(long = "emit-readinfo")]
    emit_readinfo: bool,

    /// In-memory sort buffer per file (K/M/G suffix, or plain bytes).
    #[arg(long = "sort-mem", value_name = "SIZE", default_value = "1G")]
    sort_mem: String,

    /// Temp directory for temp out of memory sort files (default: --outdir).
    #[arg(long = "sort-tmp-dir", value_name = "DIR")]
    tmp_dir: Option<String>,

    /// Number of sort threads (default: 1).
    #[arg(long = "sort-threads", value_name = "N", default_value_t = 1)]
    sort_threads: usize,

    /// Keep the intermediate sorted PAFs instead of deleting them at the end.
    #[arg(long = "keep-sorted-paf")]
    keep_sorted_paf: bool,

    /// Skip the differing-reads + genomic region tables that `compare` writes
    /// by default.
    #[arg(long = "skip-find-aln-diff")]
    skip_find_aln_diff: bool,
}

pub fn run(args: &CompareArgs) -> Result<()> {
    // ── Step 0: setup ─────────────────────────────────────────────────────────
    // Distinct labels are required: they name the per-set output files, so equal
    // labels would silently overwrite A's alninfo/readinfo with B's.
    if args.label_a == args.label_b {
        bail!(
            "--label-a and --label-b are both '{}' — they must differ (they name \
             the per-set output files)",
            args.label_a
        );
    }
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

    let path = |name: String| outdir.join(name).to_string_lossy().into_owned();
    let a_sorted = path(format!("{}.{}.sorted.paf.gz", args.prefix, args.label_a));
    let b_sorted = path(format!("{}.{}.sorted.paf.gz", args.prefix, args.label_b));
    let a_alninfo = path(format!("{}.{}.alninfo.tsv.gz", args.prefix, args.label_a));
    let b_alninfo = path(format!("{}.{}.alninfo.tsv.gz", args.prefix, args.label_b));
    let a_readinfo = path(format!("{}.{}.readinfo.tsv.gz", args.prefix, args.label_a));
    let b_readinfo = path(format!("{}.{}.readinfo.tsv.gz", args.prefix, args.label_b));
    let compare_tsv = args
        .format
        .writes_tsv()
        .then(|| path(format!("{}.compare.tsv.gz", args.prefix)));
    let compare_parquet = args
        .format
        .writes_parquet()
        .then(|| path(format!("{}.compare.parquet", args.prefix)));
    let summary_out = path(format!("{}.compare.summary.tsv", args.prefix));
    // `find-aln-diff`'s default-mode (`--space query --compare-by all`) output
    // paths, fused into this same pass unless `--skip-find-aln-diff`. Same
    // naming convention as standalone `find-aln-diff`'s default run, so a
    // later standalone re-run against this table never collides.
    let diff_reads_out = path(format!("{}.query_diff_reads.tsv.gz", args.prefix));
    let diff_regions_a_out = path(format!("{}.query_diff_regions.A.bed.gz", args.prefix));
    let diff_regions_b_out = path(format!("{}.query_diff_regions.B.bed.gz", args.prefix));

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
        eprintln!("[INFO] comparing in one pass ({})...", describe_outputs(&compare_tsv, &compare_parquet));
    } else {
        eprintln!(
            "[INFO] Step 3/3 — comparing in one pass ({})...",
            describe_outputs(&compare_tsv, &compare_parquet)
        );
    }
    let mut summary = CompareSummary::default();
    let diff_acc = if args.skip_find_aln_diff {
        None
    } else {
        Some(AlnDiffAccumulator::new(&diff_reads_out, DiffSpace::Query, CompareBy::All)?)
    };
    let result = compare_sorted_pafs(
        &a_in,
        &b_in,
        &args.label_a,
        &args.label_b,
        compare_tsv.as_deref(),
        if args.emit_readinfo { Some(&a_readinfo) } else { None },
        if args.emit_readinfo { Some(&b_readinfo) } else { None },
        if args.emit_alninfo { Some(&a_alninfo) } else { None },
        if args.emit_alninfo { Some(&b_alninfo) } else { None },
        compare_parquet.as_deref(),
        args.allow_id_mismatch,
        &mut summary,
        diff_acc,
        &diff_regions_a_out,
        &diff_regions_b_out,
        &describe_outputs(&compare_tsv, &compare_parquet),
    );
    let (counts, diff_stats) = match result {
        Ok(v) => v,
        Err(e) => {
            // The compare pass can fail partway (e.g. --presorted inputs that are
            // not actually in the same order), having already written part of the
            // output. Remove the partial artifacts so the failure leaves nothing
            // half-written, then surface the error (with a hint under --presorted).
            for p in [compare_tsv.as_deref(), compare_parquet.as_deref()].into_iter().flatten() {
                let _ = fs::remove_file(p);
            }
            if args.emit_alninfo {
                let _ = fs::remove_file(&a_alninfo);
                let _ = fs::remove_file(&b_alninfo);
            }
            if args.emit_readinfo {
                let _ = fs::remove_file(&a_readinfo);
                let _ = fs::remove_file(&b_readinfo);
            }
            if !args.skip_find_aln_diff {
                let _ = fs::remove_file(&diff_reads_out);
                let _ = fs::remove_file(&diff_regions_a_out);
                let _ = fs::remove_file(&diff_regions_b_out);
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
    };
    let _ = counts;

    // ── cleanup + summary ─────────────────────────────────────────────────────
    if !args.presorted && !args.keep_sorted_paf {
        let _ = fs::remove_file(&a_sorted);
        let _ = fs::remove_file(&b_sorted);
    }
    // Aggregate summary statistics → sidecar TSV + stderr block.
    summary.write_tsv(&summary_out, &args.label_a, &args.label_b, &[])?;
    summary.render_stderr(&args.label_a, &args.label_b, &[]);
    eprintln!("Outputs in {}:", args.outdir);
    if args.emit_alninfo {
        eprintln!("  {a_alninfo}");
        eprintln!("  {b_alninfo}");
    }
    if args.emit_readinfo {
        eprintln!("  {a_readinfo}");
        eprintln!("  {b_readinfo}");
    }
    for p in [compare_tsv.as_deref(), compare_parquet.as_deref()].into_iter().flatten() {
        eprintln!("  {p}");
    }
    eprintln!("  {summary_out}");
    if let Some(stats) = &diff_stats {
        eprintln!("  {diff_reads_out}");
        eprintln!("  {diff_regions_a_out}");
        eprintln!("  {diff_regions_b_out}");
        eprintln!(
            "  ({} differing reads found (space=query, compare-by=all); \
             see compare.summary.tsv for the full tally)",
            stats.n_diff_rows
        );
        if stats.n_bad_interval > 0 {
            eprintln!(
                "  ({} intervals skipped: unparseable or degenerate coordinates)",
                stats.n_bad_interval
            );
        }
    }
    if args.keep_sorted_paf {
        eprintln!("  {a_sorted}");
        eprintln!("  {b_sorted}");
    }

    // The fused default (space=query, compare-by=all) output above covers the
    // most common case; point at the standalone command for the modes it
    // doesn't cover (or, under --skip-find-aln-diff, for that default mode too).
    if let Some(tsv) = &compare_tsv {
        eprintln!();
        if args.skip_find_aln_diff {
            eprintln!("For the differing reads and the genomic regions where they cluster:");
            eprintln!(
                "  maligno compare-toolkit find-aln-diff -i {tsv} --outdir {} --prefix {}",
                args.outdir, args.prefix
            );
        } else {
            eprintln!("For reference-space or junctions-based differences, run standalone:");
            eprintln!(
                "  maligno compare-toolkit find-aln-diff -i {tsv} --space reference --outdir {} --prefix {}",
                args.outdir, args.prefix
            );
        }
    }

    Ok(())
}

/// Render the comparison output path(s) for the progress line.
fn describe_outputs(tsv: &Option<String>, parquet: &Option<String>) -> String {
    let mut v: Vec<&str> = Vec::new();
    if let Some(p) = tsv {
        v.push(p);
    }
    if let Some(p) = parquet {
        v.push(p);
    }
    v.join(" + ")
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
/// per-set alninfo / readinfo tables. When `diff_acc` is `Some`, also drives
/// `find-aln-diff`'s default-mode core inline (differing reads + region
/// tables at `regions_a_out`/`regions_b_out`) — no re-read of the just-written
/// comparison table. Returns ((matched, a_only, b_only), diff stats if run).
#[allow(clippy::too_many_arguments)]
fn compare_sorted_pafs(
    a_sorted: &str,
    b_sorted: &str,
    label_a: &str,
    label_b: &str,
    compare_tsv: Option<&str>,
    readinfo_a: Option<&str>,
    readinfo_b: Option<&str>,
    alninfo_a: Option<&str>,
    alninfo_b: Option<&str>,
    parquet_out: Option<&str>,
    allow_id_mismatch: bool,
    summary: &mut CompareSummary,
    mut diff_acc: Option<AlnDiffAccumulator>,
    regions_a_out: &str,
    regions_b_out: &str,
    diff_source_desc: &str,
) -> Result<((u64, u64, u64), Option<AlnDiffStats>)> {
    // Comparison TSV + header. When the TSV is not requested the rows go to
    // `io::sink()`, the same way suppressed alninfo/readinfo outputs do, so no
    // empty file is created and `run_merge` needs no extra plumbing.
    let mut out: Box<dyn Write> = match compare_tsv {
        Some(p) => {
            let mut w = open_output(Some(p))?;
            write_compare_header(&mut w)?;
            w
        }
        None => Box::new(io::sink()),
    };

    // Optional Parquet output. Written to a plain file, never gzipped — Parquet
    // compresses per column internally.
    let mut parquet = match parquet_out {
        Some(p) => Some(ComparisonParquetWriter::new(
            fs::File::create(p).with_context(|| format!("cannot create '{p}'"))?,
        )?),
        None => None,
    };

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
        parquet.as_mut(),
        allow_id_mismatch,
        label_a,
        label_b,
        summary,
        &mut diff_acc,
    )?;

    let diff_stats = match diff_acc {
        Some(acc) => Some(acc.finish(regions_a_out, regions_b_out, label_a, label_b, diff_source_desc)?),
        None => None,
    };

    if let Some(pq) = parquet {
        pq.finish()?;
    }

    out.flush()?;
    ri_a.flush()?;
    ri_b.flush()?;
    al_a.flush()?;
    al_b.flush()?;
    Ok((counts, diff_stats))
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
    mut parquet: Option<&mut ComparisonParquetWriter<fs::File>>,
    allow_id_mismatch: bool,
    label_a: &str,
    label_b: &str,
    summary: &mut CompareSummary,
    diff_acc: &mut Option<AlnDiffAccumulator>,
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

                    // One pair of accessor closures, shared by the summary classifier,
                    // the row emitter, and (when active) the fused find-aln-diff core
                    // (closures are Copy — they capture &map_*).
                    let get_a = |c: &str| *map_a.get(c).unwrap_or(&"");
                    let get_b = |c: &str| *map_b.get(c).unwrap_or(&"");
                    let base = classify(&get_a, &get_b);
                    summary.observe(&base);
                    if let Some(acc) = diff_acc.as_mut() {
                        acc.observe_row(&ra.read_name, get_a, get_b, &base)?;
                    }

                    // One construction, both writers — so the TSV and the Parquet
                    // can never disagree about a row.
                    let row = ComparisonRow::build(
                        &ra.read_name, ra.read_len, label_a, label_b, get_a, get_b,
                    );
                    row.write_tsv_row(out)?;
                    if let Some(pq) = parquet.as_deref_mut() {
                        pq.append(&row)?;
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
