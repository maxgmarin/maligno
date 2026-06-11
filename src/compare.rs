//! The primary `compare` command: an on-rails pipeline that takes two PAFs and
//! produces, in one invocation, the per-set alninfo + readinfo tables **and** the
//! per-read comparison table.
//!
//! It owns its preconditions rather than trusting the user:
//!   1. it **sorts** both inputs by `Query_Name` (identical deterministic rule),
//!      so grouping and matching order are guaranteed, then
//!   2. it **verifies** the two PAFs carry the same `Query_Name` set (O(1) check),
//!      erroring by default if they differ, then
//!   3. builds alninfo + readinfo for each set, and
//!   4. runs the `compare-readinfo` merge-join to emit the comparison table.
//!
//! This is the porcelain over the plumbing subcommands (`paf2tables`,
//! `compare-readinfo`, …); it's equivalent to running them by hand on sorted
//! inputs, but as one safe command.
//!
//! Precondition (documented, not enforced): a `Query_Name` uniquely identifies a
//! single read/sequence — so sorting by name alone (no `Read_Len` secondary key)
//! is sufficient for the downstream `(Read_Name, Read_Len)` merge-join.

use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::compare_streaming::{self, CompareMode, CompareReadinfoArgs};
use crate::external_sort::{parse_mem, read_id_set_check, sort_paf_to_file};
use crate::io_utils::{open_input, open_output};
use crate::paf2tables::stream_readinfo;

#[derive(clap::Args, Debug)]
pub struct CompareArgs {
    /// PAF for dataset A. Use '-' for stdin; '.gz' auto-decompressed.
    #[arg(short = 'a', long = "paf-a", value_name = "a.paf")]
    paf_a: String,

    /// PAF for dataset B. Use '-' for stdin; '.gz' auto-decompressed.
    #[arg(short = 'b', long = "paf-b", value_name = "b.paf")]
    paf_b: String,

    /// Label for dataset A (names the per-set outputs and the column suffix).
    #[arg(long = "label-a", value_name = "LABEL", default_value = "SetA")]
    label_a: String,

    /// Label for dataset B.
    #[arg(long = "label-b", value_name = "LABEL", default_value = "SetB")]
    label_b: String,

    /// Output directory (created if missing). All artifacts are written here.
    #[arg(long = "outdir", value_name = "DIR")]
    outdir: String,

    /// Filename prefix for all outputs.
    #[arg(short = 'p', long = "prefix", value_name = "NAME")]
    prefix: String,

    /// Comparison view: `full` (94 cols) or `junctions` (47-col splice view).
    #[arg(long = "mode", value_enum, default_value_t = CompareMode::Full)]
    mode: CompareMode,

    /// In-memory sort buffer per file (K/M/G suffix, or plain bytes).
    #[arg(long = "mem", value_name = "SIZE", default_value = "1G")]
    mem: String,

    /// Temp directory for sort spill files (default: --outdir).
    #[arg(long = "tmp-dir", value_name = "DIR")]
    tmp_dir: Option<String>,

    /// Number of sort threads (default: ext-sort's default).
    #[arg(long = "threads", value_name = "N")]
    threads: Option<usize>,

    /// Compare the shared intersection instead of erroring when the two PAFs do
    /// not carry the exact same Query_Name set.
    #[arg(long = "allow-id-mismatch")]
    allow_id_mismatch: bool,

    /// Keep the intermediate sorted PAFs instead of deleting them at the end.
    #[arg(long = "keep-sorted")]
    keep_sorted: bool,
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
    let mem = parse_mem(&args.mem)?;
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

    // ── Step 1: sort both PAFs by Query_Name (consistent rule) ────────────────
    eprintln!(
        "[INFO] Step 1/4 — sorting both PAFs by Query_Name (mem={} bytes, tmp={})",
        mem,
        tmp_dir.display()
    );
    sort_paf_to_file(&args.paf_a, &a_sorted, mem, &tmp_dir, args.threads)
        .with_context(|| format!("sorting PAF A ({})", args.paf_a))?;
    sort_paf_to_file(&args.paf_b, &b_sorted, mem, &tmp_dir, args.threads)
        .with_context(|| format!("sorting PAF B ({})", args.paf_b))?;

    // ── Step 2: read-ID set-equality check (O(1) memory) ──────────────────────
    eprintln!("[INFO] Step 2/4 — verifying the two PAFs share the same read-ID set...");
    let chk = read_id_set_check(&a_sorted, &b_sorted, 5)?;
    eprintln!(
        "  shared: {}   only in {}: {}   only in {}: {}",
        chk.shared, args.label_a, chk.only_a, args.label_b, chk.only_b
    );
    if chk.only_a > 0 || chk.only_b > 0 {
        if !args.allow_id_mismatch {
            if !args.keep_sorted {
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

    // ── Step 3: per-set alninfo + readinfo tables ─────────────────────────────
    eprintln!("[INFO] Step 3/4 — building alninfo + readinfo for each set...");
    for (sorted, alninfo, readinfo) in [
        (&a_sorted, &a_alninfo, &a_readinfo),
        (&b_sorted, &b_alninfo, &b_readinfo),
    ] {
        let reader = open_input(sorted)?;
        let mut alninfo_w = open_output(Some(alninfo))?;
        let mut readinfo_w = open_output(Some(readinfo))?;
        stream_readinfo(
            reader,
            &mut *readinfo_w,
            Some(alninfo_w.as_mut()),
            /* strict_grouping = */ false,
        )?;
        alninfo_w.flush()?;
        readinfo_w.flush()?;
    }

    // ── Step 4: compare the two readinfo tables ───────────────────────────────
    eprintln!("[INFO] Step 4/4 — comparing readinfo tables ({})...", compare_out);
    let cmp_args = CompareReadinfoArgs {
        readinfo_a: a_readinfo.clone(),
        readinfo_b: b_readinfo.clone(),
        label_a: args.label_a.clone(),
        label_b: args.label_b.clone(),
        output: compare_out.clone(),
        ignore_row_mismatch: args.allow_id_mismatch,
        mode: args.mode.clone(),
    };
    compare_streaming::run(&cmp_args)?;

    // ── Step 5: cleanup + summary ─────────────────────────────────────────────
    if !args.keep_sorted {
        let _ = fs::remove_file(&a_sorted);
        let _ = fs::remove_file(&b_sorted);
    }
    eprintln!("[INFO] Done. Outputs in {}:", args.outdir);
    eprintln!("  {a_alninfo}");
    eprintln!("  {b_alninfo}");
    eprintln!("  {a_readinfo}");
    eprintln!("  {b_readinfo}");
    eprintln!("  {compare_out}");
    if args.keep_sorted {
        eprintln!("  {a_sorted}");
        eprintln!("  {b_sorted}");
    }
    Ok(())
}
