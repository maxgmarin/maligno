//! `query-junction-diff` — from a `compare` / `compare-toolkit merge-readinfo`
//! table, reconstruct each differing read's splice junctions on both sides
//! (A/B), paired correctly in query and genomic coordinate space, and report
//! which junctions are unsupported by the other side.
//!
//! **Why this re-derives junctions from each side's `cs` tag instead of
//! trusting the already-stored `junctions`/`genomic_junctions` columns**: in
//! `record.rs`'s `AlnInfo::from_paf`, a `-`-strand alignment's `junctions` is
//! flipped and sorted into ascending query order, but `genomic_junctions`
//! (built from the same cs-walk data) is never reordered to match. Every
//! individual value in both stored columns is correct — only the positional
//! correspondence between the two lists breaks for `-`-strand alignments.
//! Nothing else in the crate reads those two columns index-paired (every
//! existing consumer treats each as an independent, order-insensitive set),
//! so this has never affected any existing output. This command is the first
//! consumer that needs the pairing, so it reconstructs both lists itself,
//! from `cs`, with a single stable sort permutation applied to the query and
//! genomic vectors together — self-contained here, no change to `record.rs`
//! or any existing table's schema.
//!
//! **Why "query-junction-diff"**: read selection is anchored at the
//! **query** coordinate space — a read counts as differing based on its
//! query-space junction sets (`N_Junctions_OnlyA`/`N_Junctions_OnlyB`), not a
//! genomic-space comparison. The genomic coordinates reported alongside each
//! junction are the reconstructed pairing, not an independent selection
//! criterion.
//!
//! A read counts as "differing" (and gets full junction reconstruction) iff:
//!   - mapped on both sides, and its query-space junction sets differ
//!     (`N_Junctions_OnlyA`/`N_Junctions_OnlyB` not both `0`); or
//!   - mapped on only one side, and that side has at least one splice
//!     junction (`JuncCount > 0` — a junction with nothing to compare against
//!     is, by definition, unsupported by the other side).
//! Reads unmapped on both sides, both-mapped reads with identical query
//! junctions, and only-one-side reads with zero junctions are excluded (but
//! still tallied in the summary).
//!
//! Outputs (to `--outdir`, `--prefix`-named):
//!   1. `{prefix}.query_junction_diff.summary.tsv`                    — parse-time
//!      funnel counts (always uncompressed).
//!   2. `{prefix}.per_read_query_junction_diff.summary.tsv.gz`        — one row per
//!      reconstructed junction per side per differing read (always gzipped).
//!   3. `{prefix}.query_junction_diff_unmatched.A.tsv.gz`             — distinct
//!      genomic junctions called in A that were never matched in B, with how
//!      many reads support each, plus how many reads *total* (across the
//!      whole table) carry that junction on side A (always gzipped).
//!   4. `{prefix}.query_junction_diff_unmatched.B.tsv.gz`             — same,
//!      for B.
//!
//! **Two-pass design.** Pass 1 (above) reconstructs junctions only for the
//! small "differing" subset, identifying which junctions are unsupported.
//! Pass 2 re-scans the *entire* comparison table a second time, narrowly
//! column-projected, tallying total occurrences of exactly those
//! already-flagged junctions — a cheap pass needing only genomic-space set
//! membership (`genomic_junctions_A/B`), no `cs`-tag reconstruction, since it
//! never needs to know a junction's query position. This is why `--input`
//! must be Parquet: pass 2 re-reads the same file, which needs a real,
//! seekable, column-project-able source (not `stdin`, and much cheaper than
//! a second full TSV parse).

mod reconstruct;
mod rollup;
mod summary;

use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::compare_summary::{labels_from_row, require_ab_schema};
use crate::io_utils::open_output;
use crate::table_input::{open_table, InputFormat};

use reconstruct::{build_junction_records, build_side_junctions, side_position_sets, write_per_read_row};
use rollup::{accumulate_side_total, accumulate_unmatched, write_unmatched_table, TotalsAcc, UnmatchedAcc};
use summary::Summary;

/// A `TargetChr` value indicates an unmapped read when it is empty or `"*"`
/// — same convention as `compare_summary::is_mapped`, duplicated here (a
/// one-line private copy) rather than changing that function's visibility.
#[inline]
fn is_mapped(target_chr: &str) -> bool {
    !(target_chr.is_empty() || target_chr == "*")
}

/// Reconstruct one read's differing splice junctions on both sides. Comparison
/// table → `compare.tsv`/`.parquet` schema — not the original PAFs.
#[derive(clap::Args, Debug)]
pub struct QueryJunctionDiffArgs {
    /// Comparison table from `compare` / `compare-toolkit merge-readinfo` —
    /// Parquet only (`.parquet`). This command runs a second pass over the
    /// same file (see the module doc comment), so there is no `--input-format`
    /// choice and no TSV/stdin support.
    #[arg(short = 'i', long = "input", value_name = "compare.parquet")]
    input: String,

    /// Output directory (created if it does not exist).
    #[arg(long = "outdir", value_name = "DIR")]
    outdir: String,

    /// Filename prefix for all outputs.
    #[arg(long = "prefix", value_name = "STR")]
    prefix: String,
}

const NEEDED: [&str; 7] = ["cs", "Strand", "Query_Start", "Query_End", "Target_Start", "TargetChr", "JuncCount"];

pub fn run(args: &QueryJunctionDiffArgs) -> Result<()> {
    if !args.input.ends_with(".parquet") {
        bail!(
            "--input must be a Parquet file ('.parquet') — got '{}'. \
             query-junction-diff doesn't support TSV input.",
            args.input
        );
    }

    let outdir = Path::new(&args.outdir);
    fs::create_dir_all(outdir).with_context(|| format!("cannot create --outdir '{}'", args.outdir))?;

    let path = |name: String| outdir.join(name).to_string_lossy().into_owned();
    let summary_out = path(format!("{}.query_junction_diff.summary.tsv", args.prefix));
    let per_read_out = path(format!("{}.per_read_query_junction_diff.summary.tsv.gz", args.prefix));
    let unmatched_a_out = path(format!("{}.query_junction_diff_unmatched.A.tsv.gz", args.prefix));
    let unmatched_b_out = path(format!("{}.query_junction_diff_unmatched.B.tsv.gz", args.prefix));

    let wanted: Vec<String> = ["Read_Name", "Read_Len", "Label_A", "Label_B", "N_Junctions_OnlyA", "N_Junctions_OnlyB"]
        .into_iter()
        .map(String::from)
        .chain(["A", "B"].iter().flat_map(|side| NEEDED.iter().map(move |b| format!("{b}_{side}"))))
        .collect();
    let wanted_refs: Vec<&str> = wanted.iter().map(String::as_str).collect();
    let mut source = open_table(&args.input, InputFormat::Parquet, Some(&wanted_refs))?;

    let cols_owned = source.header().with_context(|| format!("reading comparison table '{}'", args.input))?;
    let cols: Vec<&str> = cols_owned.iter().map(String::as_str).collect();
    let col_index: HashMap<&str, usize> = cols.iter().copied().enumerate().map(|(i, c)| (c, i)).collect();
    require_ab_schema(&cols)?;

    let read_name_idx = *col_index.get("Read_Name").context("comparison table is missing column 'Read_Name'")?;
    let read_len_idx = *col_index.get("Read_Len").context("comparison table is missing column 'Read_Len'")?;
    let n_only_a_idx =
        *col_index.get("N_Junctions_OnlyA").context("comparison table is missing column 'N_Junctions_OnlyA'")?;
    let n_only_b_idx =
        *col_index.get("N_Junctions_OnlyB").context("comparison table is missing column 'N_Junctions_OnlyB'")?;

    let resolve = |side: &str| -> Result<HashMap<&'static str, usize>> {
        let mut m = HashMap::new();
        for base in NEEDED {
            let name = format!("{base}_{side}");
            let idx =
                *col_index.get(name.as_str()).with_context(|| format!("comparison table is missing column '{name}'"))?;
            m.insert(base, idx);
        }
        Ok(m)
    };
    let idx_a = resolve("A")?;
    let idx_b = resolve("B")?;

    let mut label_a = "A".to_string();
    let mut label_b = "B".to_string();
    let mut seen_row = false;

    let mut summary = Summary::default();
    let mut per_read_w = open_output(Some(&per_read_out))?;
    writeln!(
        per_read_w,
        "Read_Name\tside\tjunction_index\tstrand\tquery_pos\tchrom\tgenomic_start\tgenomic_end\t\
         matched_in_query\tmatched_in_genomic\tother_side_aligned"
    )?;
    let mut unmatched_acc_a: UnmatchedAcc = HashMap::new();
    let mut unmatched_acc_b: UnmatchedAcc = HashMap::new();

    while let Some(fields_owned) = source.next_row()? {
        let fields: Vec<&str> = fields_owned.iter().map(String::as_str).collect();
        if !seen_row {
            (label_a, label_b) = labels_from_row(&col_index, &fields);
            seen_row = true;
        }
        let get_a = |c: &str| -> &str { idx_a.get(c).and_then(|&i| fields.get(i)).copied().unwrap_or("") };
        let get_b = |c: &str| -> &str { idx_b.get(c).and_then(|&i| fields.get(i)).copied().unwrap_or("") };
        let read_name = fields.get(read_name_idx).copied().unwrap_or("");
        let read_len: u64 = fields.get(read_len_idx).copied().unwrap_or("0").parse().unwrap_or(0);
        let n_only_a: u64 = fields.get(n_only_a_idx).copied().unwrap_or("0").parse().unwrap_or(0);
        let n_only_b: u64 = fields.get(n_only_b_idx).copied().unwrap_or("0").parse().unwrap_or(0);

        let mapped_a = is_mapped(get_a("TargetChr"));
        let mapped_b = is_mapped(get_b("TargetChr"));
        let junc_count_a: u64 = get_a("JuncCount").parse().unwrap_or(0);
        let junc_count_b: u64 = get_b("JuncCount").parse().unwrap_or(0);

        let differing =
            summary.observe_row(mapped_a, mapped_b, n_only_a, n_only_b, junc_count_a, junc_count_b);
        if !differing {
            continue;
        }

        let pairs_a = if mapped_a {
            build_side_junctions(
                get_a("cs"),
                get_a("Strand").chars().next().unwrap_or('.'),
                get_a("Query_Start").parse().unwrap_or(0),
                get_a("Query_End").parse().unwrap_or(0),
                get_a("Target_Start").parse().unwrap_or(0),
                get_a("TargetChr"),
                read_len,
            )
        } else {
            Vec::new()
        };
        let pairs_b = if mapped_b {
            build_side_junctions(
                get_b("cs"),
                get_b("Strand").chars().next().unwrap_or('.'),
                get_b("Query_Start").parse().unwrap_or(0),
                get_b("Query_End").parse().unwrap_or(0),
                get_b("Target_Start").parse().unwrap_or(0),
                get_b("TargetChr"),
                read_len,
            )
        } else {
            Vec::new()
        };

        let (qset_a, gset_a) = side_position_sets(&pairs_a);
        let (qset_b, gset_b) = side_position_sets(&pairs_b);

        let strand_a = get_a("Strand").chars().next().unwrap_or('.');
        let strand_b = get_b("Strand").chars().next().unwrap_or('.');
        let recs_a = build_junction_records(pairs_a, strand_a, &qset_b, &gset_b, mapped_b);
        let recs_b = build_junction_records(pairs_b, strand_b, &qset_a, &gset_a, mapped_a);

        for r in &recs_a {
            write_per_read_row(&mut per_read_w, read_name, 'A', r)?;
        }
        for r in &recs_b {
            write_per_read_row(&mut per_read_w, read_name, 'B', r)?;
        }

        accumulate_unmatched(&mut unmatched_acc_a, &recs_a);
        accumulate_unmatched(&mut unmatched_acc_b, &recs_b);
    }
    per_read_w.flush()?;

    // ── Pass 2: total occurrence counts ─────────────────────────────────────
    // Seed the "junctions of interest" from pass 1's own unmatched rollups —
    // no need to write/re-read anything, they're already in memory. Then
    // re-scan the whole table once more, narrowly column-projected, tallying
    // how many reads *total* carry each already-flagged junction.
    let mut totals_a: TotalsAcc = unmatched_acc_a.keys().map(|k| (k.clone(), 0u64)).collect();
    let mut totals_b: TotalsAcc = unmatched_acc_b.keys().map(|k| (k.clone(), 0u64)).collect();

    const TOTALS_WANTED: [&str; 6] =
        ["TargetChr_A", "TargetChr_B", "genomic_junctions_A", "genomic_junctions_B", "Strand_A", "Strand_B"];
    let mut totals_source = open_table(&args.input, InputFormat::Parquet, Some(&TOTALS_WANTED))
        .with_context(|| format!("re-reading comparison table '{}' for pass 2", args.input))?;
    let totals_cols_owned = totals_source.header()?;
    let totals_col_index: HashMap<&str, usize> =
        totals_cols_owned.iter().map(String::as_str).enumerate().map(|(i, c)| (c, i)).collect();
    let totals_idx = |name: &str| -> Result<usize> {
        totals_col_index
            .get(name)
            .copied()
            .with_context(|| format!("comparison table is missing column '{name}'"))
    };
    let chrom_a_idx = totals_idx("TargetChr_A")?;
    let chrom_b_idx = totals_idx("TargetChr_B")?;
    let gj_a_idx = totals_idx("genomic_junctions_A")?;
    let gj_b_idx = totals_idx("genomic_junctions_B")?;
    let strand_a_idx = totals_idx("Strand_A")?;
    let strand_b_idx = totals_idx("Strand_B")?;

    while let Some(fields) = totals_source.next_row()? {
        let get = |i: usize| -> &str { fields.get(i).map(String::as_str).unwrap_or("") };

        let strand_a = get(strand_a_idx).chars().next().unwrap_or('.');
        accumulate_side_total(get(chrom_a_idx), strand_a, get(gj_a_idx), &mut totals_a);

        let strand_b = get(strand_b_idx).chars().next().unwrap_or('.');
        accumulate_side_total(get(chrom_b_idx), strand_b, get(gj_b_idx), &mut totals_b);
    }

    write_unmatched_table(&unmatched_a_out, &unmatched_acc_a, &totals_a)?;
    write_unmatched_table(&unmatched_b_out, &unmatched_acc_b, &totals_b)?;
    summary.write_tsv(&summary_out, &label_a, &label_b)?;

    eprintln!(
        "Query junction diff (A={label_a}, B={label_b}): {} of {} reads had differing query junctions",
        summary.n_query_junctions_different(),
        summary.n_total()
    );
    eprintln!("Outputs in {}:", args.outdir);
    eprintln!("  {summary_out}");
    eprintln!("  {per_read_out}");
    eprintln!("  {unmatched_a_out}");
    eprintln!("  {unmatched_b_out}");
    Ok(())
}
