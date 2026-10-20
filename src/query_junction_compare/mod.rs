//! `query-junction-compare` — from a `compare` / `toolkit compare-readinfo`
//! table, reconstruct each differing read's splice junctions on both sides
//! (A/B), paired correctly in query and genomic coordinate space, and report
//! which junctions are unsupported by the other side.
//!
//! **Why this re-derives junctions from each side's `cs` tag instead of
//! trusting the already-stored `query_junctions`/`genomic_junctions` columns**: in
//! `record.rs`'s `AlnInfo::from_paf`, a `-`-strand alignment's `query_junctions` is
//! flipped and sorted into ascending query order, but `genomic_junctions`
//! (built from the same cs-walk data) is never reordered to match. Every
//! individual value in both stored columns is correct — only the positional
//! correspondence between the two lists breaks for `-`-strand alignments.
//! Every other consumer treats each column as an independent, order-insensitive
//! set. This command needs the pairing, so it reconstructs both lists itself,
//! from `cs`, with a single stable sort permutation applied to the query and
//! genomic vectors together.
//!
//! **Read selection is in query space**: a read counts as differing based on
//! its query-space junction sets (`n_query_junctions_only_A`/`n_query_junctions_only_B`),
//! not a genomic-space comparison. The genomic coordinates reported alongside
//! each junction are the reconstructed pairing, not an independent selection
//! criterion. Unlike `genomic-junction-compare`, which writes a row for every
//! read with a junction on either side, only the differing reads are
//! reconstructed.
//!
//! **Junction identity in reference coordinates** is the same as in
//! `genomic-junction-compare`: `(target_name, intron_start, intron_end,
//! intron_strand)`, where `intron_strand` is the side's `genomic_intron_strand`
//! (`*` when the aligner gave no transcript strand). It is used for
//! `matched_in_genomic`, the unmatched tables, their whole-table totals and the
//! summary's distinct-junction counts. `matched_in_query` compares read
//! positions only. Each junction also carries its splice-site motif from the
//! `cs` tag, written donor-acceptor in the junction's orientation (`gt-ag`);
//! the alignment's own strand is reported separately as `strand`.
//!
//! A read counts as "differing" (and gets full junction reconstruction) iff:
//!   - mapped on both sides, and its query-space junction sets differ
//!     (`n_query_junctions_only_A`/`n_query_junctions_only_B` not both `0`); or
//!   - mapped on only one side, and that side has at least one splice
//!     junction (`n_intron_events > 0` — a junction with nothing to compare against
//!     is, by definition, unsupported by the other side).
//! Reads unmapped on both sides, both-mapped reads with identical query
//! junctions, and only-one-side reads with zero junctions are excluded (but
//! still tallied in the summary).
//!
//! Outputs (to `--outdir`, `--prefix`-named):
//!   1. `{prefix}.query_junction_compare.summary.tsv`                    — read-selection
//!      funnel counts, then distinct-junction counts per side (unique, shared,
//!      only-A/only-B, unmatched) from pass 2 (always uncompressed).
//!   2. `{prefix}.query_junction_compare.per_read_per_junc_info.tsv.gz` — one row per
//!      reconstructed junction per side per differing read, with its intron
//!      strand and motif (always gzipped).
//!   3. `{prefix}.query_junction_compare.unmatched_junctions.A.tsv.gz`  — distinct
//!      genomic junctions called in A that were never matched in B, with their
//!      motif and how many reads support each, plus how many reads *total*
//!      (across the whole table) carry that junction on side A (always gzipped).
//!   4. `{prefix}.query_junction_compare.unmatched_junctions.B.tsv.gz`  — same,
//!      for B.
//!   5. With `--emit-all-junctions`, `{prefix}.query_junction_compare.all_junctions.{A,B}.tsv.gz`
//!      — every distinct junction on that side across the whole table, in the
//!      unmatched tables' columns (the unmatched table is its subset with
//!      `n_reads_query_junctions_different > 0`).
//!
//! **Two-pass design.** Pass 1 (above) reconstructs junctions only for the
//! small "differing" subset, identifying which junctions are unsupported.
//! Pass 2 re-scans the *entire* comparison table a second time, narrowly
//! column-projected, and adds every junction of every read (from `cs`, with
//! the same `junction::side_genomic_junctions` as `genomic-junction-compare`)
//! to one `JunctionCatalog` per side: whole-table read counts and motifs. It
//! never needs a junction's query position. This is why `--input` must be
//! Parquet: pass 2 re-reads the same file, which needs a real, seekable,
//! column-project-able source (not `stdin`, and much cheaper than a second
//! full TSV parse). The unmatched tables, the all-junctions tables and the
//! summary's distinct-junction counts all come from the two catalogs.

mod reconstruct;
mod rollup;
mod summary;

use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::columns::{
    header_line, Col, CS, GENOMIC_INTRON_STRAND, N_INTRON_EVENTS, LABEL_A, LABEL_B,
    N_QUERY_JUNCTIONS_ONLY_A, N_QUERY_JUNCTIONS_ONLY_B, QJC_PER_JUNCTION_COLS, QUERY_END,
    QUERY_START, QUERY_LEN, QUERY_NAME, SIDE_A, SIDE_B, STRAND, TARGET_NAME, TARGET_START,
};
use crate::compare_summary::{
    labels_from_row, require_ab_schema, require_col, wanted_columns, SideIndex,
};
use crate::io_utils::open_output;
use crate::junction::{normalize_intron_strand, side_genomic_junctions};
use crate::table_input::{open_table, InputFormat};

use reconstruct::{build_junction_records, build_side_junctions, side_position_sets, write_per_read_row};
use rollup::{unique_junction_counts, JunctionCatalog};
use summary::Summary;

/// A `target_name` value indicates an unmapped read when it is empty or `"*"`
/// — same convention as `compare_summary::is_mapped`, duplicated here (a
/// one-line private copy) rather than changing that function's visibility.
#[inline]
fn is_mapped(target_chr: &str) -> bool {
    !(target_chr.is_empty() || target_chr == "*")
}

/// Reconstruct each differing read's splice junctions on both sides, from a
/// comparison table (Parquet), not the original PAFs.
#[derive(clap::Args, Debug)]
pub struct QueryJunctionCompareArgs {
    /// Alignment comparison table in Parquet format
    #[arg(short = 'i', long = "input", value_name = "compare.parquet")]
    input: String,

    /// Output directory (created if it does not exist).
    #[arg(short = 'o', long = "outdir", value_name = "DIR")]
    outdir: String,

    /// Filename prefix for all outputs.
    #[arg(short = 'p', long = "prefix", value_name = "STR")]
    prefix: String,

    /// Also write every distinct junction of each side, across the whole table,
    /// to `{prefix}.query_junction_compare.all_junctions.{A,B}.tsv.gz`.
    #[arg(long = "emit-all-junctions")]
    emit_all_junctions: bool,
}

/// The per-side columns pass 2 reads: what `side_genomic_junctions` needs.
const PASS2_SIDE_COLS: [Col; 4] = [TARGET_NAME, TARGET_START, GENOMIC_INTRON_STRAND, CS];

/// The per-side columns pass 1 reads, resolved for each side up front.
const SIDE_COLS: [Col; 8] = [
    CS,
    STRAND,
    GENOMIC_INTRON_STRAND,
    QUERY_START,
    QUERY_END,
    TARGET_START,
    TARGET_NAME,
    N_INTRON_EVENTS,
];

pub fn run(args: &QueryJunctionCompareArgs) -> Result<()> {
    if !args.input.ends_with(".parquet") {
        bail!(
            "--input must be a Parquet file ('.parquet') — got '{}'. \
             query-junction-compare doesn't support TSV input.",
            args.input
        );
    }

    let outdir = Path::new(&args.outdir);
    fs::create_dir_all(outdir).with_context(|| format!("cannot create --outdir '{}'", args.outdir))?;

    let path = |name: String| outdir.join(name).to_string_lossy().into_owned();
    let summary_out = path(format!("{}.query_junction_compare.summary.tsv", args.prefix));
    let per_read_out = path(format!("{}.query_junction_compare.per_read_per_junc_info.tsv.gz", args.prefix));
    let unmatched_a_out = path(format!("{}.query_junction_compare.unmatched_junctions.A.tsv.gz", args.prefix));
    let unmatched_b_out = path(format!("{}.query_junction_compare.unmatched_junctions.B.tsv.gz", args.prefix));
    let all_a_out = path(format!("{}.query_junction_compare.all_junctions.A.tsv.gz", args.prefix));
    let all_b_out = path(format!("{}.query_junction_compare.all_junctions.B.tsv.gz", args.prefix));

    let wanted = wanted_columns(
        &[QUERY_NAME, QUERY_LEN, LABEL_A, LABEL_B, N_QUERY_JUNCTIONS_ONLY_A, N_QUERY_JUNCTIONS_ONLY_B],
        &SIDE_COLS,
    );
    let wanted_refs: Vec<&str> = wanted.iter().map(String::as_str).collect();
    let mut source = open_table(&args.input, InputFormat::Parquet, Some(&wanted_refs))?;

    let cols_owned = source.header().with_context(|| format!("reading comparison table '{}'", args.input))?;
    let cols: Vec<&str> = cols_owned.iter().map(String::as_str).collect();
    let col_index: HashMap<&str, usize> = cols.iter().copied().enumerate().map(|(i, c)| (c, i)).collect();
    require_ab_schema(&cols)?;

    let read_name_idx = require_col(&col_index, QUERY_NAME.name)?;
    let read_len_idx = require_col(&col_index, QUERY_LEN.name)?;
    let n_only_a_idx = require_col(&col_index, N_QUERY_JUNCTIONS_ONLY_A.name)?;
    let n_only_b_idx = require_col(&col_index, N_QUERY_JUNCTIONS_ONLY_B.name)?;
    let idx_a = SideIndex::resolve(&col_index, &SIDE_COLS, SIDE_A)?;
    let idx_b = SideIndex::resolve(&col_index, &SIDE_COLS, SIDE_B)?;

    let mut label_a = "A".to_string();
    let mut label_b = "B".to_string();
    let mut seen_row = false;

    let mut summary = Summary::default();
    let mut per_read_w = open_output(Some(&per_read_out))?;
    writeln!(per_read_w, "{}", header_line(&QJC_PER_JUNCTION_COLS))?;
    let mut catalog_a = JunctionCatalog::default();
    let mut catalog_b = JunctionCatalog::default();

    while let Some(fields_owned) = source.next_row()? {
        let fields: Vec<&str> = fields_owned.iter().map(String::as_str).collect();
        if !seen_row {
            (label_a, label_b) = labels_from_row(&col_index, &fields);
            seen_row = true;
        }
        let get_a = |c: &str| idx_a.get(&fields, c);
        let get_b = |c: &str| idx_b.get(&fields, c);
        let read_name = fields.get(read_name_idx).copied().unwrap_or("");
        let read_len: u64 = fields.get(read_len_idx).copied().unwrap_or("0").parse().unwrap_or(0);
        let n_only_a: u64 = fields.get(n_only_a_idx).copied().unwrap_or("0").parse().unwrap_or(0);
        let n_only_b: u64 = fields.get(n_only_b_idx).copied().unwrap_or("0").parse().unwrap_or(0);

        let mapped_a = is_mapped(get_a(TARGET_NAME.name));
        let mapped_b = is_mapped(get_b(TARGET_NAME.name));
        let junc_count_a: u64 = get_a(N_INTRON_EVENTS.name).parse().unwrap_or(0);
        let junc_count_b: u64 = get_b(N_INTRON_EVENTS.name).parse().unwrap_or(0);

        let differing =
            summary.observe_row(mapped_a, mapped_b, n_only_a, n_only_b, junc_count_a, junc_count_b);
        if !differing {
            continue;
        }

        let junctions_a = if mapped_a {
            build_side_junctions(
                get_a(CS.name),
                get_a(STRAND.name).chars().next().unwrap_or('.'),
                get_a(QUERY_START.name).parse().unwrap_or(0),
                get_a(QUERY_END.name).parse().unwrap_or(0),
                get_a(TARGET_START.name).parse().unwrap_or(0),
                get_a(TARGET_NAME.name),
                read_len,
            )
        } else {
            Vec::new()
        };
        let junctions_b = if mapped_b {
            build_side_junctions(
                get_b(CS.name),
                get_b(STRAND.name).chars().next().unwrap_or('.'),
                get_b(QUERY_START.name).parse().unwrap_or(0),
                get_b(QUERY_END.name).parse().unwrap_or(0),
                get_b(TARGET_START.name).parse().unwrap_or(0),
                get_b(TARGET_NAME.name),
                read_len,
            )
        } else {
            Vec::new()
        };

        let intron_strand_a = normalize_intron_strand(get_a(GENOMIC_INTRON_STRAND.name));
        let intron_strand_b = normalize_intron_strand(get_b(GENOMIC_INTRON_STRAND.name));
        let (qset_a, gset_a) = side_position_sets(&junctions_a, intron_strand_a);
        let (qset_b, gset_b) = side_position_sets(&junctions_b, intron_strand_b);

        let strand_a = get_a(STRAND.name).chars().next().unwrap_or('.');
        let strand_b = get_b(STRAND.name).chars().next().unwrap_or('.');
        let recs_a =
            build_junction_records(junctions_a, intron_strand_a, strand_a, &qset_b, &gset_b, mapped_b);
        let recs_b =
            build_junction_records(junctions_b, intron_strand_b, strand_b, &qset_a, &gset_a, mapped_a);

        for r in &recs_a {
            write_per_read_row(&mut per_read_w, read_name, 'A', r)?;
        }
        for r in &recs_b {
            write_per_read_row(&mut per_read_w, read_name, 'B', r)?;
        }

        catalog_a.add_unmatched(&recs_a);
        catalog_b.add_unmatched(&recs_b);
    }
    per_read_w.flush()?;

    // ── Pass 2: every junction of every read, per side ──────────────────────
    // Re-scan the whole table, narrowly column-projected, adding each side's
    // junctions (whole-table read counts and motifs) to its catalog.
    let pass2_wanted = wanted_columns(&[], &PASS2_SIDE_COLS);
    let pass2_wanted_refs: Vec<&str> = pass2_wanted.iter().map(String::as_str).collect();
    let mut pass2_source = open_table(&args.input, InputFormat::Parquet, Some(&pass2_wanted_refs))
        .with_context(|| format!("re-reading comparison table '{}' for pass 2", args.input))?;
    let pass2_cols_owned = pass2_source.header()?;
    let pass2_col_index: HashMap<&str, usize> =
        pass2_cols_owned.iter().map(String::as_str).enumerate().map(|(i, c)| (c, i)).collect();
    let p2_a = SideIndex::resolve(&pass2_col_index, &PASS2_SIDE_COLS, SIDE_A)?;
    let p2_b = SideIndex::resolve(&pass2_col_index, &PASS2_SIDE_COLS, SIDE_B)?;

    while let Some(fields_owned) = pass2_source.next_row()? {
        let fields: Vec<&str> = fields_owned.iter().map(String::as_str).collect();
        for (idx, catalog) in [(&p2_a, &mut catalog_a), (&p2_b, &mut catalog_b)] {
            let get = |c: &str| idx.get(&fields, c);
            catalog.add_side(&side_genomic_junctions(
                get(TARGET_NAME.name),
                get(TARGET_START.name),
                get(GENOMIC_INTRON_STRAND.name),
                get(CS.name),
            ));
        }
    }

    catalog_a.write(&unmatched_a_out, true)?;
    catalog_b.write(&unmatched_b_out, true)?;
    if args.emit_all_junctions {
        catalog_a.write(&all_a_out, false)?;
        catalog_b.write(&all_b_out, false)?;
    }
    let unique = unique_junction_counts(&catalog_a, &catalog_b);
    summary.set_junction_counts(&unique, catalog_a.n_unmatched(), catalog_b.n_unmatched());
    summary.write_tsv(&summary_out, &label_a, &label_b)?;

    eprintln!(
        "Query junction compare (A={label_a}, B={label_b}): {} of {} reads had differing query junctions",
        summary.n_query_junctions_different(),
        summary.n_total()
    );
    eprintln!(
        "  unique junctions: A={}, B={} ({} shared, {} only in A, {} only in B)",
        unique.n_a, unique.n_b, unique.n_shared, unique.n_only_a, unique.n_only_b
    );
    eprintln!("Outputs in {}:", args.outdir);
    eprintln!("  {summary_out}");
    eprintln!("  {per_read_out}");
    eprintln!("  {unmatched_a_out}");
    eprintln!("  {unmatched_b_out}");
    if args.emit_all_junctions {
        eprintln!("  {all_a_out}");
        eprintln!("  {all_b_out}");
    }
    Ok(())
}
