//! Splice-junction-focused comparison view (the `--mode junctions` output).
//!
//! This module is a **library**: it provides the header + row emitters for the
//! 47-column junction-focused comparison, reused by both `compare --mode
//! junctions` (`compare_streaming::run`) and `pafcompare --mode junctions`
//! (`pafcompare::run`). It no longer owns a subcommand of its own.
//!
//! The view is a slim subset of the full `compare` output: a 15-column per-side
//! block (target/strand/MQ + junction/coverage fields, dropping the AS/ms scores
//! and cs-derived event counts) plus junction set-overlap metrics. The
//! genomic-junction set metrics are **always** emitted here (no flag), because
//! this view is entirely about junction differences.

use std::io::Write;

use crate::io_utils::{escape_tsv_field, fmt_float};
use crate::junction::{
    format_genomic_junction_tuple, format_junction_tuple, genomic_junction_set_diffs,
    genomic_junction_set_stats, junction_set_diffs, junction_set_stats,
    parse_genomic_junction_str, parse_junction_str,
};

// ── Column schemas ──────────────────────────────────────────────────────────

/// Per-side data columns (each appears once for A and once for B, suffixed).
const PER_SIDE_COLS: &[&str] = &[
    "TargetChr",
    "Strand",
    "MQ_Best",
    "Num_Aln",
    "Num_Aln_MaxScore",
    "JuncCount",
    "seqid_Max",
    "Query_Aln_Cov_Max",
    "junctions",
    "genomic_junctions",
    "cs",
    "Query_Start",
    "Query_End",
    "Target_Start",
    "Target_End",
];

/// Comparison metric column names (written as-is, no suffix).
const COMPARISON_COLS: &[&str] = &[
    "Strand_Match",
    "seqid_Diff",
    "QueryAlnCov_Diff",
    "N_Matched_Junctions",
    "N_Unmatched_Junctions",
    "N_Junctions_OnlyA",
    "N_Junctions_OnlyB",
    "Genomic_N_Matched_Junctions",
    "Genomic_N_Unmatched_Junctions",
    "Genomic_N_Junctions_OnlyA",
    "Genomic_N_Junctions_OnlyB",
    // Object lists at the very end — actual non-overlapping junctions, not counts.
    "Junctions_OnlyA",
    "Junctions_OnlyB",
    "Genomic_Junctions_OnlyA",
    "Genomic_Junctions_OnlyB",
];

// ── Reusable header + row emitters (shared with `pafcompare`) ────────────────

/// Write the `compare-junctions` header: keys, per-side data columns (suffixed
/// with each label), then the fixed comparison/object columns.
pub(crate) fn write_compare_junctions_header<W: Write>(
    out: &mut W,
    label_a: &str,
    label_b: &str,
) -> std::io::Result<()> {
    write!(out, "Read_Name\tRead_Len")?;
    for col in PER_SIDE_COLS {
        write!(out, "\t{col}_{label_a}")?;
    }
    for col in PER_SIDE_COLS {
        write!(out, "\t{col}_{label_b}")?;
    }
    for col in COMPARISON_COLS {
        write!(out, "\t{col}")?;
    }
    writeln!(out)
}

/// Emit one junction-comparison row given by-name column accessors for each
/// side. Single source of truth shared by `compare-junctions::run` (reading
/// from `ReadInfoReader`) and `pafcompare --junctions` (reading from in-memory
/// `ReadInfoRow`s serialized to readinfo lines).
pub(crate) fn emit_compare_junctions_row<'r, W, FA, FB>(
    out: &mut W,
    name: &str,
    len: u64,
    get_a: FA,
    get_b: FB,
) -> std::io::Result<()>
where
    W: Write,
    FA: Fn(&str) -> &'r str,
    FB: Fn(&str) -> &'r str,
{
    // Pull per-side raw values.
    let a_raw: Vec<&str> = PER_SIDE_COLS.iter().map(|c| get_a(c)).collect();
    let b_raw: Vec<&str> = PER_SIDE_COLS.iter().map(|c| get_b(c)).collect();

    // seqid and query-coverage diffs (B − A); NaN-safe via f64 parse.
    let a_seqid: f64 = get_a("seqid_Max").parse().unwrap_or(f64::NAN);
    let b_seqid: f64 = get_b("seqid_Max").parse().unwrap_or(f64::NAN);
    let a_cov: f64 = get_a("Query_Aln_Cov_Max").parse().unwrap_or(f64::NAN);
    let b_cov: f64 = get_b("Query_Aln_Cov_Max").parse().unwrap_or(f64::NAN);

    let seqid_diff = b_seqid - a_seqid;
    let qac_diff = b_cov - a_cov;

    // Strand match (boolean rendered as "true"/"false").
    let a_strand = get_a("Strand");
    let b_strand = get_b("Strand");
    let strand_match = a_strand == b_strand;

    // Query-junction set stats (same definition as `compare`).
    let juncs_a = parse_junction_str(get_a("junctions"));
    let juncs_b = parse_junction_str(get_b("junctions"));
    let (n_matched, n_only_a, n_only_b) = junction_set_stats(&juncs_a, &juncs_b);
    let n_unmatched = n_only_a + n_only_b;

    // Genomic-junction set stats. Cross-chrom safety is preserved by
    // reconstructing (chrom, start, end) tuples from the parsed (start, end)
    // pairs + the per-side TargetChr column.
    let chrom_a = get_a("TargetChr").to_string();
    let chrom_b = get_b("TargetChr").to_string();
    let pairs_a = parse_genomic_junction_str(get_a("genomic_junctions"));
    let pairs_b = parse_genomic_junction_str(get_b("genomic_junctions"));
    let gj_a: Vec<(String, u64, u64)> = pairs_a
        .into_iter()
        .map(|(s, e)| (chrom_a.clone(), s, e))
        .collect();
    let gj_b: Vec<(String, u64, u64)> = pairs_b
        .into_iter()
        .map(|(s, e)| (chrom_b.clone(), s, e))
        .collect();
    let (g_matched, g_only_a, g_only_b) = genomic_junction_set_stats(&gj_a, &gj_b);
    let g_unmatched = g_only_a + g_only_b;

    // Object diff lists (actual A-only / B-only junctions for both coord systems).
    let (j_only_a_vec, j_only_b_vec) = junction_set_diffs(&juncs_a, &juncs_b);
    let (gj_only_a_vec, gj_only_b_vec) = genomic_junction_set_diffs(&gj_a, &gj_b);
    let j_only_a_str = format_junction_tuple(&j_only_a_vec);
    let j_only_b_str = format_junction_tuple(&j_only_b_vec);
    let g_only_a_str = format_genomic_junction_tuple(&gj_only_a_vec);
    let g_only_b_str = format_genomic_junction_tuple(&gj_only_b_vec);

    // Write the row.
    write!(out, "{name}\t{len}")?;
    for f in &a_raw {
        write!(out, "\t{}", escape_tsv_field(f))?;
    }
    for f in &b_raw {
        write!(out, "\t{}", escape_tsv_field(f))?;
    }
    write!(
        out,
        "\t{strand_match}\t{}\t{}\t{n_matched}\t{n_unmatched}\t{n_only_a}\t{n_only_b}\
         \t{g_matched}\t{g_unmatched}\t{g_only_a}\t{g_only_b}",
        fmt_float(seqid_diff),
        fmt_float(qac_diff),
    )?;
    write!(
        out,
        "\t{}\t{}\t{}\t{}",
        escape_tsv_field(&j_only_a_str),
        escape_tsv_field(&j_only_b_str),
        escape_tsv_field(&g_only_a_str),
        escape_tsv_field(&g_only_b_str),
    )?;
    writeln!(out)
}
