//! Streamlined splice-junction comparison subcommand (`compare-junctions`).
//!
//! A slim variant of `compare` focused on per-read splice-junction differences.
//! Emits 28 columns total:
//!
//!   - 2 join keys: `Read_Name`, `Read_Len`
//!   - 8 per-side data columns × 2 sides (suffixed with `--label-a` / `--label-b`):
//!     `TargetRef_1st, Num_Aln, JuncCount, seqid_Max, Query_Aln_Cov_Max,
//!      junctions, genomic_junctions, cs`
//!   - 10 comparison metrics (no suffix):
//!     `seqid_Diff, QueryAlnCov_Diff,
//!      N_Matched_Junctions, N_Unmatched_Junctions, N_Junctions_OnlyA, N_Junctions_OnlyB,
//!      Genomic_N_Matched_Junctions, Genomic_N_Unmatched_Junctions,
//!      Genomic_N_Junctions_OnlyA, Genomic_N_Junctions_OnlyB`
//!
//! Same streaming merge-join algorithm and inner-join semantics as `compare`. The
//! genomic-junction set metrics are always emitted here (no flag) because this whole
//! subcommand is about junction differences.

use std::io::Write;

use anyhow::{bail, Context, Result};

use crate::compare_streaming::{ReadInfoReader, ReadKey};
use crate::io_utils::{escape_tsv_field, fmt_float, open_output};
use crate::junction::{
    format_genomic_junction_tuple, format_junction_tuple, genomic_junction_set_diffs,
    genomic_junction_set_stats, junction_set_diffs, junction_set_stats,
    parse_genomic_junction_str, parse_junction_str,
};

// ── CLI args ─────────────────────────────────────────────────────────────────

#[derive(clap::Args, Debug)]
pub struct CompareJunctionsArgs {
    /// Readinfo TSV A (must be sorted by Read_Name, Read_Len)
    #[arg(short = 'a', long = "readinfo-a", value_name = "readinfo_a.tsv")]
    pub readinfo_a: String,

    /// Readinfo TSV B (must be sorted by Read_Name, Read_Len)
    #[arg(short = 'b', long = "readinfo-b", value_name = "readinfo_b.tsv")]
    pub readinfo_b: String,

    /// Label for dataset A (used as per-side column suffix)
    #[arg(long = "label-a", value_name = "LABEL", default_value = "SetA")]
    pub label_a: String,

    /// Label for dataset B (used as per-side column suffix)
    #[arg(long = "label-b", value_name = "LABEL", default_value = "SetB")]
    pub label_b: String,

    /// Output comparison TSV file ('.gz' for gzip)
    #[arg(short = 'o', long = "output", value_name = "compare_junctions.tsv[.gz]")]
    pub output: String,

    /// Skip reads that appear in only one file instead of stopping with an error.
    /// Both readinfo files must still be lex-sorted by Read_Name for the skip
    /// heuristic to work correctly. Unmatched reads are counted in the summary.
    #[arg(long = "ignore-row-mismatch")]
    pub ignore_row_mismatch: bool,
}

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

// ── Main streaming comparison function ──────────────────────────────────────

pub fn run(args: &CompareJunctionsArgs) -> Result<()> {
    eprintln!("[INFO] Opening readinfo files...");
    eprintln!("  A: {}", args.readinfo_a);
    eprintln!("  B: {}", args.readinfo_b);

    let mut reader_a = ReadInfoReader::new(&args.readinfo_a)?;
    let mut reader_b = ReadInfoReader::new(&args.readinfo_b)?;

    let idx_name_a = reader_a
        .get_col_idx("Read_Name")
        .context("missing Read_Name in A")?;
    let idx_len_a = reader_a
        .get_col_idx("Read_Len")
        .context("missing Read_Len in A")?;
    let idx_name_b = reader_b
        .get_col_idx("Read_Name")
        .context("missing Read_Name in B")?;
    let idx_len_b = reader_b
        .get_col_idx("Read_Len")
        .context("missing Read_Len in B")?;

    let mut out = open_output(Some(&args.output))?;

    // Write header.
    write_compare_junctions_header(&mut out, &args.label_a, &args.label_b)?;

    eprintln!("[INFO] Starting merge-join comparison...");

    let mut n_a_total: u64 = 0;
    let mut n_b_total: u64 = 0;
    let mut n_merged: u64 = 0;

    while let (Some(a_fields), Some(b_fields)) = (reader_a.current(), reader_b.current()) {
        let key_a = ReadKey::from_fields(a_fields, idx_name_a, idx_len_a)?;
        let key_b = ReadKey::from_fields(b_fields, idx_name_b, idx_len_b)?;

        if key_a == key_b {
            n_a_total += 1;
            n_b_total += 1;

            emit_compare_junctions_row(
                &mut out,
                &key_a.name,
                key_a.len,
                |c| reader_a.get_col(c).unwrap_or(""),
                |c| reader_b.get_col(c).unwrap_or(""),
            )?;

            n_merged += 1;
            if n_merged % 100_000 == 0 {
                eprintln!("[INFO] Processed {} matched records...", n_merged);
            }

            reader_a.advance()?;
            reader_b.advance()?;
        } else if key_a < key_b {
            if !args.ignore_row_mismatch {
                bail!(
                    "read-name mismatch: A has {:?} but B has {:?} \
                     (A row #{}, B row #{}). Both readinfo files must list reads \
                     in the same order. Use --ignore-row-mismatch to skip \
                     unmatched reads instead of stopping.",
                    key_a.name, key_b.name, n_a_total + 1, n_b_total + 1
                );
            }
            n_a_total += 1;
            reader_a.advance()?;
        } else {
            if !args.ignore_row_mismatch {
                bail!(
                    "read-name mismatch: B has {:?} but A has {:?} \
                     (A row #{}, B row #{}). Both readinfo files must list reads \
                     in the same order. Use --ignore-row-mismatch to skip \
                     unmatched reads instead of stopping.",
                    key_b.name, key_a.name, n_a_total + 1, n_b_total + 1
                );
            }
            n_b_total += 1;
            reader_b.advance()?;
        }
    }

    // Drain remaining unmatched on either side.
    while reader_a.current().is_some() {
        n_a_total += 1;
        reader_a.advance()?;
    }
    while reader_b.current().is_some() {
        n_b_total += 1;
        reader_b.advance()?;
    }

    out.flush()?;

    // End-of-run summary on stderr.
    let a_only = n_a_total - n_merged;
    let b_only = n_b_total - n_merged;
    eprintln!("Junction comparison summary:");
    eprintln!("  Label A: {}", args.label_a);
    eprintln!("  Label B: {}", args.label_b);
    eprintln!("  rows in A (readinfo-a):     {n_a_total}");
    eprintln!("  rows in B (readinfo-b):     {n_b_total}");
    eprintln!("  matched (in both, written): {n_merged}");
    eprintln!("  A-only (dropped, not in B): {a_only}");
    eprintln!("  B-only (dropped, not in A): {b_only}");

    Ok(())
}
