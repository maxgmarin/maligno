/// Streaming merge-join comparison of two readinfo TSVs.
///
/// Memory usage: O(1) [constant, independent of file size]
/// Time: O(|A| + |B|) [after sorting both files by (Read_Name, Read_Len)]
///
/// Algorithm: Two-pointer merge-join
/// - Stream through both files simultaneously
/// - When keys match, output comparison row
/// - When keys don't match, advance the pointer with the smaller key
///
/// REQUIREMENT: Both input files must be sorted by (Read_Name, Read_Len)

use std::io::{BufRead, Write};

use anyhow::{bail, Context, Result};

use crate::compare_junctions::{emit_compare_junctions_row, write_compare_junctions_header};
use crate::io_utils::{escape_tsv_field, fmt_float, open_input, open_output};
use crate::junction::{
    format_genomic_junction_tuple, format_junction_tuple, genomic_junction_set_diffs,
    genomic_junction_set_stats, junction_distance, junction_set_diffs, junction_set_stats,
    parse_genomic_junction_str, parse_junction_str,
};

// ── Comparison mode ───────────────────────────────────────────────────────────

/// Which comparison view to emit. Shared by `compare` and `pafcompare` so the
/// two commands expose an identical `--mode` surface.
#[derive(Clone, Debug, Default, clap::ValueEnum)]
pub(crate) enum CompareMode {
    /// All per-read metrics, including genomic-junction comparison (94 cols).
    #[default]
    Full,
    /// Splice-junction-focused view (47 cols).
    Junctions,
}

// ── CLI args ─────────────────────────────────────────────────────────────────

#[derive(clap::Args, Debug)]
pub struct CompareStreamingArgs {
    /// Readinfo TSV A (must be sorted by Read_Name, Read_Len)
    #[arg(short = 'a', long = "readinfo-a", value_name = "readinfo_a.tsv")]
    pub readinfo_a: String,

    /// Readinfo TSV B (must be sorted by Read_Name, Read_Len)
    #[arg(short = 'b', long = "readinfo-b", value_name = "readinfo_b.tsv")]
    pub readinfo_b: String,

    /// Label for dataset A (used as column suffix)
    #[arg(long = "label-a", value_name = "LABEL", default_value = "SetA")]
    pub label_a: String,

    /// Label for dataset B (used as column suffix)
    #[arg(long = "label-b", value_name = "LABEL", default_value = "SetB")]
    pub label_b: String,

    /// Output comparison TSV file ('.gz' for gzip)
    #[arg(short = 'o', long = "output", value_name = "compare.tsv[.gz]")]
    pub output: String,

    /// Skip reads that appear in only one file instead of stopping with an error.
    /// Both readinfo files must still be lex-sorted by Read_Name for the skip
    /// heuristic to work correctly. Unmatched reads are counted in the summary.
    #[arg(long = "ignore-row-mismatch")]
    pub ignore_row_mismatch: bool,

    /// Comparison view: `full` (all per-read metrics incl. genomic-junction
    /// comparison, 94 cols) or `junctions` (splice-junction-focused, 47 cols).
    #[arg(long = "mode", value_enum, default_value_t = CompareMode::Full)]
    pub mode: CompareMode,
}

// ── ReadInfo column indices (must match readinfo.rs) ──────────────────────────

const READINFO_DATA_COLS: &[&str] = &[
    "TargetChr",
    "Strand",
    "MQ_Best",
    "AS_Max",
    "ms_Max",
    "Query_Aln_Cov_Max",
    "Query_Aln_Len_Max",
    "seqid_Max",
    "junctions",
    "Num_Aln",
    "Num_Aln_MaxScore",
    "JuncCount",
    "N_Match_Events",
    "N_Match_Bases",
    "N_Substitution_Events",
    "N_Substitution_Bases",
    "N_Insertion_Events",
    "N_Insertion_Bases",
    "N_Deletion_Events",
    "N_Deletion_Bases",
    "N_Splice_Junction_Events",
    "N_Splice_Junction_Bases",
    "N_SoftClipped_Bases_Start",
    "N_SoftClipped_Bases_End",
    "N_SoftClipped_Events",
    "cs",
    "genomic_junctions",
    "Query_Start",
    "Query_End",
    "Target_Start",
    "Target_End",
];

fn comparison_col_names() -> Vec<&'static str> {
    vec![
        "Strand_Match",
        "AS_Diff",
        "ms_Diff",
        "AS_Ratio",
        "ms_Ratio",
        "seqid_Diff",
        "QueryAlnLen_Diff",
        "QueryAlnCov_Diff",
        "N_Insertion_Bases_Diff",
        "N_Insertion_Bases_Ratio",
        "N_Deletion_Bases_Diff",
        "N_Deletion_Bases_Ratio",
        "N_Substitution_Bases_Diff",
        "N_Substitution_Bases_Ratio",
        "N_SoftClipped_Bases_Start_Diff",
        "N_SoftClipped_Bases_End_Diff",
        "Junction_Distance",
        "N_Unmatched_Junctions",
        "Junc_Dist_V2",
        "N_Matched_Junctions",
        "N_Junctions_OnlyA",
        "N_Junctions_OnlyB",
        "Genomic_N_Matched_Junctions",
        "Genomic_N_Unmatched_Junctions",
        "Genomic_N_Junctions_OnlyA",
        "Genomic_N_Junctions_OnlyB",
        // Object columns at the very end (the actual non-overlapping junctions,
        // as opposed to just their counts above).
        "Junctions_OnlyA",
        "Junctions_OnlyB",
        "Genomic_Junctions_OnlyA",
        "Genomic_Junctions_OnlyB",
    ]
}

// ── Reusable header + row emitters (shared with `pafcompare`) ────────────────

/// Write the `compare` output header: `Read_Name`, `Read_Len`, the per-side data
/// columns (suffixed with each label), then the comparison/object columns.
pub(crate) fn write_compare_header<W: Write>(
    out: &mut W,
    label_a: &str,
    label_b: &str,
) -> std::io::Result<()> {
    write!(out, "Read_Name\tRead_Len")?;
    for col in READINFO_DATA_COLS {
        write!(out, "\t{col}_{label_a}")?;
    }
    for col in READINFO_DATA_COLS {
        write!(out, "\t{col}_{label_b}")?;
    }
    for col in comparison_col_names() {
        write!(out, "\t{col}")?;
    }
    writeln!(out)
}

/// Emit one comparison row given by-name column accessors for each side.
///
/// `get_a` / `get_b` return the readinfo column value (or `""` if absent) — the
/// existing `.parse().unwrap_or(default)` calls below preserve identical
/// defaults to the old `reader.get_col(c).unwrap_or(default)` form. This is the
/// single source of truth for the `compare` row, shared by `compare::run`
/// (reading from `ReadInfoReader`) and `pafcompare` (reading from in-memory
/// `ReadInfoRow`s serialized to readinfo lines).
pub(crate) fn emit_compare_row<'r, W, FA, FB>(
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
    // Extract per-side passthrough fields.
    let a_raw_fields: Vec<&str> = READINFO_DATA_COLS.iter().map(|c| get_a(c)).collect();
    let b_raw_fields: Vec<&str> = READINFO_DATA_COLS.iter().map(|c| get_b(c)).collect();

    // Parse typed values from A.
    let a_as_max: i64 = get_a("AS_Max").parse().unwrap_or(0);
    let a_ms_max: i64 = get_a("ms_Max").parse().unwrap_or(0);
    let a_seqid_max: f64 = get_a("seqid_Max").parse().unwrap_or(f64::NAN);
    let a_aln_len_max: u64 = get_a("Query_Aln_Len_Max").parse().unwrap_or(0);
    let a_cov_max: f64 = get_a("Query_Aln_Cov_Max").parse().unwrap_or(f64::NAN);
    let a_n_ins_bases: u64 = get_a("N_Insertion_Bases").parse().unwrap_or(0);
    let a_n_del_bases: u64 = get_a("N_Deletion_Bases").parse().unwrap_or(0);
    let a_n_sub_bases: u64 = get_a("N_Substitution_Bases").parse().unwrap_or(0);
    let a_n_sc_start: u64 = get_a("N_SoftClipped_Bases_Start").parse().unwrap_or(0);
    let a_n_sc_end: u64 = get_a("N_SoftClipped_Bases_End").parse().unwrap_or(0);
    let a_junc_count: usize = get_a("JuncCount").parse().unwrap_or(0);
    let a_junctions = get_a("junctions");
    let a_strand = get_a("Strand");

    // Parse typed values from B.
    let b_as_max: i64 = get_b("AS_Max").parse().unwrap_or(0);
    let b_ms_max: i64 = get_b("ms_Max").parse().unwrap_or(0);
    let b_seqid_max: f64 = get_b("seqid_Max").parse().unwrap_or(f64::NAN);
    let b_aln_len_max: u64 = get_b("Query_Aln_Len_Max").parse().unwrap_or(0);
    let b_cov_max: f64 = get_b("Query_Aln_Cov_Max").parse().unwrap_or(f64::NAN);
    let b_n_ins_bases: u64 = get_b("N_Insertion_Bases").parse().unwrap_or(0);
    let b_n_del_bases: u64 = get_b("N_Deletion_Bases").parse().unwrap_or(0);
    let b_n_sub_bases: u64 = get_b("N_Substitution_Bases").parse().unwrap_or(0);
    let b_n_sc_start: u64 = get_b("N_SoftClipped_Bases_Start").parse().unwrap_or(0);
    let b_n_sc_end: u64 = get_b("N_SoftClipped_Bases_End").parse().unwrap_or(0);
    let b_junc_count: usize = get_b("JuncCount").parse().unwrap_or(0);
    let b_junctions = get_b("junctions");
    let b_strand = get_b("Strand");

    // Compute metrics.
    let strand_match = a_strand == b_strand;
    let as_diff = b_as_max - a_as_max;
    let ms_diff = b_ms_max - a_ms_max;
    let as_ratio = safe_ratio_i64(b_as_max, a_as_max);
    let ms_ratio = safe_ratio_i64(b_ms_max, a_ms_max);
    let seqid_diff = b_seqid_max - a_seqid_max;
    let qal_diff = b_aln_len_max as i64 - a_aln_len_max as i64;
    let qac_diff = b_cov_max - a_cov_max;
    let n_ins_diff = b_n_ins_bases as i64 - a_n_ins_bases as i64;
    let n_ins_ratio = safe_ratio_u64(b_n_ins_bases, a_n_ins_bases);
    let n_del_diff = b_n_del_bases as i64 - a_n_del_bases as i64;
    let n_del_ratio = safe_ratio_u64(b_n_del_bases, a_n_del_bases);
    let n_sub_diff = b_n_sub_bases as i64 - a_n_sub_bases as i64;
    let n_sub_ratio = safe_ratio_u64(b_n_sub_bases, a_n_sub_bases);
    let n_sc_start_diff = b_n_sc_start as i64 - a_n_sc_start as i64;
    let n_sc_end_diff = b_n_sc_end as i64 - a_n_sc_end as i64;

    // Junction metrics.
    let juncs_a = parse_junction_str(a_junctions);
    let juncs_b = parse_junction_str(b_junctions);
    let junction_distance_val = junction_distance(&juncs_a, &juncs_b);

    let n_junc_count_diff = (a_junc_count as i64 - b_junc_count as i64).unsigned_abs();
    let junc_dist_v2 = 50 * n_junc_count_diff;

    let (n_matched, n_only_a, n_only_b) = junction_set_stats(&juncs_a, &juncs_b);
    let n_unmatched = n_only_a + n_only_b;

    let (j_only_a_vec, j_only_b_vec) = junction_set_diffs(&juncs_a, &juncs_b);
    let j_only_a_str = format_junction_tuple(&j_only_a_vec);
    let j_only_b_str = format_junction_tuple(&j_only_b_vec);

    // Genomic-junction set comparison (always emitted in `full` mode).
    let (g_matched, g_only_a, g_only_b, g_unmatched, g_only_a_str, g_only_b_str) = {
        let a_genomic = get_a("genomic_junctions");
        let b_genomic = get_b("genomic_junctions");
        let chrom_a = get_a("TargetChr").to_string();
        let chrom_b = get_b("TargetChr").to_string();
        let pairs_a = parse_genomic_junction_str(a_genomic);
        let pairs_b = parse_genomic_junction_str(b_genomic);
        let gj_a: Vec<(String, u64, u64)> = pairs_a
            .into_iter()
            .map(|(s, e)| (chrom_a.clone(), s, e))
            .collect();
        let gj_b: Vec<(String, u64, u64)> = pairs_b
            .into_iter()
            .map(|(s, e)| (chrom_b.clone(), s, e))
            .collect();
        let (m, oa, ob) = genomic_junction_set_stats(&gj_a, &gj_b);
        let (gj_only_a_vec, gj_only_b_vec) = genomic_junction_set_diffs(&gj_a, &gj_b);
        let oa_str = format_genomic_junction_tuple(&gj_only_a_vec);
        let ob_str = format_genomic_junction_tuple(&gj_only_b_vec);
        (m, oa, ob, oa + ob, oa_str, ob_str)
    };

    // Write output row.
    write!(out, "{name}\t{len}")?;
    for f in &a_raw_fields {
        write!(out, "\t{}", escape_tsv_field(f))?;
    }
    for f in &b_raw_fields {
        write!(out, "\t{}", escape_tsv_field(f))?;
    }
    write!(
        out,
        "\t{strand_match}\t\
         {as_diff}\t{ms_diff}\t{as_ratio}\t{ms_ratio}\t{}\t{qal_diff}\t{}\t\
         {n_ins_diff}\t{n_ins_ratio}\t\
         {n_del_diff}\t{n_del_ratio}\t\
         {n_sub_diff}\t{n_sub_ratio}\t\
         {n_sc_start_diff}\t{n_sc_end_diff}\t\
         {junction_distance_val}\t{n_unmatched}\t{junc_dist_v2}\t{n_matched}\
         \t{n_only_a}\t{n_only_b}",
        fmt_float(seqid_diff),
        fmt_float(qac_diff),
    )?;
    write!(out, "\t{g_matched}\t{g_unmatched}\t{g_only_a}\t{g_only_b}")?;
    write!(
        out,
        "\t{}\t{}",
        escape_tsv_field(&j_only_a_str),
        escape_tsv_field(&j_only_b_str),
    )?;
    write!(
        out,
        "\t{}\t{}",
        escape_tsv_field(&g_only_a_str),
        escape_tsv_field(&g_only_b_str),
    )?;
    writeln!(out)
}

// ── ReadKey for sorting/comparison ──────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct ReadKey {
    pub(crate) name: String,
    pub(crate) len: u64,
}

impl ReadKey {
    pub(crate) fn from_fields(fields: &[String], idx_name: usize, idx_len: usize) -> Result<Self> {
        let name = fields[idx_name].clone();
        let len: u64 = fields[idx_len]
            .parse()
            .context("failed to parse Read_Len")?;
        Ok(ReadKey { name, len })
    }
}

// ── Helper functions for ratio calculations ─────────────────────────────────

fn safe_ratio_i64(a: i64, b: i64) -> String {
    if b == 0 {
        "NaN".to_string()
    } else {
        fmt_float(a as f64 / b as f64)
    }
}

fn safe_ratio_u64(a: u64, b: u64) -> String {
    if b == 0 {
        "NaN".to_string()
    } else {
        fmt_float(a as f64 / b as f64)
    }
}

// ── Streaming line iterator with buffering ─────────────────────────────────

pub(crate) struct ReadInfoReader {
    reader: Box<dyn BufRead>,
    header: Vec<String>,
    col_map: std::collections::HashMap<String, usize>,
    current_line: Option<Vec<String>>,
    done: bool,
}

impl ReadInfoReader {
    pub(crate) fn new(path: &str) -> Result<Self> {
        let reader = open_input(path)?;

        let mut r = ReadInfoReader {
            reader,
            header: Vec::new(),
            col_map: std::collections::HashMap::new(),
            current_line: None,
            done: false,
        };

        r.read_header()?;
        r.advance()?;
        Ok(r)
    }

    fn read_header(&mut self) -> Result<()> {
        let mut line = String::new();
        self.reader.read_line(&mut line)?;
        self.header = line
            .trim_end()
            .split('\t')
            .map(|s| s.to_string())
            .collect();

        for (i, col) in self.header.iter().enumerate() {
            self.col_map.insert(col.clone(), i);
        }

        Ok(())
    }

    pub(crate) fn advance(&mut self) -> Result<()> {
        if self.done {
            self.current_line = None;
            return Ok(());
        }

        let mut line = String::new();
        match self.reader.read_line(&mut line) {
            Ok(0) => {
                self.done = true;
                self.current_line = None;
            }
            Ok(_) => {
                if line.trim().is_empty() {
                    return self.advance(); // Skip empty lines
                }
                let fields: Vec<String> = line
                    .trim_end()
                    .split('\t')
                    .map(|s| s.to_string())
                    .collect();
                self.current_line = Some(fields);
            }
            Err(e) => return Err(e.into()),
        }

        Ok(())
    }

    pub(crate) fn current(&self) -> Option<&[String]> {
        self.current_line.as_ref().map(|v| v.as_slice())
    }

    pub(crate) fn get_col(&self, col_name: &str) -> Option<&str> {
        self.col_map
            .get(col_name)
            .and_then(|&i| {
                self.current_line
                    .as_ref()
                    .and_then(|fields| fields.get(i).map(|s| s.as_str()))
            })
    }

    pub(crate) fn get_col_idx(&self, col_name: &str) -> Option<usize> {
        self.col_map.get(col_name).copied()
    }
}

// ── Main streaming comparison function ──────────────────────────────────────

pub fn run(args: &CompareStreamingArgs) -> Result<()> {
    eprintln!("[INFO] Opening readinfo files...");
    eprintln!("  A: {}", args.readinfo_a);
    eprintln!("  B: {}", args.readinfo_b);

    let mut reader_a = ReadInfoReader::new(&args.readinfo_a)?;
    let mut reader_b = ReadInfoReader::new(&args.readinfo_b)?;

    // Get column indices
    let idx_name_a = reader_a.get_col_idx("Read_Name").context("missing Read_Name in A")?;
    let idx_len_a = reader_a.get_col_idx("Read_Len").context("missing Read_Len in A")?;
    let idx_name_b = reader_b.get_col_idx("Read_Name").context("missing Read_Name in B")?;
    let idx_len_b = reader_b.get_col_idx("Read_Len").context("missing Read_Len in B")?;

    // Open output
    let mut out = open_output(Some(&args.output))?;

    // `junctions` mode emits the splice-focused 47-col view (genomic-junction
    // metrics always on); `full` mode emits the full per-read comparison.
    let junctions_mode = matches!(args.mode, CompareMode::Junctions);

    // Write header
    if junctions_mode {
        write_compare_junctions_header(&mut out, &args.label_a, &args.label_b)?;
    } else {
        write_compare_header(&mut out, &args.label_a, &args.label_b)?;
    }

    eprintln!("[INFO] Starting comparison...");

    let mut n_a_total: u64 = 0;
    let mut n_b_total: u64 = 0;
    let mut n_merged: u64 = 0;

    while let (Some(a_fields), Some(b_fields)) = (reader_a.current(), reader_b.current()) {
        let key_a = ReadKey::from_fields(a_fields, idx_name_a, idx_len_a)?;
        let key_b = ReadKey::from_fields(b_fields, idx_name_b, idx_len_b)?;

        if key_a == key_b {
            n_a_total += 1;
            n_b_total += 1;

            if junctions_mode {
                emit_compare_junctions_row(
                    &mut out,
                    &key_a.name,
                    key_a.len,
                    |c| reader_a.get_col(c).unwrap_or(""),
                    |c| reader_b.get_col(c).unwrap_or(""),
                )?;
            } else {
                emit_compare_row(
                    &mut out,
                    &key_a.name,
                    key_a.len,
                    |c| reader_a.get_col(c).unwrap_or(""),
                    |c| reader_b.get_col(c).unwrap_or(""),
                )?;
            }

            n_merged += 1;
            if n_merged % 100000 == 0 {
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

    // Count remaining records
    while reader_a.current().is_some() {
        n_a_total += 1;
        reader_a.advance()?;
    }
    while reader_b.current().is_some() {
        n_b_total += 1;
        reader_b.advance()?;
    }

    out.flush()?;

    // ── End-of-run summary ────────────────────────────────────────────────
    let a_only = n_a_total - n_merged; // in A, absent from B
    let b_only = n_b_total - n_merged; // in B, absent from A
    eprintln!("Read comparison summary:");
    eprintln!("  Label A: {}", args.label_a);
    eprintln!("  Label B: {}", args.label_b);
    eprintln!("  rows in A (readinfo-a):     {n_a_total}");
    eprintln!("  rows in B (readinfo-b):     {n_b_total}");
    eprintln!("  matched (in both, written): {n_merged}");
    eprintln!("  A-only (dropped, not in B): {a_only}");
    eprintln!("  B-only (dropped, not in A): {b_only}");

    Ok(())
}
