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

use anyhow::{Context, Result};

use crate::io_utils::{escape_tsv_field, fmt_float, open_input, open_output};
use crate::junction::{
    genomic_junction_set_stats, junction_distance, junction_set_stats, parse_genomic_junction_str,
    parse_junction_str,
};

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
    #[arg(long = "label-a", value_name = "LABEL", default_value = "RefA")]
    pub label_a: String,

    /// Label for dataset B (used as column suffix)
    #[arg(long = "label-b", value_name = "LABEL", default_value = "RefB")]
    pub label_b: String,

    /// Output comparison TSV file ('.gz' for gzip)
    #[arg(short = 'o', long = "output", value_name = "compare.tsv[.gz]")]
    pub output: String,

    /// Also emit set-based comparison metrics for genome-coordinate junctions.
    /// Adds 4 columns at the end: Genomic_N_Matched_Junctions, Genomic_N_Unmatched_Junctions,
    /// Genomic_N_Junctions_OnlyA, Genomic_N_Junctions_OnlyB.
    /// Cross-chromosome compares are automatically disjoint because chrom is embedded
    /// in each genomic-junction tuple.
    #[arg(long = "compare-genomic-junctions")]
    pub compare_genomic_junctions: bool,
}

// ── ReadInfo column indices (must match readinfo.rs) ──────────────────────────

const READINFO_DATA_COLS: &[&str] = &[
    "TargetRef_1st",
    "AS_Max",
    "AS_Min",
    "ms_Max",
    "ms_Min",
    "Query_Aln_Cov_Max",
    "Query_Aln_Len_Max",
    "seqid_Max",
    "junctions",
    "Num_Aln",
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
    "num_bp_inserted",
    "cs",
    "genomic_junctions",
];

fn comparison_col_names(with_genomic: bool) -> Vec<&'static str> {
    let mut cols = vec![
        "AS_Diff",
        "ms_Diff",
        "AS_Ratio",
        "ms_Ratio",
        "seqid_Diff",
        "QueryAlnLen_Diff",
        "QueryAlnCov_Diff",
        "num_bp_inserted_Diff",
        "num_bp_inserted_Ratio",
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
    ];
    if with_genomic {
        cols.extend_from_slice(&[
            "Genomic_N_Matched_Junctions",
            "Genomic_N_Unmatched_Junctions",
            "Genomic_N_Junctions_OnlyA",
            "Genomic_N_Junctions_OnlyB",
        ]);
    }
    cols
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

    // Write header
    write!(out, "Read_Name\tRead_Len")?;
    for col in READINFO_DATA_COLS {
        write!(out, "\t{col}_{}", args.label_a)?;
    }
    for col in READINFO_DATA_COLS {
        write!(out, "\t{col}_{}", args.label_b)?;
    }
    for col in comparison_col_names(args.compare_genomic_junctions) {
        write!(out, "\t{col}")?;
    }
    writeln!(out)?;

    // Two-pointer merge-join
    eprintln!("[INFO] Starting merge-join comparison...");

    let mut n_a_total: u64 = 0;
    let mut n_b_total: u64 = 0;
    let mut n_merged: u64 = 0;

    while let (Some(a_fields), Some(b_fields)) = (reader_a.current(), reader_b.current()) {
        let key_a = ReadKey::from_fields(a_fields, idx_name_a, idx_len_a)?;
        let key_b = ReadKey::from_fields(b_fields, idx_name_b, idx_len_b)?;

        if key_a == key_b {
            // MATCH: compute and output comparison
            n_a_total += 1;
            n_b_total += 1;

            // Extract A fields
            let a_raw_fields: Vec<String> = READINFO_DATA_COLS
                .iter()
                .map(|col_name| {
                    reader_a
                        .get_col(col_name)
                        .unwrap_or("")
                        .to_string()
                })
                .collect();

            // Extract B fields
            let b_raw_fields: Vec<String> = READINFO_DATA_COLS
                .iter()
                .map(|col_name| {
                    reader_b
                        .get_col(col_name)
                        .unwrap_or("")
                        .to_string()
                })
                .collect();

            // Parse typed values from A
            let a_as_max: i64 = reader_a.get_col("AS_Max").unwrap_or("0").parse().unwrap_or(0);
            let a_ms_max: i64 = reader_a.get_col("ms_Max").unwrap_or("0").parse().unwrap_or(0);
            let a_seqid_max: f64 = reader_a
                .get_col("seqid_Max")
                .unwrap_or("NaN")
                .parse()
                .unwrap_or(f64::NAN);
            let a_aln_len_max: u64 = reader_a
                .get_col("Query_Aln_Len_Max")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let a_cov_max: f64 = reader_a
                .get_col("Query_Aln_Cov_Max")
                .unwrap_or("NaN")
                .parse()
                .unwrap_or(f64::NAN);
            let a_num_bp_inserted: u64 = reader_a
                .get_col("num_bp_inserted")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let a_n_ins_bases: u64 = reader_a
                .get_col("N_Insertion_Bases")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let a_n_del_bases: u64 = reader_a
                .get_col("N_Deletion_Bases")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let a_n_sub_bases: u64 = reader_a
                .get_col("N_Substitution_Bases")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let a_n_sc_start: u64 = reader_a
                .get_col("N_SoftClipped_Bases_Start")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let a_n_sc_end: u64 = reader_a
                .get_col("N_SoftClipped_Bases_End")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let a_junc_count: usize = reader_a
                .get_col("JuncCount")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let a_junctions = reader_a.get_col("junctions").unwrap_or("");

            // Parse typed values from B
            let b_as_max: i64 = reader_b.get_col("AS_Max").unwrap_or("0").parse().unwrap_or(0);
            let b_ms_max: i64 = reader_b.get_col("ms_Max").unwrap_or("0").parse().unwrap_or(0);
            let b_seqid_max: f64 = reader_b
                .get_col("seqid_Max")
                .unwrap_or("NaN")
                .parse()
                .unwrap_or(f64::NAN);
            let b_aln_len_max: u64 = reader_b
                .get_col("Query_Aln_Len_Max")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let b_cov_max: f64 = reader_b
                .get_col("Query_Aln_Cov_Max")
                .unwrap_or("NaN")
                .parse()
                .unwrap_or(f64::NAN);
            let b_num_bp_inserted: u64 = reader_b
                .get_col("num_bp_inserted")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let b_n_ins_bases: u64 = reader_b
                .get_col("N_Insertion_Bases")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let b_n_del_bases: u64 = reader_b
                .get_col("N_Deletion_Bases")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let b_n_sub_bases: u64 = reader_b
                .get_col("N_Substitution_Bases")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let b_n_sc_start: u64 = reader_b
                .get_col("N_SoftClipped_Bases_Start")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let b_n_sc_end: u64 = reader_b
                .get_col("N_SoftClipped_Bases_End")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let b_junc_count: usize = reader_b
                .get_col("JuncCount")
                .unwrap_or("0")
                .parse()
                .unwrap_or(0);
            let b_junctions = reader_b.get_col("junctions").unwrap_or("");

            // Compute metrics
            let as_diff = b_as_max - a_as_max;
            let ms_diff = b_ms_max - a_ms_max;
            let as_ratio = safe_ratio_i64(b_as_max, a_as_max);
            let ms_ratio = safe_ratio_i64(b_ms_max, a_ms_max);
            let seqid_diff = b_seqid_max - a_seqid_max;
            let qal_diff = b_aln_len_max as i64 - a_aln_len_max as i64;
            let qac_diff = b_cov_max - a_cov_max;
            let bp_ins_diff = b_num_bp_inserted as i64 - a_num_bp_inserted as i64;
            let bp_ins_ratio = safe_ratio_u64(b_num_bp_inserted, a_num_bp_inserted);
            let n_ins_diff = b_n_ins_bases as i64 - a_n_ins_bases as i64;
            let n_ins_ratio = safe_ratio_u64(b_n_ins_bases, a_n_ins_bases);
            let n_del_diff = b_n_del_bases as i64 - a_n_del_bases as i64;
            let n_del_ratio = safe_ratio_u64(b_n_del_bases, a_n_del_bases);
            let n_sub_diff = b_n_sub_bases as i64 - a_n_sub_bases as i64;
            let n_sub_ratio = safe_ratio_u64(b_n_sub_bases, a_n_sub_bases);
            let n_sc_start_diff = b_n_sc_start as i64 - a_n_sc_start as i64;
            let n_sc_end_diff = b_n_sc_end as i64 - a_n_sc_end as i64;

            // Junction metrics
            let juncs_a = parse_junction_str(a_junctions);
            let juncs_b = parse_junction_str(b_junctions);
            let junction_distance_val = junction_distance(&juncs_a, &juncs_b); // col 72 (unchanged)

            // Count-difference kept ONLY to preserve Junc_Dist_V2's original values.
            let n_junc_count_diff = (a_junc_count as i64 - b_junc_count as i64).unsigned_abs();
            let junc_dist_v2 = 50 * n_junc_count_diff; // col 74 (unchanged)

            // Set-based overlap stats (deduped both sides → internally consistent).
            let (n_matched, n_only_a, n_only_b) = junction_set_stats(&juncs_a, &juncs_b); // 75/76/77
            let n_unmatched = n_only_a + n_only_b; // col 73 (REDEFINED: set symmetric difference)

            // Genomic-junction set comparison (flag-gated).
            // Cross-chromosome safety is automatic: chrom is part of each tuple element.
            let (g_matched, g_only_a, g_only_b, g_unmatched) = if args.compare_genomic_junctions {
                let a_genomic = reader_a.get_col("genomic_junctions").unwrap_or("");
                let b_genomic = reader_b.get_col("genomic_junctions").unwrap_or("");
                let gj_a = parse_genomic_junction_str(a_genomic);
                let gj_b = parse_genomic_junction_str(b_genomic);
                let (m, oa, ob) = genomic_junction_set_stats(&gj_a, &gj_b);
                (m, oa, ob, oa + ob)
            } else {
                (0, 0, 0, 0)
            };

            // Write output row
            write!(out, "{}\t{}", key_a.name, key_a.len)?;

            // A fields (with escaping for special characters)
            for f in &a_raw_fields {
                let escaped = escape_tsv_field(f);
                write!(out, "\t{escaped}")?;
            }
            // B fields (with escaping for special characters)
            for f in &b_raw_fields {
                let escaped = escape_tsv_field(f);
                write!(out, "\t{escaped}")?;
            }

            // Comparison metrics
            write!(
                out,
                "\t{as_diff}\t{ms_diff}\t{as_ratio}\t{ms_ratio}\t{}\t{qal_diff}\t{}\t\
                 {bp_ins_diff}\t{bp_ins_ratio}\t\
                 {n_ins_diff}\t{n_ins_ratio}\t\
                 {n_del_diff}\t{n_del_ratio}\t\
                 {n_sub_diff}\t{n_sub_ratio}\t\
                 {n_sc_start_diff}\t{n_sc_end_diff}\t\
                 {junction_distance_val}\t{n_unmatched}\t{junc_dist_v2}\t{n_matched}\
                 \t{n_only_a}\t{n_only_b}",
                fmt_float(seqid_diff),
                fmt_float(qac_diff),
            )?;
            if args.compare_genomic_junctions {
                write!(
                    out,
                    "\t{g_matched}\t{g_unmatched}\t{g_only_a}\t{g_only_b}",
                )?;
            }
            writeln!(out)?;

            n_merged += 1;

            // Print progress
            if n_merged % 100000 == 0 {
                eprintln!("[INFO] Processed {} matched records...", n_merged);
            }

            // Advance both
            reader_a.advance()?;
            reader_b.advance()?;
        } else if key_a < key_b {
            // A is behind: advance A
            n_a_total += 1;
            reader_a.advance()?;
        } else {
            // B is behind: advance B
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
    // Inner join: only reads matched on (Read_Name, Read_Len) are written.
    // Reads present in only one file are dropped here, but reported so the
    // user knows how many rows were lost from each side.
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
