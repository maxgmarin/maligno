//! Provides the `compare-pipeline merge-readinfo` command (two readinfo TSVs →
//! comparison) plus
//! the merge-join machinery it shares with the primary `compare` command
//! (`compare.rs`): `ReadKey` and `ReadInfoReader`.
//!
//! The comparison **table schema** — the column lists, the row type and the two
//! TSV writers — lives in `comparison_row.rs`, not here. This module only decides
//! *which* pairs of readinfo rows get compared; `comparison_row` decides what a
//! comparison row contains.
//!
//! Algorithm: streaming two-pointer merge-join, O(1) memory, O(|A|+|B|) time.
//! - Stream through both files simultaneously; on key match, emit a row.
//! - On mismatch: error by default, or (with `--ignore-row-mismatch`) advance the
//!   pointer with the smaller key.
//!
//! REQUIREMENT: both input readinfo files must be sorted by (Read_Name, Read_Len).

use std::io::{BufRead, Write};

use anyhow::{bail, Context, Result};

use crate::comparison_row::{write_compare_header, ComparisonRow};
use crate::parquet_out::{is_parquet_path, ComparisonParquetWriter};
use crate::io_utils::{open_input, open_output};

// ── CLI args ─────────────────────────────────────────────────────────────────

/// Validate a `--label-a` / `--label-b` value.
///
/// Since v0.13.0 the label is written as a **data value** (the `Label_A` /
/// `Label_B` columns) on every comparison row, and it is also interpolated into
/// per-set output filenames, so a tab/newline would corrupt the TSV and a path
/// separator would redirect output. Underscores are fine — side identity comes
/// from the fixed `_A` / `_B` column suffixes, not from parsing the label.
pub(crate) fn validate_set_label(s: &str) -> Result<String, String> {
    if s.is_empty() {
        return Err("label must not be empty".to_string());
    }
    if let Some(bad) = s.chars().find(|c| matches!(c, '\t' | '\n' | '\r' | '/' | '\\')) {
        return Err(format!(
            "label must not contain {bad:?} (tabs/newlines would corrupt the output \
             table; path separators would redirect the per-set output files)"
        ));
    }
    Ok(s.to_string())
}

#[derive(clap::Args, Debug)]
pub struct MergeReadinfoArgs {
    /// Readinfo TSV A (must be sorted by Read_Name, Read_Len)
    #[arg(short = 'a', long = "readinfo-a", value_name = "readinfo_a.tsv")]
    pub readinfo_a: String,

    /// Readinfo TSV B (must be sorted by Read_Name, Read_Len)
    #[arg(short = 'b', long = "readinfo-b", value_name = "readinfo_b.tsv")]
    pub readinfo_b: String,

    /// Name for dataset A, recorded in the output's `Label_A` column
    #[arg(long = "label-a", value_name = "LABEL", default_value = "SetA",
          value_parser = validate_set_label)]
    pub label_a: String,

    /// Name for dataset B, recorded in the output's `Label_B` column
    #[arg(long = "label-b", value_name = "LABEL", default_value = "SetB",
          value_parser = validate_set_label)]
    pub label_b: String,

    /// Output comparison TSV file ('.gz' for gzip)
    #[arg(short = 'o', long = "output", value_name = "compare.tsv[.gz]")]
    pub output: String,

    /// Skip reads that appear in only one file instead of stopping with an error.
    /// Both readinfo files must still be lex-sorted by Read_Name for the skip
    /// heuristic to work correctly. Unmatched reads are counted in the summary.
    #[arg(long = "ignore-row-mismatch")]
    pub ignore_row_mismatch: bool,
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

pub fn run(args: &MergeReadinfoArgs) -> Result<()> {
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

    // Output format follows the extension, matching how `open_output` already
    // dispatches on `.gz`: `-o x.parquet` writes Parquet, anything else (including
    // `-` for stdout) writes the TSV. Exactly one of the two is produced here —
    // unlike `compare`, which owns its filenames and can write both.
    let mut parquet = if is_parquet_path(&args.output) {
        eprintln!("[INFO] Output format: Parquet (from the '.parquet' extension)");
        Some(ComparisonParquetWriter::new(
            std::fs::File::create(&args.output)
                .with_context(|| format!("cannot create '{}'", args.output))?,
        )?)
    } else {
        None
    };

    // The TSV writer is only opened when Parquet was not selected, so we never
    // create a stray empty .tsv alongside a .parquet.
    let mut out = match parquet {
        Some(_) => None,
        None => {
            let mut w = open_output(Some(&args.output))?;
            write_compare_header(&mut w)?;
            Some(w)
        }
    };

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

            // One construction, whichever writer is active.
            let row = ComparisonRow::build(
                &key_a.name,
                key_a.len,
                &args.label_a,
                &args.label_b,
                |c| reader_a.get_col(c).unwrap_or(""),
                |c| reader_b.get_col(c).unwrap_or(""),
            );
            match (&mut out, &mut parquet) {
                (Some(w), _) => row.write_tsv_row(w)?,
                (None, Some(pq)) => pq.append(&row)?,
                (None, None) => unreachable!("one writer is always active"),
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

    // Close whichever writer is active. Parquet must be finished explicitly — its
    // footer (schema + row-group index + metadata) is written on close.
    if let Some(mut w) = out {
        w.flush()?;
    }
    if let Some(pq) = parquet {
        let n = pq.finish()?;
        eprintln!("[INFO] Wrote {n} rows to {}", args.output);
    }

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
