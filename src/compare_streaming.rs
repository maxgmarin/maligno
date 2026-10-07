//! `toolkit merge-readinfo`: two readinfo tables → the alignment comparison table.
//!
//! The pairing, classification and row writing are `compare`'s own
//! ([`crate::compare::run_merge`]); this module only supplies a readinfo-TSV
//! [`ReadSource`], so both commands produce the same table and summary from the
//! same per-read rows. Reads are matched on `Read_Name`; both inputs must be
//! sorted by `Read_Name` in byte order (`LC_ALL=C`).

use std::io::{self, BufRead, Write};

use anyhow::{Context, Result};

use crate::compare::{run_merge, ReadSource, SourceRead};
use crate::compare_summary::CompareSummary;
use crate::comparison_row::write_compare_header;
use crate::io_utils::{open_input, open_output};
use crate::parquet_out::{is_parquet_path, ComparisonParquetWriter};

// ── CLI args ─────────────────────────────────────────────────────────────────

/// Validate a `--label-a` / `--label-b` value.
///
/// The label is written as a **data value** (the `Label_A` /
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
    /// Readinfo TSV A (must be sorted by Read_Name)
    #[arg(short = 'a', long = "readinfo-a", value_name = "readinfo_a.tsv[.gz]")]
    pub readinfo_a: String,

    /// Readinfo TSV B (must be sorted by Read_Name)
    #[arg(short = 'b', long = "readinfo-b", value_name = "readinfo_b.tsv[.gz]")]
    pub readinfo_b: String,

    /// Name for dataset A, recorded in the output's `Label_A` column
    #[arg(long = "label-a", value_name = "LABEL", default_value = "SetA",
          value_parser = validate_set_label)]
    pub label_a: String,

    /// Name for dataset B, recorded in the output's `Label_B` column
    #[arg(long = "label-b", value_name = "LABEL", default_value = "SetB",
          value_parser = validate_set_label)]
    pub label_b: String,

    /// Output alignment comparison table (.parquet or .tsv[.gz])
    #[arg(short = 'o', long = "output", value_name = "compare.parquet|compare.tsv[.gz]")]
    pub output: String,

    /// Compare the shared intersection of reads instead of erroring when the two
    /// inputs do not carry the exact same "Read_Name" set.
    #[arg(long = "allow-id-mismatch")]
    pub allow_id_mismatch: bool,
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

    pub(crate) fn get_col_idx(&self, col_name: &str) -> Option<usize> {
        self.col_map.get(col_name).copied()
    }
}

// ── Readinfo-TSV read source ────────────────────────────────────────────────

/// A readinfo TSV as a [`ReadSource`]: one row per read, columns taken from the
/// file's own header, so column order and extra columns don't matter.
struct TsvReadSource {
    reader: ReadInfoReader,
    idx_name: usize,
    idx_len: usize,
}

impl TsvReadSource {
    fn open(path: &str, side: &str) -> Result<Self> {
        let reader = ReadInfoReader::new(path)?;
        let idx_name = reader
            .get_col_idx("Read_Name")
            .with_context(|| format!("missing Read_Name in {side}"))?;
        let idx_len = reader
            .get_col_idx("Read_Len")
            .with_context(|| format!("missing Read_Len in {side}"))?;
        Ok(TsvReadSource { reader, idx_name, idx_len })
    }
}

impl ReadSource for TsvReadSource {
    fn columns(&self) -> Vec<String> {
        self.reader.header.clone()
    }

    fn next_read(&mut self) -> Result<Option<SourceRead>> {
        let Some(fields) = self.reader.current() else {
            return Ok(None);
        };
        let name = fields
            .get(self.idx_name)
            .context("readinfo row is missing its Read_Name field")?
            .clone();
        let len: u64 = fields
            .get(self.idx_len)
            .context("readinfo row is missing its Read_Len field")?
            .parse()
            .context("failed to parse Read_Len")?;
        let line = fields.join("\t");
        self.reader.advance()?;
        Ok(Some(SourceRead { name, len, line }))
    }
}

// ── Command ─────────────────────────────────────────────────────────────────

pub fn run(args: &MergeReadinfoArgs) -> Result<()> {
    eprintln!("[INFO] Opening readinfo files...");
    eprintln!("  A: {}", args.readinfo_a);
    eprintln!("  B: {}", args.readinfo_b);

    let mut src_a = TsvReadSource::open(&args.readinfo_a, "A")?;
    let mut src_b = TsvReadSource::open(&args.readinfo_b, "B")?;

    // Output format follows the extension, matching how `open_output` already
    // dispatches on `.gz`: `-o x.parquet` writes Parquet, anything else (including
    // `-` for stdout) writes the TSV. Exactly one of the two is produced here —
    // unlike `compare`, which owns its filenames and can write both. The unused
    // side is an `io::sink()`, so no stray empty file is created.
    let mut parquet = if is_parquet_path(&args.output) {
        eprintln!("[INFO] Output format: Parquet (from the '.parquet' extension)");
        Some(ComparisonParquetWriter::new(
            std::fs::File::create(&args.output)
                .with_context(|| format!("cannot create '{}'", args.output))?,
        )?)
    } else {
        None
    };
    let mut out: Box<dyn Write> = match parquet {
        Some(_) => Box::new(io::sink()),
        None => {
            let mut w = open_output(Some(&args.output))?;
            write_compare_header(&mut w)?;
            w
        }
    };

    eprintln!("[INFO] Starting comparison...");
    let mut summary = CompareSummary::default();
    let (n_matched, n_a_only, n_b_only) = run_merge(
        &mut src_a,
        &mut src_b,
        &mut out,
        parquet.as_mut(),
        args.allow_id_mismatch,
        " Re-run with --allow-id-mismatch to compare only the reads present in both files.",
        &args.label_a,
        &args.label_b,
        &mut summary,
        &mut None,
    )?;

    // Parquet must be finished explicitly — its footer (schema + row-group index
    // + metadata) is written on close.
    out.flush()?;
    if let Some(pq) = parquet {
        let n = pq.finish()?;
        eprintln!("[INFO] Wrote {n} rows to {}", args.output);
    }

    eprintln!(
        "  shared: {n_matched}   only in {}: {n_a_only}   only in {}: {n_b_only}",
        args.label_a, args.label_b
    );
    summary.render_stderr_brief(&args.label_a, &args.label_b);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Write a small readinfo-style TSV and return its path.
    fn write_tsv(name: &str, rows: &[&str]) -> String {
        let path = std::env::temp_dir().join(format!(
            "maligno_test_merge_readinfo_{}_{name}.tsv",
            std::process::id()
        ));
        let mut body = String::from("Read_Name\tRead_Len\tTargetChr\n");
        for r in rows {
            body.push_str(r);
            body.push('\n');
        }
        std::fs::write(&path, body).unwrap();
        path.to_string_lossy().into_owned()
    }

    fn args(a: &str, b: &str, out: &str, allow_id_mismatch: bool) -> MergeReadinfoArgs {
        MergeReadinfoArgs {
            readinfo_a: a.to_string(),
            readinfo_b: b.to_string(),
            label_a: "A".to_string(),
            label_b: "B".to_string(),
            output: out.to_string(),
            allow_id_mismatch,
        }
    }

    /// (Read_Name, Read_Len) of every data row in a comparison TSV.
    fn rows(out: &str) -> Vec<(String, String)> {
        std::fs::read_to_string(out)
            .unwrap()
            .lines()
            .skip(1)
            .map(|l| {
                let mut f = l.split('\t');
                (f.next().unwrap().to_string(), f.next().unwrap().to_string())
            })
            .collect()
    }

    #[test]
    fn matches_on_read_name_only_like_compare() {
        // r1's Read_Len differs between the sides: it is still compared, and the
        // output Read_Len is side A's, as in `compare`.
        let a = write_tsv("len_a", &["r1\t100\tchr1", "r2\t50\tchr2"]);
        let b = write_tsv("len_b", &["r1\t101\tchr1", "r2\t50\t*"]);
        let out = write_tsv("len_out", &[]);
        run(&args(&a, &b, &out, false)).unwrap();
        assert_eq!(
            rows(&out),
            vec![("r1".into(), "100".into()), ("r2".into(), "50".into())]
        );
    }

    #[test]
    fn read_id_mismatch_errors_unless_allowed() {
        let a = write_tsv("idm_a", &["r1\t100\tchr1", "r2\t50\tchr2"]);
        let b = write_tsv("idm_b", &["r1\t100\tchr1", "r3\t70\tchr3"]);
        let out = write_tsv("idm_out", &[]);

        let err = run(&args(&a, &b, &out, false)).unwrap_err().to_string();
        assert!(err.contains("--allow-id-mismatch"), "unexpected error: {err}");

        run(&args(&a, &b, &out, true)).unwrap();
        assert_eq!(rows(&out), vec![("r1".into(), "100".into())]);
    }
}
