//! TSV-vs-Parquet dispatch for reading a two-sided comparison table.
//!
//! Shared by `find-aln-diff` and `compare-pipeline summary`: both stream an
//! existing `compare` / `compare-pipeline merge-readinfo` table row-by-row,
//! needing only a subset of its columns by name, insensitive to which format
//! produced the file.

use std::io::BufRead;

use anyhow::{bail, Context, Result};

use crate::io_utils::open_input;
use crate::parquet_out::{is_parquet_path, ParquetRowReader};

/// Serialization of a comparison-table `--input`. `auto` (default) picks
/// Parquet for a `.parquet`-named path and TSV otherwise — the same
/// extension convention `compare`/`compare-pipeline merge-readinfo` already
/// use on the output side.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum InputFormat {
    #[default]
    Auto,
    Tsv,
    Parquet,
}

/// One data row's cells, in header/column order, from either input format.
pub enum RowSource {
    Tsv(Box<dyn BufRead>),
    Parquet(ParquetRowReader),
}

impl RowSource {
    /// The header/column names, in file order.
    pub fn header(&mut self) -> Result<Vec<String>> {
        match self {
            RowSource::Tsv(r) => {
                let mut header = String::new();
                if r.read_line(&mut header)? == 0 {
                    bail!("comparison table is empty");
                }
                Ok(header
                    .trim_end_matches(['\n', '\r'])
                    .split('\t')
                    .map(String::from)
                    .collect())
            }
            RowSource::Parquet(pr) => Ok(pr.columns().to_vec()),
        }
    }

    /// The next data row as owned cells, in header order. `None` at EOF. Blank
    /// TSV lines are skipped, matching the pre-Parquet behavior.
    pub fn next_row(&mut self) -> Result<Option<Vec<String>>> {
        match self {
            RowSource::Tsv(r) => loop {
                let mut line = String::new();
                if r.read_line(&mut line)? == 0 {
                    return Ok(None);
                }
                let line = line.trim_end_matches(['\n', '\r']);
                if line.is_empty() {
                    continue;
                }
                return Ok(Some(line.split('\t').map(String::from).collect()));
            },
            RowSource::Parquet(pr) => pr.next_row(),
        }
    }
}

/// Open a two-sided comparison table for row-by-row reading. `format`
/// resolves `Auto` by sniffing `path`'s `.parquet` extension (via
/// `is_parquet_path`); `-` (stdin) is rejected up front for Parquet, since
/// `ParquetRecordBatchReaderBuilder` needs a seekable file. `columns` is
/// forwarded to `ParquetRowReader::open` (see its doc for what `None` vs
/// `Some` means) and ignored for TSV, which is always read in full.
pub fn open_table(path: &str, format: InputFormat, columns: Option<&[&str]>) -> Result<RowSource> {
    let want_parquet = match format {
        InputFormat::Parquet => true,
        InputFormat::Tsv => false,
        InputFormat::Auto => path != "-" && is_parquet_path(path),
    };
    if want_parquet && path == "-" {
        bail!("--input-format parquet cannot read from stdin ('-'); pass a Parquet file path");
    }
    if want_parquet {
        Ok(RowSource::Parquet(
            ParquetRowReader::open(path, columns)
                .with_context(|| format!("opening comparison table '{path}'"))?,
        ))
    } else {
        Ok(RowSource::Tsv(
            open_input(path).with_context(|| format!("opening comparison table '{path}'"))?,
        ))
    }
}
