//! Parquet output for the comparison table.
//!
//! Selected by `compare --format tsv|parquet|both` (default `both`), and by
//! `compare-readinfo -o x.parquet`, which picks Parquet from the extension.
//!
//! Why it exists: the comparison table is written once and re-read many times —
//! by the notebooks, and by `find-query-diff` re-run standalone with
//! `--compare-by junctions`. Parquet stores each column separately, so a reader
//! touching a few of the 96 columns never pays for the rest. Note the cost being
//! avoided is **splitting every row into 96 fields**, not decompression:
//! gunzipping the whole 408 MB table takes only 0.22 s.
//!
//! Measured on the 507,365-row Splice-vs-SpliceHQ table (DuckDB 1.5.5, default CSV
//! sampling, best of three): the 16 columns `find-query-diff` needs take 2.01 s
//! from `tsv.gz` versus 0.55 s from Parquet (3.6×); a two-column aggregate drops
//! from 0.85 s to 0.02 s (48×); a full 96-column scan gains least, 3.97 s to
//! 2.31 s (1.7×). Writing is faster too — 4.5 s versus 10.8 s, measured with
//! maligno itself. The cost is size: Parquet+zstd is ~31% larger here (85 MB vs
//! 65 MB), because six long-string columns dominate the table and gzip compresses
//! across the whole row stream where Parquet compresses each column alone.
//!
//! (Benchmark those numbers with DuckDB's *default* sampling. `sample_size=-1`
//! forces a full-file type-inference scan before any row is read, which maligno
//! never pays and which inflates the TSV side 2–4×.)
//!
//! **The schema is derived, never restated.** Column names and order come from
//! `comparison_row`'s `READINFO_DATA_COLS` and `comparison_col_names()`, the same
//! lists that drive the TSV header, so the two cannot drift. This module adds only
//! the name → type mapping, and a test asserts that mapping is exhaustive.
//!
//! Layout is **flat** — 96 top-level columns named exactly as the TSV header —
//! rather than nested `aln_a`/`aln_b`/`diff` structs, because the consumers are
//! pandas/polars/DuckDB, where `read_parquet` should be a drop-in for `read_csv`.

use std::fs::File;
use std::io::Write;
use std::sync::Arc;

use anyhow::{bail, Context, Result};
use arrow_array::builder::{
    BooleanBuilder, Float64Builder, Int64Builder, StringBuilder, UInt64Builder,
};
use arrow_array::{
    Array, ArrayRef, BooleanArray, Float64Array, Int64Array, RecordBatch, StringArray,
    UInt64Array,
};
use arrow_schema::{DataType, Field, Schema, SchemaRef};
use parquet::arrow::arrow_reader::{ParquetRecordBatchReader, ParquetRecordBatchReaderBuilder};
use parquet::arrow::{ArrowWriter, ProjectionMask};
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::metadata::KeyValue;
use parquet::file::properties::WriterProperties;

use crate::comparison_row::{
    comparison_col_names, ComparisonRow, DiffValue, READINFO_DATA_COLS,
};
use crate::io_utils::fmt_float;

/// Bumped when the Parquet column set, types or null semantics change. Recorded in
/// the file footer so a reader can detect a schema it does not understand.
const SCHEMA_VERSION: &str = "1";

/// Rows buffered before a row group is flushed. Bounds peak memory at
/// O(row group) rather than O(reads).
///
/// 20k was chosen by measurement on the 507,365-row table: row-group size is a
/// pure memory dial here, with no cost in speed or size. 20k vs 100k rows gave
/// 4.53s vs 4.57s to write (noise) and 86.3 MB either way, but 306 MB vs 490 MB
/// peak RSS. There is a further ~270 MB of fixed overhead from parquet's 96
/// per-column writers and their zstd contexts that no row-group choice affects.
const ROW_GROUP_ROWS: usize = 20_000;

// ── Column types ─────────────────────────────────────────────────────────────

/// The Arrow type of one comparison-table column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ColType {
    Str,
    U64,
    I64,
    F64,
    Bool,
}

impl ColType {
    fn arrow(self) -> DataType {
        match self {
            ColType::Str => DataType::Utf8,
            ColType::U64 => DataType::UInt64,
            ColType::I64 => DataType::Int64,
            ColType::F64 => DataType::Float64,
            ColType::Bool => DataType::Boolean,
        }
    }
}

/// Type of a per-side column, by its unsuffixed name. `None` means the name is not
/// a known per-side column — `per_side_types_are_exhaustive` makes that a test
/// failure rather than a runtime surprise.
fn per_side_type(name: &str) -> Option<ColType> {
    Some(match name {
        // Locus & span.
        "TargetChr" | "Strand" => ColType::Str,
        "Target_Start" | "Target_End" | "Query_Start" | "Query_End" => ColType::U64,
        // Alignment selection & score.
        "MQ_Best" | "Num_Aln" | "Num_Aln_MaxScore" => ColType::U64,
        "AS_Max" | "ms_Max" => ColType::I64,
        // Identity & coverage.
        "seqid_Max" | "Query_Aln_Cov_Max" => ColType::F64,
        "Query_Aln_Len_Max" => ColType::U64,
        // Junction and cs-derived counts.
        "JuncCount"
        | "N_Splice_Junction_Events"
        | "N_Splice_Junction_Bases"
        | "N_Match_Events"
        | "N_Match_Bases"
        | "N_Substitution_Events"
        | "N_Substitution_Bases"
        | "N_Insertion_Events"
        | "N_Insertion_Bases"
        | "N_Deletion_Events"
        | "N_Deletion_Bases"
        | "N_SoftClipped_Events"
        | "N_SoftClipped_Bases_Start"
        | "N_SoftClipped_Bases_End" => ColType::U64,
        // Long strings.
        "junctions" | "genomic_junctions" | "cs" => ColType::Str,
        _ => return None,
    })
}

/// Type of a comparison-block column, by name.
fn comparison_type(name: &str) -> Option<ColType> {
    Some(match name {
        "Strand_Match" => ColType::Bool,
        "seqid_Diff" | "QueryAlnCov_Diff" => ColType::F64,
        "AS_Ratio" | "ms_Ratio" => ColType::F64,
        "N_Substitution_Bases_Ratio" | "N_Insertion_Bases_Ratio" | "N_Deletion_Bases_Ratio" => {
            ColType::F64
        }
        "QueryAlnLen_Diff"
        | "AS_Diff"
        | "ms_Diff"
        | "N_Substitution_Bases_Diff"
        | "N_Insertion_Bases_Diff"
        | "N_Deletion_Bases_Diff"
        | "N_SoftClipped_Bases_Start_Diff"
        | "N_SoftClipped_Bases_End_Diff" => ColType::I64,
        "N_Matched_Junctions"
        | "N_Unmatched_Junctions"
        | "N_Junctions_OnlyA"
        | "N_Junctions_OnlyB"
        | "Junction_Distance"
        | "Junc_Dist_V2"
        | "Genomic_N_Matched_Junctions"
        | "Genomic_N_Unmatched_Junctions"
        | "Genomic_N_Junctions_OnlyA"
        | "Genomic_N_Junctions_OnlyB" => ColType::U64,
        "Junctions_OnlyA" | "Junctions_OnlyB" | "Genomic_Junctions_OnlyA"
        | "Genomic_Junctions_OnlyB" => ColType::Str,
        _ => return None,
    })
}

/// Build the flat 96-column schema, deriving names and order from
/// `comparison_row`'s lists.
///
/// **Nullability** encodes what a null *means*: "undefined", nothing else.
/// - Per-side columns are nullable: they are transported as strings, so a value
///   can be absent (a readinfo table missing a column) or unparseable. Writing a
///   fabricated `0` there is exactly the corruption this format is meant to avoid.
/// - Float columns are nullable: `NaN` means "undefined" (identity with no
///   alignment, a ratio over a zero denominator) and is stored as a null.
/// - The keys, `Strand_Match` and the computed integer comparison columns are
///   **not** nullable — they are always derived and always present.
///
/// Note `*` (`TargetChr`/`Strand`) and `()` (empty junction set) are **kept as
/// strings**, not nulled: they are meaningful values, and `*` is the mapping
/// indicator downstream code tests against.
pub(crate) fn build_schema() -> SchemaRef {
    let mut fields = Vec::with_capacity(96);
    fields.push(Field::new("Read_Name", DataType::Utf8, false));
    fields.push(Field::new("Read_Len", DataType::UInt64, false));
    fields.push(Field::new("Label_A", DataType::Utf8, false));
    fields.push(Field::new("Label_B", DataType::Utf8, false));
    for side in ["A", "B"] {
        for col in READINFO_DATA_COLS {
            let ty = per_side_type(col)
                .unwrap_or_else(|| panic!("no Arrow type declared for per-side column {col:?}"));
            fields.push(Field::new(format!("{col}_{side}"), ty.arrow(), true));
        }
    }
    for col in comparison_col_names() {
        let ty = comparison_type(col)
            .unwrap_or_else(|| panic!("no Arrow type declared for comparison column {col:?}"));
        // Only the float columns can be undefined; everything else is computed.
        let nullable = ty == ColType::F64;
        fields.push(Field::new(col, ty.arrow(), nullable));
    }
    Arc::new(Schema::new(fields))
}

// ── Builders ─────────────────────────────────────────────────────────────────

/// One column's Arrow builder, tagged so the appenders stay type-checked.
enum ColBuilder {
    Str(StringBuilder),
    U64(UInt64Builder),
    I64(Int64Builder),
    F64(Float64Builder),
    Bool(BooleanBuilder),
}

impl ColBuilder {
    fn new(ty: ColType) -> Self {
        match ty {
            ColType::Str => ColBuilder::Str(StringBuilder::new()),
            ColType::U64 => ColBuilder::U64(UInt64Builder::new()),
            ColType::I64 => ColBuilder::I64(Int64Builder::new()),
            ColType::F64 => ColBuilder::F64(Float64Builder::new()),
            ColType::Bool => ColBuilder::Bool(BooleanBuilder::new()),
        }
    }

    /// Append a per-side transport value, which arrives as a `&str`.
    ///
    /// An empty string becomes a null for every type: for `cs` that is the
    /// unmapped case, and for a numeric column it means the value was absent.
    /// An unparseable numeric also becomes a null rather than a fabricated zero.
    /// `NaN` in a float column becomes a null, since it means "undefined".
    fn append_raw(&mut self, s: &str) {
        match self {
            ColBuilder::Str(b) => {
                if s.is_empty() {
                    b.append_null()
                } else {
                    b.append_value(s)
                }
            }
            ColBuilder::U64(b) => b.append_option(s.parse::<u64>().ok()),
            ColBuilder::I64(b) => b.append_option(s.parse::<i64>().ok()),
            ColBuilder::F64(b) => {
                b.append_option(s.parse::<f64>().ok().filter(|v| !v.is_nan()))
            }
            ColBuilder::Bool(b) => b.append_option(s.parse::<bool>().ok()),
        }
    }

    /// Append an already-typed comparison-block value.
    fn append_diff(&mut self, v: DiffValue<'_>) -> Result<()> {
        match (self, v) {
            (ColBuilder::Bool(b), DiffValue::Bool(x)) => b.append_value(x),
            (ColBuilder::I64(b), DiffValue::I64(x)) => b.append_value(x),
            (ColBuilder::U64(b), DiffValue::U64(x)) => b.append_value(x),
            // NaN means "undefined" — store it as a null, not a float NaN.
            (ColBuilder::F64(b), DiffValue::F64(x)) => {
                b.append_option(if x.is_nan() { None } else { Some(x) })
            }
            (ColBuilder::Str(b), DiffValue::Str(x)) => b.append_value(x),
            (_, v) => anyhow::bail!(
                "comparison column type does not match the value produced for it ({v:?}) \
                 — the type table in parquet_out.rs and AlignmentDiff::values() disagree"
            ),
        }
        Ok(())
    }

    fn finish(&mut self) -> ArrayRef {
        match self {
            ColBuilder::Str(b) => Arc::new(b.finish()),
            ColBuilder::U64(b) => Arc::new(b.finish()),
            ColBuilder::I64(b) => Arc::new(b.finish()),
            ColBuilder::F64(b) => Arc::new(b.finish()),
            ColBuilder::Bool(b) => Arc::new(b.finish()),
        }
    }
}

// ── The writer ───────────────────────────────────────────────────────────────

/// Streams `ComparisonRow`s into a Parquet file, one row group per
/// `ROW_GROUP_ROWS` rows.
///
/// Rows are appended into per-column builders as they stream past, so no row is
/// retained beyond its row group and peak memory does not scale with read count.
pub(crate) struct ComparisonParquetWriter<W: Write + Send> {
    writer: ArrowWriter<W>,
    schema: SchemaRef,
    cols: Vec<ColBuilder>,
    buffered: usize,
    rows_written: u64,
}

impl<W: Write + Send> ComparisonParquetWriter<W> {
    pub(crate) fn new(sink: W) -> Result<Self> {
        let schema = build_schema();
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(
                ZstdLevel::try_new(3).context("invalid zstd level")?,
            ))
            .set_key_value_metadata(Some(vec![KeyValue::new(
                "maligno_comparison_schema_version".to_string(),
                SCHEMA_VERSION.to_string(),
            )]))
            .build();
        let writer = ArrowWriter::try_new(sink, schema.clone(), Some(props))
            .context("creating the Parquet writer")?;
        let cols = Self::fresh_builders(&schema);
        Ok(Self { writer, schema, cols, buffered: 0, rows_written: 0 })
    }

    fn fresh_builders(schema: &SchemaRef) -> Vec<ColBuilder> {
        schema
            .fields()
            .iter()
            .map(|f| {
                ColBuilder::new(match f.data_type() {
                    DataType::Utf8 => ColType::Str,
                    DataType::UInt64 => ColType::U64,
                    DataType::Int64 => ColType::I64,
                    DataType::Float64 => ColType::F64,
                    DataType::Boolean => ColType::Bool,
                    other => unreachable!("build_schema never emits {other:?}"),
                })
            })
            .collect()
    }

    /// Append one row, flushing a row group when the buffer is full.
    pub(crate) fn append(&mut self, row: &ComparisonRow<'_>) -> Result<()> {
        let n = READINFO_DATA_COLS.len();
        self.cols[0].append_raw(row.read_name);
        match &mut self.cols[1] {
            ColBuilder::U64(b) => b.append_value(row.read_len),
            _ => unreachable!("Read_Len is UInt64"),
        }
        self.cols[2].append_raw(row.label_a);
        self.cols[3].append_raw(row.label_b);
        for (i, v) in row.aln_a.raw().iter().enumerate() {
            self.cols[4 + i].append_raw(v);
        }
        for (i, v) in row.aln_b.raw().iter().enumerate() {
            self.cols[4 + n + i].append_raw(v);
        }
        for (i, v) in row.diff.values().into_iter().enumerate() {
            self.cols[4 + 2 * n + i].append_diff(v)?;
        }
        self.buffered += 1;
        self.rows_written += 1;
        if self.buffered >= ROW_GROUP_ROWS {
            self.flush_row_group()?;
        }
        Ok(())
    }

    /// Turn the buffered rows into a row group and reset the builders.
    fn flush_row_group(&mut self) -> Result<()> {
        if self.buffered == 0 {
            return Ok(());
        }
        let arrays: Vec<ArrayRef> = self.cols.iter_mut().map(|c| c.finish()).collect();
        let batch = RecordBatch::try_new(self.schema.clone(), arrays)
            .context("assembling a Parquet row group")?;
        self.writer.write(&batch).context("writing a Parquet row group")?;
        // `finish()` already reset each builder; drop the buffered count to match.
        self.buffered = 0;
        Ok(())
    }

    /// Flush the partial row group and close the file. Returns the row count.
    pub(crate) fn finish(mut self) -> Result<u64> {
        self.flush_row_group()?;
        self.writer.close().context("closing the Parquet file")?;
        Ok(self.rows_written)
    }
}

/// Which serialization(s) of the comparison table to write.
///
/// Declared as a clap `ValueEnum`, so an unrecognised value is rejected before the
/// run starts, with the accepted values listed:
/// `error: invalid value 'xml' for '--format <FORMAT>' [possible values: tsv,
/// parquet, both]`.
///
/// The default is `both` for backward compatibility: existing scripts expect
/// the gzipped TSV. `find-query-diff` and `compare-summary` can read either
/// format (`--input-format`), so a Parquet-only run works standalone with
/// both of them.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum OutputFormat {
    /// Gzipped TSV only.
    Tsv,
    /// Parquet only.
    Parquet,
    /// Both TSV and Parquet.
    #[default]
    Both,
}

impl OutputFormat {
    pub(crate) fn writes_tsv(self) -> bool {
        matches!(self, OutputFormat::Tsv | OutputFormat::Both)
    }
    pub(crate) fn writes_parquet(self) -> bool {
        matches!(self, OutputFormat::Parquet | OutputFormat::Both)
    }
}

/// True when `path` should be written as Parquet.
///
/// Extension-based, matching the repo's existing convention: `open_input` /
/// `open_output` dispatch on `.gz` the same way, and nothing in maligno sniffs
/// magic bytes. `-` (stdout) has no `.parquet` suffix, so it keeps writing TSV.
pub(crate) fn is_parquet_path(path: &str) -> bool {
    path.ends_with(".parquet")
}

/// Render one Arrow cell back to the string a TSV reader would see for the
/// same value: a null becomes `NaN` for a float column (the "undefined"
/// convention from `append_raw`/`append_diff`) and an empty cell otherwise.
/// Used by the round-trip test below and by `find-query-diff`'s Parquet
/// reader (`ParquetRowReader`), so both input formats hand identical strings
/// to the same downstream parsing/comparison code.
pub(crate) fn arrow_cell_to_string(a: &dyn Array, row: usize) -> String {
    if a.is_null(row) {
        return match a.data_type() {
            DataType::Float64 => "NaN".to_string(),
            _ => String::new(),
        };
    }
    match a.data_type() {
        DataType::Utf8 => a.as_any().downcast_ref::<StringArray>().unwrap().value(row).to_string(),
        DataType::UInt64 => a.as_any().downcast_ref::<UInt64Array>().unwrap().value(row).to_string(),
        DataType::Int64 => a.as_any().downcast_ref::<Int64Array>().unwrap().value(row).to_string(),
        DataType::Float64 => fmt_float(a.as_any().downcast_ref::<Float64Array>().unwrap().value(row)),
        DataType::Boolean => a.as_any().downcast_ref::<BooleanArray>().unwrap().value(row).to_string(),
        other => panic!("unexpected type {other:?}"),
    }
}

/// Row-by-row Parquet reader for the comparison table's **input** side.
/// Cells come out already stringified via `arrow_cell_to_string`, so a
/// consumer written against tab-split TSV rows (column-name → index lookup,
/// then string comparisons/parses) works unchanged against either format.
///
/// Requires a seekable file, so `-` (stdin) is not supported — callers must
/// reject that combination before opening.
pub(crate) struct ParquetRowReader {
    reader: ParquetRecordBatchReader,
    columns: Vec<String>,
    current: Option<RecordBatch>,
    row_in_batch: usize,
}

impl ParquetRowReader {
    /// `columns`: `None` reads every column, in the file's schema order.
    /// `Some(names)` projects to just those columns — the reader never
    /// decodes the rest, which is where Parquet's per-column storage
    /// actually pays off (see the benchmark in this module's doc comment) —
    /// and bails with a clear error if any requested name isn't in the
    /// file's schema. Row order always matches the file's schema order
    /// restricted to the requested set, not the order `columns` was given
    /// in (this is `ProjectionMask`'s documented behavior).
    pub(crate) fn open(path: &str, columns: Option<&[&str]>) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("cannot open '{path}'"))?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .with_context(|| format!("'{path}' is not a valid Parquet file"))?;
        let all_columns: Vec<String> =
            builder.schema().fields().iter().map(|f| f.name().clone()).collect();

        let (columns, builder) = match columns {
            None => (all_columns, builder),
            Some(wanted) => {
                let missing: Vec<&str> = wanted
                    .iter()
                    .copied()
                    .filter(|w| !all_columns.iter().any(|c| c == w))
                    .collect();
                if !missing.is_empty() {
                    bail!(
                        "'{path}' is missing column(s) needed for this command: {}",
                        missing.join(", ")
                    );
                }
                let projected: Vec<String> =
                    all_columns.into_iter().filter(|c| wanted.contains(&c.as_str())).collect();
                let mask = ProjectionMask::columns(builder.parquet_schema(), wanted.iter().copied());
                (projected, builder.with_projection(mask))
            }
        };

        let reader = builder
            .build()
            .with_context(|| format!("reading Parquet file '{path}'"))?;
        Ok(Self { reader, columns, current: None, row_in_batch: 0 })
    }

    pub(crate) fn columns(&self) -> &[String] {
        &self.columns
    }

    /// The next data row as owned string cells, in column order. `None` at EOF.
    pub(crate) fn next_row(&mut self) -> Result<Option<Vec<String>>> {
        loop {
            if let Some(batch) = &self.current {
                if self.row_in_batch < batch.num_rows() {
                    let row: Vec<String> = (0..batch.num_columns())
                        .map(|c| arrow_cell_to_string(batch.column(c).as_ref(), self.row_in_batch))
                        .collect();
                    self.row_in_batch += 1;
                    return Ok(Some(row));
                }
            }
            match self.reader.next() {
                Some(batch) => {
                    self.current = Some(batch.context("reading a Parquet row group")?);
                    self.row_in_batch = 0;
                }
                None => return Ok(None),
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io_utils::escape_tsv_field;

    fn lookup(pairs: &[(&str, &'static str)], col: &str) -> &'static str {
        pairs.iter().find(|(k, _)| *k == col).map(|(_, v)| *v).unwrap_or("")
    }

    /// Write rows to a temp Parquet file, read it back as one batch, and also
    /// render the same rows to TSV — so the two can be compared cell by cell.
    fn write_and_read(
        rows: &[(&[(&str, &'static str)], &[(&str, &'static str)])],
    ) -> (RecordBatch, Vec<Vec<String>>) {
        let path = std::env::temp_dir().join(format!(
            "maligno_pq_test_{}_{:?}.parquet",
            std::process::id(),
            std::thread::current().id()
        ));
        let mut tsv_rows = Vec::new();
        {
            let mut w = ComparisonParquetWriter::new(std::fs::File::create(&path).unwrap())
                .expect("writer");
            for (a, b) in rows {
                let row = ComparisonRow::build(
                    "read1",
                    100,
                    "SetA",
                    "SetB",
                    |c| lookup(a, c),
                    |c| lookup(b, c),
                );
                let mut buf: Vec<u8> = Vec::new();
                row.write_tsv_row(&mut buf).unwrap();
                let line = String::from_utf8(buf).unwrap();
                tsv_rows.push(
                    line.trim_end_matches('\n').split('\t').map(String::from).collect(),
                );
                w.append(&row).expect("append");
            }
            assert_eq!(w.finish().unwrap(), rows.len() as u64);
        }
        let file = std::fs::File::open(&path).unwrap();
        let mut reader = ParquetRecordBatchReaderBuilder::try_new(file)
            .unwrap()
            .build()
            .unwrap();
        let batch = reader.next().expect("one batch").unwrap();
        let _ = std::fs::remove_file(&path);
        (batch, tsv_rows)
    }

    /// The documented inverse: render a Parquet cell back to its TSV form.
    fn cell_to_tsv(batch: &RecordBatch, col: usize, row: usize) -> String {
        arrow_cell_to_string(batch.column(col).as_ref(), row)
    }

    // ── the schema is derived, not restated ──────────────────────────────────

    #[test]
    fn every_column_has_a_declared_arrow_type() {
        for c in READINFO_DATA_COLS {
            assert!(per_side_type(c).is_some(), "per-side column {c:?} has no Arrow type");
        }
        for c in comparison_col_names() {
            assert!(comparison_type(c).is_some(), "comparison column {c:?} has no Arrow type");
        }
    }

    #[test]
    fn type_table_has_no_entries_for_columns_that_do_not_exist() {
        // Guards the other direction: a renamed column must not leave a stale entry
        // behind. Counted rather than enumerated, so it cannot drift out of date.
        let per_side_named = [
            "TargetChr", "Strand", "Target_Start", "Target_End", "Query_Start", "Query_End",
            "MQ_Best", "Num_Aln", "Num_Aln_MaxScore", "AS_Max", "ms_Max", "seqid_Max",
            "Query_Aln_Cov_Max", "Query_Aln_Len_Max", "JuncCount",
            "N_Splice_Junction_Events", "N_Splice_Junction_Bases", "N_Match_Events",
            "N_Match_Bases", "N_Substitution_Events", "N_Substitution_Bases",
            "N_Insertion_Events", "N_Insertion_Bases", "N_Deletion_Events",
            "N_Deletion_Bases", "N_SoftClipped_Events", "N_SoftClipped_Bases_Start",
            "N_SoftClipped_Bases_End", "junctions", "genomic_junctions", "cs",
        ];
        assert_eq!(per_side_named.len(), READINFO_DATA_COLS.len());
        for c in per_side_named {
            assert!(READINFO_DATA_COLS.contains(&c), "stale type-table entry: {c:?}");
        }
    }

    #[test]
    fn schema_field_names_match_the_tsv_header_exactly() {
        let mut buf: Vec<u8> = Vec::new();
        crate::comparison_row::write_compare_header(&mut buf).unwrap();
        let header = String::from_utf8(buf).unwrap();
        let tsv: Vec<&str> = header.trim_end_matches('\n').split('\t').collect();
        let schema = build_schema();
        let arrow: Vec<&str> = schema.fields().iter().map(|f| f.name().as_str()).collect();
        assert_eq!(arrow.len(), 96);
        assert_eq!(arrow, tsv, "Arrow schema and TSV header have drifted apart");
    }

    #[test]
    fn diff_values_agree_with_the_declared_comparison_types() {
        let (batch, _) = write_and_read(&[(&[], &[])]);
        let names = comparison_col_names();
        let base = 4 + 2 * READINFO_DATA_COLS.len();
        for (i, name) in names.iter().enumerate() {
            let want = comparison_type(name).unwrap().arrow();
            let got = batch.column(base + i).data_type();
            assert_eq!(got, &want, "column {name:?} type mismatch");
        }
    }

    // ── the round trip ───────────────────────────────────────────────────────

    #[test]
    fn parquet_round_trips_to_the_tsv_row() {
        // A mapped/mapped pair and an unmapped/mapped pair, so nulls, "*", "()" and
        // real zeros all appear.
        let mapped_a = [
            ("TargetChr", "chr1"), ("Strand", "+"), ("cs", ":100"),
            ("Query_Start", "0"), ("Query_End", "100"),
            ("Target_Start", "500"), ("Target_End", "600"),
            ("seqid_Max", "0.99"), ("Query_Aln_Cov_Max", "1.0"),
            ("AS_Max", "90"), ("ms_Max", "100"), ("junctions", "(10, 20)"),
            ("genomic_junctions", "((510, 520),)"), ("N_Insertion_Bases", "2"),
        ];
        let mapped_b = [
            ("TargetChr", "chr1"), ("Strand", "+"), ("cs", ":50*at:49"),
            ("Query_Start", "0"), ("Query_End", "100"),
            ("Target_Start", "500"), ("Target_End", "600"),
            ("seqid_Max", "0.98"), ("Query_Aln_Cov_Max", "1.0"),
            ("AS_Max", "88"), ("ms_Max", "99"), ("junctions", "(10, 30)"),
            ("genomic_junctions", "((510, 530),)"), ("N_Insertion_Bases", "4"),
        ];
        let unmapped = [
            ("TargetChr", "*"), ("Strand", "*"), ("cs", ""),
            ("Query_Start", "0"), ("Query_End", "0"),
            ("Target_Start", "0"), ("Target_End", "0"),
            ("seqid_Max", "NaN"), ("Query_Aln_Cov_Max", "0.0"),
            ("AS_Max", "0"), ("ms_Max", "0"), ("junctions", "()"),
            ("genomic_junctions", "()"),
        ];
        let (batch, tsv) = write_and_read(&[
            (&mapped_a, &mapped_b),
            (&unmapped, &mapped_b),
            (&unmapped, &unmapped),
        ]);
        let schema = build_schema();
        let n_side = READINFO_DATA_COLS.len();
        // Which columns are escaped on the way into the TSV: both per-side blocks
        // and the four trailing object lists. Not the keys or the labels.
        let escaped = |c: usize| (4..4 + 2 * n_side).contains(&c) || c >= 92;

        for (r, tsv_row) in tsv.iter().enumerate() {
            assert_eq!(tsv_row.len(), 96);
            for c in 0..96 {
                let recovered = cell_to_tsv(&batch, c, r);
                let expect = if escaped(c) { escape_tsv_field(&recovered) } else { recovered };
                assert_eq!(
                    expect, tsv_row[c],
                    "row {r}, column {:?} does not round-trip", schema.field(c).name()
                );
            }
        }
    }

    // ── the null rule ────────────────────────────────────────────────────────

    #[test]
    fn nulls_mean_undefined_and_nothing_else() {
        let unmapped = [
            ("TargetChr", "*"), ("Strand", "*"), ("cs", ""),
            ("seqid_Max", "NaN"), ("Query_Aln_Cov_Max", "0.0"),
            ("junctions", "()"), ("genomic_junctions", "()"), ("AS_Max", "0"),
        ];
        let (batch, _) = write_and_read(&[(&unmapped, &unmapped)]);
        let idx = |name: &str| {
            build_schema().fields().iter().position(|f| f.name() == name).unwrap()
        };
        let is_null = |name: &str| batch.column(idx(name)).is_null(0);

        // Undefined -> null.
        assert!(is_null("cs_A"), "an unmapped cs is undefined, so null");
        assert!(is_null("seqid_Max_A"), "NaN identity is undefined, so null");
        assert!(is_null("AS_Ratio"), "a zero denominator is undefined, so null");
        assert!(is_null("seqid_Diff"), "a diff against NaN is undefined, so null");

        // Meaningful values -> preserved, NOT nulled.
        assert!(!is_null("TargetChr_A"));
        assert_eq!(
            batch.column(idx("TargetChr_A")).as_any().downcast_ref::<StringArray>().unwrap().value(0),
            "*", "'*' is the mapping indicator and must survive"
        );
        assert!(!is_null("junctions_A"));
        assert_eq!(
            batch.column(idx("junctions_A")).as_any().downcast_ref::<StringArray>().unwrap().value(0),
            "()", "an empty junction set is a value, not a null"
        );
        assert!(!is_null("Query_Aln_Cov_Max_A"), "0.0 coverage is a real zero");
        assert_eq!(
            batch.column(idx("Query_Aln_Cov_Max_A")).as_any().downcast_ref::<Float64Array>().unwrap().value(0),
            0.0
        );
        assert!(!is_null("AS_Max_A"), "a zero score is a real zero, not a null");
    }

    #[test]
    fn a_missing_numeric_column_becomes_null_not_a_fabricated_zero() {
        // Both accessors return "" for an absent column. Writing 0 there would
        // invent data; null says "we do not know".
        let (batch, _) = write_and_read(&[(&[], &[])]);
        let idx = |name: &str| {
            build_schema().fields().iter().position(|f| f.name() == name).unwrap()
        };
        for c in ["Target_Start_A", "MQ_Best_A", "N_Match_Events_B", "AS_Max_A"] {
            assert!(batch.column(idx(c)).is_null(0), "{c} should be null, not 0");
        }
    }

    // ── the escaping asymmetry (no dataset exercises this) ───────────────────

    #[test]
    fn parquet_stores_the_unescaped_value_while_the_tsv_escapes_it_again() {
        // Per-side values arrive already escaped once (readinfo escaped them on the
        // way out) and the TSV writer escapes again. Parquet stores what it was
        // handed, so the TSV cell is one further escaping of the Parquet cell.
        let a = [("cs", "a\\b"), ("TargetChr", "chr1")];
        let (batch, tsv) = write_and_read(&[(&a, &[])]);
        let idx = build_schema().fields().iter().position(|f| f.name() == "cs_A").unwrap();
        let pq = batch.column(idx).as_any().downcast_ref::<StringArray>().unwrap().value(0);
        assert_eq!(pq, "a\\b", "Parquet holds the value as received");
        assert_eq!(tsv[0][idx], "a\\\\b", "the TSV holds it escaped once more");
        assert_eq!(escape_tsv_field(pq), tsv[0][idx], "tsv == escape(parquet)");
    }

    // ── metadata ─────────────────────────────────────────────────────────────

    #[test]
    fn footer_carries_the_schema_version() {
        let path = std::env::temp_dir()
            .join(format!("maligno_pq_meta_{}.parquet", std::process::id()));
        {
            let w = ComparisonParquetWriter::new(std::fs::File::create(&path).unwrap()).unwrap();
            w.finish().unwrap();
        }
        let file = std::fs::File::open(&path).unwrap();
        let b = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
        let kv = b.metadata().file_metadata().key_value_metadata().cloned().unwrap_or_default();
        let found = kv.iter().find(|k| k.key == "maligno_comparison_schema_version");
        assert_eq!(found.and_then(|k| k.value.clone()).as_deref(), Some(SCHEMA_VERSION));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn output_format_selects_the_right_writers() {
        assert!(OutputFormat::Tsv.writes_tsv());
        assert!(!OutputFormat::Tsv.writes_parquet());
        assert!(!OutputFormat::Parquet.writes_tsv());
        assert!(OutputFormat::Parquet.writes_parquet());
        assert!(OutputFormat::Both.writes_tsv());
        assert!(OutputFormat::Both.writes_parquet());
        // `both` is the default while find-query-diff still reads the TSV; see the
        // enum's doc comment.
        assert_eq!(OutputFormat::default(), OutputFormat::Both);
    }

    #[test]
    fn extension_selects_the_format() {
        assert!(is_parquet_path("x.parquet"));
        assert!(is_parquet_path("/a/b/AvsB.compare.parquet"));
        assert!(!is_parquet_path("x.tsv.gz"));
        assert!(!is_parquet_path("x.tsv"));
        assert!(!is_parquet_path("-"), "stdout keeps writing TSV");
        assert!(!is_parquet_path("x.parquet.gz"), "Parquet compresses internally");
    }

    // ── ParquetRowReader ─────────────────────────────────────────────────────

    /// Write one row to a temp Parquet file (reusing the writer above) and
    /// return its path. Caller is responsible for removing the file.
    fn write_temp_parquet(tag: &str) -> String {
        let path = std::env::temp_dir()
            .join(format!("maligno_pq_row_reader_{tag}_{}.parquet", std::process::id()));
        let a = [("TargetChr", "chr1"), ("Strand", "+"), ("cs", ":100"),
            ("Query_Start", "0"), ("Query_End", "100"),
            ("Target_Start", "500"), ("Target_End", "600")];
        let b = [("TargetChr", "chr1"), ("Strand", "+"), ("cs", ":50*at:49"),
            ("Query_Start", "0"), ("Query_End", "100"),
            ("Target_Start", "500"), ("Target_End", "600")];
        let mut w = ComparisonParquetWriter::new(std::fs::File::create(&path).unwrap())
            .expect("writer");
        let row = ComparisonRow::build("read1", 100, "SetA", "SetB", |c| lookup(&a, c), |c| lookup(&b, c));
        w.append(&row).expect("append");
        assert_eq!(w.finish().unwrap(), 1);
        path.to_string_lossy().into_owned()
    }

    #[test]
    fn unprojected_reader_yields_every_column_in_schema_order() {
        let path = write_temp_parquet("full");
        let mut r = ParquetRowReader::open(&path, None).unwrap();
        assert_eq!(r.columns().len(), 96, "no projection: every column present");
        assert_eq!(r.columns(), build_schema().fields().iter().map(|f| f.name().clone()).collect::<Vec<_>>().as_slice());
        let row = r.next_row().unwrap().expect("one row");
        assert_eq!(row.len(), 96);
        assert!(r.next_row().unwrap().is_none(), "only one row was written");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn projected_reader_yields_only_the_requested_columns_in_schema_order() {
        let path = write_temp_parquet("proj");
        // Deliberately out of schema order and with a duplicate-ish adjacent
        // pair, to confirm the reader restores file/schema order regardless.
        let wanted = ["Target_End_B", "TargetChr_A", "Strand_A", "Read_Name"];
        let mut r = ParquetRowReader::open(&path, Some(&wanted)).unwrap();
        // Schema order (Read_Name is column 0; per-side block is A-then-B,
        // READINFO_DATA_COLS order): Read_Name, TargetChr_A, Strand_A, ..., Target_End_B.
        assert_eq!(r.columns(), &["Read_Name", "TargetChr_A", "Strand_A", "Target_End_B"]);
        let row = r.next_row().unwrap().expect("one row");
        assert_eq!(row, vec!["read1", "chr1", "+", "600"]);
        assert!(r.next_row().unwrap().is_none());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn projecting_a_nonexistent_column_is_a_clear_error() {
        let path = write_temp_parquet("missing");
        let wanted = ["TargetChr_A", "Not_A_Real_Column"];
        let err = match ParquetRowReader::open(&path, Some(&wanted)) {
            Ok(_) => panic!("expected an error for a nonexistent column"),
            Err(e) => e.to_string(),
        };
        assert!(err.contains("Not_A_Real_Column"), "unexpected error: {err}");
        let _ = std::fs::remove_file(&path);
    }
}
