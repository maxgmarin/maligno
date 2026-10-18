//! The per-read comparison row: one typed value per output line.
//!
//! This module owns the comparison table's **rows**: the row type, how each
//! A-vs-B value is computed, and the TSV writers. The column names, order and
//! types live in `columns.rs` (`COMPARE_KEY_COLS`, `COMPARE_SIDE_COLS`,
//! `COMPARE_DIFF_COLS`).
//!
//! Splitting the row into typed parts lets the diff logic be tested on its own,
//! and means the per-side parsing is written once and used for both sides.
//!
//! The split here is deliberate, and it is about what each layer is *for*:
//!
//! - [`AlignmentMetrics`] — the 10 values per side that the comparison math needs,
//!   parsed once. No lifetime, `Copy`; this is what [`AlignmentDiff::compute`] is
//!   tested against.
//! - [`AlignmentRow`] — one side's 35 columns exactly as read (`&str`), plus the
//!   parsed view. The raw strings are what the TSV writer emits.
//! - [`AlignmentDiff`] — the computed A-vs-B block, owned.
//!
//! **The per-side columns are transported as `&str`, not re-rendered from parsed
//! values.** That is load-bearing for two reasons. First, the writer loops over
//! `COMPARE_SIDE_COLS` exactly as the header does, so the two cannot drift apart.
//! Second, `toolkit compare-readinfo` accepts user-written TSVs: re-rendering would turn
//! `1e-05` into `0.00001`, `+5` into `5`, and would fabricate `0` where an empty
//! cell passes through. Most of these columns feed no computation at all;
//! parsing them would be a lossy operation performed for looks.

use std::io::Write;

use crate::columns::{
    compare_header, index_of, AS_SCORE, COMPARE_SIDE_COLS, QUERY_JUNCTIONS, MS_SCORE,
    N_DELETION_BASES, N_INSERTION_BASES, N_SOFTCLIP_BASES_END, N_SOFTCLIP_BASES_START,
    N_SUBSTITUTION_BASES, QUERY_ALN_COV, QUERY_ALN_LEN, ALN_IDENTITY, STRAND,
};
use crate::io_utils::{escape_tsv_field, fmt_float};
use crate::junction::{
    format_junction_tuple, junction_set_diffs, junction_set_stats,
    parse_junction_str,
};

// ── Column indices ───────────────────────────────────────────────────────────

/// Number of per-side data columns. `AlignmentRow::raw` is sized from this, so the
/// array and the header loop are provably the same length.
const N_SIDE_COLS: usize = COMPARE_SIDE_COLS.len();

const I_STRAND: usize = index_of(&COMPARE_SIDE_COLS, STRAND);
const I_JUNCTIONS: usize = index_of(&COMPARE_SIDE_COLS, QUERY_JUNCTIONS);

// ── Ratio helpers ────────────────────────────────────────────────────────────

/// `num / den`, or `NaN` when the denominator is zero.
///
/// The zero test is on the **integer**, before the cast. That matters: `5.0 / 0.0` is `inf` in IEEE
/// arithmetic, so a bare float division would render `inf` for every read whose
/// denominator side scored 0, where this renders `NaN`.
fn ratio_i64(num: i64, den: i64) -> f64 {
    if den == 0 {
        f64::NAN
    } else {
        num as f64 / den as f64
    }
}

// ── Per-side types ───────────────────────────────────────────────────────────

/// The values a side contributes to the comparison math, parsed once.
///
/// Unparseable values default to `0` for the integer fields and `f64::NAN` for
/// the two float fields.
#[derive(Debug, Clone, Copy)]
pub(crate) struct AlignmentMetrics {
    pub(crate) aln_score: i64,
    pub(crate) ms: i64,
    pub(crate) seqid: f64,
    pub(crate) query_aln_cov: f64,
    pub(crate) query_aln_len: u64,
    pub(crate) n_ins_bases: u64,
    pub(crate) n_del_bases: u64,
    pub(crate) n_sub_bases: u64,
    pub(crate) n_sc_start: u64,
    pub(crate) n_sc_end: u64,
}

/// One side of the comparison: the raw transport columns plus the parsed view.
pub(crate) struct AlignmentRow<'r> {
    /// The 35 per-side columns in `COMPARE_SIDE_COLS` order, exactly as read.
    raw: [&'r str; N_SIDE_COLS],
    pub(crate) metrics: AlignmentMetrics,
}

impl<'r> AlignmentRow<'r> {
    /// Read one side out of a readinfo row via a by-name accessor.
    ///
    /// The explicit `'r` is **required**: `impl Fn(&str) -> &str` in argument
    /// position elides to `for<'x> Fn(&'x str) -> &'x str`, which would tie the
    /// returned string to the *key*. Neither call site satisfies that — both return
    /// data borrowed from the row or reader, not from the column name.
    pub(crate) fn from_accessor(get: impl Fn(&str) -> &'r str) -> Self {
        let raw = std::array::from_fn(|i| get(COMPARE_SIDE_COLS[i].name));
        let metrics = AlignmentMetrics {
            aln_score: get(AS_SCORE.name).parse().unwrap_or(0),
            ms: get(MS_SCORE.name).parse().unwrap_or(0),
            seqid: get(ALN_IDENTITY.name).parse().unwrap_or(f64::NAN),
            query_aln_cov: get(QUERY_ALN_COV.name).parse().unwrap_or(f64::NAN),
            query_aln_len: get(QUERY_ALN_LEN.name).parse().unwrap_or(0),
            n_ins_bases: get(N_INSERTION_BASES.name).parse().unwrap_or(0),
            n_del_bases: get(N_DELETION_BASES.name).parse().unwrap_or(0),
            n_sub_bases: get(N_SUBSTITUTION_BASES.name).parse().unwrap_or(0),
            n_sc_start: get(N_SOFTCLIP_BASES_START.name).parse().unwrap_or(0),
            n_sc_end: get(N_SOFTCLIP_BASES_END.name).parse().unwrap_or(0),
        };
        Self { raw, metrics }
    }

    /// The 35 per-side columns in `COMPARE_SIDE_COLS` order, exactly as read.
    /// The Parquet writer walks this in step with the column list.
    pub(crate) fn raw(&self) -> &[&'r str; N_SIDE_COLS] {
        &self.raw
    }

    fn strand(&self) -> &'r str {
        self.raw[I_STRAND]
    }
    fn junctions(&self) -> &'r str {
        self.raw[I_JUNCTIONS]
    }
}

// ── The computed comparison block ────────────────────────────────────────────

/// The A-vs-B metrics. Field order matches `COMPARE_DIFF_COLS`.
///
/// Sign conventions: every
/// `*_diff` is **B − A**, every ratio is **B / A** (so the denominator is the A
/// side, and A is what the zero-guard tests).
#[derive(Debug, Clone)]
pub(crate) struct AlignmentDiff {
    pub(crate) strand_match: bool,
    pub(crate) seqid_diff: f64,
    pub(crate) query_aln_cov_diff: f64,
    pub(crate) query_aln_len_diff: i64,
    pub(crate) as_diff: i64,
    pub(crate) ms_diff: i64,
    pub(crate) as_ratio: f64,
    pub(crate) ms_ratio: f64,
    pub(crate) n_sub_bases_diff: i64,
    pub(crate) n_ins_bases_diff: i64,
    pub(crate) n_del_bases_diff: i64,
    pub(crate) n_sc_start_diff: i64,
    pub(crate) n_sc_end_diff: i64,
    pub(crate) n_matched_junctions: u64,
    pub(crate) n_unmatched_junctions: u64,
    pub(crate) n_junctions_only_a: u64,
    pub(crate) n_junctions_only_b: u64,
    pub(crate) junctions_only_a: String,
    pub(crate) junctions_only_b: String,
}

impl AlignmentDiff {
    pub(crate) fn compute(a: &AlignmentRow, b: &AlignmentRow) -> Self {
        let (ma, mb) = (&a.metrics, &b.metrics);

        // Query-space junction sets.
        let juncs_a = parse_junction_str(a.junctions());
        let juncs_b = parse_junction_str(b.junctions());
        let (n_matched, n_only_a, n_only_b) = junction_set_stats(&juncs_a, &juncs_b);
        let (j_only_a_vec, j_only_b_vec) = junction_set_diffs(&juncs_a, &juncs_b);

        Self {
            strand_match: a.strand() == b.strand(),
            seqid_diff: mb.seqid - ma.seqid,
            query_aln_cov_diff: mb.query_aln_cov - ma.query_aln_cov,
            query_aln_len_diff: mb.query_aln_len as i64 - ma.query_aln_len as i64,
            as_diff: mb.aln_score - ma.aln_score,
            ms_diff: mb.ms - ma.ms,
            as_ratio: ratio_i64(mb.aln_score, ma.aln_score),
            ms_ratio: ratio_i64(mb.ms, ma.ms),
            n_sub_bases_diff: mb.n_sub_bases as i64 - ma.n_sub_bases as i64,
            n_ins_bases_diff: mb.n_ins_bases as i64 - ma.n_ins_bases as i64,
            n_del_bases_diff: mb.n_del_bases as i64 - ma.n_del_bases as i64,
            n_sc_start_diff: mb.n_sc_start as i64 - ma.n_sc_start as i64,
            n_sc_end_diff: mb.n_sc_end as i64 - ma.n_sc_end as i64,
            n_matched_junctions: n_matched,
            n_unmatched_junctions: n_only_a + n_only_b,
            n_junctions_only_a: n_only_a,
            n_junctions_only_b: n_only_b,
            junctions_only_a: format_junction_tuple(&j_only_a_vec),
            junctions_only_b: format_junction_tuple(&j_only_b_vec),
        }
    }
}

/// One comparison-block value, tagged with its type.
///
/// [`AlignmentDiff::values`] yields these in `COMPARE_DIFF_COLS` order so a
/// second writer (Parquet) can consume the block without restating the order. The
/// TSV writer keeps its own explicit `write!` sequence — the two are cross-checked
/// by `parquet_matches_tsv_for_the_comparison_block` in `parquet_out`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum DiffValue<'a> {
    Bool(bool),
    I64(i64),
    U64(u64),
    /// `NaN` here means "undefined"; the Parquet writer stores it as a null.
    F64(f64),
    Str(&'a str),
}

impl AlignmentDiff {
    /// The 19 comparison values in `COMPARE_DIFF_COLS` order.
    pub(crate) fn values(&self) -> Vec<DiffValue<'_>> {
        use DiffValue::*;
        vec![
            // Orientation & identity.
            Bool(self.strand_match),
            F64(self.seqid_diff),
            F64(self.query_aln_cov_diff),
            I64(self.query_aln_len_diff),
            // Score.
            I64(self.as_diff),
            I64(self.ms_diff),
            F64(self.as_ratio),
            F64(self.ms_ratio),
            // Event-count diffs.
            I64(self.n_sub_bases_diff),
            I64(self.n_ins_bases_diff),
            I64(self.n_del_bases_diff),
            I64(self.n_sc_start_diff),
            I64(self.n_sc_end_diff),
            // Query-space junction set comparison.
            U64(self.n_matched_junctions),
            U64(self.n_unmatched_junctions),
            U64(self.n_junctions_only_a),
            U64(self.n_junctions_only_b),
            // Object lists last.
            Str(&self.junctions_only_a),
            Str(&self.junctions_only_b),
        ]
    }
}

// ── The row, and its TSV serialization ───────────────────────────────────────

/// One output line: the join keys, the two set labels, both sides, and the diff.
pub(crate) struct ComparisonRow<'r> {
    pub(crate) read_name: &'r str,
    pub(crate) read_len: u64,
    pub(crate) label_a: &'r str,
    pub(crate) label_b: &'r str,
    pub(crate) aln_a: AlignmentRow<'r>,
    pub(crate) aln_b: AlignmentRow<'r>,
    pub(crate) diff: AlignmentDiff,
}

impl<'r> ComparisonRow<'r> {
    /// Build a row from one by-name accessor per side.
    pub(crate) fn build(
        read_name: &'r str,
        read_len: u64,
        label_a: &'r str,
        label_b: &'r str,
        get_a: impl Fn(&str) -> &'r str,
        get_b: impl Fn(&str) -> &'r str,
    ) -> Self {
        let aln_a = AlignmentRow::from_accessor(get_a);
        let aln_b = AlignmentRow::from_accessor(get_b);
        let diff = AlignmentDiff::compute(&aln_a, &aln_b);
        Self { read_name, read_len, label_a, label_b, aln_a, aln_b, diff }
    }

    /// Write the row. Field order here is checked against the header by
    /// `header_and_row_field_counts_agree`.
    ///
    /// Escaping: the 70 per-side fields and
    /// the 2 object-list fields are escaped; `query_name` and the two labels are
    /// **not** (labels are validated by `validate_set_label` instead). Note the
    /// per-side values were already escaped once when the readinfo row was written,
    /// so this is a second pass over them — `escape_tsv_field` is not idempotent,
    /// and that double escaping is intentional; changing it would change the TSV output.
    pub(crate) fn write_tsv_row<W: Write>(&self, out: &mut W) -> std::io::Result<()> {
        let d = &self.diff;
        write!(
            out,
            "{}\t{}\t{}\t{}",
            self.read_name, self.read_len, self.label_a, self.label_b
        )?;
        for f in &self.aln_a.raw {
            write!(out, "\t{}", escape_tsv_field(f))?;
        }
        for f in &self.aln_b.raw {
            write!(out, "\t{}", escape_tsv_field(f))?;
        }
        // Orientation & identity, then score.
        write!(
            out,
            "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            d.strand_match,
            fmt_float(d.seqid_diff),
            fmt_float(d.query_aln_cov_diff),
            d.query_aln_len_diff,
            d.as_diff,
            d.ms_diff,
            fmt_float(d.as_ratio),
            fmt_float(d.ms_ratio),
        )?;
        // Event-count diffs.
        write!(
            out,
            "\t{}\t{}\t{}\t{}\t{}",
            d.n_sub_bases_diff,
            d.n_ins_bases_diff,
            d.n_del_bases_diff,
            d.n_sc_start_diff,
            d.n_sc_end_diff,
        )?;
        // Query-space junction set comparison.
        write!(
            out,
            "\t{}\t{}\t{}\t{}",
            d.n_matched_junctions,
            d.n_unmatched_junctions,
            d.n_junctions_only_a,
            d.n_junctions_only_b,
        )?;
        // Object lists last.
        write!(
            out,
            "\t{}\t{}",
            escape_tsv_field(&d.junctions_only_a),
            escape_tsv_field(&d.junctions_only_b),
        )?;
        writeln!(out)
    }
}

/// Write the comparison header: the join keys, the two set-label columns, the
/// per-side data columns (suffixed `_A` / `_B`), then the comparison columns.
///
/// Side suffixes are **fixed** (`_A` / `_B`), never the user's label — the
/// human-readable labels are carried as the `label_A` / `label_B` data columns
/// instead, so column names are stable across datasets and unambiguous even when a
/// label itself contains an underscore.
pub(crate) fn write_compare_header<W: Write>(out: &mut W) -> std::io::Result<()> {
    writeln!(out, "{}", compare_header().join("\t"))
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::columns::COMPARE_DIFF_COLS;

    /// Look a column up in a `(name, value)` slice, mimicking the accessors used at
    /// both call sites: an absent column yields `""`.
    fn lookup(pairs: &[(&str, &'static str)], col: &str) -> &'static str {
        pairs
            .iter()
            .find(|(k, _)| *k == col)
            .map(|(_, v)| *v)
            .unwrap_or("")
    }

    fn side(pairs: &[(&str, &'static str)]) -> AlignmentRow<'static> {
        AlignmentRow::from_accessor(|c| lookup(pairs, c))
    }

    fn diff_of(
        a: &[(&str, &'static str)],
        b: &[(&str, &'static str)],
    ) -> AlignmentDiff {
        AlignmentDiff::compute(&side(a), &side(b))
    }

    /// Render one full row and split it into its 93 fields.
    fn render(a: &[(&str, &'static str)], b: &[(&str, &'static str)]) -> Vec<String> {
        let mut buf: Vec<u8> = Vec::new();
        ComparisonRow::build("read1", 100, "SetA", "SetB", |c| lookup(a, c), |c| lookup(b, c))
            .write_tsv_row(&mut buf)
            .expect("write to Vec cannot fail");
        let s = String::from_utf8(buf).expect("output is utf-8");
        assert!(s.ends_with('\n'), "row must end with a newline");
        assert!(!s.trim_end().ends_with('\t'), "row must not end with a tab");
        s.trim_end_matches('\n').split('\t').map(String::from).collect()
    }

    fn header_fields() -> Vec<String> {
        let mut buf: Vec<u8> = Vec::new();
        write_compare_header(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        s.trim_end_matches('\n').split('\t').map(String::from).collect()
    }

    /// Index of an output column in the 93-column header.
    fn col(name: &str) -> usize {
        header_fields()
            .iter()
            .position(|c| c == name)
            .unwrap_or_else(|| panic!("no output column named {name}"))
    }

    // ── schema shape ─────────────────────────────────────────────────────────

    #[test]
    fn header_and_row_field_counts_agree() {
        let h = header_fields();
        let r = render(&[], &[]);
        assert_eq!(h.len(), 93, "header should have 93 columns");
        assert_eq!(
            r.len(),
            h.len(),
            "row emits {} fields but the header names {} — the two lists have drifted",
            r.len(),
            h.len()
        );
        assert_eq!(
            4 + 2 * COMPARE_SIDE_COLS.len() + COMPARE_DIFF_COLS.len(),
            93
        );
    }

    #[test]
    fn header_layout_is_keys_then_a_then_b_then_comparison() {
        let h = header_fields();
        assert_eq!(&h[..4], &["query_name", "query_len", "label_A", "label_B"]);
        assert_eq!(h[4], "target_name_A");
        assert_eq!(h[38], "genomic_junctions_A");
        assert_eq!(h[39], "target_name_B");
        assert_eq!(h[73], "genomic_junctions_B");
        assert_eq!(h[74], "strand_match");
        assert_eq!(h[92], "query_junctions_only_B");
        for side in ["A", "B"] {
            // The representative's tp:A, then the alignment-type counts.
            let m = col(&format!("mapq_{side}"));
            assert_eq!(h[m + 1], format!("tp_tag_{side}"));
            let n = col(&format!("n_alns_{side}"));
            assert_eq!(h[n + 1], format!("n_alns_primsupp_{side}"));
            assert_eq!(h[n + 2], format!("n_alns_secondary_{side}"));
            // cs follows query_aln_len, then the two strands.
            let l = col(&format!("query_aln_len_{side}"));
            assert_eq!(h[l + 1], format!("cs_{side}"));
            assert_eq!(h[l + 2], format!("query_transcript_strand_{side}"));
            assert_eq!(h[l + 3], format!("genomic_intron_strand_{side}"));
        }
    }

    // ── escaping (no test dataset exercises this: real data has no \, tab, CR or LF) ──

    #[test]
    fn per_side_values_are_escaped_and_escaping_is_applied_again() {
        // Per-side values reach us already escaped once (readinfo escaped them on the
        // way out), and this writer escapes again. A single backslash must therefore
        // come out as two, not one.
        let a = [("cs", "a\\b")];
        let r = render(&a, &[]);
        assert_eq!(r[col("cs_A")], "a\\\\b");
    }

    #[test]
    fn tabs_in_a_per_side_value_are_escaped_so_the_row_keeps_its_shape() {
        let a = [("cs", "a\tb")];
        let r = render(&a, &[]);
        assert_eq!(r[col("cs_A")], "a\\tb");
        assert_eq!(r.len(), 93, "an embedded tab must not add a field");
    }

    #[test]
    fn keys_and_labels_are_not_escaped() {
        // Labels are validated by `validate_set_label` rather than escaped, and
        // query_name is passed through as-is. Backslashes must survive untouched.
        let mut buf: Vec<u8> = Vec::new();
        ComparisonRow::build("read\\1", 100, "Set\\A", "Set\\B", |_| "", |_| "")
            .write_tsv_row(&mut buf)
            .unwrap();
        let s = String::from_utf8(buf).unwrap();
        let f: Vec<&str> = s.trim_end_matches('\n').split('\t').collect();
        assert_eq!(f[0], "read\\1");
        assert_eq!(f[2], "Set\\A");
        assert_eq!(f[3], "Set\\B");
    }

    // ── the zero-denominator guard ───────────────────────────────────────────

    #[test]
    fn zero_denominator_yields_nan_not_infinity() {
        // The guard tests the integer before the cast. A bare float division would
        // give inf here, changing as_score_ratio for every read whose A side scored 0.
        assert!(ratio_i64(5, 0).is_nan(), "5/0 must be NaN, not inf");
        assert!(ratio_i64(-5, 0).is_nan(), "-5/0 must be NaN, not -inf");
        assert!(ratio_i64(0, 0).is_nan());
        assert_eq!(fmt_float(ratio_i64(5, 0)), "NaN");
        assert_eq!(fmt_float(ratio_i64(-5, 0)), "NaN");
    }

    #[test]
    fn both_ratio_columns_render_nan_when_the_a_side_is_zero() {
        let a = [("as_score", "0"), ("ms_score", "0")];
        let b = [("as_score", "50"), ("ms_score", "60")];
        let r = render(&a, &b);
        for c in ["as_score_ratio", "ms_score_ratio"] {
            assert_eq!(r[col(c)], "NaN", "{c} should be NaN when the A side is 0");
        }
    }

    #[test]
    fn ratios_are_b_over_a() {
        let a = [("as_score", "10"), ("ms_score", "10")];
        let b = [("as_score", "20"), ("ms_score", "5")];
        let d = diff_of(&a, &b);
        assert_eq!(d.as_ratio, 2.0, "as_score_ratio is B/A");
        assert_eq!(d.ms_ratio, 0.5, "ms_score_ratio is B/A");
    }

    // ── sign conventions ─────────────────────────────────────────────────────

    #[test]
    fn every_diff_is_b_minus_a() {
        let a = [
            ("as_score", "10"),
            ("ms_score", "10"),
            ("query_aln_len", "100"),
            ("n_substitution_bases", "1"),
            ("n_insertion_bases", "2"),
            ("n_deletion_bases", "3"),
            ("n_softclip_bases_start", "4"),
            ("n_softclip_bases_end", "5"),
            ("aln_identity", "0.9"),
            ("query_aln_cov", "0.8"),
        ];
        let b = [
            ("as_score", "13"),
            ("ms_score", "14"),
            ("query_aln_len", "110"),
            ("n_substitution_bases", "2"),
            ("n_insertion_bases", "4"),
            ("n_deletion_bases", "6"),
            ("n_softclip_bases_start", "9"),
            ("n_softclip_bases_end", "1"),
            ("aln_identity", "1.0"),
            ("query_aln_cov", "0.5"),
        ];
        let d = diff_of(&a, &b);
        assert_eq!(d.as_diff, 3);
        assert_eq!(d.ms_diff, 4);
        assert_eq!(d.query_aln_len_diff, 10);
        assert_eq!(d.n_sub_bases_diff, 1);
        assert_eq!(d.n_ins_bases_diff, 2);
        assert_eq!(d.n_del_bases_diff, 3);
        assert_eq!(d.n_sc_start_diff, 5);
        assert_eq!(d.n_sc_end_diff, -4, "diffs may be negative");
        assert!((d.seqid_diff - 0.1).abs() < 1e-12);
        assert!((d.query_aln_cov_diff + 0.3).abs() < 1e-12);
    }

    // ── strand ───────────────────────────────────────────────────────────────

    #[test]
    fn strand_match_compares_the_raw_strand_strings() {
        assert!(diff_of(&[("strand", "+")], &[("strand", "+")]).strand_match);
        assert!(!diff_of(&[("strand", "+")], &[("strand", "-")]).strand_match);
        // Two unmapped sides agree with each other.
        assert!(diff_of(&[("strand", "*")], &[("strand", "*")]).strand_match);
        assert!(!diff_of(&[("strand", "*")], &[("strand", "+")]).strand_match);
    }

    #[test]
    fn strand_match_renders_lowercase() {
        let r = render(&[("strand", "+")], &[("strand", "+")]);
        assert_eq!(r[col("strand_match")], "true");
        let r = render(&[("strand", "+")], &[("strand", "-")]);
        assert_eq!(r[col("strand_match")], "false");
    }

    // ── parse defaults and the unmapped case ─────────────────────────────────

    #[test]
    fn missing_columns_parse_to_the_documented_defaults() {
        let m = side(&[]).metrics;
        assert_eq!(m.aln_score, 0);
        assert_eq!(m.ms, 0);
        assert_eq!(m.query_aln_len, 0);
        assert_eq!(m.n_ins_bases, 0);
        assert_eq!(m.n_del_bases, 0);
        assert_eq!(m.n_sub_bases, 0);
        assert_eq!(m.n_sc_start, 0);
        assert_eq!(m.n_sc_end, 0);
        assert!(m.seqid.is_nan(), "aln_identity defaults to NaN, not 0");
        assert!(m.query_aln_cov.is_nan());
    }

    #[test]
    fn a_missing_column_stays_an_empty_cell_and_is_not_rendered_as_zero() {
        // The per-side columns are transported as &str precisely so this holds. If
        // they were re-rendered from parsed values, these would become "0".
        let r = render(&[], &[]);
        for c in ["n_match_events_A", "target_start_A", "cs_B", "mapq_B"] {
            assert_eq!(r[col(c)], "", "{c} should stay empty, not become 0");
        }
    }

    #[test]
    fn unmapped_side_propagates_nan_through_the_float_diffs() {
        // What an unmapped read looks like coming out of readinfo: identity/coverage are NaN.
        let unmapped = [("target_name", "*"), ("strand", "*"), ("aln_identity", "NaN"),
                        ("query_aln_cov", "NaN"), ("cs", ""), ("query_junctions", "()")];
        let mapped = [("target_name", "chr1"), ("strand", "+"), ("aln_identity", "0.99"),
                      ("query_aln_cov", "1.0"), ("cs", ":100"), ("query_junctions", "()")];
        let r = render(&unmapped, &mapped);
        assert_eq!(r[col("aln_identity_diff")], "NaN");
        assert_eq!(r[col("query_aln_cov_diff")], "NaN");
        assert_eq!(r[col("target_name_A")], "*");
        assert_eq!(r[col("cs_A")], "", "an unmapped cs is an empty cell");
    }

    // ── junction set comparison wiring ───────────────────────────────────────

    #[test]
    fn junction_counts_and_object_lists_line_up() {
        let a = [("query_junctions", "(10, 20, 30)")];
        let b = [("query_junctions", "(20, 30, 40)")];
        let d = diff_of(&a, &b);
        assert_eq!(d.n_matched_junctions, 2);
        assert_eq!(d.n_junctions_only_a, 1);
        assert_eq!(d.n_junctions_only_b, 1);
        assert_eq!(
            d.n_unmatched_junctions,
            d.n_junctions_only_a + d.n_junctions_only_b
        );
        assert_eq!(d.junctions_only_a, "(10,)");
        assert_eq!(d.junctions_only_b, "(40,)");
    }
}
