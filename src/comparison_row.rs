//! The per-read comparison row: one typed value per output line.
//!
//! This module owns the **entire comparison-table schema** — the per-side column
//! list, the comparison column list, the row type, and the two TSV writers. It is
//! the single place to look when asking "what columns does `compare` emit, in what
//! order, and how is each one computed?".
//!
//! Why it exists: until v0.14.1 the row was assembled and serialized in one
//! ~145-line function that pulled values out of a `HashMap<&str, &str>`, parsed
//! them inline, computed 30 metrics as flat locals, and wrote straight to the
//! output. The A side, the B side and the diff block never existed as values, so
//! the diff logic could not be tested without running the whole pipeline, and the
//! per-side parse block had to be written twice — once per side, mirror-imaged.
//!
//! The split here is deliberate, and it is about what each layer is *for*:
//!
//! - [`AlignmentMetrics`] — the 11 values per side that the comparison math needs,
//!   parsed once. No lifetime, `Copy`; this is what [`AlignmentDiff::compute`] is
//!   tested against.
//! - [`AlignmentRow`] — one side's 31 columns exactly as read (`&str`), plus the
//!   parsed view. The raw strings are what the TSV writer emits.
//! - [`AlignmentDiff`] — the computed A-vs-B block, owned.
//!
//! **The per-side columns are transported as `&str`, not re-rendered from parsed
//! values.** That is load-bearing for two reasons. First, the writer loops over
//! `READINFO_DATA_COLS` exactly as the header does, so the two cannot drift apart.
//! Second, `compare-readinfo` accepts arbitrary user TSVs: re-rendering would turn
//! `1e-05` into `0.00001`, `+5` into `5`, and — because a missing column yields
//! `""` from the accessors — would fabricate `0` where today an empty cell passes
//! through. 14 of these columns feed no computation at all; parsing them would be a
//! lossy operation performed for looks.

use std::io::Write;

use crate::io_utils::{escape_tsv_field, fmt_float};
use crate::junction::{
    format_genomic_junction_tuple, format_junction_tuple, genomic_junction_set_diffs,
    genomic_junction_set_stats, junction_distance, junction_set_diffs, junction_set_stats,
    parse_genomic_junction_str, parse_junction_str,
};

// ── Column schema ────────────────────────────────────────────────────────────

/// The per-side data columns, each emitted once suffixed `_A` and once `_B`.
///
/// These are read out of the readinfo table **by name**, so this order is free to
/// differ from `readinfo.rs`'s `READINFO_HEADER`. Before v0.14.0 it was that header
/// verbatim; it is now grouped by topic — locus, alignment selection,
/// identity/coverage, junction counts, cs-derived event counts, and the three long
/// strings last — so the 96-column table is readable via
/// `head -1 | tr '\t' '\n' | nl`. **This divergence from `READINFO_HEADER` is
/// deliberate; do not "resync" the two.**
pub(crate) const READINFO_DATA_COLS: &[&str] = &[
    // Locus & span.
    "TargetChr",
    "Strand",
    "Target_Start",
    "Target_End",
    "Query_Start",
    "Query_End",
    // Alignment selection & score.
    "MQ_Best",
    "Num_Aln",
    "Num_Aln_MaxScore",
    "AS_Max",
    "ms_Max",
    // Identity & coverage.
    "seqid_Max",
    "Query_Aln_Cov_Max",
    "Query_Aln_Len_Max",
    // Junction counts.
    "JuncCount",
    "N_Splice_Junction_Events",
    "N_Splice_Junction_Bases",
    // cs-derived event counts.
    "N_Match_Events",
    "N_Match_Bases",
    "N_Substitution_Events",
    "N_Substitution_Bases",
    "N_Insertion_Events",
    "N_Insertion_Bases",
    "N_Deletion_Events",
    "N_Deletion_Bases",
    "N_SoftClipped_Events",
    "N_SoftClipped_Bases_Start",
    "N_SoftClipped_Bases_End",
    // Long strings last — these dominate the table's on-disk size.
    "junctions",
    "genomic_junctions",
    "cs",
];

/// Number of per-side data columns. `AlignmentRow::raw` is sized from this, so the
/// array and the header loop are provably the same length.
const N_SIDE_COLS: usize = READINFO_DATA_COLS.len();

/// The A-vs-B comparison columns, emitted after both per-side blocks. Grouped to
/// mirror `READINFO_DATA_COLS`: orientation/identity, score, event diffs,
/// query-space junctions, genomic-space junctions, then the object lists.
///
/// `AlignmentDiff`'s field order and `write_tsv_row`'s emission order must both
/// match this list; `header_and_row_field_counts_agree` pins the length.
pub(crate) fn comparison_col_names() -> Vec<&'static str> {
    vec![
        // Orientation & identity.
        "Strand_Match",
        "seqid_Diff",
        "QueryAlnCov_Diff",
        "QueryAlnLen_Diff",
        // Score.
        "AS_Diff",
        "ms_Diff",
        "AS_Ratio",
        "ms_Ratio",
        // Event-count diffs.
        "N_Substitution_Bases_Diff",
        "N_Substitution_Bases_Ratio",
        "N_Insertion_Bases_Diff",
        "N_Insertion_Bases_Ratio",
        "N_Deletion_Bases_Diff",
        "N_Deletion_Bases_Ratio",
        "N_SoftClipped_Bases_Start_Diff",
        "N_SoftClipped_Bases_End_Diff",
        // Query-space junction set comparison.
        "N_Matched_Junctions",
        "N_Unmatched_Junctions",
        "N_Junctions_OnlyA",
        "N_Junctions_OnlyB",
        "Junction_Distance",
        "Junc_Dist_V2",
        // Genomic-space junction set comparison.
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

// ── Compile-time column indices ──────────────────────────────────────────────

/// `str` equality usable in a `const fn` (`==` on `&str` is not const yet).
const fn str_eq(a: &str, b: &str) -> bool {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    if a.len() != b.len() {
        return false;
    }
    let mut i = 0;
    while i < a.len() {
        if a[i] != b[i] {
            return false;
        }
        i += 1;
    }
    true
}

/// Index of a column in `READINFO_DATA_COLS`, resolved at compile time. A typo in
/// a column name is therefore a build failure, not a silently-empty field.
const fn col_idx(name: &str) -> usize {
    let mut i = 0;
    while i < READINFO_DATA_COLS.len() {
        if str_eq(READINFO_DATA_COLS[i], name) {
            return i;
        }
        i += 1;
    }
    panic!("column name is not in READINFO_DATA_COLS")
}

const I_TARGET_CHR: usize = col_idx("TargetChr");
const I_STRAND: usize = col_idx("Strand");
const I_JUNCTIONS: usize = col_idx("junctions");
const I_GENOMIC_JUNCTIONS: usize = col_idx("genomic_junctions");

// ── Ratio helpers ────────────────────────────────────────────────────────────

/// `num / den`, or `NaN` when the denominator is zero.
///
/// The zero test is on the **integer**, before the cast, which is what the
/// pre-v0.14.1 `safe_ratio_i64` did. That matters: `5.0 / 0.0` is `inf` in IEEE
/// arithmetic, so a bare float division would render `inf` for every read whose
/// denominator side scored 0, where this renders `NaN`.
fn ratio_i64(num: i64, den: i64) -> f64 {
    if den == 0 {
        f64::NAN
    } else {
        num as f64 / den as f64
    }
}

/// Unsigned counterpart of [`ratio_i64`], with the same zero-denominator rule.
fn ratio_u64(num: u64, den: u64) -> f64 {
    if den == 0 {
        f64::NAN
    } else {
        num as f64 / den as f64
    }
}

// ── Per-side types ───────────────────────────────────────────────────────────

/// The values a side contributes to the comparison math, parsed once.
///
/// Defaults match the pre-v0.14.1 inline `.parse().unwrap_or(..)` calls exactly:
/// `0` for the integer fields, `f64::NAN` for the two float fields.
#[derive(Debug, Clone, Copy)]
pub(crate) struct AlignmentMetrics {
    pub(crate) as_max: i64,
    pub(crate) ms_max: i64,
    pub(crate) seqid_max: f64,
    pub(crate) query_aln_cov_max: f64,
    pub(crate) query_aln_len_max: u64,
    pub(crate) n_ins_bases: u64,
    pub(crate) n_del_bases: u64,
    pub(crate) n_sub_bases: u64,
    pub(crate) n_sc_start: u64,
    pub(crate) n_sc_end: u64,
    /// `usize`, as before v0.14.1 — feeds `Junc_Dist_V2`.
    pub(crate) junc_count: usize,
}

/// One side of the comparison: the raw transport columns plus the parsed view.
pub(crate) struct AlignmentRow<'r> {
    /// The 31 per-side columns in `READINFO_DATA_COLS` order, exactly as read.
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
        let raw = std::array::from_fn(|i| get(READINFO_DATA_COLS[i]));
        let metrics = AlignmentMetrics {
            as_max: get("AS_Max").parse().unwrap_or(0),
            ms_max: get("ms_Max").parse().unwrap_or(0),
            seqid_max: get("seqid_Max").parse().unwrap_or(f64::NAN),
            query_aln_cov_max: get("Query_Aln_Cov_Max").parse().unwrap_or(f64::NAN),
            query_aln_len_max: get("Query_Aln_Len_Max").parse().unwrap_or(0),
            n_ins_bases: get("N_Insertion_Bases").parse().unwrap_or(0),
            n_del_bases: get("N_Deletion_Bases").parse().unwrap_or(0),
            n_sub_bases: get("N_Substitution_Bases").parse().unwrap_or(0),
            n_sc_start: get("N_SoftClipped_Bases_Start").parse().unwrap_or(0),
            n_sc_end: get("N_SoftClipped_Bases_End").parse().unwrap_or(0),
            junc_count: get("JuncCount").parse().unwrap_or(0),
        };
        Self { raw, metrics }
    }

    fn target_chr(&self) -> &'r str {
        self.raw[I_TARGET_CHR]
    }
    fn strand(&self) -> &'r str {
        self.raw[I_STRAND]
    }
    fn junctions(&self) -> &'r str {
        self.raw[I_JUNCTIONS]
    }
    fn genomic_junctions(&self) -> &'r str {
        self.raw[I_GENOMIC_JUNCTIONS]
    }
}

// ── The computed comparison block ────────────────────────────────────────────

/// The A-vs-B metrics. Field order matches [`comparison_col_names`].
///
/// Sign conventions, preserved verbatim from the pre-v0.14.1 emitter: every
/// `*_diff` is **B − A**, every ratio is **B / A** (so the denominator is the A
/// side, and A is what the zero-guard tests) — except `junc_dist_v2`, whose inner
/// subtraction is **A − B**. `unsigned_abs()` makes that last one numerically
/// moot, but it is left as it was rather than silently normalized.
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
    pub(crate) n_sub_bases_ratio: f64,
    pub(crate) n_ins_bases_diff: i64,
    pub(crate) n_ins_bases_ratio: f64,
    pub(crate) n_del_bases_diff: i64,
    pub(crate) n_del_bases_ratio: f64,
    pub(crate) n_sc_start_diff: i64,
    pub(crate) n_sc_end_diff: i64,
    pub(crate) n_matched_junctions: u64,
    pub(crate) n_unmatched_junctions: u64,
    pub(crate) n_junctions_only_a: u64,
    pub(crate) n_junctions_only_b: u64,
    pub(crate) junction_distance: u64,
    pub(crate) junc_dist_v2: u64,
    pub(crate) genomic_n_matched_junctions: u64,
    pub(crate) genomic_n_unmatched_junctions: u64,
    pub(crate) genomic_n_junctions_only_a: u64,
    pub(crate) genomic_n_junctions_only_b: u64,
    pub(crate) junctions_only_a: String,
    pub(crate) junctions_only_b: String,
    pub(crate) genomic_junctions_only_a: String,
    pub(crate) genomic_junctions_only_b: String,
}

impl AlignmentDiff {
    pub(crate) fn compute(a: &AlignmentRow, b: &AlignmentRow) -> Self {
        let (ma, mb) = (&a.metrics, &b.metrics);

        // Query-space junction sets.
        let juncs_a = parse_junction_str(a.junctions());
        let juncs_b = parse_junction_str(b.junctions());
        let (n_matched, n_only_a, n_only_b) = junction_set_stats(&juncs_a, &juncs_b);
        let (j_only_a_vec, j_only_b_vec) = junction_set_diffs(&juncs_a, &juncs_b);

        // Genomic-space junction sets. The chrom is reattached to every pair (it was
        // dropped from the stored tuples in v0.2.3) so that junctions on different
        // contigs cannot match; `format_genomic_junction_tuple` drops it again on the
        // way out. Do not unify those two behaviours.
        let chrom_a = a.target_chr().to_string();
        let chrom_b = b.target_chr().to_string();
        let gj_a: Vec<(String, u64, u64)> = parse_genomic_junction_str(a.genomic_junctions())
            .into_iter()
            .map(|(s, e)| (chrom_a.clone(), s, e))
            .collect();
        let gj_b: Vec<(String, u64, u64)> = parse_genomic_junction_str(b.genomic_junctions())
            .into_iter()
            .map(|(s, e)| (chrom_b.clone(), s, e))
            .collect();
        let (g_matched, g_only_a, g_only_b) = genomic_junction_set_stats(&gj_a, &gj_b);
        let (gj_only_a_vec, gj_only_b_vec) = genomic_junction_set_diffs(&gj_a, &gj_b);

        Self {
            strand_match: a.strand() == b.strand(),
            seqid_diff: mb.seqid_max - ma.seqid_max,
            query_aln_cov_diff: mb.query_aln_cov_max - ma.query_aln_cov_max,
            query_aln_len_diff: mb.query_aln_len_max as i64 - ma.query_aln_len_max as i64,
            as_diff: mb.as_max - ma.as_max,
            ms_diff: mb.ms_max - ma.ms_max,
            as_ratio: ratio_i64(mb.as_max, ma.as_max),
            ms_ratio: ratio_i64(mb.ms_max, ma.ms_max),
            n_sub_bases_diff: mb.n_sub_bases as i64 - ma.n_sub_bases as i64,
            n_sub_bases_ratio: ratio_u64(mb.n_sub_bases, ma.n_sub_bases),
            n_ins_bases_diff: mb.n_ins_bases as i64 - ma.n_ins_bases as i64,
            n_ins_bases_ratio: ratio_u64(mb.n_ins_bases, ma.n_ins_bases),
            n_del_bases_diff: mb.n_del_bases as i64 - ma.n_del_bases as i64,
            n_del_bases_ratio: ratio_u64(mb.n_del_bases, ma.n_del_bases),
            n_sc_start_diff: mb.n_sc_start as i64 - ma.n_sc_start as i64,
            n_sc_end_diff: mb.n_sc_end as i64 - ma.n_sc_end as i64,
            n_matched_junctions: n_matched,
            n_unmatched_junctions: n_only_a + n_only_b,
            n_junctions_only_a: n_only_a,
            n_junctions_only_b: n_only_b,
            junction_distance: junction_distance(&juncs_a, &juncs_b),
            // NOTE: A − B here, unlike every other diff. Preserved as-is.
            junc_dist_v2: 50 * (ma.junc_count as i64 - mb.junc_count as i64).unsigned_abs(),
            genomic_n_matched_junctions: g_matched,
            genomic_n_unmatched_junctions: g_only_a + g_only_b,
            genomic_n_junctions_only_a: g_only_a,
            genomic_n_junctions_only_b: g_only_b,
            junctions_only_a: format_junction_tuple(&j_only_a_vec),
            junctions_only_b: format_junction_tuple(&j_only_b_vec),
            genomic_junctions_only_a: format_genomic_junction_tuple(&gj_only_a_vec),
            genomic_junctions_only_b: format_genomic_junction_tuple(&gj_only_b_vec),
        }
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
    /// Escaping matches the pre-v0.14.1 emitter exactly: the 62 per-side fields and
    /// the 4 object-list fields are escaped; `Read_Name` and the two labels are
    /// **not** (labels are validated by `validate_set_label` instead). Note the
    /// per-side values were already escaped once when the readinfo row was written,
    /// so this is a second pass over them — `escape_tsv_field` is not idempotent,
    /// and that double escaping is preserved behaviour.
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
            "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            d.n_sub_bases_diff,
            fmt_float(d.n_sub_bases_ratio),
            d.n_ins_bases_diff,
            fmt_float(d.n_ins_bases_ratio),
            d.n_del_bases_diff,
            fmt_float(d.n_del_bases_ratio),
            d.n_sc_start_diff,
            d.n_sc_end_diff,
        )?;
        // Query-space junction set comparison.
        write!(
            out,
            "\t{}\t{}\t{}\t{}\t{}\t{}",
            d.n_matched_junctions,
            d.n_unmatched_junctions,
            d.n_junctions_only_a,
            d.n_junctions_only_b,
            d.junction_distance,
            d.junc_dist_v2,
        )?;
        // Genomic-space junction set comparison.
        write!(
            out,
            "\t{}\t{}\t{}\t{}",
            d.genomic_n_matched_junctions,
            d.genomic_n_unmatched_junctions,
            d.genomic_n_junctions_only_a,
            d.genomic_n_junctions_only_b,
        )?;
        // Object lists last.
        write!(
            out,
            "\t{}\t{}\t{}\t{}",
            escape_tsv_field(&d.junctions_only_a),
            escape_tsv_field(&d.junctions_only_b),
            escape_tsv_field(&d.genomic_junctions_only_a),
            escape_tsv_field(&d.genomic_junctions_only_b),
        )?;
        writeln!(out)
    }
}

/// Write the comparison header: the join keys, the two set-label columns, the
/// per-side data columns (suffixed `_A` / `_B`), then the comparison columns.
///
/// Side suffixes are **fixed** (`_A` / `_B`), never the user's label — the
/// human-readable labels are carried as the `Label_A` / `Label_B` data columns
/// instead, so column names are stable across datasets and unambiguous even when a
/// label itself contains an underscore.
pub(crate) fn write_compare_header<W: Write>(out: &mut W) -> std::io::Result<()> {
    write!(out, "Read_Name\tRead_Len\tLabel_A\tLabel_B")?;
    for col in READINFO_DATA_COLS {
        write!(out, "\t{col}_A")?;
    }
    for col in READINFO_DATA_COLS {
        write!(out, "\t{col}_B")?;
    }
    for col in comparison_col_names() {
        write!(out, "\t{col}")?;
    }
    writeln!(out)
}

/// Emit one comparison row given a by-name column accessor for each side.
///
/// Single source of truth for the `compare` row, shared by `compare-readinfo`
/// (reading from a `ReadInfoReader`) and the fused `compare` pass (reading from
/// in-memory `ReadInfoRow`s serialized to readinfo lines). Accessors return `""`
/// for an absent column, which flows through to an empty output cell.
pub(crate) fn emit_compare_row<'r, W, FA, FB>(
    out: &mut W,
    name: &'r str,
    len: u64,
    label_a: &'r str,
    label_b: &'r str,
    get_a: FA,
    get_b: FB,
) -> std::io::Result<()>
where
    W: Write,
    FA: Fn(&str) -> &'r str,
    FB: Fn(&str) -> &'r str,
{
    ComparisonRow::build(name, len, label_a, label_b, get_a, get_b).write_tsv_row(out)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

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

    /// Render one full row and split it into its 96 fields.
    fn render(a: &[(&str, &'static str)], b: &[(&str, &'static str)]) -> Vec<String> {
        let mut buf: Vec<u8> = Vec::new();
        emit_compare_row(
            &mut buf,
            "read1",
            100,
            "SetA",
            "SetB",
            |c| lookup(a, c),
            |c| lookup(b, c),
        )
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

    /// Index of an output column in the 96-column header.
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
        assert_eq!(h.len(), 96, "header should have 96 columns");
        assert_eq!(
            r.len(),
            h.len(),
            "row emits {} fields but the header names {} — the two lists have drifted",
            r.len(),
            h.len()
        );
        assert_eq!(
            4 + 2 * READINFO_DATA_COLS.len() + comparison_col_names().len(),
            96
        );
    }

    #[test]
    fn header_layout_is_keys_then_a_then_b_then_comparison() {
        let h = header_fields();
        assert_eq!(&h[..4], &["Read_Name", "Read_Len", "Label_A", "Label_B"]);
        assert_eq!(h[4], "TargetChr_A");
        assert_eq!(h[34], "cs_A");
        assert_eq!(h[35], "TargetChr_B");
        assert_eq!(h[65], "cs_B");
        assert_eq!(h[66], "Strand_Match");
        assert_eq!(h[95], "Genomic_Junctions_OnlyB");
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
        assert_eq!(r.len(), 96, "an embedded tab must not add a field");
    }

    #[test]
    fn keys_and_labels_are_not_escaped() {
        // Labels are validated by `validate_set_label` rather than escaped, and
        // Read_Name is passed through as-is. Backslashes must survive untouched.
        let mut buf: Vec<u8> = Vec::new();
        emit_compare_row(&mut buf, "read\\1", 100, "Set\\A", "Set\\B", |_| "", |_| "").unwrap();
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
        // give inf here, changing AS_Ratio for every read whose A side scored 0.
        assert!(ratio_i64(5, 0).is_nan(), "5/0 must be NaN, not inf");
        assert!(ratio_i64(-5, 0).is_nan(), "-5/0 must be NaN, not -inf");
        assert!(ratio_i64(0, 0).is_nan());
        assert!(ratio_u64(5, 0).is_nan());
        assert!(ratio_u64(0, 0).is_nan());
        assert_eq!(fmt_float(ratio_i64(5, 0)), "NaN");
        assert_eq!(fmt_float(ratio_i64(-5, 0)), "NaN");
    }

    #[test]
    fn all_four_ratio_columns_render_nan_when_the_a_side_is_zero() {
        let a = [
            ("AS_Max", "0"),
            ("ms_Max", "0"),
            ("N_Substitution_Bases", "0"),
            ("N_Insertion_Bases", "0"),
            ("N_Deletion_Bases", "0"),
        ];
        let b = [
            ("AS_Max", "50"),
            ("ms_Max", "60"),
            ("N_Substitution_Bases", "7"),
            ("N_Insertion_Bases", "8"),
            ("N_Deletion_Bases", "9"),
        ];
        let r = render(&a, &b);
        for c in [
            "AS_Ratio",
            "ms_Ratio",
            "N_Substitution_Bases_Ratio",
            "N_Insertion_Bases_Ratio",
            "N_Deletion_Bases_Ratio",
        ] {
            assert_eq!(r[col(c)], "NaN", "{c} should be NaN when the A side is 0");
        }
    }

    #[test]
    fn ratios_are_b_over_a() {
        let a = [("AS_Max", "10"), ("ms_Max", "10")];
        let b = [("AS_Max", "20"), ("ms_Max", "5")];
        let d = diff_of(&a, &b);
        assert_eq!(d.as_ratio, 2.0, "AS_Ratio is B/A");
        assert_eq!(d.ms_ratio, 0.5, "ms_Ratio is B/A");
    }

    // ── sign conventions ─────────────────────────────────────────────────────

    #[test]
    fn every_diff_is_b_minus_a() {
        let a = [
            ("AS_Max", "10"),
            ("ms_Max", "10"),
            ("Query_Aln_Len_Max", "100"),
            ("N_Substitution_Bases", "1"),
            ("N_Insertion_Bases", "2"),
            ("N_Deletion_Bases", "3"),
            ("N_SoftClipped_Bases_Start", "4"),
            ("N_SoftClipped_Bases_End", "5"),
            ("seqid_Max", "0.9"),
            ("Query_Aln_Cov_Max", "0.8"),
        ];
        let b = [
            ("AS_Max", "13"),
            ("ms_Max", "14"),
            ("Query_Aln_Len_Max", "110"),
            ("N_Substitution_Bases", "2"),
            ("N_Insertion_Bases", "4"),
            ("N_Deletion_Bases", "6"),
            ("N_SoftClipped_Bases_Start", "9"),
            ("N_SoftClipped_Bases_End", "1"),
            ("seqid_Max", "1.0"),
            ("Query_Aln_Cov_Max", "0.5"),
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

    #[test]
    fn junc_dist_v2_uses_the_absolute_count_difference() {
        // Its inner subtraction is A - B, unlike every other diff; unsigned_abs()
        // makes the result symmetric, which is what this pins.
        let fwd = diff_of(&[("JuncCount", "2")], &[("JuncCount", "5")]);
        let rev = diff_of(&[("JuncCount", "5")], &[("JuncCount", "2")]);
        assert_eq!(fwd.junc_dist_v2, 150);
        assert_eq!(rev.junc_dist_v2, 150);
        assert_eq!(
            diff_of(&[("JuncCount", "3")], &[("JuncCount", "3")]).junc_dist_v2,
            0
        );
    }

    // ── strand ───────────────────────────────────────────────────────────────

    #[test]
    fn strand_match_compares_the_raw_strand_strings() {
        assert!(diff_of(&[("Strand", "+")], &[("Strand", "+")]).strand_match);
        assert!(!diff_of(&[("Strand", "+")], &[("Strand", "-")]).strand_match);
        // Two unmapped sides agree with each other.
        assert!(diff_of(&[("Strand", "*")], &[("Strand", "*")]).strand_match);
        assert!(!diff_of(&[("Strand", "*")], &[("Strand", "+")]).strand_match);
    }

    #[test]
    fn strand_match_renders_lowercase() {
        let r = render(&[("Strand", "+")], &[("Strand", "+")]);
        assert_eq!(r[col("Strand_Match")], "true");
        let r = render(&[("Strand", "+")], &[("Strand", "-")]);
        assert_eq!(r[col("Strand_Match")], "false");
    }

    // ── parse defaults and the unmapped case ─────────────────────────────────

    #[test]
    fn missing_columns_parse_to_the_documented_defaults() {
        let m = side(&[]).metrics;
        assert_eq!(m.as_max, 0);
        assert_eq!(m.ms_max, 0);
        assert_eq!(m.query_aln_len_max, 0);
        assert_eq!(m.n_ins_bases, 0);
        assert_eq!(m.n_del_bases, 0);
        assert_eq!(m.n_sub_bases, 0);
        assert_eq!(m.n_sc_start, 0);
        assert_eq!(m.n_sc_end, 0);
        assert_eq!(m.junc_count, 0);
        assert!(m.seqid_max.is_nan(), "seqid_Max defaults to NaN, not 0");
        assert!(m.query_aln_cov_max.is_nan());
    }

    #[test]
    fn a_missing_column_stays_an_empty_cell_and_is_not_rendered_as_zero() {
        // The per-side columns are transported as &str precisely so this holds. If
        // they were re-rendered from parsed values, these would become "0".
        let r = render(&[], &[]);
        for c in ["N_Match_Events_A", "Target_Start_A", "cs_B", "MQ_Best_B"] {
            assert_eq!(r[col(c)], "", "{c} should stay empty, not become 0");
        }
    }

    #[test]
    fn unmapped_side_propagates_nan_through_the_float_diffs() {
        // What an unmapped read looks like coming out of readinfo: seqid/cov are NaN.
        let unmapped = [("TargetChr", "*"), ("Strand", "*"), ("seqid_Max", "NaN"),
                        ("Query_Aln_Cov_Max", "NaN"), ("cs", ""), ("junctions", "()")];
        let mapped = [("TargetChr", "chr1"), ("Strand", "+"), ("seqid_Max", "0.99"),
                      ("Query_Aln_Cov_Max", "1.0"), ("cs", ":100"), ("junctions", "()")];
        let r = render(&unmapped, &mapped);
        assert_eq!(r[col("seqid_Diff")], "NaN");
        assert_eq!(r[col("QueryAlnCov_Diff")], "NaN");
        assert_eq!(r[col("TargetChr_A")], "*");
        assert_eq!(r[col("cs_A")], "", "an unmapped cs is an empty cell");
    }

    // ── junction set comparison wiring ───────────────────────────────────────

    #[test]
    fn junction_counts_and_object_lists_line_up() {
        let a = [("junctions", "(10, 20, 30)")];
        let b = [("junctions", "(20, 30, 40)")];
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

    #[test]
    fn genomic_junctions_on_different_contigs_never_match() {
        // The chrom is reattached before the set comparison, so identical coordinate
        // pairs on different contigs must not overlap.
        let a = [("TargetChr", "chr1"), ("genomic_junctions", "((100, 200),)")];
        let b = [("TargetChr", "chr2"), ("genomic_junctions", "((100, 200),)")];
        let d = diff_of(&a, &b);
        assert_eq!(d.genomic_n_matched_junctions, 0);
        assert_eq!(d.genomic_n_junctions_only_a, 1);
        assert_eq!(d.genomic_n_junctions_only_b, 1);
        // ...but the rendered tuples drop the chrom again.
        assert_eq!(d.genomic_junctions_only_a, "((100, 200),)");
    }

    #[test]
    fn same_contig_genomic_junctions_match() {
        let a = [("TargetChr", "chr1"), ("genomic_junctions", "((100, 200),)")];
        let b = [("TargetChr", "chr1"), ("genomic_junctions", "((100, 200),)")];
        let d = diff_of(&a, &b);
        assert_eq!(d.genomic_n_matched_junctions, 1);
        assert_eq!(d.genomic_n_unmatched_junctions, 0);
        assert_eq!(d.genomic_junctions_only_a, "()");
    }
}
