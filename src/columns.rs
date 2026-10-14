//! Every column name of the alninfo, readinfo and comparison tables, defined once.
//!
//! Each column is a [`Col`] constant: its header name plus its type. The ordered
//! lists below ([`ALNINFO_COLS`], [`READINFO_COLS`], [`COMPARE_KEY_COLS`],
//! [`COMPARE_SIDE_COLS`], [`COMPARE_DIFF_COLS`]) are what the writers build their
//! headers from, what the Parquet schema is derived from, and what the readers
//! resolve column indices against. Code outside this module refers to a column
//! through its constant, never through a string literal, so renaming a column is
//! a one-line edit here and a misspelled constant is a build error.
//!
//! Naming: lowercase snake_case; counts start with `n_`; per-side comparison
//! columns end in `_A` / `_B`; a comparison column is `{col}_diff` (B − A) or
//! `{col}_ratio` (B / A) of the per-side column `{col}`. A value that appears in
//! several tables has one name (and one constant) everywhere.

/// The type of a column's values. Drives the Parquet schema for the comparison
/// table; for alninfo/readinfo it documents what the TSV cell holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ColType {
    Str,
    U64,
    I64,
    F64,
    Bool,
}

/// One output column: its header name and value type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Col {
    pub(crate) name: &'static str,
    pub(crate) ty: ColType,
}

const fn col(name: &'static str, ty: ColType) -> Col {
    Col { name, ty }
}

use ColType::{Bool, Str, F64, I64, U64};

// ── Side suffixes ────────────────────────────────────────────────────────────

/// Side identifier for dataset A: the suffix of every A-side comparison column.
pub(crate) const SIDE_A: &str = "A";
/// Side identifier for dataset B.
pub(crate) const SIDE_B: &str = "B";
/// Both sides, in output order.
pub(crate) const SIDES: [&str; 2] = [SIDE_A, SIDE_B];

/// A per-side column's name in the comparison table: `{name}_{side}`.
pub(crate) fn sided(col: Col, side: &str) -> String {
    format!("{}_{side}", col.name)
}

// ── Per-alignment columns (alninfo, and readinfo for the representative) ────
//
// readinfo describes each read by its representative alignment, so every
// per-alignment column it carries has the same name (and constant) as in alninfo.

pub(crate) const QUERY_NAME: Col = col("query_name", Str);
pub(crate) const QUERY_LEN: Col = col("query_len", U64);
pub(crate) const QUERY_START: Col = col("query_start", U64);
pub(crate) const QUERY_END: Col = col("query_end", U64);
pub(crate) const STRAND: Col = col("strand", Str);
pub(crate) const TARGET_NAME: Col = col("target_name", Str);
pub(crate) const TARGET_START: Col = col("target_start", U64);
pub(crate) const TARGET_END: Col = col("target_end", U64);
pub(crate) const MAPQ: Col = col("mapq", U64);
pub(crate) const MS_SCORE: Col = col("ms_score", I64);
pub(crate) const AS_SCORE: Col = col("as_score", I64);
pub(crate) const CS: Col = col("cs", Str);
pub(crate) const N_MATCH_EVENTS: Col = col("n_match_events", U64);
pub(crate) const N_MATCH_BASES: Col = col("n_match_bases", U64);
pub(crate) const N_SUBSTITUTION_EVENTS: Col = col("n_substitution_events", U64);
pub(crate) const N_SUBSTITUTION_BASES: Col = col("n_substitution_bases", U64);
pub(crate) const N_INSERTION_EVENTS: Col = col("n_insertion_events", U64);
pub(crate) const N_INSERTION_BASES: Col = col("n_insertion_bases", U64);
pub(crate) const N_DELETION_EVENTS: Col = col("n_deletion_events", U64);
pub(crate) const N_DELETION_BASES: Col = col("n_deletion_bases", U64);
pub(crate) const N_INTRON_EVENTS: Col = col("n_intron_events", U64);
pub(crate) const N_INTRON_BASES: Col = col("n_intron_bases", U64);
pub(crate) const N_SOFTCLIP_BASES_START: Col = col("n_softclip_bases_start", U64);
pub(crate) const N_SOFTCLIP_BASES_END: Col = col("n_softclip_bases_end", U64);
pub(crate) const N_SOFTCLIP_EVENTS: Col = col("n_softclip_events", U64);
pub(crate) const QUERY_JUNCTIONS: Col = col("query_junctions", Str);
pub(crate) const N_QUERY_JUNCTIONS: Col = col("n_query_junctions", U64);
pub(crate) const GENOMIC_JUNCTIONS: Col = col("genomic_junctions", Str);
pub(crate) const ALN_IDENTITY: Col = col("aln_identity", F64);
pub(crate) const QUERY_ALN_LEN: Col = col("query_aln_len", U64);
pub(crate) const QUERY_ALN_COV: Col = col("query_aln_cov", F64);
pub(crate) const TP_TAG: Col = col("tp_tag", Str);

// ── alninfo-only columns ─────────────────────────────────────────────────────

pub(crate) const TARGET_LEN: Col = col("target_len", U64);
pub(crate) const N_RESIDUE_MATCHES: Col = col("n_residue_matches", U64);
pub(crate) const ALN_BLOCK_LEN: Col = col("aln_block_len", U64);
pub(crate) const TARGET_START_1BASED: Col = col("target_start_1based", U64);

// ── readinfo-only columns: counts over all of a read's alignments ────────────

pub(crate) const N_ALNS: Col = col("n_alns", U64);
pub(crate) const N_ALNS_TIED_TOP_SCORE: Col = col("n_alns_tied_top_score", U64);
pub(crate) const N_ALNS_TP_P: Col = col("n_alns_tp_p", U64);
pub(crate) const N_ALNS_TP_S: Col = col("n_alns_tp_s", U64);

// ── Table layouts ────────────────────────────────────────────────────────────

/// The alninfo table, in output order. `AlnInfo::write_row` emits its fields in
/// this order (pinned by a test in `record.rs`).
pub(crate) const ALNINFO_COLS: [Col; 36] = [
    QUERY_NAME,
    QUERY_LEN,
    QUERY_START,
    QUERY_END,
    STRAND,
    TARGET_NAME,
    TARGET_LEN,
    TARGET_START,
    TARGET_END,
    N_RESIDUE_MATCHES,
    ALN_BLOCK_LEN,
    MAPQ,
    MS_SCORE,
    AS_SCORE,
    CS,
    N_MATCH_EVENTS,
    N_MATCH_BASES,
    N_SUBSTITUTION_EVENTS,
    N_SUBSTITUTION_BASES,
    N_INSERTION_EVENTS,
    N_INSERTION_BASES,
    N_DELETION_EVENTS,
    N_DELETION_BASES,
    N_INTRON_EVENTS,
    N_INTRON_BASES,
    N_SOFTCLIP_BASES_START,
    N_SOFTCLIP_BASES_END,
    N_SOFTCLIP_EVENTS,
    QUERY_JUNCTIONS,
    N_QUERY_JUNCTIONS,
    TARGET_START_1BASED,
    ALN_IDENTITY,
    QUERY_ALN_LEN,
    QUERY_ALN_COV,
    GENOMIC_JUNCTIONS,
    TP_TAG,
];

/// The readinfo table, in output order. `collapse_group` builds its fields in
/// this order (pinned by a test in `readinfo.rs`).
pub(crate) const READINFO_COLS: [Col; 36] = [
    QUERY_NAME,
    QUERY_LEN,
    TARGET_NAME,
    STRAND,
    MAPQ,
    AS_SCORE,
    MS_SCORE,
    QUERY_ALN_COV,
    QUERY_ALN_LEN,
    ALN_IDENTITY,
    QUERY_JUNCTIONS,
    N_ALNS,
    N_ALNS_TIED_TOP_SCORE,
    N_QUERY_JUNCTIONS,
    N_MATCH_EVENTS,
    N_MATCH_BASES,
    N_SUBSTITUTION_EVENTS,
    N_SUBSTITUTION_BASES,
    N_INSERTION_EVENTS,
    N_INSERTION_BASES,
    N_DELETION_EVENTS,
    N_DELETION_BASES,
    N_INTRON_EVENTS,
    N_INTRON_BASES,
    N_SOFTCLIP_BASES_START,
    N_SOFTCLIP_BASES_END,
    N_SOFTCLIP_EVENTS,
    CS,
    GENOMIC_JUNCTIONS,
    QUERY_START,
    QUERY_END,
    TARGET_START,
    TARGET_END,
    TP_TAG,
    N_ALNS_TP_P,
    N_ALNS_TP_S,
];

// ── Comparison table ─────────────────────────────────────────────────────────

pub(crate) const LABEL_A: Col = col("label_A", Str);
pub(crate) const LABEL_B: Col = col("label_B", Str);

/// The comparison table's leading columns: the join keys and the two set labels.
pub(crate) const COMPARE_KEY_COLS: [Col; 4] = [QUERY_NAME, QUERY_LEN, LABEL_A, LABEL_B];

/// The per-side columns, each emitted once suffixed `_A` and once `_B`, read out
/// of each side's readinfo row by name.
///
/// This order is deliberately different from [`READINFO_COLS`]. It is grouped by
/// topic (locus, alignment selection, identity/coverage, junction counts,
/// cs-derived event counts, and the three long strings last) so the comparison
/// table is readable via `head -1 | tr '\t' '\n' | nl`. Do not "resync" the two.
pub(crate) const COMPARE_SIDE_COLS: [Col; 33] = [
    // Locus & span.
    TARGET_NAME,
    STRAND,
    TARGET_START,
    TARGET_END,
    QUERY_START,
    QUERY_END,
    // Alignment selection & score.
    MAPQ,
    N_ALNS,
    N_ALNS_TP_P,
    N_ALNS_TP_S,
    N_ALNS_TIED_TOP_SCORE,
    AS_SCORE,
    MS_SCORE,
    // Identity & coverage.
    ALN_IDENTITY,
    QUERY_ALN_COV,
    QUERY_ALN_LEN,
    // Junction counts.
    N_QUERY_JUNCTIONS,
    N_INTRON_EVENTS,
    N_INTRON_BASES,
    // cs-derived event counts.
    N_MATCH_EVENTS,
    N_MATCH_BASES,
    N_SUBSTITUTION_EVENTS,
    N_SUBSTITUTION_BASES,
    N_INSERTION_EVENTS,
    N_INSERTION_BASES,
    N_DELETION_EVENTS,
    N_DELETION_BASES,
    N_SOFTCLIP_EVENTS,
    N_SOFTCLIP_BASES_START,
    N_SOFTCLIP_BASES_END,
    // Long strings last — these dominate the table's on-disk size.
    QUERY_JUNCTIONS,
    GENOMIC_JUNCTIONS,
    CS,
];

pub(crate) const STRAND_MATCH: Col = col("strand_match", Bool);
pub(crate) const ALN_IDENTITY_DIFF: Col = col("aln_identity_diff", F64);
pub(crate) const QUERY_ALN_COV_DIFF: Col = col("query_aln_cov_diff", F64);
pub(crate) const QUERY_ALN_LEN_DIFF: Col = col("query_aln_len_diff", I64);
pub(crate) const AS_SCORE_DIFF: Col = col("as_score_diff", I64);
pub(crate) const MS_SCORE_DIFF: Col = col("ms_score_diff", I64);
pub(crate) const AS_SCORE_RATIO: Col = col("as_score_ratio", F64);
pub(crate) const MS_SCORE_RATIO: Col = col("ms_score_ratio", F64);
pub(crate) const N_SUBSTITUTION_BASES_DIFF: Col = col("n_substitution_bases_diff", I64);
pub(crate) const N_SUBSTITUTION_BASES_RATIO: Col = col("n_substitution_bases_ratio", F64);
pub(crate) const N_INSERTION_BASES_DIFF: Col = col("n_insertion_bases_diff", I64);
pub(crate) const N_INSERTION_BASES_RATIO: Col = col("n_insertion_bases_ratio", F64);
pub(crate) const N_DELETION_BASES_DIFF: Col = col("n_deletion_bases_diff", I64);
pub(crate) const N_DELETION_BASES_RATIO: Col = col("n_deletion_bases_ratio", F64);
pub(crate) const N_SOFTCLIP_BASES_START_DIFF: Col = col("n_softclip_bases_start_diff", I64);
pub(crate) const N_SOFTCLIP_BASES_END_DIFF: Col = col("n_softclip_bases_end_diff", I64);
pub(crate) const N_QUERY_JUNCTIONS_SHARED: Col = col("n_query_junctions_shared", U64);
pub(crate) const N_QUERY_JUNCTIONS_NOT_SHARED: Col = col("n_query_junctions_not_shared", U64);
pub(crate) const N_QUERY_JUNCTIONS_ONLY_A: Col = col("n_query_junctions_only_A", U64);
pub(crate) const N_QUERY_JUNCTIONS_ONLY_B: Col = col("n_query_junctions_only_B", U64);
pub(crate) const QUERY_JUNCTION_POS_DISTANCE: Col = col("query_junction_pos_distance", U64);
pub(crate) const N_GENOMIC_JUNCTIONS_SHARED: Col = col("n_genomic_junctions_shared", U64);
pub(crate) const N_GENOMIC_JUNCTIONS_NOT_SHARED: Col = col("n_genomic_junctions_not_shared", U64);
pub(crate) const N_GENOMIC_JUNCTIONS_ONLY_A: Col = col("n_genomic_junctions_only_A", U64);
pub(crate) const N_GENOMIC_JUNCTIONS_ONLY_B: Col = col("n_genomic_junctions_only_B", U64);
pub(crate) const QUERY_JUNCTIONS_ONLY_A: Col = col("query_junctions_only_A", Str);
pub(crate) const QUERY_JUNCTIONS_ONLY_B: Col = col("query_junctions_only_B", Str);
pub(crate) const GENOMIC_JUNCTIONS_ONLY_A: Col = col("genomic_junctions_only_A", Str);
pub(crate) const GENOMIC_JUNCTIONS_ONLY_B: Col = col("genomic_junctions_only_B", Str);

/// The A-vs-B comparison columns, emitted after both per-side blocks. Grouped to
/// mirror [`COMPARE_SIDE_COLS`]: orientation/identity, score, event diffs,
/// query-space junctions, genomic-space junctions, then the object lists.
///
/// `AlignmentDiff::values()` and `write_tsv_row` both emit in this order; the
/// Parquet round-trip test in `parquet_out.rs` checks them against each other.
pub(crate) const COMPARE_DIFF_COLS: [Col; 29] = [
    // Orientation & identity.
    STRAND_MATCH,
    ALN_IDENTITY_DIFF,
    QUERY_ALN_COV_DIFF,
    QUERY_ALN_LEN_DIFF,
    // Score.
    AS_SCORE_DIFF,
    MS_SCORE_DIFF,
    AS_SCORE_RATIO,
    MS_SCORE_RATIO,
    // Event-count diffs.
    N_SUBSTITUTION_BASES_DIFF,
    N_SUBSTITUTION_BASES_RATIO,
    N_INSERTION_BASES_DIFF,
    N_INSERTION_BASES_RATIO,
    N_DELETION_BASES_DIFF,
    N_DELETION_BASES_RATIO,
    N_SOFTCLIP_BASES_START_DIFF,
    N_SOFTCLIP_BASES_END_DIFF,
    // Query-space junction set comparison.
    N_QUERY_JUNCTIONS_SHARED,
    N_QUERY_JUNCTIONS_NOT_SHARED,
    N_QUERY_JUNCTIONS_ONLY_A,
    N_QUERY_JUNCTIONS_ONLY_B,
    QUERY_JUNCTION_POS_DISTANCE,
    // Genomic-space junction set comparison.
    N_GENOMIC_JUNCTIONS_SHARED,
    N_GENOMIC_JUNCTIONS_NOT_SHARED,
    N_GENOMIC_JUNCTIONS_ONLY_A,
    N_GENOMIC_JUNCTIONS_ONLY_B,
    // Object columns at the very end (the actual non-overlapping junctions, as
    // opposed to just their counts above).
    QUERY_JUNCTIONS_ONLY_A,
    QUERY_JUNCTIONS_ONLY_B,
    GENOMIC_JUNCTIONS_ONLY_A,
    GENOMIC_JUNCTIONS_ONLY_B,
];

/// Number of columns in the comparison table.
pub(crate) const N_COMPARE_COLS: usize =
    COMPARE_KEY_COLS.len() + 2 * COMPARE_SIDE_COLS.len() + COMPARE_DIFF_COLS.len();

/// The comparison table's full header, in output order.
pub(crate) fn compare_header() -> Vec<String> {
    let mut h: Vec<String> = COMPARE_KEY_COLS.iter().map(|c| c.name.to_string()).collect();
    for side in SIDES {
        h.extend(COMPARE_SIDE_COLS.iter().map(|&c| sided(c, side)));
    }
    h.extend(COMPARE_DIFF_COLS.iter().map(|c| c.name.to_string()));
    h
}

/// The per-side columns the read classifier (`compare_summary::classify`) and
/// `find-aln-diff`'s interval builder read. `toolkit summary` and
/// `find-aln-diff` resolve exactly these, per side, from a comparison table.
pub(crate) const CLASSIFY_SIDE_COLS: [Col; 9] = [
    TARGET_NAME,
    STRAND,
    CS,
    QUERY_START,
    QUERY_END,
    TARGET_START,
    TARGET_END,
    QUERY_JUNCTIONS,
    GENOMIC_JUNCTIONS,
];

// ── Classification columns ───────────────────────────────────────────────────

pub(crate) const QUERY_IDENTICAL_SAME_STRAND: &str = "query_identical_same_strand";
pub(crate) const QUERY_IDENTICAL_REVCOMP: &str = "query_identical_revcomp";
pub(crate) const QUERY_JUNCTIONS_IDENTICAL: &str = "query_junctions_identical";
pub(crate) const REF_SAME_POSITION_SAME_ALN: &str = "ref_same_position_same_aln";
pub(crate) const REF_SAME_POSITION_DIFF_ALN: &str = "ref_same_position_diff_aln";
pub(crate) const REF_DIFF_POSITION_SAME_ALN: &str = "ref_diff_position_same_aln";
pub(crate) const REF_DIFF_POSITION_DIFF_ALN: &str = "ref_diff_position_diff_aln";
pub(crate) const REF_SAME_POSITION_SAME_JUNCTIONS: &str = "ref_same_position_same_junctions";

/// The 8 classification booleans: columns of `find-aln-diff`'s reads tables, and
/// (under the same names) rows of the summary tables.
pub(crate) const CLASSIFICATION_COLS: [&str; 8] = [
    QUERY_IDENTICAL_SAME_STRAND,
    QUERY_IDENTICAL_REVCOMP,
    QUERY_JUNCTIONS_IDENTICAL,
    REF_SAME_POSITION_SAME_ALN,
    REF_SAME_POSITION_DIFF_ALN,
    REF_DIFF_POSITION_SAME_ALN,
    REF_DIFF_POSITION_DIFF_ALN,
    REF_SAME_POSITION_SAME_JUNCTIONS,
];

// ── Summary tables ───────────────────────────────────────────────────────────

/// Header of the two-column summary tables (`compare`, `toolkit summary`,
/// `find-aln-diff`, `query-junction-diff`).
pub(crate) const SUMMARY_HEADER: [&str; 2] = ["category", "count"];

// ── Helpers ──────────────────────────────────────────────────────────────────

/// Tab-join a list of columns into a header line (no trailing newline).
pub(crate) fn header_line(cols: &[Col]) -> String {
    cols.iter().map(|c| c.name).collect::<Vec<_>>().join("\t")
}

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

/// Index of `target` in `cols`, resolved at compile time when used in a `const`.
/// A column that is not in the list is a build failure, not a wrong index.
pub(crate) const fn index_of(cols: &[Col], target: Col) -> usize {
    let mut i = 0;
    while i < cols.len() {
        if str_eq(cols[i].name, target.name) {
            return i;
        }
        i += 1;
    }
    panic!("column is not in this table")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_unique(cols: &[Col], table: &str) {
        for (i, a) in cols.iter().enumerate() {
            for b in &cols[i + 1..] {
                assert_ne!(a.name, b.name, "{table} lists column {:?} twice", a.name);
            }
        }
    }

    #[test]
    fn no_table_lists_a_column_twice() {
        assert_unique(&ALNINFO_COLS, "alninfo");
        assert_unique(&READINFO_COLS, "readinfo");
        assert_unique(&COMPARE_SIDE_COLS, "compare per-side block");
        let h = compare_header();
        for (i, a) in h.iter().enumerate() {
            assert!(!h[i + 1..].contains(a), "comparison table lists {a:?} twice");
        }
    }

    #[test]
    fn every_compare_side_column_is_a_readinfo_column() {
        // The per-side block is read out of readinfo rows by name.
        for c in COMPARE_SIDE_COLS {
            assert!(READINFO_COLS.contains(&c), "{:?} is not a readinfo column", c.name);
        }
        for c in CLASSIFY_SIDE_COLS {
            assert!(COMPARE_SIDE_COLS.contains(&c), "{:?} is not a comparison column", c.name);
        }
    }

    #[test]
    fn compare_side_block_is_readinfo_minus_the_join_keys_and_tp_tag() {
        // Every readinfo column except query_name/query_len (the join keys) and
        // tp_tag appears once per side.
        assert_eq!(COMPARE_SIDE_COLS.len(), READINFO_COLS.len() - 3);
        assert!(!COMPARE_SIDE_COLS.contains(&TP_TAG));
    }

    #[test]
    fn index_of_finds_columns_by_name() {
        assert_eq!(index_of(&ALNINFO_COLS, QUERY_NAME), 0);
        assert_eq!(index_of(&ALNINFO_COLS, TP_TAG), 35);
        assert_eq!(index_of(&READINFO_COLS, N_ALNS_TP_S), 35);
    }

    #[test]
    fn compare_header_is_keys_then_a_then_b_then_diff() {
        let h = compare_header();
        assert_eq!(h.len(), N_COMPARE_COLS);
        assert_eq!(h.len(), 99);
        assert_eq!(h[4], sided(TARGET_NAME, SIDE_A));
        assert_eq!(h[4 + COMPARE_SIDE_COLS.len()], sided(TARGET_NAME, SIDE_B));
    }
}
