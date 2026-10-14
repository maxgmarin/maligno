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
//! Where alninfo and readinfo carry a column under the same name, both lists
//! share one constant.

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

// ── Columns shared by alninfo and readinfo ───────────────────────────────────

pub(crate) const STRAND: Col = col("Strand", Str);
pub(crate) const QUERY_START: Col = col("Query_Start", U64);
pub(crate) const QUERY_END: Col = col("Query_End", U64);
pub(crate) const TARGET_START: Col = col("Target_Start", U64);
pub(crate) const TARGET_END: Col = col("Target_End", U64);
pub(crate) const CS: Col = col("cs", Str);
pub(crate) const N_MATCH_EVENTS: Col = col("N_Match_Events", U64);
pub(crate) const N_MATCH_BASES: Col = col("N_Match_Bases", U64);
pub(crate) const N_SUBSTITUTION_EVENTS: Col = col("N_Substitution_Events", U64);
pub(crate) const N_SUBSTITUTION_BASES: Col = col("N_Substitution_Bases", U64);
pub(crate) const N_INSERTION_EVENTS: Col = col("N_Insertion_Events", U64);
pub(crate) const N_INSERTION_BASES: Col = col("N_Insertion_Bases", U64);
pub(crate) const N_DELETION_EVENTS: Col = col("N_Deletion_Events", U64);
pub(crate) const N_DELETION_BASES: Col = col("N_Deletion_Bases", U64);
pub(crate) const N_SPLICE_JUNCTION_EVENTS: Col = col("N_Splice_Junction_Events", U64);
pub(crate) const N_SPLICE_JUNCTION_BASES: Col = col("N_Splice_Junction_Bases", U64);
pub(crate) const N_SOFTCLIPPED_BASES_START: Col = col("N_SoftClipped_Bases_Start", U64);
pub(crate) const N_SOFTCLIPPED_BASES_END: Col = col("N_SoftClipped_Bases_End", U64);
pub(crate) const N_SOFTCLIPPED_EVENTS: Col = col("N_SoftClipped_Events", U64);
pub(crate) const JUNCTIONS: Col = col("junctions", Str);
pub(crate) const GENOMIC_JUNCTIONS: Col = col("genomic_junctions", Str);
pub(crate) const TP_TAG: Col = col("tp_tag", Str);

// ── alninfo-only columns ─────────────────────────────────────────────────────

pub(crate) const QUERY_NAME: Col = col("Query_Name", Str);
pub(crate) const QUERY_LEN: Col = col("Query_Len", U64);
pub(crate) const TARGET_NAME: Col = col("Target_Name", Str);
pub(crate) const TARGET_LEN: Col = col("Target_Len", U64);
pub(crate) const NUM_RESIDUE_MATCHES: Col = col("Num_Residue_Matches", U64);
pub(crate) const ALN_BLOCK_LEN: Col = col("Aln_Block_Len", U64);
pub(crate) const MQ: Col = col("MQ", U64);
pub(crate) const MS: Col = col("ms", I64);
pub(crate) const AS: Col = col("AS", I64);
pub(crate) const SPLICE_JUNCTION_COUNT: Col = col("splice_junction_count", U64);
pub(crate) const TARGET_START_1BASED: Col = col("Target_Start_1based", U64);
pub(crate) const SEQID: Col = col("seqid", F64);
pub(crate) const QUERY_ALN_LEN: Col = col("Query_Aln_Len", U64);
pub(crate) const QUERY_ALN_COV: Col = col("Query_Aln_Cov", F64);

// ── readinfo-only columns ────────────────────────────────────────────────────

pub(crate) const READ_NAME: Col = col("Read_Name", Str);
pub(crate) const READ_LEN: Col = col("Read_Len", U64);
pub(crate) const TARGET_CHR: Col = col("TargetChr", Str);
pub(crate) const MQ_BEST: Col = col("MQ_Best", U64);
pub(crate) const AS_MAX: Col = col("AS_Max", I64);
pub(crate) const MS_MAX: Col = col("ms_Max", I64);
pub(crate) const QUERY_ALN_COV_MAX: Col = col("Query_Aln_Cov_Max", F64);
pub(crate) const QUERY_ALN_LEN_MAX: Col = col("Query_Aln_Len_Max", U64);
pub(crate) const SEQID_MAX: Col = col("seqid_Max", F64);
pub(crate) const NUM_ALN: Col = col("Num_Aln", U64);
pub(crate) const NUM_ALN_MAXSCORE: Col = col("Num_Aln_MaxScore", U64);
pub(crate) const JUNC_COUNT: Col = col("JuncCount", U64);
pub(crate) const NUM_ALN_TPP: Col = col("Num_Aln_tpP", U64);
pub(crate) const NUM_ALN_TPS: Col = col("Num_Aln_tpS", U64);

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
    NUM_RESIDUE_MATCHES,
    ALN_BLOCK_LEN,
    MQ,
    MS,
    AS,
    CS,
    N_MATCH_EVENTS,
    N_MATCH_BASES,
    N_SUBSTITUTION_EVENTS,
    N_SUBSTITUTION_BASES,
    N_INSERTION_EVENTS,
    N_INSERTION_BASES,
    N_DELETION_EVENTS,
    N_DELETION_BASES,
    N_SPLICE_JUNCTION_EVENTS,
    N_SPLICE_JUNCTION_BASES,
    N_SOFTCLIPPED_BASES_START,
    N_SOFTCLIPPED_BASES_END,
    N_SOFTCLIPPED_EVENTS,
    JUNCTIONS,
    SPLICE_JUNCTION_COUNT,
    TARGET_START_1BASED,
    SEQID,
    QUERY_ALN_LEN,
    QUERY_ALN_COV,
    GENOMIC_JUNCTIONS,
    TP_TAG,
];

/// The readinfo table, in output order. `collapse_group` builds its fields in
/// this order (pinned by a test in `readinfo.rs`).
pub(crate) const READINFO_COLS: [Col; 36] = [
    READ_NAME,
    READ_LEN,
    TARGET_CHR,
    STRAND,
    MQ_BEST,
    AS_MAX,
    MS_MAX,
    QUERY_ALN_COV_MAX,
    QUERY_ALN_LEN_MAX,
    SEQID_MAX,
    JUNCTIONS,
    NUM_ALN,
    NUM_ALN_MAXSCORE,
    JUNC_COUNT,
    N_MATCH_EVENTS,
    N_MATCH_BASES,
    N_SUBSTITUTION_EVENTS,
    N_SUBSTITUTION_BASES,
    N_INSERTION_EVENTS,
    N_INSERTION_BASES,
    N_DELETION_EVENTS,
    N_DELETION_BASES,
    N_SPLICE_JUNCTION_EVENTS,
    N_SPLICE_JUNCTION_BASES,
    N_SOFTCLIPPED_BASES_START,
    N_SOFTCLIPPED_BASES_END,
    N_SOFTCLIPPED_EVENTS,
    CS,
    GENOMIC_JUNCTIONS,
    QUERY_START,
    QUERY_END,
    TARGET_START,
    TARGET_END,
    TP_TAG,
    NUM_ALN_TPP,
    NUM_ALN_TPS,
];

// ── Comparison table ─────────────────────────────────────────────────────────

pub(crate) const LABEL_A: Col = col("Label_A", Str);
pub(crate) const LABEL_B: Col = col("Label_B", Str);

/// The comparison table's leading columns: the join keys and the two set labels.
pub(crate) const COMPARE_KEY_COLS: [Col; 4] = [READ_NAME, READ_LEN, LABEL_A, LABEL_B];

/// The per-side columns, each emitted once suffixed `_A` and once `_B`, read out
/// of each side's readinfo row by name.
///
/// This order is deliberately different from [`READINFO_COLS`]. It is grouped by
/// topic (locus, alignment selection, identity/coverage, junction counts,
/// cs-derived event counts, and the three long strings last) so the comparison
/// table is readable via `head -1 | tr '\t' '\n' | nl`. Do not "resync" the two.
pub(crate) const COMPARE_SIDE_COLS: [Col; 33] = [
    // Locus & span.
    TARGET_CHR,
    STRAND,
    TARGET_START,
    TARGET_END,
    QUERY_START,
    QUERY_END,
    // Alignment selection & score.
    MQ_BEST,
    NUM_ALN,
    NUM_ALN_TPP,
    NUM_ALN_TPS,
    NUM_ALN_MAXSCORE,
    AS_MAX,
    MS_MAX,
    // Identity & coverage.
    SEQID_MAX,
    QUERY_ALN_COV_MAX,
    QUERY_ALN_LEN_MAX,
    // Junction counts.
    JUNC_COUNT,
    N_SPLICE_JUNCTION_EVENTS,
    N_SPLICE_JUNCTION_BASES,
    // cs-derived event counts.
    N_MATCH_EVENTS,
    N_MATCH_BASES,
    N_SUBSTITUTION_EVENTS,
    N_SUBSTITUTION_BASES,
    N_INSERTION_EVENTS,
    N_INSERTION_BASES,
    N_DELETION_EVENTS,
    N_DELETION_BASES,
    N_SOFTCLIPPED_EVENTS,
    N_SOFTCLIPPED_BASES_START,
    N_SOFTCLIPPED_BASES_END,
    // Long strings last — these dominate the table's on-disk size.
    JUNCTIONS,
    GENOMIC_JUNCTIONS,
    CS,
];

pub(crate) const STRAND_MATCH: Col = col("Strand_Match", Bool);
pub(crate) const SEQID_DIFF: Col = col("seqid_Diff", F64);
pub(crate) const QUERY_ALN_COV_DIFF: Col = col("QueryAlnCov_Diff", F64);
pub(crate) const QUERY_ALN_LEN_DIFF: Col = col("QueryAlnLen_Diff", I64);
pub(crate) const AS_DIFF: Col = col("AS_Diff", I64);
pub(crate) const MS_DIFF: Col = col("ms_Diff", I64);
pub(crate) const AS_RATIO: Col = col("AS_Ratio", F64);
pub(crate) const MS_RATIO: Col = col("ms_Ratio", F64);
pub(crate) const N_SUBSTITUTION_BASES_DIFF: Col = col("N_Substitution_Bases_Diff", I64);
pub(crate) const N_SUBSTITUTION_BASES_RATIO: Col = col("N_Substitution_Bases_Ratio", F64);
pub(crate) const N_INSERTION_BASES_DIFF: Col = col("N_Insertion_Bases_Diff", I64);
pub(crate) const N_INSERTION_BASES_RATIO: Col = col("N_Insertion_Bases_Ratio", F64);
pub(crate) const N_DELETION_BASES_DIFF: Col = col("N_Deletion_Bases_Diff", I64);
pub(crate) const N_DELETION_BASES_RATIO: Col = col("N_Deletion_Bases_Ratio", F64);
pub(crate) const N_SOFTCLIPPED_BASES_START_DIFF: Col = col("N_SoftClipped_Bases_Start_Diff", I64);
pub(crate) const N_SOFTCLIPPED_BASES_END_DIFF: Col = col("N_SoftClipped_Bases_End_Diff", I64);
pub(crate) const N_MATCHED_JUNCTIONS: Col = col("N_Matched_Junctions", U64);
pub(crate) const N_UNMATCHED_JUNCTIONS: Col = col("N_Unmatched_Junctions", U64);
pub(crate) const N_JUNCTIONS_ONLY_A: Col = col("N_Junctions_OnlyA", U64);
pub(crate) const N_JUNCTIONS_ONLY_B: Col = col("N_Junctions_OnlyB", U64);
pub(crate) const JUNCTION_DISTANCE: Col = col("Junction_Distance", U64);
pub(crate) const JUNC_DIST_V2: Col = col("Junc_Dist_V2", U64);
pub(crate) const GENOMIC_N_MATCHED_JUNCTIONS: Col = col("Genomic_N_Matched_Junctions", U64);
pub(crate) const GENOMIC_N_UNMATCHED_JUNCTIONS: Col = col("Genomic_N_Unmatched_Junctions", U64);
pub(crate) const GENOMIC_N_JUNCTIONS_ONLY_A: Col = col("Genomic_N_Junctions_OnlyA", U64);
pub(crate) const GENOMIC_N_JUNCTIONS_ONLY_B: Col = col("Genomic_N_Junctions_OnlyB", U64);
pub(crate) const JUNCTIONS_ONLY_A: Col = col("Junctions_OnlyA", Str);
pub(crate) const JUNCTIONS_ONLY_B: Col = col("Junctions_OnlyB", Str);
pub(crate) const GENOMIC_JUNCTIONS_ONLY_A: Col = col("Genomic_Junctions_OnlyA", Str);
pub(crate) const GENOMIC_JUNCTIONS_ONLY_B: Col = col("Genomic_Junctions_OnlyB", Str);

/// The A-vs-B comparison columns, emitted after both per-side blocks. Grouped to
/// mirror [`COMPARE_SIDE_COLS`]: orientation/identity, score, event diffs,
/// query-space junctions, genomic-space junctions, then the object lists.
///
/// `AlignmentDiff::values()` and `write_tsv_row` both emit in this order; the
/// Parquet round-trip test in `parquet_out.rs` checks them against each other.
pub(crate) const COMPARE_DIFF_COLS: [Col; 30] = [
    // Orientation & identity.
    STRAND_MATCH,
    SEQID_DIFF,
    QUERY_ALN_COV_DIFF,
    QUERY_ALN_LEN_DIFF,
    // Score.
    AS_DIFF,
    MS_DIFF,
    AS_RATIO,
    MS_RATIO,
    // Event-count diffs.
    N_SUBSTITUTION_BASES_DIFF,
    N_SUBSTITUTION_BASES_RATIO,
    N_INSERTION_BASES_DIFF,
    N_INSERTION_BASES_RATIO,
    N_DELETION_BASES_DIFF,
    N_DELETION_BASES_RATIO,
    N_SOFTCLIPPED_BASES_START_DIFF,
    N_SOFTCLIPPED_BASES_END_DIFF,
    // Query-space junction set comparison.
    N_MATCHED_JUNCTIONS,
    N_UNMATCHED_JUNCTIONS,
    N_JUNCTIONS_ONLY_A,
    N_JUNCTIONS_ONLY_B,
    JUNCTION_DISTANCE,
    JUNC_DIST_V2,
    // Genomic-space junction set comparison.
    GENOMIC_N_MATCHED_JUNCTIONS,
    GENOMIC_N_UNMATCHED_JUNCTIONS,
    GENOMIC_N_JUNCTIONS_ONLY_A,
    GENOMIC_N_JUNCTIONS_ONLY_B,
    // Object columns at the very end (the actual non-overlapping junctions, as
    // opposed to just their counts above).
    JUNCTIONS_ONLY_A,
    JUNCTIONS_ONLY_B,
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
    TARGET_CHR,
    STRAND,
    CS,
    QUERY_START,
    QUERY_END,
    TARGET_START,
    TARGET_END,
    JUNCTIONS,
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
pub(crate) const SUMMARY_HEADER: [&str; 2] = ["Category", "Count"];

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
        // Every readinfo column except Read_Name/Read_Len (the join keys) and
        // tp_tag appears once per side.
        assert_eq!(COMPARE_SIDE_COLS.len(), READINFO_COLS.len() - 3);
        assert!(!COMPARE_SIDE_COLS.contains(&TP_TAG));
    }

    #[test]
    fn index_of_finds_columns_by_name() {
        assert_eq!(index_of(&ALNINFO_COLS, QUERY_NAME), 0);
        assert_eq!(index_of(&ALNINFO_COLS, TP_TAG), 35);
        assert_eq!(index_of(&READINFO_COLS, NUM_ALN_TPS), 35);
    }

    #[test]
    fn compare_header_is_keys_then_a_then_b_then_diff() {
        let h = compare_header();
        assert_eq!(h.len(), N_COMPARE_COLS);
        assert_eq!(h.len(), 100);
        assert_eq!(h[4], sided(TARGET_CHR, SIDE_A));
        assert_eq!(h[4 + COMPARE_SIDE_COLS.len()], sided(TARGET_CHR, SIDE_B));
    }
}
