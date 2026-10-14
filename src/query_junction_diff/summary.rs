//! Summary counts for `query-junction-diff`: the read-selection funnel (every
//! row read from the input table lands in exactly one bucket), then pass 2's
//! distinct-junction counts.

use anyhow::Result;

use crate::columns::SUMMARY_HEADER;
use crate::io_utils::open_output;

use super::rollup::UniqueJunctionCounts;

#[derive(Default)]
pub(super) struct Summary {
    n_total: u64,
    n_neither_aligned: u64,
    n_both_aligned_query_junctions_identical: u64,
    n_aligned_only_no_junctions: u64,
    n_query_junctions_different_aligned_both: u64,
    n_query_junctions_different_aligned_only_a: u64,
    n_query_junctions_different_aligned_only_b: u64,
    n_unique_junctions_a: u64,
    n_unique_junctions_b: u64,
    n_unique_junctions_shared: u64,
    n_unique_junctions_only_a: u64,
    n_unique_junctions_only_b: u64,
    n_unmatched_junctions_a: u64,
    n_unmatched_junctions_b: u64,
}

impl Summary {
    /// Returns whether this row is "differing" (i.e. should get full junction
    /// reconstruction), tallying it into the right bucket either way.
    pub(super) fn observe_row(
        &mut self,
        mapped_a: bool,
        mapped_b: bool,
        n_only_a: u64,
        n_only_b: u64,
        junc_count_a: u64,
        junc_count_b: u64,
    ) -> bool {
        self.n_total += 1;
        match (mapped_a, mapped_b) {
            (false, false) => {
                self.n_neither_aligned += 1;
                false
            }
            (true, true) => {
                if n_only_a == 0 && n_only_b == 0 {
                    self.n_both_aligned_query_junctions_identical += 1;
                    false
                } else {
                    self.n_query_junctions_different_aligned_both += 1;
                    true
                }
            }
            (true, false) => {
                if junc_count_a > 0 {
                    self.n_query_junctions_different_aligned_only_a += 1;
                    true
                } else {
                    self.n_aligned_only_no_junctions += 1;
                    false
                }
            }
            (false, true) => {
                if junc_count_b > 0 {
                    self.n_query_junctions_different_aligned_only_b += 1;
                    true
                } else {
                    self.n_aligned_only_no_junctions += 1;
                    false
                }
            }
        }
    }

    /// Record pass 2's distinct-junction counts and the row counts of the two
    /// unmatched-junction tables.
    pub(super) fn set_junction_counts(&mut self, unique: &UniqueJunctionCounts, n_unmatched_a: u64, n_unmatched_b: u64) {
        self.n_unique_junctions_a = unique.n_a;
        self.n_unique_junctions_b = unique.n_b;
        self.n_unique_junctions_shared = unique.n_shared;
        self.n_unique_junctions_only_a = unique.n_only_a;
        self.n_unique_junctions_only_b = unique.n_only_b;
        self.n_unmatched_junctions_a = n_unmatched_a;
        self.n_unmatched_junctions_b = n_unmatched_b;
    }

    pub(super) fn n_query_junctions_different(&self) -> u64 {
        self.n_query_junctions_different_aligned_both
            + self.n_query_junctions_different_aligned_only_a
            + self.n_query_junctions_different_aligned_only_b
    }

    pub(super) fn n_total(&self) -> u64 {
        self.n_total
    }

    fn rows(&self) -> Vec<(String, u64)> {
        vec![
            ("n_total".to_string(), self.n_total),
            ("n_neither_aligned".to_string(), self.n_neither_aligned),
            (
                "n_both_aligned_query_junctions_identical".to_string(),
                self.n_both_aligned_query_junctions_identical,
            ),
            ("n_aligned_only_no_junctions".to_string(), self.n_aligned_only_no_junctions),
            ("n_query_junctions_different".to_string(), self.n_query_junctions_different()),
            (
                "n_query_junctions_different_aligned_both".to_string(),
                self.n_query_junctions_different_aligned_both,
            ),
            (
                "n_query_junctions_different_aligned_only_A".to_string(),
                self.n_query_junctions_different_aligned_only_a,
            ),
            (
                "n_query_junctions_different_aligned_only_B".to_string(),
                self.n_query_junctions_different_aligned_only_b,
            ),
            ("n_unique_junctions_A".to_string(), self.n_unique_junctions_a),
            ("n_unique_junctions_B".to_string(), self.n_unique_junctions_b),
            ("n_unique_junctions_shared".to_string(), self.n_unique_junctions_shared),
            ("n_unique_junctions_only_A".to_string(), self.n_unique_junctions_only_a),
            ("n_unique_junctions_only_B".to_string(), self.n_unique_junctions_only_b),
            ("n_unmatched_junctions_A".to_string(), self.n_unmatched_junctions_a),
            ("n_unmatched_junctions_B".to_string(), self.n_unmatched_junctions_b),
        ]
    }

    /// Same `category`/`count` + `label_A`/`label_B` provenance-row
    /// convention as `CompareSummary::write_tsv` (compare_summary.rs) — a
    /// dedicated writer for this dedicated struct, not a change to
    /// `CompareSummary` itself.
    pub(super) fn write_tsv(&self, path: &str, label_a: &str, label_b: &str) -> Result<()> {
        use std::io::Write;
        let mut w = open_output(Some(path))?;
        writeln!(w, "{}", SUMMARY_HEADER.join("\t"))?;
        writeln!(w, "label_A\t{label_a}")?;
        writeln!(w, "label_B\t{label_b}")?;
        for (k, v) in self.rows() {
            writeln!(w, "{k}\t{v}")?;
        }
        w.flush()?;
        Ok(())
    }
}
