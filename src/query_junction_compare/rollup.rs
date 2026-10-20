//! Per-side catalog of genomic junctions, keyed by
//! `(target_name, intron_start, intron_end, intron_strand)`.
//!
//! Pass 1 (see `mod.rs`) adds, for each differing read, the junctions the other
//! side doesn't support; pass 2 adds every junction of every read. From one
//! catalog per side come the unmatched table (the junctions with an unmatched
//! occurrence), the all-junctions table (`--emit-all-junctions`), and the
//! summary's unique / shared / only-A / only-B counts.

use std::collections::{HashMap, HashSet};
use std::io::Write;

use anyhow::Result;

use crate::columns::{header_line, QJC_JUNCTION_COLS};
use crate::cs_parser::{note_motif, oriented_motif};
use crate::io_utils::open_output;
use crate::junction::SideJunctions;

use super::reconstruct::{GenomicKey, JunctionRecord};

/// One junction's counts and motif on one side.
#[derive(Debug, Default, Clone, PartialEq)]
pub(super) struct JunctionTally {
    /// Reads in the whole table whose alignment on this side has the junction.
    pub n_reads_total: u64,
    /// Differing reads where the junction is unmatched on the other side.
    pub n_reads_unmatched: u64,
    /// First motif among the unmatched occurrences, preferring real letters.
    pub motif_unmatched: Option<[u8; 4]>,
    /// First motif among all occurrences, preferring real letters.
    pub motif_any: Option<[u8; 4]>,
}

impl JunctionTally {
    /// The reported motif: from the unmatched occurrences when there are any,
    /// otherwise from all occurrences.
    fn motif(&self) -> Option<[u8; 4]> {
        self.motif_unmatched.or(self.motif_any)
    }
}

/// Every junction seen on one side, with its [`JunctionTally`].
///
/// Two-level (target_name → `(start, end, intron_strand)` → tally) rather than
/// one map keyed on the full tuple: the name is looked up once per row instead
/// of allocated once per junction. Memory scales with the number of *distinct*
/// junctions, not reads.
#[derive(Default)]
pub(super) struct JunctionCatalog(HashMap<String, HashMap<(u64, u64, char), JunctionTally>>);

impl JunctionCatalog {
    fn tally(&mut self, target_name: &str) -> &mut HashMap<(u64, u64, char), JunctionTally> {
        if !self.0.contains_key(target_name) {
            self.0.insert(target_name.to_string(), HashMap::new());
        }
        self.0.get_mut(target_name).expect("inserted above")
    }

    /// Pass 1: one differing read's junction records on this side. Each
    /// junction unmatched on the other side counts once per read, however many
    /// times the read reports it.
    pub(super) fn add_unmatched(&mut self, records: &[JunctionRecord]) {
        let mut seen_this_read: HashSet<GenomicKey> = HashSet::new();
        for r in records {
            if r.matched_in_genomic {
                continue;
            }
            let key = r.genomic_key();
            if seen_this_read.insert(key.clone()) {
                let (target_name, start, end, strand) = key;
                let t = self.tally(&target_name).entry((start, end, strand)).or_default();
                t.n_reads_unmatched += 1;
                note_motif(&mut t.motif_unmatched, r.junction.motif);
            }
        }
    }

    /// Pass 2: every junction of one read's alignment on this side.
    pub(super) fn add_side(&mut self, side: &SideJunctions) {
        let Some((target_name, ..)) = side.junctions.first() else {
            return;
        };
        let by_pos = self.tally(target_name);
        for (&(_, start, end, strand), &motif) in side.junctions.iter().zip(&side.motifs) {
            let t = by_pos.entry((start, end, strand)).or_default();
            t.n_reads_total += 1;
            note_motif(&mut t.motif_any, motif);
        }
    }

    /// Distinct junctions on this side.
    pub(super) fn len(&self) -> u64 {
        self.0.values().map(|m| m.len() as u64).sum()
    }

    /// Distinct junctions with at least one unmatched occurrence.
    pub(super) fn n_unmatched(&self) -> u64 {
        self.0.values().flat_map(|m| m.values()).filter(|t| t.n_reads_unmatched > 0).count() as u64
    }

    /// How many of this side's distinct junctions also appear in `other`.
    fn n_shared_with(&self, other: &JunctionCatalog) -> u64 {
        let mut n = 0;
        for (target_name, by_pos) in &self.0 {
            if let Some(other_by_pos) = other.0.get(target_name) {
                n += by_pos.keys().filter(|k| other_by_pos.contains_key(k)).count() as u64;
            }
        }
        n
    }

    /// Write the junction table (`QJC_JUNCTION_COLS`), sorted by
    /// `(target_name, intron_start, intron_end, intron_strand)`: every
    /// junction, or with `unmatched_only` just those with an unmatched
    /// occurrence (the unmatched table).
    pub(super) fn write(&self, path: &str, unmatched_only: bool) -> Result<()> {
        let mut names: Vec<&String> = self.0.keys().collect();
        names.sort();
        let mut w = open_output(Some(path))?;
        writeln!(w, "{}", header_line(&QJC_JUNCTION_COLS))?;
        for target_name in names {
            let by_pos = &self.0[target_name];
            let mut keys: Vec<&(u64, u64, char)> = by_pos.keys().collect();
            keys.sort();
            for key @ (start, end, strand) in keys {
                let t = &by_pos[key];
                if unmatched_only && t.n_reads_unmatched == 0 {
                    continue;
                }
                let motif = t.motif().map_or(".".to_string(), |m| oriented_motif(m, *strand));
                writeln!(
                    w,
                    "{target_name}\t{start}\t{end}\t{strand}\t{motif}\t{}\t{}",
                    t.n_reads_unmatched, t.n_reads_total
                )?;
            }
        }
        w.flush()?;
        Ok(())
    }
}

/// Distinct-junction counts for the summary: per side, and split into shared
/// vs only-one-side across the whole table.
pub(super) struct UniqueJunctionCounts {
    pub n_a: u64,
    pub n_b: u64,
    pub n_shared: u64,
    pub n_only_a: u64,
    pub n_only_b: u64,
}

pub(super) fn unique_junction_counts(a: &JunctionCatalog, b: &JunctionCatalog) -> UniqueJunctionCounts {
    let n_a = a.len();
    let n_b = b.len();
    let n_shared = a.n_shared_with(b);
    UniqueJunctionCounts { n_a, n_b, n_shared, n_only_a: n_a - n_shared, n_only_b: n_b - n_shared }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::reconstruct::{build_junction_records, build_side_junctions};
    use crate::junction::side_genomic_junctions;

    fn side(target: &str, strand: &str, cs: &str) -> SideJunctions {
        side_genomic_junctions(target, "0", strand, cs)
    }

    /// One read's records with every junction unmatched (other side unaligned).
    fn unmatched_records(cs: &str) -> Vec<JunctionRecord> {
        let none = (HashSet::new(), HashSet::new());
        let j = build_side_junctions(cs, '+', 0, 20, 0, "chr1", 20);
        build_junction_records(j, '+', '+', &none.0, &none.1, false)
    }

    /// Write `c` and read it back. `name` keeps the temp file unique per test,
    /// since tests run in parallel.
    fn written(c: &JunctionCatalog, unmatched_only: bool, name: &str) -> String {
        let path = std::env::temp_dir()
            .join(format!("maligno_test_catalog_{}_{name}_{unmatched_only}.tsv", std::process::id()));
        let p = path.to_str().unwrap();
        c.write(p, unmatched_only).unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        let _ = std::fs::remove_file(&path);
        content
    }

    #[test]
    fn add_side_counts_reads_per_junction() {
        let mut c = JunctionCatalog::default();
        c.add_side(&side("chr1", "+", ":10~gt100ag:5~gt50ag:5"));
        c.add_side(&side("chr1", "+", ":10~gt100ag:5"));
        c.add_side(&side("*", "*", ""));
        assert_eq!(c.len(), 2);
        assert_eq!(c.0["chr1"][&(10, 110, '+')].n_reads_total, 2);
        assert_eq!(c.0["chr1"][&(115, 165, '+')].n_reads_total, 1);
    }

    #[test]
    fn junctions_are_keyed_on_target_name_and_intron_strand() {
        let mut c = JunctionCatalog::default();
        c.add_side(&side("chr1", "+", ":10~gt100ag:5"));
        c.add_side(&side("chr1", "-", ":10~ct100ac:5"));
        c.add_side(&side("chr2", "+", ":10~gt100ag:5"));
        assert_eq!(c.len(), 3);
    }

    #[test]
    fn unique_junction_counts_splits_shared_and_only() {
        let mut a = JunctionCatalog::default();
        a.add_side(&side("chr1", "+", ":100~gt150ag:150~gt400ag:5"));
        a.add_side(&side("chr2", "-", ":10~ct10ac:5"));
        let mut b = JunctionCatalog::default();
        b.add_side(&side("chr1", "+", ":100~gt150ag:5"));
        b.add_side(&side("chr1", "+", ":900~gt50ag:5"));
        b.add_side(&side("chr1", "-", ":400~ct400ac:5")); // same coords as A's 2nd, other strand
        b.add_side(&side("chr3", "+", ":5~gt1ag:5"));
        let c = unique_junction_counts(&a, &b);
        assert_eq!((c.n_a, c.n_b), (3, 4));
        assert_eq!(c.n_shared, 1);
        assert_eq!((c.n_only_a, c.n_only_b), (2, 3));
    }

    #[test]
    fn add_unmatched_counts_reads_and_keeps_a_real_motif() {
        let mut c = JunctionCatalog::default();
        c.add_unmatched(&unmatched_records(":10~nn100nn:10"));
        c.add_unmatched(&unmatched_records(":10~gt100ag:10"));
        let t = &c.0["chr1"][&(10, 110, '+')];
        assert_eq!((t.n_reads_unmatched, t.motif_unmatched), (2, Some(*b"gtag")));
        assert_eq!(c.n_unmatched(), 1);
    }

    #[test]
    fn the_motif_comes_from_unmatched_reads_when_there_are_any() {
        let mut c = JunctionCatalog::default();
        // Junction 1: unmatched only with `nn`, but real letters in another read.
        c.add_unmatched(&unmatched_records(":10~nn100nn:10"));
        c.add_side(&side("chr1", "+", ":10~nn100nn:10"));
        c.add_side(&side("chr1", "+", ":10~gt100ag:10"));
        // Junction 2: never unmatched.
        c.add_side(&side("chr1", "+", ":300~gc100ag:10"));
        let all = written(&c, false, "motif_source");
        let lines: Vec<&str> = all.lines().collect();
        assert_eq!(lines[1], "chr1\t10\t110\t+\tnn-nn\t1\t2");
        assert_eq!(lines[2], "chr1\t300\t400\t+\tgc-ag\t0\t1");
    }

    #[test]
    fn the_unmatched_table_is_the_all_junctions_table_filtered() {
        let mut c = JunctionCatalog::default();
        c.add_unmatched(&unmatched_records(":10~gt100ag:10"));
        c.add_side(&side("chr1", "+", ":10~gt100ag:10"));
        c.add_side(&side("chr1", "-", ":300~ct100ac:10"));
        c.add_side(&side("chr2", "+", ":5~gt20ag:5"));
        let all = written(&c, false, "filtered");
        let unmatched = written(&c, true, "filtered");
        let filtered: Vec<&str> = all
            .lines()
            .enumerate()
            .filter(|(i, l)| *i == 0 || !l.split('\t').nth(5).is_some_and(|n| n == "0"))
            .map(|(_, l)| l)
            .collect();
        assert_eq!(unmatched.lines().collect::<Vec<_>>(), filtered);
        assert_eq!(
            all.lines().next().unwrap(),
            "target_name\tintron_start\tintron_end\tintron_strand\tmotif\t\
             n_reads_query_junctions_different\tn_reads_with_junction_total"
        );
        assert_eq!(all.lines().count(), 4);
        assert_eq!(all.lines().nth(2).unwrap(), "chr1\t300\t400\t-\tgt-ag\t0\t1");
    }
}
