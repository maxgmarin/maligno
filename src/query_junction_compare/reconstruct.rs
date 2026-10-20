//! Per-read, per-side splice-junction reconstruction and cross-matching.
//!
//! See the parent module's doc comment for why this re-derives junctions
//! from each side's `cs` tag instead of trusting the comparison table's own
//! `query_junctions`/`genomic_junctions` columns.

use std::collections::HashSet;
use std::io::Write;

use anyhow::Result;

use crate::cs_parser::{oriented_motif, parse_cs};

/// One junction of one side's alignment: its query position, paired with its
/// reference coordinates and splice-site motif.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct SideJunction {
    pub query_pos: i64,
    pub target_name: String,
    pub start: u64,
    pub end: u64,
    /// Reference-forward motif letters (`gtag`, `ctac`, `nnnn`).
    pub motif: [u8; 4],
}

/// A junction's identity in reference coordinates:
/// `(target_name, intron_start, intron_end, intron_strand)`.
pub(super) type GenomicKey = (String, u64, u64, char);

/// One reconstructed splice junction on one side of one read, paired and
/// cross-matched against the other side.
pub(super) struct JunctionRecord {
    pub junction_index: u32,
    pub junction: SideJunction,
    /// The intron strand on the reference (the side's `genomic_intron_strand`).
    pub intron_strand: char,
    /// The alignment's strand.
    pub strand: char,
    pub matched_in_query: bool,
    pub matched_in_genomic: bool,
    pub other_side_aligned: bool,
}

impl JunctionRecord {
    pub(super) fn genomic_key(&self) -> GenomicKey {
        let j = &self.junction;
        (j.target_name.clone(), j.start, j.end, self.intron_strand)
    }
}

/// Reconstruct one side's final, query-ascending junction list from its `cs`
/// tag — applying the same offset/strand conversion as `record.rs`'s
/// `AlnInfo::from_paf`, but keeping each junction's query position, reference
/// coordinates and motif paired throughout via a single stable-sort
/// permutation (see module doc comment in `mod.rs`).
pub(super) fn build_side_junctions(
    cs: &str,
    strand: char,
    query_start: u64,
    query_end: u64,
    target_start: u64,
    target_name: &str,
    read_len: u64,
) -> Vec<SideJunction> {
    if cs.is_empty() {
        return Vec::new();
    }
    let stats = parse_cs(cs);
    let n = stats.raw_junctions.len();
    if n == 0 {
        return Vec::new();
    }

    // Same offset formula as record.rs:AlnInfo::from_paf.
    let offset: u64 = if strand == '+' { query_start } else { read_len.saturating_sub(query_end) };

    let mut query_vals: Vec<i64> = stats.raw_junctions.iter().map(|&q| (offset + q) as i64).collect();
    if strand == '-' {
        let ql = read_len as i64;
        query_vals = query_vals.into_iter().map(|j| ql - j).collect();
    }

    // Sort by one shared, stable permutation (derived from the query values)
    // rather than sorting `query_vals` alone — this is what keeps each
    // junction's query position, coordinates and motif paired after sorting,
    // unlike record.rs's independent `junctions.sort_unstable()`.
    let mut perm: Vec<usize> = (0..n).collect();
    perm.sort_by_key(|&i| query_vals[i]);

    perm.into_iter()
        .map(|i| {
            let (s, e) = stats.raw_genomic_junctions[i];
            SideJunction {
                query_pos: query_vals[i],
                target_name: target_name.to_string(),
                start: target_start + s,
                end: target_start + e,
                motif: stats.raw_junction_motifs[i],
            }
        })
        .collect()
}

/// Build this side's membership sets for O(1) per-element lookup by the other
/// side: query positions, and genomic keys with this side's intron strand.
pub(super) fn side_position_sets(
    junctions: &[SideJunction],
    intron_strand: char,
) -> (HashSet<i64>, HashSet<GenomicKey>) {
    let qset = junctions.iter().map(|j| j.query_pos).collect();
    let gset = junctions
        .iter()
        .map(|j| (j.target_name.clone(), j.start, j.end, intron_strand))
        .collect();
    (qset, gset)
}

pub(super) fn build_junction_records(
    junctions: Vec<SideJunction>,
    intron_strand: char,
    strand: char,
    other_query_set: &HashSet<i64>,
    other_genomic_set: &HashSet<GenomicKey>,
    other_side_aligned: bool,
) -> Vec<JunctionRecord> {
    junctions
        .into_iter()
        .enumerate()
        .map(|(i, j)| {
            let key = (j.target_name.clone(), j.start, j.end, intron_strand);
            JunctionRecord {
                junction_index: (i + 1) as u32,
                matched_in_query: other_side_aligned && other_query_set.contains(&j.query_pos),
                matched_in_genomic: other_side_aligned && other_genomic_set.contains(&key),
                junction: j,
                intron_strand,
                strand,
                other_side_aligned,
            }
        })
        .collect()
}

/// One row of the per-read table, in `QJC_PER_JUNCTION_COLS` order.
pub(super) fn write_per_read_row(w: &mut dyn Write, read_name: &str, side: char, rec: &JunctionRecord) -> Result<()> {
    let j = &rec.junction;
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        read_name,
        side,
        rec.junction_index,
        j.query_pos,
        j.target_name,
        j.start,
        j.end,
        rec.intron_strand,
        oriented_motif(j.motif, rec.intron_strand),
        rec.strand,
        rec.matched_in_query as u8,
        rec.matched_in_genomic as u8,
        rec.other_side_aligned as u8,
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Two introns: GT-AG at ref [10, 110) and an `nn` one at [160, 260), relative
    // to the alignment start; with target_start 100: [110, 210) and [260, 360).
    const CS: &str = ":10~gt100ag:50~nn100nn:20";

    #[test]
    fn plus_strand_junctions_keep_query_position_coordinates_and_motif_paired() {
        let j = build_side_junctions(CS, '+', 5, 85, 100, "chr1", 100);
        assert_eq!(j.len(), 2);
        assert_eq!((j[0].query_pos, j[0].start, j[0].end, j[0].motif), (15, 110, 210, *b"gtag"));
        assert_eq!((j[1].query_pos, j[1].start, j[1].end, j[1].motif), (65, 260, 360, *b"nnnn"));
    }

    #[test]
    fn minus_strand_junctions_are_sorted_by_query_position_with_their_coordinates_and_motif() {
        // On `-`, the reference-last intron comes first along the read.
        let j = build_side_junctions(CS, '-', 5, 85, 100, "chr1", 100);
        assert_eq!(j.len(), 2);
        assert_eq!((j[0].query_pos, j[0].start, j[0].motif), (25, 260, *b"nnnn"));
        assert_eq!((j[1].query_pos, j[1].start, j[1].motif), (75, 110, *b"gtag"));
    }

    #[test]
    fn genomic_matching_needs_the_same_intron_strand() {
        let a = build_side_junctions(CS, '+', 5, 85, 100, "chr1", 100);
        let b = a.clone();
        let (qset_b, gset_b) = side_position_sets(&b, '-');
        let recs = build_junction_records(a.clone(), '+', '+', &qset_b, &gset_b, true);
        assert!(recs.iter().all(|r| r.matched_in_query && !r.matched_in_genomic));

        let (qset_b, gset_b) = side_position_sets(&b, '+');
        let recs = build_junction_records(a, '+', '+', &qset_b, &gset_b, true);
        assert!(recs.iter().all(|r| r.matched_in_query && r.matched_in_genomic));
    }

    #[test]
    fn per_read_row_writes_the_motif_in_the_junction_orientation() {
        let j = build_side_junctions(":10~ct100ac:10", '+', 0, 20, 0, "chr1", 20);
        let (q, g) = (HashSet::new(), HashSet::new());
        for (strand, expect) in [('-', "gt-ag"), ('+', "ct-ac"), ('*', "ct-ac")] {
            let recs = build_junction_records(j.clone(), strand, '+', &q, &g, false);
            let mut out = Vec::new();
            write_per_read_row(&mut out, "r1", 'A', &recs[0]).unwrap();
            let row = String::from_utf8(out).unwrap();
            let cells: Vec<&str> = row.trim_end().split('\t').collect();
            assert_eq!(cells.len(), crate::columns::QJC_PER_JUNCTION_COLS.len());
            assert_eq!((cells[7], cells[8], cells[9]), (&strand.to_string()[..], expect, "+"));
        }
    }
}
