/// Parse a Python-format tuple string into a Vec<i64>.
/// "()" → []
/// "(108,)" → [108]
/// "(108, 168, 359)" → [108, 168, 359]
pub fn parse_junction_str(s: &str) -> Vec<i64> {
    let trimmed = s.trim().trim_start_matches('(').trim_end_matches(')');
    if trimmed.is_empty() {
        return Vec::new();
    }
    trimmed
        .split(',')
        .filter_map(|tok| {
            let t = tok.trim();
            if t.is_empty() {
                None
            } else {
                t.parse::<i64>().ok()
            }
        })
        .collect()
}

/// Set-overlap stats between two junction lists (deduplicated on both sides).
/// Returns `(n_overlap, n_only_a, n_only_b)` where:
///   n_overlap = |A ∩ B|, n_only_a = |A \ B|, n_only_b = |B \ A|.
///
/// Because both sides are deduplicated, the parts stay internally consistent:
///   n_overlap + n_only_a == |set(A)|  and  n_overlap + n_only_b == |set(B)|.
pub fn junction_set_stats(a: &[i64], b: &[i64]) -> (u64, u64, u64) {
    use std::collections::HashSet;
    let set_a: HashSet<i64> = a.iter().copied().collect();
    let set_b: HashSet<i64> = b.iter().copied().collect();
    let overlap = set_a.intersection(&set_b).count() as u64;
    let only_a = set_a.len() as u64 - overlap;
    let only_b = set_b.len() as u64 - overlap;
    (overlap, only_a, only_b)
}

// ───────────────────────────────────────────────────────────────────────────
// Genomic junctions: extraction from `cs` + set comparison
// ───────────────────────────────────────────────────────────────────────────

/// A genomic junction as `toolkit genomic-junction-compare` identifies it:
/// `(contig, intron_start, intron_end, intron_strand)`, 0-based half-open, with
/// `intron_strand` `+`, `-` or `*` (unknown; its own value, matching only `*`).
pub type GenomicJunction = (String, u64, u64, char);

/// A `genomic_intron_strand` value as `+`, `-` or `*` (unknown).
pub(crate) fn normalize_intron_strand(s: &str) -> char {
    match s.chars().next() {
        Some(c @ ('+' | '-')) => c,
        _ => '*',
    }
}

/// One side's genomic junctions, rebuilt from its `cs`. Shared by
/// `genomic-junction-compare` and `query-junction-compare`.
#[derive(Debug, Default, PartialEq)]
pub(crate) struct SideJunctions {
    pub intron_strand: char,
    /// In `cs` order, i.e. reference-ascending. An alignment's introns never
    /// repeat (its reference position only moves forward), so these are distinct.
    pub junctions: Vec<GenomicJunction>,
    /// Reference-forward motif letters, parallel to `junctions`.
    pub motifs: Vec<[u8; 4]>,
}

/// One side's junctions from its comparison-table columns: every `~` op of
/// `cs`, placed at `target_start`, with the side's intron strand. An unmapped
/// side (`target_name` empty or `*`) has none.
pub(crate) fn side_genomic_junctions(
    target_name: &str,
    target_start: &str,
    genomic_intron_strand: &str,
    cs: &str,
) -> SideJunctions {
    if target_name.is_empty() || target_name == "*" {
        return SideJunctions { intron_strand: '*', ..Default::default() };
    }
    let gs = normalize_intron_strand(genomic_intron_strand);
    let tstart: u64 = target_start.parse().unwrap_or(0);
    let stats = crate::cs_parser::parse_cs(cs);
    let junctions = stats
        .raw_genomic_junctions
        .iter()
        .map(|&(s, e)| (target_name.to_string(), tstart + s, tstart + e, gs))
        .collect();
    SideJunctions { intron_strand: gs, junctions, motifs: stats.raw_junction_motifs }
}

/// Set-overlap stats for genome-coordinate junctions, treating each
/// [`GenomicJunction`] as a set element. Deduplicated on both sides.
///
/// Different contigs, or different strands, are different set elements, so
/// they never match.
pub fn genomic_junction_set_stats(
    a: &[GenomicJunction],
    b: &[GenomicJunction],
) -> (u64, u64, u64) {
    use std::collections::HashSet;
    let set_a: HashSet<&GenomicJunction> = a.iter().collect();
    let set_b: HashSet<&GenomicJunction> = b.iter().collect();
    let overlap = set_a.intersection(&set_b).count() as u64;
    let only_a = set_a.len() as u64 - overlap;
    let only_b = set_b.len() as u64 - overlap;
    (overlap, only_a, only_b)
}

// ───────────────────────────────────────────────────────────────────────────
// Set-difference list helpers (return the actual non-overlapping junctions,
// not just the counts). Both sides are deduplicated, and the resulting lists
// are sorted for deterministic output.
// ───────────────────────────────────────────────────────────────────────────

/// Symmetric counterpart to `junction_set_stats`: return the actual
/// `(only_a, only_b)` lists of query-coordinate junctions, sorted ascending.
///
/// `only_a` contains every value in A that is not in B (and vice versa). Both
/// returned vecs are deduplicated, so their lengths match `n_only_a` /
/// `n_only_b` from `junction_set_stats`.
pub fn junction_set_diffs(a: &[i64], b: &[i64]) -> (Vec<i64>, Vec<i64>) {
    use std::collections::HashSet;
    let set_a: HashSet<i64> = a.iter().copied().collect();
    let set_b: HashSet<i64> = b.iter().copied().collect();
    let mut only_a: Vec<i64> = set_a.difference(&set_b).copied().collect();
    let mut only_b: Vec<i64> = set_b.difference(&set_a).copied().collect();
    only_a.sort();
    only_b.sort();
    (only_a, only_b)
}

/// Symmetric counterpart to `genomic_junction_set_stats`: return the actual
/// `(only_a, only_b)` lists of genome-coordinate junctions, sorted by
/// `(contig, start, end, strand)`.
pub fn genomic_junction_set_diffs(
    a: &[GenomicJunction],
    b: &[GenomicJunction],
) -> (Vec<GenomicJunction>, Vec<GenomicJunction>) {
    use std::collections::HashSet;
    let set_a: HashSet<&GenomicJunction> = a.iter().collect();
    let set_b: HashSet<&GenomicJunction> = b.iter().collect();
    let mut only_a: Vec<GenomicJunction> =
        set_a.difference(&set_b).map(|t| (*t).clone()).collect();
    let mut only_b: Vec<GenomicJunction> =
        set_b.difference(&set_a).map(|t| (*t).clone()).collect();
    only_a.sort();
    only_b.sort();
    (only_a, only_b)
}

// ───────────────────────────────────────────────────────────────────────────
// String formatters — Python-tuple style, matching the per-side
// `query_junctions` / `genomic_junctions` column format produced by record.rs.
// ───────────────────────────────────────────────────────────────────────────

/// Render `&[i64]` as a Python tuple-of-ints string:
/// ```text
///   []          → "()"
///   [108]       → "(108,)"          (trailing comma — Python 1-tuple form)
///   [108, 168]  → "(108, 168)"
/// ```
pub fn format_junction_tuple(juncs: &[i64]) -> String {
    match juncs.len() {
        0 => "()".to_string(),
        1 => format!("({},)", juncs[0]),
        _ => {
            let mut s = String::from("(");
            for (i, j) in juncs.iter().enumerate() {
                if i > 0 {
                    s.push_str(", ");
                }
                s.push_str(&j.to_string());
            }
            s.push(')');
            s
        }
    }
}

/// Render `&[GenomicJunction]` as a Python tuple of `(start, end)` pairs — the
/// same format as the `genomic_junctions` columns of alninfo, readinfo and the
/// comparison table (see `record.rs`):
/// ```text
///   []                                  → "()"
///   [("chr22", 100, 250, '+')]          → "((100, 250),)"
///   [(c, s1, e1, x), (c, s2, e2, x)]    → "((s1, e1), (s2, e2))"
/// ```
/// The contig and strand are dropped: one alignment's junctions share them, and
/// `toolkit genomic-junction-compare` writes them once per row (`target_name_*`,
/// `genomic_intron_strand_*`). For a self-contained form that keeps them, see
/// [`format_genomic_junction_full_tuple`].
pub fn format_genomic_junction_tuple(juncs: &[GenomicJunction]) -> String {
    let items: Vec<String> = juncs.iter().map(|(_contig, st, en, _strand)| format!("({st}, {en})")).collect();
    match items.len() {
        0 => "()".to_string(),
        1 => format!("({},)", items[0]),
        _ => format!("({})", items.join(", ")),
    }
}

/// Render `&[GenomicJunction]` as a Python tuple of **complete** junction tuples,
/// `(contig, start, end, strand)`, parseable with `ast.literal_eval`:
/// ```text
///   []                                       → "()"
///   [("chr22", 100, 250, '+')]               → "(('chr22', 100, 250, '+'),)"
///   [(c, s1, e1, '-'), (c, s2, e2, '-')]     → "(('c', s1, e1, '-'), ('c', s2, e2, '-'))"
/// ```
/// Each junction then carries its full identity, so lists from different rows
/// can be pooled or compared without the row's `target_name_*` /
/// `genomic_intron_strand_*` columns. Contig names are escaped as Python string
/// literals ([`py_str`]).
///
/// Not used by any output. Kept, with its tests, so `genomic-junction-compare`
/// can write self-contained junction lists by calling this instead of
/// [`format_genomic_junction_tuple`].
#[allow(dead_code)]
pub fn format_genomic_junction_full_tuple(juncs: &[GenomicJunction]) -> String {
    let items: Vec<String> = juncs
        .iter()
        .map(|(contig, st, en, strand)| {
            format!("({}, {st}, {en}, {})", py_str(contig), py_str(&strand.to_string()))
        })
        .collect();
    match items.len() {
        0 => "()".to_string(),
        1 => format!("({},)", items[0]),
        _ => format!("({})", items.join(", ")),
    }
}

/// A Python single-quoted string literal, escaping what would break the
/// literal or a TSV cell (backslash, quote, tab, newline, carriage return).
#[allow(dead_code)]
fn py_str(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('\'');
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\'' => out.push_str("\\'"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            _ => out.push(c),
        }
    }
    out.push('\'');
    out
}

// ───────────────────────────────────────────────────────────────────────────
// Tests
// ───────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intron_strand_is_normalized() {
        assert_eq!(normalize_intron_strand("+"), '+');
        assert_eq!(normalize_intron_strand("-"), '-');
        assert_eq!(normalize_intron_strand("*"), '*');
        assert_eq!(normalize_intron_strand(""), '*');
    }

    #[test]
    fn side_junctions_are_placed_on_the_reference_with_strand_and_motif() {
        let s = side_genomic_junctions("chr1", "1000", "+", ":10~gt100ag:5~nn50nn:5");
        assert_eq!(s.intron_strand, '+');
        assert_eq!(
            s.junctions,
            vec![("chr1".to_string(), 1010, 1110, '+'), ("chr1".to_string(), 1115, 1165, '+')]
        );
        assert_eq!(s.motifs, vec![*b"gtag", *b"nnnn"]);
        let unmapped = side_genomic_junctions("*", "0", "+", ":5~gt20ag:5");
        assert_eq!(unmapped, SideJunctions { intron_strand: '*', ..Default::default() });
    }

    // ── genomic_junction_set_stats ─────────────────────────────────────────

    #[test]
    fn genomic_set_same_chrom_partial_overlap() {
        let a = vec![
            ("chr1".into(), 100, 200, '+'),
            ("chr1".into(), 300, 400, '+'),
        ];
        let b = vec![
            ("chr1".into(), 300, 400, '+'),
            ("chr1".into(), 500, 600, '+'),
        ];
        // overlap = {(chr1, 300, 400)}; only_a = {(chr1, 100, 200)}; only_b = {(chr1, 500, 600)}
        assert_eq!(genomic_junction_set_stats(&a, &b), (1, 1, 1));
    }

    #[test]
    fn genomic_set_cross_chrom_disjoint() {
        // Same coordinates on different chromosomes must NOT match.
        let a = vec![("chr1".into(), 100, 250, '+')];
        let b = vec![("chr2".into(), 100, 250, '+')];
        assert_eq!(genomic_junction_set_stats(&a, &b), (0, 1, 1));
    }

    #[test]
    fn genomic_set_strand_is_part_of_the_junction() {
        // Same coordinates on opposite strands are different junctions, and an
        // unknown strand ('*') matches only another '*'.
        let plus = vec![("chr1".into(), 100, 250, '+')];
        let minus = vec![("chr1".into(), 100, 250, '-')];
        let unknown = vec![("chr1".into(), 100, 250, '*')];
        assert_eq!(genomic_junction_set_stats(&plus, &minus), (0, 1, 1));
        assert_eq!(genomic_junction_set_stats(&plus, &unknown), (0, 1, 1));
        assert_eq!(genomic_junction_set_stats(&unknown, &unknown), (1, 0, 0));
    }

    #[test]
    fn genomic_set_identical() {
        let v = vec![
            ("chr1".into(), 100, 200, '+'),
            ("chr1".into(), 300, 400, '+'),
        ];
        assert_eq!(genomic_junction_set_stats(&v, &v), (2, 0, 0));
    }

    #[test]
    fn genomic_set_one_side_empty() {
        let a = vec![("chr1".into(), 1, 2, '+'), ("chr1".into(), 3, 4, '+')];
        let b: Vec<GenomicJunction> = vec![];
        assert_eq!(genomic_junction_set_stats(&a, &b), (0, 2, 0));
        assert_eq!(genomic_junction_set_stats(&b, &a), (0, 0, 2));
    }

    // ── junction_set_diffs ─────────────────────────────────────────────────

    #[test]
    fn junction_diffs_partial_overlap() {
        // overlap {30, 50}; only_a {10}; only_b {70}
        let (oa, ob) = junction_set_diffs(&[10, 30, 50], &[30, 50, 70]);
        assert_eq!(oa, vec![10]);
        assert_eq!(ob, vec![70]);
    }

    #[test]
    fn junction_diffs_identical() {
        let (oa, ob) = junction_set_diffs(&[1, 2, 3], &[1, 2, 3]);
        assert!(oa.is_empty());
        assert!(ob.is_empty());
    }

    #[test]
    fn junction_diffs_disjoint_and_sorted() {
        // Unsorted inputs → sorted outputs.
        let (oa, ob) = junction_set_diffs(&[50, 10, 30], &[90, 70]);
        assert_eq!(oa, vec![10, 30, 50]);
        assert_eq!(ob, vec![70, 90]);
    }

    #[test]
    fn junction_diffs_counts_match_set_stats() {
        let a = [10, 30, 50, 70];
        let b = [30, 50];
        let (overlap, n_only_a, n_only_b) = junction_set_stats(&a, &b);
        let (oa, ob) = junction_set_diffs(&a, &b);
        assert_eq!(oa.len() as u64, n_only_a);
        assert_eq!(ob.len() as u64, n_only_b);
        assert_eq!(overlap, 2);
    }

    // ── genomic_junction_set_diffs ─────────────────────────────────────────

    #[test]
    fn genomic_diffs_cross_chrom_disjoint() {
        // Same coordinates on different chromosomes must NOT match.
        let a = vec![("chr1".into(), 100, 250, '+')];
        let b = vec![("chr2".into(), 100, 250, '+')];
        let (oa, ob) = genomic_junction_set_diffs(&a, &b);
        assert_eq!(oa, vec![("chr1".into(), 100, 250, '+')]);
        assert_eq!(ob, vec![("chr2".into(), 100, 250, '+')]);
    }

    #[test]
    fn genomic_diffs_partial_overlap() {
        let a = vec![
            ("chr1".into(), 100, 200, '+'),
            ("chr1".into(), 300, 400, '+'),
        ];
        let b = vec![
            ("chr1".into(), 300, 400, '+'),
            ("chr1".into(), 500, 600, '+'),
        ];
        let (oa, ob) = genomic_junction_set_diffs(&a, &b);
        assert_eq!(oa, vec![("chr1".into(), 100, 200, '+')]);
        assert_eq!(ob, vec![("chr1".into(), 500, 600, '+')]);
    }

    // ── format_junction_tuple ──────────────────────────────────────────────

    #[test]
    fn format_junction_empty() {
        assert_eq!(format_junction_tuple(&[]), "()");
    }

    #[test]
    fn format_junction_single_has_trailing_comma() {
        assert_eq!(format_junction_tuple(&[108]), "(108,)");
    }

    #[test]
    fn format_junction_multiple() {
        assert_eq!(format_junction_tuple(&[108, 168, 359]), "(108, 168, 359)");
    }

    // ── format_genomic_junction_tuple ──────────────────────────────────────

    #[test]
    fn format_genomic_empty() {
        let v: Vec<GenomicJunction> = vec![];
        assert_eq!(format_genomic_junction_tuple(&v), "()");
    }

    #[test]
    fn format_genomic_single_has_trailing_comma() {
        // Contig and strand are dropped: the pairs format of the genomic_junctions columns.
        let v = vec![("chr22".to_string(), 100, 250, '+')];
        assert_eq!(format_genomic_junction_tuple(&v), "((100, 250),)");
    }

    #[test]
    fn format_genomic_pairs_drop_contig_and_strand() {
        let v = vec![("chr22".to_string(), 100, 250, '-'), ("chr22".to_string(), 400, 800, '-')];
        assert_eq!(format_genomic_junction_tuple(&v), "((100, 250), (400, 800))");
    }

    #[test]
    fn format_genomic_full_empty_and_single() {
        let empty: Vec<GenomicJunction> = vec![];
        assert_eq!(format_genomic_junction_full_tuple(&empty), "()");
        let v = vec![("chr22".to_string(), 100, 250, '+')];
        assert_eq!(format_genomic_junction_full_tuple(&v), "(('chr22', 100, 250, '+'),)");
    }

    #[test]
    fn format_genomic_full_multiple() {
        let v = vec![
            ("chr22".to_string(), 100, 250, '-'),
            ("chr22".to_string(), 400, 800, '*'),
        ];
        assert_eq!(
            format_genomic_junction_full_tuple(&v),
            "(('chr22', 100, 250, '-'), ('chr22', 400, 800, '*'))"
        );
    }

    #[test]
    fn format_genomic_contig_is_escaped_as_a_python_string() {
        // A quote or tab in a contig name must neither end the literal nor add a
        // TSV field.
        let v = vec![("ab\tc'd".to_string(), 1, 2, '-')];
        assert_eq!(format_genomic_junction_full_tuple(&v), "(('ab\\tc\\'d', 1, 2, '-'),)");
    }

}
