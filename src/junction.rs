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

/// Count the number of junctions in a Python-format tuple string.
/// "()" → 0, "(108,)" → 1, "(108, 168)" → 2
pub fn junction_count_str(s: &str) -> usize {
    let trimmed = s.trim().trim_start_matches('(').trim_end_matches(')');
    if trimmed.is_empty() {
        return 0;
    }
    trimmed
        .split(',')
        .filter(|tok| !tok.trim().is_empty())
        .count()
}

/// Sum of |a_i - b_i| for sorted pairs up to min(len(a), len(b)).
pub fn junction_distance(a: &[i64], b: &[i64]) -> u64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).unsigned_abs())
        .sum()
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
// Genomic-junction (chrom, start, end) tuple parsing + set comparison
// ───────────────────────────────────────────────────────────────────────────

/// Parse a Python tuple-of-tuples string of `(start, end)` integer pairs.
///
/// As of v0.2.3 the serialized form for `genomic_junctions` drops the chromosome
/// from each inner tuple — chrom is available separately via the per-row
/// `TargetChr` (alninfo) / `TargetChr_<label>` (compare) column. The caller
/// reconstructs full `(chrom, start, end)` tuples by combining this output with
/// the relevant per-row chrom.
///
/// Accepts the format written by `write_genomic_junction_tuple` in record.rs and
/// `format_genomic_junction_tuple` in this module:
///
/// ```text
/// "()"                          → []
/// "((100, 250),)"               → [(100, 250)]
/// "((100, 250), (400, 800))"    → [(100, 250), (400, 800)]
/// ```
///
/// Defensive against malformed input — returns the pairs it could parse.
pub fn parse_genomic_junction_str(s: &str) -> Vec<(u64, u64)> {
    let bytes = s.as_bytes();
    let n = bytes.len();
    let mut out: Vec<(u64, u64)> = Vec::new();
    let mut i = 0;

    while i < n {
        // Look for "(" followed by a digit — start of an inner (start, end) pair.
        if bytes[i] == b'(' && i + 1 < n && bytes[i + 1].is_ascii_digit() {
            i += 1; // past '('

            // Parse start integer.
            let start_idx = i;
            while i < n && bytes[i].is_ascii_digit() {
                i += 1;
            }
            let start: u64 = s[start_idx..i].parse().unwrap_or(0);

            // Skip "," + optional whitespace.
            while i < n && (bytes[i] == b',' || bytes[i].is_ascii_whitespace()) {
                i += 1;
            }

            // Parse end integer.
            let end_idx = i;
            while i < n && bytes[i].is_ascii_digit() {
                i += 1;
            }
            let end: u64 = s[end_idx..i].parse().unwrap_or(0);

            // Skip to closing ')'.
            while i < n && bytes[i] != b')' {
                i += 1;
            }
            if i < n {
                i += 1;
            }

            out.push((start, end));
        } else {
            i += 1;
        }
    }

    out
}

/// Set-overlap stats for genome-coordinate junctions, treating each
/// `(chrom, start, end)` triple as a set element. Deduplicated on both sides.
///
/// Cross-chromosome safety: different `chrom` strings → different set elements →
/// junctions on different contigs cannot accidentally match.
pub fn genomic_junction_set_stats(
    a: &[(String, u64, u64)],
    b: &[(String, u64, u64)],
) -> (u64, u64, u64) {
    use std::collections::HashSet;
    let set_a: HashSet<&(String, u64, u64)> = a.iter().collect();
    let set_b: HashSet<&(String, u64, u64)> = b.iter().collect();
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
/// `(only_a, only_b)` lists of genome-coordinate junction tuples, sorted by
/// `(chrom, start, end)`.
pub fn genomic_junction_set_diffs(
    a: &[(String, u64, u64)],
    b: &[(String, u64, u64)],
) -> (Vec<(String, u64, u64)>, Vec<(String, u64, u64)>) {
    use std::collections::HashSet;
    let set_a: HashSet<&(String, u64, u64)> = a.iter().collect();
    let set_b: HashSet<&(String, u64, u64)> = b.iter().collect();
    let mut only_a: Vec<(String, u64, u64)> =
        set_a.difference(&set_b).map(|t| (*t).clone()).collect();
    let mut only_b: Vec<(String, u64, u64)> =
        set_b.difference(&set_a).map(|t| (*t).clone()).collect();
    only_a.sort();
    only_b.sort();
    (only_a, only_b)
}

// ───────────────────────────────────────────────────────────────────────────
// String formatters — Python-tuple style, matching the per-side
// `junctions` / `genomic_junctions` column format produced by record.rs.
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

/// Render `&[(String, u64, u64)]` as a Python tuple-of-tuples string.
///
/// As of v0.2.3 the chrom field is **intentionally dropped on serialization** —
/// chrom is available separately via the per-row `TargetChr` sibling column, so
/// embedding it in every tuple was redundant noise. The in-memory tuple keeps
/// its `String` chrom (used by `genomic_junction_set_stats` and friends for
/// cross-chromosome safety); only the serialized form is shorter.
/// ```text
///   []                              → "()"
///   [("chr22", 100, 250)]           → "((100, 250),)"
///   [(c1, s1, e1), (c2, s2, e2)]   → "((s1, e1), (s2, e2))"
/// ```
pub fn format_genomic_junction_tuple(juncs: &[(String, u64, u64)]) -> String {
    match juncs.len() {
        0 => "()".to_string(),
        1 => {
            let (_chrom, st, en) = &juncs[0];
            format!("(({st}, {en}),)")
        }
        _ => {
            let mut s = String::from("(");
            for (i, (_chrom, st, en)) in juncs.iter().enumerate() {
                if i > 0 {
                    s.push_str(", ");
                }
                s.push_str(&format!("({st}, {en})"));
            }
            s.push(')');
            s
        }
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Tests
// ───────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── parse_genomic_junction_str ─────────────────────────────────────────

    #[test]
    fn genomic_parse_empty() {
        assert!(parse_genomic_junction_str("()").is_empty());
        assert!(parse_genomic_junction_str("").is_empty());
    }

    #[test]
    fn genomic_parse_single() {
        assert_eq!(
            parse_genomic_junction_str("((100, 250),)"),
            vec![(100, 250)]
        );
    }

    #[test]
    fn genomic_parse_multiple() {
        assert_eq!(
            parse_genomic_junction_str("((100, 250), (400, 800))"),
            vec![(100, 250), (400, 800)]
        );
    }

    #[test]
    fn genomic_parse_extra_whitespace_tolerated() {
        // Be defensive about loose whitespace between integers / commas.
        assert_eq!(
            parse_genomic_junction_str("((100,250), (400,  800))"),
            vec![(100, 250), (400, 800)]
        );
    }

    // ── genomic_junction_set_stats ─────────────────────────────────────────

    #[test]
    fn genomic_set_same_chrom_partial_overlap() {
        let a = vec![
            ("chr1".into(), 100, 200),
            ("chr1".into(), 300, 400),
        ];
        let b = vec![
            ("chr1".into(), 300, 400),
            ("chr1".into(), 500, 600),
        ];
        // overlap = {(chr1, 300, 400)}; only_a = {(chr1, 100, 200)}; only_b = {(chr1, 500, 600)}
        assert_eq!(genomic_junction_set_stats(&a, &b), (1, 1, 1));
    }

    #[test]
    fn genomic_set_cross_chrom_disjoint() {
        // Same coordinates on different chromosomes must NOT match.
        let a = vec![("chr1".into(), 100, 250)];
        let b = vec![("chr2".into(), 100, 250)];
        assert_eq!(genomic_junction_set_stats(&a, &b), (0, 1, 1));
    }

    #[test]
    fn genomic_set_identical() {
        let v = vec![
            ("chr1".into(), 100, 200),
            ("chr1".into(), 300, 400),
        ];
        assert_eq!(genomic_junction_set_stats(&v, &v), (2, 0, 0));
    }

    #[test]
    fn genomic_set_one_side_empty() {
        let a = vec![("chr1".into(), 1, 2), ("chr1".into(), 3, 4)];
        let b: Vec<(String, u64, u64)> = vec![];
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
        let a = vec![("chr1".into(), 100, 250)];
        let b = vec![("chr2".into(), 100, 250)];
        let (oa, ob) = genomic_junction_set_diffs(&a, &b);
        assert_eq!(oa, vec![("chr1".into(), 100, 250)]);
        assert_eq!(ob, vec![("chr2".into(), 100, 250)]);
    }

    #[test]
    fn genomic_diffs_partial_overlap() {
        let a = vec![
            ("chr1".into(), 100, 200),
            ("chr1".into(), 300, 400),
        ];
        let b = vec![
            ("chr1".into(), 300, 400),
            ("chr1".into(), 500, 600),
        ];
        let (oa, ob) = genomic_junction_set_diffs(&a, &b);
        assert_eq!(oa, vec![("chr1".into(), 100, 200)]);
        assert_eq!(ob, vec![("chr1".into(), 500, 600)]);
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
        let v: Vec<(String, u64, u64)> = vec![];
        assert_eq!(format_genomic_junction_tuple(&v), "()");
    }

    #[test]
    fn format_genomic_single_has_trailing_comma() {
        // Chrom field is dropped on serialization (v0.2.3+).
        let v = vec![("chr22".to_string(), 100, 250)];
        assert_eq!(format_genomic_junction_tuple(&v), "((100, 250),)");
    }

    #[test]
    fn format_genomic_multiple() {
        let v = vec![
            ("chr22".to_string(), 100, 250),
            ("chr22".to_string(), 400, 800),
        ];
        assert_eq!(
            format_genomic_junction_tuple(&v),
            "((100, 250), (400, 800))"
        );
    }

    #[test]
    fn format_genomic_chrom_value_is_ignored() {
        // Chrom is dropped on serialization regardless of its value (no escaping needed).
        let v = vec![("ab\tc'd".to_string(), 1, 2)];
        assert_eq!(format_genomic_junction_tuple(&v), "((1, 2),)");
    }

    // ── round-trip: format → parse loses chrom (by design) ─────────────────

    #[test]
    fn format_then_parse_genomic_roundtrip_pairs_only() {
        // After v0.2.3 the serialized form drops chrom, so the round-trip
        // returns just (start, end) pairs. The caller reconstructs full tuples
        // by combining the per-row TargetChr with these pairs.
        let v = vec![
            ("chr22".to_string(), 100, 250),
            ("chrX".to_string(), 999, 1500),
        ];
        let s = format_genomic_junction_tuple(&v);
        let parsed = parse_genomic_junction_str(&s);
        assert_eq!(parsed, vec![(100, 250), (999, 1500)]);
    }
}
