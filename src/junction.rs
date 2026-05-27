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
