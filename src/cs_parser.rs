/// Statistics accumulated from a single cs-tag string.
#[derive(Debug, Default)]
pub struct CsStats {
    pub n_match_events:           u32,
    pub n_match_bases:            u64,
    pub n_substitution_events:    u32,
    pub n_substitution_bases:     u64,
    pub n_insertion_events:       u32,
    pub n_insertion_bases:        u64,
    pub n_deletion_events:        u32,
    pub n_deletion_bases:         u64,
    pub n_splice_junction_events: u32,
    pub n_splice_junction_bases:  u64,
    /// Query position (0-based, in aligned-query space) of every `~` op.
    /// Used to derive junction coordinates after applying the strand offset.
    pub raw_junctions: Vec<u64>,
    /// Ref-relative `(intron_start, intron_end)` pairs (0-based half-open) for each `~` op.
    /// To get absolute reference coordinates, add the alignment's `target_start`.
    /// Always reference-ascending (the cs walk traverses the ref 5' → 3').
    pub raw_genomic_junctions: Vec<(u64, u64)>,
}

/// Parse a minimap2 cs tag and return accumulated statistics.
///
/// Handles both short-form matches (`:N`) and long-form matches (`=ACGT`),
/// as well as substitutions (`*xy`), insertions (`+seq`), deletions (`-seq`),
/// and splice junctions (`~nnNNNnn`).
///
/// The cs tag is expected to be lowercase (minimap2 default), but uppercase
/// letters are also accepted and treated equivalently.
pub fn parse_cs(cs: &str) -> CsStats {
    let bytes = cs.as_bytes();
    let len   = bytes.len();
    let mut pos: usize = 0;
    let mut q_pos: u64 = 0;
    // Ref position (0-based, relative to the start of the aligned region).
    // Advances on ref-consuming ops: matches (`:N`, `=ACGT`), substitution (`*xy`),
    // deletion (`-seq`), and splice (`~NNN`). Does NOT advance on insertion (`+seq`).
    let mut r_pos: u64 = 0;
    let mut stats = CsStats::default();

    while pos < len {
        let op = bytes[pos];
        pos += 1;

        match op {
            // ── Short-form match  :N ──────────────────────────────────────
            b':' => {
                let n = parse_uint(bytes, &mut pos);
                stats.n_match_events  += 1;
                stats.n_match_bases   += n;
                q_pos                 += n;
                r_pos                 += n;
            }

            // ── Long-form match  =ACGT... ─────────────────────────────────
            b'=' => {
                let start = pos;
                while pos < len && bytes[pos].is_ascii_alphabetic() { pos += 1; }
                let n = (pos - start) as u64;
                stats.n_match_events  += 1;
                stats.n_match_bases   += n;
                q_pos                 += n;
                r_pos                 += n;
            }

            // ── Substitution  *xy ─────────────────────────────────────────
            // two chars: ref base, query base
            b'*' => {
                pos += 2; // skip both bases
                stats.n_substitution_events += 1;
                stats.n_substitution_bases  += 1;
                q_pos                       += 1;
                r_pos                       += 1;
            }

            // ── Insertion  +seq ───────────────────────────────────────────
            // Insertions advance query only (no ref bases consumed).
            b'+' => {
                let start = pos;
                while pos < len && bytes[pos].is_ascii_alphabetic() { pos += 1; }
                let n = (pos - start) as u64;
                stats.n_insertion_events += 1;
                stats.n_insertion_bases  += n;
                q_pos                    += n;
                // r_pos unchanged
            }

            // ── Deletion  -seq ────────────────────────────────────────────
            // Deletions advance ref only (no query bases consumed).
            b'-' => {
                let start = pos;
                while pos < len && bytes[pos].is_ascii_alphabetic() { pos += 1; }
                let n = (pos - start) as u64;
                stats.n_deletion_events += 1;
                stats.n_deletion_bases  += n;
                // q_pos unchanged
                r_pos                   += n;
            }

            // ── Splice junction  ~nn[intron_len]nn ────────────────────────
            // format: two alpha (motif prefix) + digits (intron len) + two alpha (suffix)
            // Splice ops advance ref only (the intron span on the reference).
            b'~' => {
                // skip 2-char motif prefix (e.g. "ct" or "gt")
                pos = pos.saturating_add(2).min(len);
                let intron_len = parse_uint(bytes, &mut pos);
                // skip 2-char motif suffix (e.g. "ac" or "ag")
                pos = pos.saturating_add(2).min(len);
                stats.n_splice_junction_events += 1;
                stats.n_splice_junction_bases  += intron_len;
                stats.raw_junctions.push(q_pos);
                // Record (r_pos, r_pos + intron_len) BEFORE advancing.
                stats.raw_genomic_junctions.push((r_pos, r_pos + intron_len));
                // q_pos unchanged; r_pos advances by intron length.
                r_pos += intron_len;
            }

            // skip unrecognised bytes (defensive)
            _ => {}
        }
    }

    stats
}

/// Parse an ASCII unsigned integer starting at `bytes[*pos]`, advancing `*pos`
/// past the digits.  Returns 0 if no digits are present.
#[inline]
fn parse_uint(bytes: &[u8], pos: &mut usize) -> u64 {
    let mut n: u64 = 0;
    while *pos < bytes.len() && bytes[*pos].is_ascii_digit() {
        n = n * 10 + (bytes[*pos] - b'0') as u64;
        *pos += 1;
    }
    n
}

// ───────────────────────────────────────────────────────────────────────────
// Tests
// ───────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn juncs(cs: &str) -> Vec<(u64, u64)> {
        parse_cs(cs).raw_genomic_junctions
    }

    #[test]
    fn empty_cs() {
        assert!(juncs("").is_empty());
    }

    #[test]
    fn no_splice_ops() {
        // Matches + substitution but no ~
        assert!(juncs(":50*ag:30").is_empty());
    }

    #[test]
    fn single_intron_after_match() {
        // 100 bases of match, then a 50 bp intron
        assert_eq!(juncs(":100~ct50ag:30"), vec![(100, 150)]);
    }

    #[test]
    fn deletion_advances_ref_before_intron() {
        // 50 match + 3 bp deletion (ref consumed) + 100 bp intron
        // r_pos at splice = 50 + 3 = 53
        assert_eq!(juncs(":50-acg~ct100ag:20"), vec![(53, 153)]);
    }

    #[test]
    fn insertion_does_not_advance_ref_before_intron() {
        // 50 match + 2 bp insertion (NO ref consumed) + 100 bp intron
        // r_pos at splice = 50 (insertion ignored on ref)
        assert_eq!(juncs(":50+ac~ct100ag:20"), vec![(50, 150)]);
    }

    #[test]
    fn substitution_advances_ref_before_intron() {
        // 50 match + 1 substitution (1 ref base) + 100 bp intron
        // r_pos at splice = 51
        assert_eq!(juncs(":50*ag~ct100ag:20"), vec![(51, 151)]);
    }

    #[test]
    fn multiple_introns_ascending() {
        // 50 match, 100 bp intron, 30 match, 75 bp intron, 20 match
        // 1st intron starts at r=50 → (50, 150)
        // After 30 match: r = 150 + 30 = 180
        // 2nd intron: (180, 255)
        assert_eq!(
            juncs(":50~ct100ag:30~gt75ac:20"),
            vec![(50, 150), (180, 255)]
        );
    }

    #[test]
    fn long_form_match_advances_ref() {
        // =ACGT (4 bases) + intron of length 50
        assert_eq!(juncs("=ACGT~ct50ag:20"), vec![(4, 54)]);
    }

    #[test]
    fn q_pos_unchanged_at_splice() {
        // After: 50 match, +5 ins, ~100 intron, 20 match
        // q_pos at splice = 50 + 5 (insertion) = 55 (not advanced by splice)
        let stats = parse_cs(":50+aaaaa~ct100ag:20");
        assert_eq!(stats.raw_junctions, vec![55]);
        assert_eq!(stats.raw_genomic_junctions, vec![(50, 150)]);
    }
}
