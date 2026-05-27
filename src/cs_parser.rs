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
            }

            // ── Long-form match  =ACGT... ─────────────────────────────────
            b'=' => {
                let start = pos;
                while pos < len && bytes[pos].is_ascii_alphabetic() { pos += 1; }
                let n = (pos - start) as u64;
                stats.n_match_events  += 1;
                stats.n_match_bases   += n;
                q_pos                 += n;
            }

            // ── Substitution  *xy ─────────────────────────────────────────
            // two chars: ref base, query base
            b'*' => {
                pos += 2; // skip both bases
                stats.n_substitution_events += 1;
                stats.n_substitution_bases  += 1;
                q_pos                       += 1;
            }

            // ── Insertion  +seq ───────────────────────────────────────────
            b'+' => {
                let start = pos;
                while pos < len && bytes[pos].is_ascii_alphabetic() { pos += 1; }
                let n = (pos - start) as u64;
                stats.n_insertion_events += 1;
                stats.n_insertion_bases  += n;
                q_pos                    += n;
            }

            // ── Deletion  -seq ────────────────────────────────────────────
            // Deletions do NOT advance q_pos (they remove ref bases, not query)
            b'-' => {
                let start = pos;
                while pos < len && bytes[pos].is_ascii_alphabetic() { pos += 1; }
                let n = (pos - start) as u64;
                stats.n_deletion_events += 1;
                stats.n_deletion_bases  += n;
                // q_pos unchanged
            }

            // ── Splice junction  ~nn[intron_len]nn ────────────────────────
            // format: two alpha (motif prefix) + digits (intron len) + two alpha (suffix)
            // Splice junctions do NOT advance q_pos.
            b'~' => {
                // skip 2-char motif prefix (e.g. "ct" or "gt")
                pos = pos.saturating_add(2).min(len);
                let intron_len = parse_uint(bytes, &mut pos);
                // skip 2-char motif suffix (e.g. "ac" or "ag")
                pos = pos.saturating_add(2).min(len);
                stats.n_splice_junction_events += 1;
                stats.n_splice_junction_bases  += intron_len;
                stats.raw_junctions.push(q_pos);
                // q_pos unchanged
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
