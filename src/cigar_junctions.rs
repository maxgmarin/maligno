//! CIGAR-based intron (splice junction) extractor.
//!
//! Standalone utility — **not yet wired into the pipeline.** The pipeline currently
//! derives genomic junctions from the cs tag (see [`crate::cs_parser`]). This module
//! provides the parallel CIGAR-walk implementation for future use (e.g. SAM/BAM
//! input paths that don't carry a cs tag).
//!
//! Output is identical in shape to the cs-tag version: a `Vec<(u64, u64)>` of
//! ref-relative `(intron_start, intron_end)` pairs, 0-based half-open. To get
//! absolute reference coordinates, the caller adds the alignment's `target_start`.
//!
//! Reference-consuming CIGAR ops: `M`, `=`, `X`, `D`, `N`.
//! Non-ref-consuming ops:         `I`, `S`, `H`, `P`.
//!
//! Matches intron-prospector's CIGAR walk (`junctions_extractor.cc`) — coordinates
//! are always reference-ascending because the walk strictly advances `r_pos`.

/// Extract intron coordinates from a SAM/PAF-style CIGAR string.
///
/// Returns ref-relative `(intron_start, intron_end)` pairs (0-based half-open),
/// in reference-ascending order. Add `target_start` for absolute coordinates.
///
/// Defensive against malformed input: unknown op bytes are silently skipped;
/// a trailing digit run with no op character is dropped.
pub fn genomic_junctions_from_cigar(cigar: &str) -> Vec<(u64, u64)> {
    let bytes = cigar.as_bytes();
    let len = bytes.len();
    let mut pos: usize = 0;
    let mut r_pos: u64 = 0;
    let mut out: Vec<(u64, u64)> = Vec::new();

    while pos < len {
        // Parse the integer run length.
        let start = pos;
        while pos < len && bytes[pos].is_ascii_digit() {
            pos += 1;
        }
        if pos == start || pos >= len {
            // No digits found, or digits with no following op — malformed; stop.
            break;
        }
        // SAFETY: the slice contains only ASCII digits → valid UTF-8.
        let n: u64 = unsafe { std::str::from_utf8_unchecked(&bytes[start..pos]) }
            .parse()
            .unwrap_or(0);
        let op = bytes[pos];
        pos += 1;

        match op {
            b'N' => {
                // Skipped region on reference = intron.
                out.push((r_pos, r_pos + n));
                r_pos += n;
            }
            // Ref-consuming alignment ops.
            b'M' | b'=' | b'X' | b'D' => {
                r_pos += n;
            }
            // Non-ref-consuming ops: I (insertion), S (soft clip), H (hard clip), P (padding).
            b'I' | b'S' | b'H' | b'P' => {
                // r_pos unchanged
            }
            // Unknown op — skip silently.
            _ => {}
        }
    }

    out
}

// ───────────────────────────────────────────────────────────────────────────
// Tests
// ───────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_cigar() {
        assert!(genomic_junctions_from_cigar("").is_empty());
    }

    #[test]
    fn no_n_op() {
        assert!(genomic_junctions_from_cigar("100M").is_empty());
        assert!(genomic_junctions_from_cigar("50M5I45M").is_empty());
    }

    #[test]
    fn single_intron() {
        assert_eq!(
            genomic_junctions_from_cigar("100M50N100M"),
            vec![(100, 150)]
        );
    }

    #[test]
    fn two_introns_ascending() {
        // 100 match → r=100
        // 50 N intron → (100, 150), r=150
        // 100 match → r=250
        // 200 N intron → (250, 450), r=450
        // 50 match
        assert_eq!(
            genomic_junctions_from_cigar("100M50N100M200N50M"),
            vec![(100, 150), (250, 450)]
        );
    }

    #[test]
    fn deletion_before_intron_advances_ref() {
        // 100M + 5D = r at 105 before the 45N intron → (105, 150)
        assert_eq!(
            genomic_junctions_from_cigar("100M5D45N50M"),
            vec![(105, 150)]
        );
    }

    #[test]
    fn insertion_before_intron_does_not_advance_ref() {
        // 50M + 10I = r still at 50 before the 50N intron → (50, 100)
        assert_eq!(
            genomic_junctions_from_cigar("50M10I50N50M"),
            vec![(50, 100)]
        );
    }

    #[test]
    fn soft_clips_do_not_advance_ref() {
        // 10S is clipped query, doesn't touch ref
        assert_eq!(
            genomic_junctions_from_cigar("10S100M50N100M10S"),
            vec![(100, 150)]
        );
    }

    #[test]
    fn hard_clips_do_not_advance_ref() {
        assert_eq!(
            genomic_junctions_from_cigar("10H100M50N100M10H"),
            vec![(100, 150)]
        );
    }

    #[test]
    fn extended_cigar_equals_and_x() {
        // 50= + 50N + 50X → (50, 100); X is a ref-consuming op (mismatch)
        assert_eq!(
            genomic_junctions_from_cigar("50=50N50X"),
            vec![(50, 100)]
        );
    }

    #[test]
    fn padding_does_not_advance_ref() {
        // P (padding) is not a ref-consuming op
        assert_eq!(
            genomic_junctions_from_cigar("50M10P50N50M"),
            vec![(50, 100)]
        );
    }

    #[test]
    fn malformed_trailing_digits_dropped() {
        // 50M followed by digits with no op — drop the trailing garbage gracefully.
        assert_eq!(
            genomic_junctions_from_cigar("50M50N50M99"),
            vec![(50, 100)]
        );
    }
}
