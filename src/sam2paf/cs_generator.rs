/// Generate a cs tag string from a SAM MD string, the CIGAR merged-op list,
/// and the query SEQ.
///
/// This is a direct port of the MD→cs logic in `paf_sam2paf` from
/// `paftools.MGM.js`, using the same two-pointer walk:
///
///   - The MD string drives the outer loop (matches, mismatches, deletions).
///   - The merged CIGAR list drives the inner loop to interleave insertions and
///     soft-clips that are invisible to MD but must appear in cs.
///
/// `long_cs = true` emits `=ACGT` style; `false` emits `:N` style for matches.
///
/// Returns `Err` on an inconsistency between the MD and CIGAR.
pub fn generate_cs(
    merged_cigar: &[(u32, u8)],
    md: &str,
    seq: &str,
    long_cs: bool,
    lineno: u64,
) -> anyhow::Result<String> {
    use super::md::{MdIter, MdToken};

    let seq_bytes = seq.as_bytes();
    let mut cs = String::with_capacity(md.len() * 2);

    let mut k: usize = 0;  // index into merged_cigar
    let mut cy: usize = 0; // query position tracked by the cigar walk
    let mut my: usize = 0; // query position tracked by the MD walk

    for token in MdIter::new(md) {
        match token {
            MdToken::Deletion(bases) => {
                // A deletion in MD corresponds to a 'D' op in CIGAR.
                // Emit the deleted reference bases and advance past the D op.
                // Lowercased to match minimap2's cs convention (deletion bases
                // are lowercase in BOTH short and long forms; only `=ACGT`
                // matches are uppercase in long mode).
                cs.push('-');
                push_lower(&mut cs, bases);
                // The D op in merged_cigar sits at index k; skip it.
                k += 1;
            }

            token @ (MdToken::Match(_) | MdToken::Mismatch(_)) => {
                // Determine whether this is a mismatch (ml=1) or a run of
                // matches (ml=count).  For mismatches we need the ref base.
                let (is_mm, ref_base, mut ml) = match token {
                    MdToken::Mismatch(b) => (true,  b,    1u32),
                    MdToken::Match(n)    => (false, 0u8,  n),
                    _ => unreachable!(),
                };

                // Walk the merged CIGAR until we have consumed `ml` query
                // bases from M blocks, interleaving + and S ops as we go.
                loop {
                    if k >= merged_cigar.len() {
                        if ml != 0 {
                            return Err(anyhow::anyhow!(
                                "line {lineno}: MD tag is inconsistent with CIGAR \
                                 (ml={ml} remains but CIGAR is exhausted)"
                            ));
                        }
                        break;
                    }

                    let (cl, op) = merged_cigar[k];
                    // Stop at the next D; it belongs to a future MD deletion token.
                    if op == b'D' { break; }

                    let cy_end = cy + cl as usize;

                    match op {
                        b'I' => {
                            // Insertion: emit +bases, advance query, next op.
                            // Inserted bases are lowercased to match minimap2's
                            // cs convention (insertion bases lowercase in BOTH
                            // short and long forms).
                            cs.push('+');
                            push_lower(&mut cs, &seq[cy..cy_end]);
                            cy = cy_end;
                            my = cy_end;
                            k += 1;
                        }
                        b'S' => {
                            // Soft-clip: advance query but do NOT emit anything.
                            cy = cy_end;
                            my = cy_end;
                            k += 1;
                        }
                        b'M' => {
                            if my + (ml as usize) < cy_end {
                                // ml fits strictly inside this M block.
                                if ml > 0 {
                                    if is_mm {
                                        // *<ref_base><query_base>. ref_base is
                                        // already lowercased by MdIter::Mismatch;
                                        // lowercase the query base too so the
                                        // emission matches minimap2's `*xy`
                                        // all-lowercase convention (true even in
                                        // long-cs mode — only `=ACGT` is upper).
                                        cs.push('*');
                                        cs.push(ref_base as char);
                                        cs.push((seq_bytes[my] as char).to_ascii_lowercase());
                                    } else if long_cs {
                                        cs.push('=');
                                        cs.push_str(&seq[my..my + ml as usize]);
                                    } else {
                                        push_match_short(&mut cs, ml);
                                    }
                                }
                                my += ml as usize;
                                ml = 0;
                                break; // stay on this M op for the next MD token
                            } else {
                                // ml reaches past (or exactly to) the end of
                                // this M block; consume the rest of the block.
                                let dl = cy_end - my;
                                if dl > 0 {
                                    if is_mm {
                                        // Pre-existing latent bug fix: if the MD
                                        // token is a Mismatch and it lands at the
                                        // very end of an M block (so `my + 1 == cy_end`,
                                        // putting us in the `else` branch), the
                                        // original code emitted a match `:1` / `=X`
                                        // instead of the substitution `*xy`. By
                                        // construction, Mismatch implies ml == 1
                                        // and therefore dl == 1, so emit the cs sub.
                                        cs.push('*');
                                        cs.push(ref_base as char);
                                        cs.push((seq_bytes[my] as char).to_ascii_lowercase());
                                    } else if long_cs {
                                        cs.push('=');
                                        cs.push_str(&seq[my..my + dl]);
                                    } else {
                                        push_match_short(&mut cs, dl as u32);
                                    }
                                }
                                cy = cy_end;
                                my += dl;
                                ml -= dl as u32;
                                k += 1;
                                // ml may now be 0; if so, the next M-op check
                                // will find `my < cy_end` trivially true and
                                // push nothing, then break.
                            }
                        }
                        b'N' => {
                            // Splice junction / reference skip (intron).
                            // N consumes the reference only — neither the query nor
                            // MD characters. Emit a placeholder cs splice op
                            // `~nn<len>nn`: the donor/acceptor motif bytes can't be
                            // recovered from MD + CIGAR + SEQ alone (they'd need the
                            // ref sequence), but maligno's cs parser only inspects
                            // the 2-byte motif slots structurally, so the `nn`
                            // placeholders preserve downstream junction extraction.
                            use std::fmt::Write;
                            cs.push('~');
                            cs.push_str("nn");
                            write!(cs, "{cl}").unwrap();
                            cs.push_str("nn");
                            // cy and my unchanged (N consumes neither query nor MD);
                            // advance k past the N op and continue the inner loop
                            // on the next CIGAR op for the current MD match token.
                            k += 1;
                        }
                        _ => {
                            return Err(anyhow::anyhow!(
                                "line {lineno}: unexpected CIGAR op '{}' during MD walk",
                                op as char
                            ));
                        }
                    }
                }

                if ml != 0 {
                    return Err(anyhow::anyhow!(
                        "line {lineno}: MD tag is inconsistent with CIGAR \
                         (ml={ml} remains after cigar walk)"
                    ));
                }
            }
        }
    }

    Ok(cs)
}

/// Append a short-form match token `:N` to `s`.
#[inline]
fn push_match_short(s: &mut String, n: u32) {
    use std::fmt::Write;
    write!(s, ":{n}").unwrap(); // infallible for String
}

/// Append `src` to `s` with each ASCII letter lowercased.
/// Used for cs op bases (`*xy`, `+seq`, `-seq`) to match minimap2's lowercase
/// convention — both short and long forms use lowercase for non-`=` ops.
#[inline]
fn push_lower(s: &mut String, src: &str) {
    for b in src.bytes() {
        s.push(b.to_ascii_lowercase() as char);
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Tests
// ───────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::cigar::build_merged_cigar;

    #[test]
    fn splice_junction_emits_tilde_op_short_cs() {
        // CIGAR: 5M3N5M — five aligned bases, 3 bp intron, five more aligned bases.
        // MD: "10" — 10 perfect matches across both M blocks (MD has no entries
        // for intron skips).
        // SEQ: 10 query bases.
        let merged = build_merged_cigar("5M3N5M");
        let cs = generate_cs(&merged, "10", "AAAAAGGGGG", false, 1).unwrap();
        assert_eq!(cs, ":5~nn3nn:5");
    }

    #[test]
    fn splice_junction_emits_tilde_op_long_cs() {
        // Same setup but --long-cs format: matches become `=ACGT` instead of `:N`.
        let merged = build_merged_cigar("5M3N5M");
        let cs = generate_cs(&merged, "10", "AAAAAGGGGG", true, 1).unwrap();
        assert_eq!(cs, "=AAAAA~nn3nn=GGGGG");
    }

    #[test]
    fn multiple_introns_in_one_alignment() {
        // CIGAR: 3M50N3M100N4M — two introns inside a 10-base aligned region.
        // MD: 10
        let merged = build_merged_cigar("3M50N3M100N4M");
        let cs = generate_cs(&merged, "10", "AAACCCGGGT", false, 1).unwrap();
        assert_eq!(cs, ":3~nn50nn:3~nn100nn:4");
    }

    #[test]
    fn intron_with_mismatch_after() {
        // CIGAR: 3M50N3M — 6 aligned bases with a mismatch.
        // MD: "4A1" — 4 matches, ref=A mismatch at aligned pos 4, 1 match.
        // (Aligned positions 0..5 in alignment order; intron between 2 and 3.)
        // SEQ has the query bases at the mismatch position; ref is A there.
        let merged = build_merged_cigar("3M50N3M");
        // SEQ at MD-pos 4 (5th aligned base) is G (a substitution: ref A → query G).
        let cs = generate_cs(&merged, "4A1", "AAACGG", false, 1).unwrap();
        // Expected walk:
        //   MD token Match(4): consume 3M (emit :3), N (emit ~nn50nn),
        //     then 1 base of next M (emit :1). ml = 0.
        //   MD token Mismatch(b'a'): emit *ag (both bases lowercased per minimap2).
        //   MD token Match(1): emit :1.
        assert_eq!(cs, ":3~nn50nn:1*ag:1");
    }

    // ─── Lowercase normalization (matches minimap2's cs convention) ────────

    #[test]
    fn lowercase_substitution_in_short_cs() {
        // 1M with ref=A, query=G → expect "*ag" (both lowercased).
        let merged = build_merged_cigar("1M");
        let cs = generate_cs(&merged, "A", "G", false, 1).unwrap();
        assert_eq!(cs, "*ag");
    }

    #[test]
    fn lowercase_insertion_in_short_cs() {
        // 1M2I1M with MD "2" → :1+gg:1 (inserted GG lowercased).
        let merged = build_merged_cigar("1M2I1M");
        let cs = generate_cs(&merged, "2", "AGGC", false, 1).unwrap();
        assert_eq!(cs, ":1+gg:1");
    }

    #[test]
    fn lowercase_deletion() {
        // 1M2D1M with MD "1^AG1" → :1-ag:1 (deleted AG lowercased).
        let merged = build_merged_cigar("1M2D1M");
        let cs = generate_cs(&merged, "1^AG1", "AC", false, 1).unwrap();
        assert_eq!(cs, ":1-ag:1");
    }

    #[test]
    fn long_cs_equals_stays_uppercase_substitution_lowercased() {
        // 11M with MD "5A5", SEQ AAAAAGCCCCC. In long mode the matches use
        // `=ACGT...` (uppercase), but the substitution `*` keeps lowercase
        // per minimap2's convention (only `=` matches are uppercase).
        let merged = build_merged_cigar("11M");
        let cs = generate_cs(&merged, "5A5", "AAAAAGCCCCC", true, 1).unwrap();
        assert_eq!(cs, "=AAAAA*ag=CCCCC");
    }

    #[test]
    fn long_cs_with_intron_and_substitution() {
        // 5M50N5M, MD "4A5", SEQ AAAAGCCCCC (a mismatch right before the intron).
        // Expect long match (uppercase) + substitution (lowercase) + intron + long match.
        let merged = build_merged_cigar("5M50N5M");
        let cs = generate_cs(&merged, "4A5", "AAAAGCCCCC", true, 1).unwrap();
        assert_eq!(cs, "=AAAA*ag~nn50nn=CCCCC");
    }
}
