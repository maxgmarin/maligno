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
                cs.push('-');
                cs.push_str(bases);
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
                            cs.push('+');
                            cs.push_str(&seq[cy..cy_end]);
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
                                        // *<ref_base><query_base>
                                        cs.push('*');
                                        cs.push(ref_base as char);
                                        cs.push(seq_bytes[my] as char);
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
                                if long_cs {
                                    cs.push('=');
                                    cs.push_str(&seq[my..my + dl]);
                                } else {
                                    push_match_short(&mut cs, dl as u32);
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
