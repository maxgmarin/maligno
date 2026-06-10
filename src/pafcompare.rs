//! Subcommand `pafcompare`: fused paired-PAF → comparison TSV in one pass.
//!
//! Equivalent to the full discrete pipeline
//! (`paf2alninfo | readinfo` on each side, then `compare`), but without
//! materializing any intermediate files. Both PAFs are read in **lock-step**:
//! one read's group is pulled from each side, collapsed to a `ReadInfoRow`, and
//! emitted as a comparison row.
//!
//! ## Ordering contract
//!
//! By default `pafcompare` is strict: it requires that both PAFs list the
//! **same `Query_Names` in the same order**. On the first read-name mismatch
//! — or if one side runs out of reads before the other — it prints an ERROR
//! to stderr and stops with a non-zero exit.
//!
//! With `--ignore-row-mismatch`, reads that appear in only one PAF are skipped
//! using a lex-name comparison to decide which side to advance. This requires
//! both PAFs to be lex-sorted by Query_Name. Unmatched reads are counted in
//! the end-of-run summary.
//!
//! Output is byte-identical to the discrete pipeline because each side's
//! `ReadInfoRow` is serialized to the exact readinfo-TSV line `readinfo` would
//! produce, then fed into the same `emit_compare_row` /
//! `emit_compare_junctions_row` used by `compare` (`--mode full` / `junctions`).

use std::collections::HashMap;

use anyhow::{bail, Result};

use crate::compare_junctions::{emit_compare_junctions_row, write_compare_junctions_header};
use crate::compare_streaming::{emit_compare_row, write_compare_header, CompareMode};
use crate::io_utils::open_input;
use crate::io_utils::open_output;
use crate::paf_groups::PafGroups;
use crate::readinfo::{collapse_group, AlnRow, ReadInfoRow, READINFO_HEADER};

#[derive(clap::Args, Debug)]
pub struct PafCompareArgs {
    /// PAF for dataset A. Use '-' for stdin; '.gz' auto-decompressed.
    #[arg(short = 'a', long = "paf-a", value_name = "a.paf")]
    paf_a: String,

    /// PAF for dataset B. Use '-' for stdin; '.gz' auto-decompressed.
    #[arg(short = 'b', long = "paf-b", value_name = "b.paf")]
    paf_b: String,

    /// Label for dataset A (used as per-side column suffix).
    #[arg(long = "label-a", value_name = "LABEL", default_value = "SetA")]
    label_a: String,

    /// Label for dataset B (used as per-side column suffix).
    #[arg(long = "label-b", value_name = "LABEL", default_value = "SetB")]
    label_b: String,

    /// Output comparison TSV file ('.gz' for gzip).
    #[arg(short = 'o', long = "output", value_name = "compare.tsv[.gz]")]
    output: String,

    /// Comparison view: `full` (all per-read metrics incl. genomic-junction
    /// comparison, 94 cols) or `junctions` (47-col splice-focused). Mirrors
    /// `compare --mode`.
    #[arg(long = "mode", value_enum, default_value_t = CompareMode::Full)]
    mode: CompareMode,

    /// Skip reads that appear in only one PAF instead of stopping with an error.
    /// Requires both PAFs to be lex-sorted by Query_Name for the skip heuristic
    /// to work correctly. Unmatched reads are counted in the summary.
    #[arg(long = "ignore-row-mismatch")]
    ignore_row_mismatch: bool,
}

fn collapse(mut group: Vec<AlnRow>) -> ReadInfoRow {
    collapse_group(&mut group)
}

/// Serialize a `ReadInfoRow` into the exact readinfo-TSV line and return it
/// (without trailing newline) so it can be parsed into a by-name column map.
fn readinfo_line(ri: &ReadInfoRow) -> Result<String> {
    let mut buf: Vec<u8> = Vec::with_capacity(256);
    ri.write(&mut buf)?;
    let s = String::from_utf8(buf).expect("readinfo serialization is valid UTF-8");
    Ok(s.trim_end_matches(|c| c == '\n' || c == '\r').to_string())
}

/// Drain remaining groups from a `PafGroups` reader, returning the count.
fn drain<R: std::io::BufRead>(groups: &mut PafGroups<R>) -> Result<u64> {
    let mut n: u64 = 0;
    while groups.next_group()?.is_some() {
        n += 1;
    }
    Ok(n)
}

pub fn run(args: &PafCompareArgs) -> Result<()> {
    eprintln!("[INFO] Opening PAF files...");
    eprintln!("  A: {}", args.paf_a);
    eprintln!("  B: {}", args.paf_b);

    let mut out = open_output(Some(&args.output))?;

    if matches!(args.mode, CompareMode::Junctions) {
        write_compare_junctions_header(&mut out, &args.label_a, &args.label_b)?;
    } else {
        write_compare_header(&mut out, &args.label_a, &args.label_b)?;
    }

    // Lock-step zip needs only identical ordering, not byte-lex sort, so the
    // group reader's lex-decrease warning is disabled unless --ignore-row-mismatch
    // is set (in which case lex order is required for the skip heuristic).
    let mut groups_a = PafGroups::new(open_input(&args.paf_a)?, false);
    let mut groups_b = PafGroups::new(open_input(&args.paf_b)?, false);

    let header_cols: Vec<&str> = READINFO_HEADER.split('\t').collect();

    // Pending-row pattern: hold one collapsed ReadInfoRow per side so that
    // after a skip (--ignore-row-mismatch), the lex-larger side's row is
    // reused in the next iteration instead of being discarded.
    let mut pending_a: Option<ReadInfoRow> =
        groups_a.next_group()?.map(collapse);
    let mut pending_b: Option<ReadInfoRow> =
        groups_b.next_group()?.map(collapse);

    let mut n_matched: u64 = 0;
    let mut n_a_only: u64 = 0;
    let mut n_b_only: u64 = 0;

    loop {
        match (pending_a.take(), pending_b.take()) {
            (None, None) => break,

            // One side exhausted — strict mode errors, lenient mode drains.
            (Some(ra), None) => {
                if !args.ignore_row_mismatch {
                    bail!(
                        "PAF A has more reads than PAF B (B exhausted after \
                         {n_matched} matched; next unmatched A read is {:?}). \
                         pafcompare requires both PAFs to list the same \
                         Query_Names in the same order. Use \
                         --ignore-row-mismatch to allow mismatches.",
                        ra.read_name
                    );
                }
                n_a_only += 1 + drain(&mut groups_a)?;
                break;
            }
            (None, Some(rb)) => {
                if !args.ignore_row_mismatch {
                    bail!(
                        "PAF B has more reads than PAF A (A exhausted after \
                         {n_matched} matched; next unmatched B read is {:?}). \
                         pafcompare requires both PAFs to list the same \
                         Query_Names in the same order. Use \
                         --ignore-row-mismatch to allow mismatches.",
                        rb.read_name
                    );
                }
                n_b_only += 1 + drain(&mut groups_b)?;
                break;
            }

            // Both sides have a row.
            (Some(ra), Some(rb)) => {
                if ra.read_name == rb.read_name {
                    // MATCH: emit and advance both.
                    let line_a = readinfo_line(&ra)?;
                    let line_b = readinfo_line(&rb)?;
                    let map_a: HashMap<&str, &str> =
                        header_cols.iter().copied().zip(line_a.split('\t')).collect();
                    let map_b: HashMap<&str, &str> =
                        header_cols.iter().copied().zip(line_b.split('\t')).collect();

                    if matches!(args.mode, CompareMode::Junctions) {
                        emit_compare_junctions_row(
                            &mut out,
                            &ra.read_name,
                            ra.read_len,
                            |c| *map_a.get(c).unwrap_or(&""),
                            |c| *map_b.get(c).unwrap_or(&""),
                        )?;
                    } else {
                        emit_compare_row(
                            &mut out,
                            &ra.read_name,
                            ra.read_len,
                            |c| *map_a.get(c).unwrap_or(&""),
                            |c| *map_b.get(c).unwrap_or(&""),
                        )?;
                    }

                    n_matched += 1;
                    if n_matched % 100_000 == 0 {
                        eprintln!("[INFO] Compared {n_matched} reads...");
                    }
                    pending_a = groups_a.next_group()?.map(collapse);
                    pending_b = groups_b.next_group()?.map(collapse);
                } else if !args.ignore_row_mismatch {
                    bail!(
                        "read-name mismatch at read #{}: A has {:?} but B has {:?}. \
                         pafcompare requires both PAFs to list the same Query_Names \
                         in the same order. Use --ignore-row-mismatch to skip \
                         unmatched reads (requires lex-sorted PAFs).",
                        n_matched + 1,
                        ra.read_name,
                        rb.read_name,
                    );
                } else if ra.read_name < rb.read_name {
                    // A is behind: skip A's read, keep B for next iteration.
                    n_a_only += 1;
                    pending_b = Some(rb);
                    pending_a = groups_a.next_group()?.map(collapse);
                } else {
                    // B is behind: skip B's read, keep A for next iteration.
                    n_b_only += 1;
                    pending_a = Some(ra);
                    pending_b = groups_b.next_group()?.map(collapse);
                }
            }
        }
    }

    out.flush()?;

    eprintln!("Paired-PAF comparison summary:");
    eprintln!("  Label A: {}", args.label_a);
    eprintln!("  Label B: {}", args.label_b);
    eprintln!("  reads matched (written):    {n_matched}");
    eprintln!("  A-only (skipped, not in B): {n_a_only}");
    eprintln!("  B-only (skipped, not in A): {n_b_only}");

    Ok(())
}
