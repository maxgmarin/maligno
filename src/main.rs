//! maligno — unified PAF alignment-comparison toolkit.
//!
//! Two ways to get from a pair of PAFs (sample A and B) to a per-read comparison
//! table — both produce **identical** comparison results:
//!
//!   1. On-rails (primary) — `compare`:
//!        `maligno compare -a A.paf -b B.paf --outdir results/ --prefix AvsB`
//!      Sorts both PAFs (consistent order), verifies they share the same read-ID
//!      set, then writes the per-set alninfo + readinfo tables AND the comparison
//!      table. Makes the key assumptions for you and does them well.
//!
//!   2. Manual building blocks (full control):
//!        `maligno paf2tables -i A.sorted.paf --alninfo A.alninfo.tsv.gz --readinfo A.readinfo.tsv.gz`
//!        (same for B), then
//!        `maligno compare-readinfo -a A.readinfo.tsv.gz -b B.readinfo.tsv.gz -o compare.tsv.gz`
//!
//! Commands:
//!
//!   1. `compare`            two PAFs → a results directory (sorted-then-compared,
//!                                      on-rails). PRIMARY. `--mode full` (default,
//!                                      94 cols) or `--mode junctions` (47-col view).
//!   2. `paf2tables`         PAF                → alninfo TSV (35 cols) and/or readinfo TSV (33 cols)
//!   3. `compare-readinfo`   two readinfo TSVs  → per-read comparison TSV (same `--mode`)
//!   4. `sam2paf`            SAM                → PAF  (utility; use before paf2tables/compare)
//!   5. `utils-readinfo`     alninfo TSV        → per-read summary TSV (low-level utility;
//!                                                most users want `paf2tables --readinfo`)
//!
//! The comparison itself is a streaming merge-join (constant memory): only reads
//! present in BOTH inputs (matched on Read_Name + Read_Len) produce an output row.
//! `compare` guarantees the inputs are sorted and share the same read-ID set
//! (erroring otherwise, unless `--allow-id-mismatch`); the plumbing commands trust
//! the caller and use `--ignore-row-mismatch` for the lenient path.

// ── Pipeline modules ──────────────────────────────────────────────────────────
mod cigar_junctions;    // CIGAR-based intron extractor (utility; not yet wired in)
mod compare_junctions;  // junction-view (47-col) header/row emitters (library; --mode junctions)
mod compare_streaming;  // `compare-readinfo` command + shared comparison core
mod cs_parser;          // cs-tag parser  (PAF → alninfo path; also extracts genomic junctions)
mod io_utils;
mod junction;
mod paf;
mod compare;            // primary `compare` command (on-rails: sort → tables → compare)
mod external_sort;      // in-process PAF sort (ext-sort) + read-ID set check
mod paf2tables;         // PAF → alninfo and/or readinfo, one pass
mod paf_groups;         // shared PAF → per-read group reader (paf2tables / compare)
mod readinfo;           // shared collapse library + utils-readinfo entry point
mod record;

// ── sam2paf utility submodule ─────────────────────────────────────────────────
mod sam2paf; // SAM → PAF converter (self-contained; owns cigar/convert/cs_generator/md)

use anyhow::Result;
use clap::{Parser, Subcommand};

use compare::CompareArgs;
use compare_streaming::CompareReadinfoArgs;
use paf2tables::Paf2TablesArgs;
use readinfo::ReadInfoArgs;
use sam2paf::Sam2pafArgs;

/// Unified PAF alignment-comparison toolkit.
#[derive(Parser, Debug)]
#[command(name = "maligno", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Two PAFs → sorted, verified, and compared into a results directory. The primary entry point.
    Compare(CompareArgs),
    /// PAF → alninfo (35 cols) and/or readinfo (33 cols).
    Paf2tables(Paf2TablesArgs),
    /// Two readinfo TSVs → per-read comparison TSV (--mode full|junctions; strict order by default).
    CompareReadinfo(CompareReadinfoArgs),
    /// SAM → PAF converter (utility; use before paf2tables/compare to start the pipeline).
    Sam2paf(Sam2pafArgs),
    /// [utility] alninfo TSV → per-read summary TSV. Most users want `paf2tables --readinfo`.
    #[command(name = "utils-readinfo")]
    UtilsReadinfo(ReadInfoArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Compare(args)          => compare::run(args),
        Commands::Paf2tables(args)       => paf2tables::run(args),
        Commands::CompareReadinfo(args)  => compare_streaming::run(args),
        Commands::Sam2paf(args)          => sam2paf::run(args),
        Commands::UtilsReadinfo(args)    => readinfo::run(args),
    }
}
