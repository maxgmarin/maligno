//! maligno — unified PAF alignment-comparison toolkit.
//!
//! Two ways to get from a pair of PAFs (sample A and B) to a per-read comparison
//! table — both produce **identical** comparison results:
//!
//!   1. On-rails (primary) — `compare`:
//!        `maligno compare -a A.paf -b B.paf --outdir results/ --prefix AvsB`
//!      Sorts both PAFs (consistent order), verifies they share the same read-ID
//!      set. Then it writes the per-set alninfo + readinfo tables AND the comparison
//!      table. 
//!
//!   2. Manual building blocks (full control):
//!        `maligno paf2tables -i A.sorted.paf --alninfo A.alninfo.tsv.gz --readinfo A.readinfo.tsv.gz`
//!        (same for B), then
//!        `maligno compare-readinfo -a A.readinfo.tsv.gz -b B.readinfo.tsv.gz -o compare.tsv.gz`
//!
//! Commands:
//!
//!   1. `compare`            end-to-end comparison of all input alignments
//!                           (PRIMARY analysis entry point). `--mode full`
//!                           (default, 94 cols) or `--mode junctions` (47-col view).
//!   2. `sam2paf`            SAM → PAF converter (utility; use before paf2tables/compare)
//!   3. `paf2tables`         PAF → alninfo TSV and/or readinfo TSV tables
//!   4. `compare-readinfo`   two readinfo TSVs → per-read comparison TSV (same `--mode`)
//!   5. `compare-summary`    comparison TSV → aggregate summary statistics
//!                           (alignment status + query/reference identity)
//!   6. `find-query-diff`   comparison TSV → query-different reads + the merged
//!                           genomic regions where they cluster
//!
//! The comparison itself is a streaming merge-join (constant memory): only reads
//! present in BOTH inputs (matched on Read_Name + Read_Len) produce an output row.
//! `compare` guarantees the inputs are sorted and share the same read-ID set.


// ── Pipeline modules ──────────────────────────────────────────────────────────
mod cigar_junctions;    // CIGAR-based intron extractor (utility; not yet wired in)
mod compare_junctions;  // junction-view (47-col) header/row emitters (library; --mode junctions)
mod compare_streaming;  // `compare-readinfo` command + shared comparison core
mod compare_summary;    // `compare-summary` command + shared classifier/accumulator
mod find_query_diff;   // `find-query-diff` command (query-different reads + regions)
mod interval_merge;     // generic sort+sweep interval merge (bedtools merge -c -o count)
mod cs_parser;          // cs-tag parser  (PAF → alninfo path; also extracts genomic junctions)
mod io_utils;
mod junction;
mod paf;
mod compare;            // primary `compare` command (on-rails: sort → tables → compare)
mod external_sort;      // in-process PAF sort (ext-sort) + read-ID set check
mod paf2tables;         // PAF → alninfo and/or readinfo, one pass
mod paf_groups;         // shared PAF → per-read group reader (paf2tables / compare)
mod readinfo;           // shared collapse library (utils-readinfo CLI unregistered; run()/ReadInfoArgs/flush_group kept for future reuse)
mod record;

// ── sam2paf utility submodule ─────────────────────────────────────────────────
mod sam2paf; // SAM → PAF converter (self-contained; owns cigar/convert/cs_generator/md)

use anyhow::Result;
use clap::{Parser, Subcommand};

use compare::CompareArgs;
use compare_streaming::CompareReadinfoArgs;
use compare_summary::CompareSummaryArgs;
use find_query_diff::FindQueryDiffArgs;
use paf2tables::Paf2TablesArgs;
use sam2paf::Sam2pafArgs;

/// Unified alignment-comparison toolkit.
#[derive(Parser, Debug)]
#[command(name = "maligno", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// end-to-end comparison of all input alignments (Primary analysis entry point).
    Compare(CompareArgs),
    /// SAM -> PAF converter (conversion utility).
    Sam2paf(Sam2pafArgs),
    /// PAF -> alninfo TSV and/or readinfo TSV tables.
    Paf2tables(Paf2TablesArgs),
    /// Two readinfo TSVs -> per-read comparison TSV (--mode full|junctions).
    CompareReadinfo(CompareReadinfoArgs),
    /// Comparison TSV → aggregate summary statistics
    CompareSummary(CompareSummaryArgs),
    /// Comparison TSV → query-different reads and the merged genomic regions where they cluster.
    #[command(name = "find-query-diff")]
    FindQueryDiff(FindQueryDiffArgs),
}




fn main() -> Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Compare(args)          => compare::run(args),
        Commands::Sam2paf(args)          => sam2paf::run(args),
        Commands::Paf2tables(args)       => paf2tables::run(args),
        Commands::CompareReadinfo(args)  => compare_streaming::run(args),
        Commands::CompareSummary(args)   => compare_summary::run(args),
        Commands::FindQueryDiff(args)    => find_query_diff::run(args),
    }
}
