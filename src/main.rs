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
//!   2. Manual building blocks (full control), grouped under `compare-pipeline`:
//!        `maligno compare-pipeline paf2tables -i A.sorted.paf --alninfo A.alninfo.tsv.gz --readinfo A.readinfo.tsv.gz`
//!        (same for B), then
//!        `maligno compare-pipeline merge-readinfo -a A.readinfo.tsv.gz -b B.readinfo.tsv.gz -o compare.tsv.gz`
//!
//! Commands:
//!
//!   1. `compare`                        end-to-end comparison of all input
//!                                       alignments (PRIMARY analysis entry
//!                                       point). Emits the single 96-column
//!                                       comparison table.
//!   2. `sam2paf`                        SAM → PAF converter (utility; use
//!                                       before compare-pipeline/compare)
//!   3. `find-aln-diff`                  comparison table → differing reads +
//!                                       the merged genomic regions where they
//!                                       cluster, in query or reference space
//!   4. `compare-pipeline paf2tables`    PAF → alninfo TSV and/or readinfo TSV
//!                                       tables — a decomposed piece of what
//!                                       `compare` does internally
//!   5. `compare-pipeline merge-readinfo` two readinfo TSVs → per-read
//!                                       comparison TSV — ditto
//!   6. `compare-pipeline summary`       comparison table → aggregate summary
//!                                       statistics (alignment status +
//!                                       query/reference identity) — the same
//!                                       thing `compare` tallies inline
//!
//! The comparison itself is a streaming merge-join (constant memory): only reads
//! present in BOTH inputs (matched on Read_Name + Read_Len) produce an output row.
//! `compare` guarantees the inputs are sorted and share the same read-ID set.


// ── Pipeline modules ──────────────────────────────────────────────────────────
mod cigar_junctions;    // CIGAR-based intron extractor (utility; not yet wired in)
mod comparison_row;     // comparison-table schema: column lists, row type, TSV writers
mod parquet_out;        // Parquet writer for the comparison table (schema derived from comparison_row)
mod compare_streaming;  // `compare-pipeline merge-readinfo` command + merge-join machinery
mod compare_summary;    // `compare-pipeline summary` command + shared classifier/accumulator
mod find_query_diff;   // `find-aln-diff` command (differing reads + regions, query or reference space)
mod interval_merge;     // generic sort+sweep interval merge (bedtools merge -c -o count)
mod cs_parser;          // cs-tag parser  (PAF → alninfo path; also extracts genomic junctions)
mod io_utils;
mod table_input;        // TSV/Parquet dispatch for reading a two-sided comparison table
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
use compare_streaming::MergeReadinfoArgs;
use compare_summary::CompareSummaryArgs;
use find_query_diff::FindAlnDiffArgs;
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
    /// Comparison table → differing reads and the merged genomic regions where
    /// they cluster, in query or reference space.
    #[command(name = "find-aln-diff")]
    FindAlnDiff(FindAlnDiffArgs),
    /// Lower-level building blocks and analysis steps used internally by `compare`.
    #[command(name = "compare-pipeline")]
    ComparePipeline {
        #[command(subcommand)]
        command: PipelineCommands,
    },
}

#[derive(Subcommand, Debug)]
enum PipelineCommands {
    /// PAF -> alninfo TSV and/or readinfo TSV tables.
    Paf2tables(Paf2TablesArgs),
    /// Two readinfo TSVs -> per-read comparison TSV.
    #[command(name = "merge-readinfo")]
    MergeReadinfo(MergeReadinfoArgs),
    /// Comparison table → aggregate summary statistics.
    Summary(CompareSummaryArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Compare(args)      => compare::run(args),
        Commands::Sam2paf(args)      => sam2paf::run(args),
        Commands::FindAlnDiff(args)  => find_query_diff::run(args),
        Commands::ComparePipeline { command } => match command {
            PipelineCommands::Paf2tables(args)    => paf2tables::run(args),
            PipelineCommands::MergeReadinfo(args) => compare_streaming::run(args),
            PipelineCommands::Summary(args)       => compare_summary::run(args),
        },
    }
}
