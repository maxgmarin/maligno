//! maligno — unified alignment-comparison toolkit (PAF, SAM and BAM input).
//!
//! Two ways to get from a pair of alignment files (sets A and B) to the per-read
//! alignment comparison table; both produce **identical** tables:
//!
//!   1. `compare` (primary entry point):
//!        `maligno compare -a A.bam -b B.bam --outdir results/ --prefix AvsB`
//!      Sorts both inputs by read name, verifies they share the same read-ID set,
//!      then writes the alignment comparison table and, by default,
//!      `find-aln-diff`'s differing-reads and region tables
//!      (`--skip-find-aln-diff` to opt out). The per-set alninfo and readinfo
//!      tables are opt-in (`--emit-alninfo` / `--emit-readinfo`).
//!
//!   2. Building blocks, grouped under `toolkit`:
//!        `maligno toolkit paf2tables -i A.sorted.paf --readinfo A.readinfo.tsv.gz`
//!        (same for B), then
//!        `maligno toolkit merge-readinfo -a A.readinfo.tsv.gz -b B.readinfo.tsv.gz -o compare.tsv.gz`
//!
//! Commands:
//!
//!   - `compare`                      two alignment files -> alignment comparison table
//!                                    (+ summary and differing reads/regions)
//!   - `sam2paf`                      SAM/BAM -> PAF converter (also run internally
//!                                    by `compare` on SAM/BAM inputs)
//!   - `toolkit paf2tables`           PAF -> alninfo and/or readinfo tables
//!   - `toolkit merge-readinfo`       two readinfo tables -> alignment comparison table
//!   - `toolkit summary`              alignment comparison table -> summary statistics
//!   - `toolkit find-aln-diff`        alignment comparison table -> differing reads and
//!                                    the regions where they cluster
//!   - `toolkit query-junction-diff`  alignment comparison table (Parquet) -> per-read,
//!                                    per-side splice-junction reconstruction and diff
//!
//! Only reads present in both inputs, matched on Read_Name, produce an output row.


// ── Pipeline modules ──────────────────────────────────────────────────────────
mod aln_input;          // PAF/SAM/BAM input detection; BAM → SAM text → sam2paf stream
mod cigar_junctions;    // CIGAR-based intron extractor (utility; not yet wired in)
mod comparison_row;     // comparison-table schema: column lists, row type, TSV writers
mod parquet_out;        // Parquet writer for the comparison table (schema derived from comparison_row)
mod compare_streaming;  // `toolkit merge-readinfo` command + merge-join machinery
mod compare_summary;    // `toolkit summary` command + shared classifier/accumulator
mod find_query_diff;   // `find-aln-diff` command (differing reads + regions, query or reference space)
mod query_junction_diff; // `query-junction-diff` command (per-read, per-side splice-junction reconstruction + diff)
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
mod readinfo;           // shared collapse library (standalone alninfo → readinfo CLI not registered; run()/ReadInfoArgs/flush_group kept for future reuse)
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
use query_junction_diff::QueryJunctionDiffArgs;
use sam2paf::Sam2pafArgs;

/// Unified alignment-comparison toolkit.
#[derive(Parser, Debug)]
#[command(name = "maligno", version, about, disable_help_subcommand = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// End-to-end comparison of two alignment files (primary entry point).
    Compare(CompareArgs),
    /// SAM/BAM -> PAF converter.
    Sam2paf(Sam2pafArgs),
    /// Building blocks, and follow-up analyses on an alignment comparison table.
    #[command(name = "toolkit")]
    Toolkit {
        #[command(subcommand)]
        command: ToolkitCommands,
    },
}

#[derive(Subcommand, Debug)]
enum ToolkitCommands {
    /// PAF -> alninfo and/or readinfo tables.
    Paf2tables(Paf2TablesArgs),
    /// Two readinfo tables -> alignment comparison table (TSV or Parquet).
    #[command(name = "merge-readinfo")]
    MergeReadinfo(MergeReadinfoArgs),
    /// Alignment comparison table -> aggregate summary statistics.
    Summary(CompareSummaryArgs),
    /// Alignment comparison table -> differing reads and the regions where they cluster.
    #[command(name = "find-aln-diff")]
    FindAlnDiff(FindAlnDiffArgs),
    /// Alignment comparison table -> per-read, per-side splice-junction reconstruction and diff.
    #[command(name = "query-junction-diff")]
    QueryJunctionDiff(QueryJunctionDiffArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Compare(args) => compare::run(args),
        Commands::Sam2paf(args) => sam2paf::run(args),
        Commands::Toolkit { command } => match command {
            ToolkitCommands::Paf2tables(args)    => paf2tables::run(args),
            ToolkitCommands::MergeReadinfo(args) => compare_streaming::run(args),
            ToolkitCommands::Summary(args)       => compare_summary::run(args),
            ToolkitCommands::FindAlnDiff(args)   => find_query_diff::run(args),
            ToolkitCommands::QueryJunctionDiff(args) => query_junction_diff::run(args),
        },
    }
}
