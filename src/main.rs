//! maligno — unified alignment-comparison toolkit (PAF, SAM and BAM input).
//!
//! Two ways to get from a pair of alignment files (sets A and B) to the per-read
//! alignment comparison table; both produce **identical** tables:
//!
//!   1. `compare` (primary entry point):
//!        `maligno compare -a A.bam -b B.bam --outdir results/ --prefix AvsB`
//!      Sorts both inputs by read name, verifies they share the same read-ID set,
//!      then writes the alignment comparison table and its summary. The per-set
//!      alninfo and readinfo tables are opt-in (`--emit-alninfo` / `--emit-readinfo`).
//!
//!   2. Building blocks, grouped under `toolkit`:
//!        `maligno toolkit paf2tables -i A.paf --readinfo A.readinfo.tsv.gz`
//!        (same for B; for SAM/BAM, pipe in `maligno sam2paf -p -U A.bam`), then
//!        `maligno toolkit compare-readinfo -a A.readinfo.tsv.gz -b B.readinfo.tsv.gz -o compare.tsv.gz`
//!
//! Commands:
//!
//!   - `compare`                      two alignment files -> alignment comparison table
//!                                    (+ summary)
//!   - `sam2paf`                      SAM/BAM -> PAF converter (also run internally
//!                                    by `compare` on SAM/BAM inputs)
//!   - `toolkit paf2tables`           PAF -> alninfo and/or readinfo tables (sorts by
//!                                    query_name first)
//!   - `toolkit compare-readinfo`     two readinfo tables -> alignment comparison table
//!   - `toolkit compare-summary`      alignment comparison table -> summary statistics
//!   - `toolkit find-aln-diff`        alignment comparison table (Parquet) -> differing
//!                                    reads and the regions where they cluster
//!   - `toolkit query-junction-compare`  alignment comparison table (Parquet) -> per-read,
//!                                    per-side splice junctions of the reads whose
//!                                    query-space junctions differ
//!   - `toolkit genomic-junction-compare`  alignment comparison table -> A-vs-B splice
//!                                    junctions in reference coordinates (strand, motifs)
//!
//! Only reads present in both inputs, matched on query_name, produce an output row.


// ── Pipeline modules ──────────────────────────────────────────────────────────
mod aln_input;          // PAF/SAM/BAM input detection; BAM → SAM text → sam2paf stream
mod cigar_junctions;    // CIGAR-based intron extractor (utility; not yet wired in)
mod columns;            // every output column name and type, defined once (alninfo, readinfo, comparison)
mod comparison_row;     // comparison row type, A-vs-B metrics, TSV writers
mod parquet_out;        // Parquet writer for the comparison table (schema derived from columns)
mod compare_streaming;  // `toolkit compare-readinfo` command + merge-join machinery
mod compare_summary;    // `toolkit compare-summary` command + shared classifier/accumulator
mod find_aln_diff;      // `find-aln-diff` command (differing reads + regions, query or reference mode)
mod query_junction_compare; // `query-junction-compare` command (per-read, per-side splice junctions of reads whose query-space junctions differ)
mod genomic_junction_compare; // `genomic-junction-compare` command (A-vs-B junctions in reference coordinates, with strand + motif)
mod interval_merge;     // generic sort+sweep interval merge (bedtools merge -c -o count)
mod cs_parser;          // cs-tag parser  (PAF → alninfo path; also extracts genomic junctions)
mod io_utils;
mod table_input;        // TSV/Parquet dispatch for reading a two-sided comparison table
mod junction;
mod paf;
mod compare;            // primary `compare` command (on-rails: sort → tables → compare)
mod external_sort;      // in-process PAF sort (ext-sort) + read-ID set check
mod paf2tables;         // PAF → sort → alninfo and/or readinfo, one pass
mod paf_groups;         // shared PAF → per-read group reader (paf2tables / compare)
mod readinfo;           // shared collapse library (standalone alninfo → readinfo CLI not registered; run()/ReadInfoArgs/flush_group kept for future reuse)
mod record;

// ── sam2paf utility submodule ─────────────────────────────────────────────────
mod sam2paf; // SAM → PAF converter (self-contained; owns cigar/convert/cs_generator/md)

use anyhow::Result;
use clap::{Parser, Subcommand};

use compare::CompareArgs;
use compare_streaming::CompareReadinfoArgs;
use compare_summary::CompareSummaryArgs;
use find_aln_diff::FindAlnDiffArgs;
use paf2tables::Paf2TablesArgs;
use query_junction_compare::QueryJunctionCompareArgs;
use genomic_junction_compare::GenomicJunctionCompareArgs;
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
    #[command(name = "compare-readinfo")]
    CompareReadinfo(CompareReadinfoArgs),
    /// Alignment comparison table -> aggregate summary statistics.
    #[command(name = "compare-summary")]
    CompareSummary(CompareSummaryArgs),
    /// Alignment comparison table -> differing reads and the regions where they cluster.
    #[command(name = "find-aln-diff")]
    FindAlnDiff(FindAlnDiffArgs),
    /// Alignment comparison table -> per-read, per-side splice junctions of the reads whose query-space junctions differ.
    #[command(name = "query-junction-compare")]
    QueryJunctionCompare(QueryJunctionCompareArgs),
    /// Alignment comparison table -> A-vs-B splice junctions in reference coordinates (strand, motifs).
    #[command(name = "genomic-junction-compare")]
    GenomicJunctionCompare(GenomicJunctionCompareArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Compare(args) => compare::run(args),
        Commands::Sam2paf(args) => sam2paf::run(args),
        Commands::Toolkit { command } => match command {
            ToolkitCommands::Paf2tables(args)    => paf2tables::run(args),
            ToolkitCommands::CompareReadinfo(args) => compare_streaming::run(args),
            ToolkitCommands::CompareSummary(args)  => compare_summary::run(args),
            ToolkitCommands::FindAlnDiff(args)   => find_aln_diff::run(args),
            ToolkitCommands::QueryJunctionCompare(args) => query_junction_compare::run(args),
            ToolkitCommands::GenomicJunctionCompare(args) => genomic_junction_compare::run(args),
        },
    }
}
