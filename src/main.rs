//! maligno — unified PAF alignment-comparison toolkit.
//!
//! Two ways to get from a pair of PAFs (sample A and B) to a per-read comparison
//! table — both produce **identical** results:
//!
//!   1. One-pass (primary):
//!        `maligno compare -a A.paf -b B.paf -o compare.tsv.gz`
//!      Does everything internally; no intermediate files.
//!
//!   2. Explicit intermediates (useful when you also want the per-file tables):
//!        `maligno paf2tables -i A.paf --readinfo A.readinfo.tsv.gz`
//!        `maligno paf2tables -i B.paf --readinfo B.readinfo.tsv.gz`
//!        `maligno compare-readinfo -a A.readinfo.tsv.gz -b B.readinfo.tsv.gz -o compare.tsv.gz`
//!
//! Commands:
//!
//!   1. `compare`            two PAFs           → per-read comparison TSV (one pass). PRIMARY.
//!                                                `--mode full` (default, 94 cols incl.
//!                                                genomic-junction comparison) or
//!                                                `--mode junctions` (47-col splice view).
//!   2. `paf2tables`         PAF                → alninfo TSV (35 cols) and/or readinfo TSV (33 cols)
//!   3. `compare-readinfo`   two readinfo TSVs  → per-read comparison TSV (same `--mode`)
//!   4. `sam2paf`            SAM                → PAF  (utility; use before paf2tables/compare)
//!   5. `utils-readinfo`     alninfo TSV        → per-read summary TSV (low-level utility;
//!                                                most users want `paf2tables --readinfo`)
//!
//! Comparison is a streaming merge-join (constant memory): only reads present in
//! BOTH inputs (matched on Read_Name + Read_Len) produce an output row. By default
//! a read-name mismatch is an error; `--ignore-row-mismatch` skips reads present
//! in only one input (counted in the end-of-run summary). Unaligned reads still
//! appear (with Num_Aln=0), so they DO get compared.

// ── Pipeline modules ──────────────────────────────────────────────────────────
mod cigar_junctions;    // CIGAR-based intron extractor (utility; not yet wired in)
mod compare_junctions;  // junction-view (47-col) header/row emitters (library; --mode junctions)
mod compare_streaming;  // `compare-readinfo` command + shared comparison core
mod cs_parser;          // cs-tag parser  (PAF → alninfo path; also extracts genomic junctions)
mod io_utils;
mod junction;
mod paf;
mod paf2tables;         // PAF → alninfo and/or readinfo, one pass
mod paf_groups;         // shared PAF → per-read group reader (paf2tables / compare)
mod pafcompare;         // primary `compare` command (paired-PAF → comparison, one pass)
mod readinfo;           // shared collapse library + utils-readinfo entry point
mod record;

// ── sam2paf utility submodule ─────────────────────────────────────────────────
mod sam2paf; // SAM → PAF converter (self-contained; owns cigar/convert/cs_generator/md)

use anyhow::Result;
use clap::{Parser, Subcommand};

use compare_streaming::CompareReadinfoArgs;
use paf2tables::Paf2TablesArgs;
use pafcompare::CompareArgs;
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
    /// Two PAFs → per-read comparison TSV in one pass. The primary entry point.
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
        Commands::Compare(args)          => pafcompare::run(args),
        Commands::Paf2tables(args)       => paf2tables::run(args),
        Commands::CompareReadinfo(args)  => compare_streaming::run(args),
        Commands::Sam2paf(args)          => sam2paf::run(args),
        Commands::UtilsReadinfo(args)    => readinfo::run(args),
    }
}
