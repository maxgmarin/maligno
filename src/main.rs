//! maligno — unified PAF alignment-comparison toolkit.
//!
//! **Start here:** `paf2tables` is the primary PAF entry point — it turns a PAF
//! into the per-alignment info table, the per-read summary table, or both in a
//! single pass:
//!
//!   1. `paf2tables`         PAF → alninfo TSV (35 cols) and/or readinfo TSV (33 cols)
//!
//! Comparison + utility subcommands:
//!
//!   2. `sam2paf`            SAM                → PAF  (utility; use before paf2tables)
//!   3. `compare`            two readinfo TSVs  → per-read comparison TSV. `--mode full`
//!                                                (default, 94 cols incl. genomic-junction
//!                                                comparison) or `--mode junctions` (47-col)
//!   4. `pafcompare`         two PAFs           → comparison TSV (one pass; same `--mode`;
//!                                                requires identical read-name order)
//!   5. `utils-readinfo`     alninfo TSV        → per-read summary TSV (low-level utility;
//!                                                most users want `paf2tables --readinfo`)
//!
//! Full pipeline (including optional SAM conversion):
//!
//!     SAM ──sam2paf──▶ PAF ──paf2tables──▶ readinfo.tsv ─┐
//!                                                         ├─compare─▶ compare.tsv
//!     SAM ──sam2paf──▶ PAF ──paf2tables──▶ readinfo.tsv ─┘
//!
//! The comparison is a streaming merge-join (constant memory): only reads present
//! in BOTH readinfo files (matched on Read_Name + Read_Len) produce an output row.
//! By default a read-name mismatch is an error; `--ignore-row-mismatch` skips
//! reads present in only one file (counted in the end-of-run summary). Unaligned
//! reads still appear in each readinfo file (with Num_Aln=0), so they DO get compared.

// ── Pipeline modules ──────────────────────────────────────────────────────────
mod cigar_junctions;    // CIGAR-based intron extractor (utility; not yet wired in)
mod compare_junctions;  // junction-view header/row emitters (library; compare --mode junctions, pafcompare)
mod compare_streaming;
mod cs_parser;          // cs-tag parser  (PAF → alninfo path; also extracts genomic junctions)
mod io_utils;
mod junction;
mod paf;
mod paf2tables;         // primary PAF entry point (alninfo and/or readinfo, one pass)
mod paf_groups;         // shared PAF → per-read group reader (paf2tables / pafcompare)
mod pafcompare;         // fused paired-PAF → comparison (one pass)
mod readinfo;           // shared collapse library + utils-readinfo entry point
mod record;

// ── sam2paf utility submodule ─────────────────────────────────────────────────
mod sam2paf; // SAM → PAF converter (self-contained; owns cigar/convert/cs_generator/md)

use anyhow::Result;
use clap::{Parser, Subcommand};

use compare_streaming::CompareStreamingArgs;
use paf2tables::Paf2TablesArgs;
use pafcompare::PafCompareArgs;
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
    /// PAF → alninfo (35 cols) and/or readinfo (33 cols). The primary PAF entry point.
    Paf2tables(Paf2TablesArgs),
    /// SAM → PAF converter (utility; use before paf2tables to start the pipeline).
    Sam2paf(Sam2pafArgs),
    /// Compare two readinfo TSVs. --mode full (default) or junctions; strict order by default.
    Compare(CompareStreamingArgs),
    /// Fused paired-PAF → comparison in one pass (--mode full|junctions; identical read-name order).
    Pafcompare(PafCompareArgs),
    /// [utility] alninfo TSV → per-read summary TSV. Most users want `paf2tables --readinfo`.
    #[command(name = "utils-readinfo")]
    UtilsReadinfo(ReadInfoArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Paf2tables(args)       => paf2tables::run(args),
        Commands::Sam2paf(args)          => sam2paf::run(args),
        Commands::Compare(args)          => compare_streaming::run(args),
        Commands::Pafcompare(args)       => pafcompare::run(args),
        Commands::UtilsReadinfo(args)    => readinfo::run(args),
    }
}
