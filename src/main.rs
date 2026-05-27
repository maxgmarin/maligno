//! maligno — unified PAF alignment-comparison toolkit.
//!
//! A single binary with four subcommands covering the full pipeline:
//!
//!   1. `paf2alninfo`  PAF                → per-alignment info TSV (35 cols)
//!   2. `readinfo`     alninfo TSV        → per-read summary TSV   (28 cols)
//!   3. `compare`      two readinfo TSVs  → per-read comparison TSV (77 cols)
//!   4. `sam2paf`      SAM                → PAF  (utility; use before paf2alninfo)
//!
//! Full pipeline (including optional SAM conversion):
//!
//!     SAM ──sam2paf──▶ PAF ──paf2alninfo──▶ alninfo.tsv ──readinfo──▶ readinfo.tsv ─┐
//!                                                                                    ├─compare─▶ compare.tsv
//!     SAM ──sam2paf──▶ PAF ──paf2alninfo──▶ alninfo.tsv ──readinfo──▶ readinfo.tsv ─┘
//!
//! The comparison is a streaming merge-join (constant memory) and an INNER
//! join: only reads present in BOTH readinfo files (matched on Read_Name +
//! Read_Len) produce an output row. Reads present in only one file are dropped
//! but counted in the end-of-run summary. Note that unaligned reads still
//! appear in each readinfo file (with Num_Aln=0), so they DO get compared.

// ── Pipeline modules ──────────────────────────────────────────────────────────
mod compare_streaming;
mod cs_parser;       // cs-tag parser  (PAF → alninfo path)
mod io_utils;
mod junction;
mod paf;
mod paf2alninfo;
mod readinfo;
mod record;

// ── sam2paf utility submodule ─────────────────────────────────────────────────
mod sam2paf; // SAM → PAF converter (self-contained; owns cigar/convert/cs_generator/md)

use anyhow::Result;
use clap::{Parser, Subcommand};

use compare_streaming::CompareStreamingArgs;
use paf2alninfo::Paf2AlnInfoArgs;
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
    /// PAF → per-alignment info TSV (keeps unaligned reads as zeroed rows).
    Paf2alninfo(Paf2AlnInfoArgs),
    /// alninfo TSV → per-read summary TSV (best alignment by ms, then AS).
    Readinfo(ReadInfoArgs),
    /// Streaming merge-join comparison of two sorted readinfo TSVs (constant memory).
    Compare(CompareStreamingArgs),
    /// SAM → PAF converter (utility; use before paf2alninfo to start the pipeline).
    Sam2paf(Sam2pafArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Commands::Paf2alninfo(args) => paf2alninfo::run(args),
        Commands::Readinfo(args)    => readinfo::run(args),
        Commands::Compare(args)     => compare_streaming::run(args),
        Commands::Sam2paf(args)     => sam2paf::run(args),
    }
}
