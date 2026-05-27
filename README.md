<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/maligno-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="docs/maligno-light.png">
  <img alt="maligno" src="docs/maligno-light.png" width="300">
</picture>

# maligno

**maligno** a unified toolkit for alignment processing and cross-file comparison. Subcommands cover the full pipeline from raw alignments (BAM/PAF) to per-read comparison statistics.

---

## Pipeline

```
  SAM/BAM ──sam2paf──▶ PAF ──paf2alninfo──▶ alninfo.tsv ──readinfo──▶ readinfo.tsv ─┐
                            (ref A)                                                 │
                                                                                    ├─ compare ─▶ compare.tsv
  SAM/BAM ──sam2paf──▶ PAF ──paf2alninfo──▶ alninfo.tsv ──readinfo──▶ readinfo.tsv ─┘            
                            (ref B)                                              
```

| Subcommand    | Input                              | Output                                  |
|---------------|------------------------------------|-----------------------------------------|
| `paf2alninfo` | PAF (`-i`, `.gz`/`-` ok)           | per-alignment info TSV (`-o`, 35 cols)  |
| `readinfo`    | alninfo TSV (`-i`, `.gz`/`-` ok)   | per-read summary TSV (`-o`, 28 cols)    |
| `compare`     | two readinfo TSVs (`-a`, `-b`)     | per-read comparison TSV (`-o`, 77 cols) |
| `sam2paf`     | SAM file or stdin (`-`)            | PAF written to stdout                   |

All inputs/outputs transparently support gzip (`.gz` suffix) and stdin/stdout (`-`).

---

## Build

```bash
# Native (macOS / Linux host)
cargo build --release
# → target/release/maligno

# Static Linux binary for HPC (no runtime deps)
cargo build --release --target x86_64-unknown-linux-musl
# → target/x86_64-unknown-linux-musl/release/maligno
```

---

## Install

Build and copy the `maligno` binary into a `bin/` directory of your choice:

```bash
# Default: installs to ~/.cargo/bin/maligno (usually already on PATH)
cargo install --path .

# Or choose the install root — cargo appends bin/ automatically
cargo install --path . --root ~/.local      # → ~/.local/bin/maligno

# Add --force to overwrite a previous install when rebuilding
cargo install --path . --root ~/.local --force
```

Ensure the target `bin/` directory is on your `PATH` (e.g. add
`export PATH="$HOME/.local/bin:$PATH"` to your shell rc) so you can run `maligno`
from anywhere. `cargo install` builds for the host machine; for the HPC static
binary use the musl cross-build above.

---

## Usage

### Full pipeline from BAM

```bash
BIN=./target/release/maligno

# 0. BAM → PAF  (requires samtools; -h preserves @SQ header for contig lengths)
samtools view -h refA.bam | $BIN sam2paf -U - | gzip > refA.paf.gz
samtools view -h refB.bam | $BIN sam2paf -U - | gzip > refB.paf.gz

# 1. PAF → alninfo  (one row per alignment, 35 cols)
$BIN paf2alninfo -i refA.paf.gz -o refA.alninfo.tsv.gz
$BIN paf2alninfo -i refB.paf.gz -o refB.alninfo.tsv.gz

# 2. alninfo → readinfo  (one row per read; best alignment chosen by ms, then AS)
$BIN readinfo -i refA.alninfo.tsv.gz -o refA.readinfo.tsv.gz
$BIN readinfo -i refB.alninfo.tsv.gz -o refB.readinfo.tsv.gz

# 3. Compare the two readinfo files (streaming merge-join, constant memory)
$BIN compare \
  -a refA.readinfo.tsv.gz --label-a RefA \
  -b refB.readinfo.tsv.gz --label-b RefB \
  -o RefA_vs_RefB.compare.tsv.gz
```

### Starting from PAF directly

```bash
BIN=./target/release/maligno

$BIN paf2alninfo -i refA.paf -o refA.alninfo.tsv.gz
$BIN readinfo    -i refA.alninfo.tsv.gz -o refA.readinfo.tsv.gz
# ... same compare step as above
```

---

## How each step works

### `paf2alninfo`

Parses each PAF record (12 mandatory fields + `ms:i`, `AS:i`, `cs:Z` tags), walks the
`cs` tag to accumulate match/substitution/insertion/deletion/splice statistics, computes
soft-clip lengths, junction coordinates (strand-aware), and derived scalars
(`seqid`, `Query_Aln_Len`, `Query_Aln_Cov`). Output is sorted by
(Query_Name, Query_Start, Query_End) unless `--no-sort`.

**Unaligned reads are kept.** A PAF record with `Target_Name == "*"` produces a full row
with zeroed alignment statistics, allowing unaligned reads to flow through the entire
pipeline.

### `readinfo`

Groups alninfo rows by `Query_Name` (contiguous in sorted input) and collapses each group
to one summary row:

- **Best alignment** is the row with the highest `ms`, ties broken by highest `AS`.
- **Aggregates over all alignments of the read:** `AS_Max`/`AS_Min`, `ms_Max`/`ms_Min`,
  `Query_Aln_Cov_Max`, `Query_Aln_Len_Max`, `seqid_Max`.
- `Num_Aln` counts only aligned rows (`Target_Name != "*"`), so a read that is present but
  entirely unaligned gets a row with `Num_Aln = 0` and zeroed stats.

### `compare`

A two-pointer **merge-join** over two sorted readinfo files, matching on
**(Read_Name, Read_Len)**. For each matched read it emits the 26 data columns from each
side (suffixed with `--label-a` / `--label-b`) plus 23 comparison metrics
(`AS_Diff`, `ms_Ratio`, `seqid_Diff`, `Junction_Distance`, `N_Matched_Junctions`, …).

**Junction set comparison.** Junctions are compared as **sets** of query coordinates
(deduplicated on both sides):

| Column | Meaning |
|--------|---------|
| `N_Matched_Junctions`    | size of the overlap, `\|A ∩ B\|` |
| `N_Junctions_OnlyA`      | junctions found only in A, `\|A \ B\|` |
| `N_Junctions_OnlyB`      | junctions found only in B, `\|B \ A\|` |
| `N_Unmatched_Junctions`  | junctions **not** in the overlap (set symmetric difference, `OnlyA + OnlyB`) |

These stay internally consistent: `N_Matched_Junctions + N_Junctions_OnlyA` equals the
junction count of A, and likewise for B. (`Junction_Distance` and `Junc_Dist_V2` are
retained positional/legacy metrics.)

**Join semantics — inner join.** Only reads present in **both** files produce an output
row. Reads present in only one file are dropped but counted in the end-of-run summary
printed to stderr:

```
Read comparison summary:
  Label A: Splice
  Label B: SpliceHQ
  rows in A (readinfo-a):     505748
  rows in B (readinfo-b):     505749
  matched (in both, written): 505742
  A-only (dropped, not in B): 6
  B-only (dropped, not in A): 7
```

**Sorting requirement.** The merge-join requires both readinfo files to be sorted by
(Read_Name, Read_Len). This is automatically satisfied when files come from this tool's
`paf2alninfo → readinfo` chain. If sorting externally:

```bash
LC_ALL=C sort -t$'\t' -k1,1 -k2,2n input.readinfo.tsv   # keep header separately
```

### `sam2paf`

Converts SAM alignments to PAF format. A high-performance port of the `sam2paf`
sub-command from paftools.js — output is byte-for-byte compatible.

Key flags:

| Flag | Meaning |
|------|---------|
| `-U` | Emit placeholder PAF records for unmapped reads (recommended for full pipeline) |
| `-p` | Primary + supplementary alignments only (skip secondary FLAG 0x100) |
| `-P` | Primary alignments only (skip secondary and supplementary) |
| `-L` | Output cs tag in long form (`=ACGT` instead of `:N`) |

> **Note:** Pass `-U` to keep unaligned reads in the pipeline (they become `Num_Aln = 0`
> rows in readinfo rather than disappearing entirely).

---


## Test data

`test_data/` contains two Chr22-scale PAFs (gzipped, ~3.3 MB each) for end-to-end testing.

```bash
BIN=./target/release/maligno

time $BIN paf2alninfo -i test_data/Splice.AlnToHG38.PriAln.paf.gz   -o /tmp/Splice.alninfo.tsv.gz
time $BIN paf2alninfo -i test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz -o /tmp/SpliceHQ.alninfo.tsv.gz

time $BIN readinfo -i /tmp/Splice.alninfo.tsv.gz   -o /tmp/Splice.readinfo.tsv.gz
time $BIN readinfo -i /tmp/SpliceHQ.alninfo.tsv.gz -o /tmp/SpliceHQ.readinfo.tsv.gz

time $BIN compare \
  -a /tmp/Splice.readinfo.tsv.gz   --label-a Splice \
  -b /tmp/SpliceHQ.readinfo.tsv.gz --label-b SpliceHQ \
  -o /tmp/Splice_vs_SpliceHQ.compare.tsv.gz

# Look at number of reads with non-concordant junction positions (query)
zcat < /tmp/Splice_vs_SpliceHQ.compare.tsv.gz | awk -F'\t' 'NR==1 || $73 > 0' | wc -l 

# Look at columns related to splice junction info (A vs B)
zcat < /tmp/Splice_vs_SpliceHQ.compare.tsv.gz | awk -F'\t' 'NR==1 || $73 > 0' | cut -f 1,2,13,39,11,37,73  | less -S


# Verify column counts (expect 35, 28, 77)
zcat < /tmp/Splice.alninfo.tsv.gz              | awk -F'\t' '{print NF}' | sort | uniq -c
zcat < /tmp/Splice.readinfo.tsv.gz             | awk -F'\t' '{print NF}' | sort | uniq -c
zcat < /tmp/Splice_vs_SpliceHQ.compare.tsv.gz | awk -F'\t' '{print NF}' | sort | uniq -c
```

---

## Source layout

```
src/
├── main.rs                 — CLI dispatcher (clap subcommands)
├── paf2alninfo.rs          — PAF → per-alignment info TSV
├── readinfo.rs             — alninfo → per-read summary TSV
├── compare_streaming.rs    — streaming merge-join comparison
├── record.rs               — AlnInfo struct + TSV serialisation
├── paf.rs                  — PAF record parser
├── cs_parser.rs            — cs-tag parser (PAF → stats)
├── io_utils.rs             — open_input / open_output (gzip transparent)
├── junction.rs             — junction coordinate utilities
└── sam2paf/
    ├── mod.rs              — sam2paf CLI args + run()
    ├── convert.rs          — SAM → PAF conversion logic
    ├── cigar.rs            — CIGAR string parser
    ├── md.rs               — MD-tag iterator
    └── cs_generator.rs     — cs-tag generator (MD + CIGAR → cs string)
```
