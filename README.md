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
| `readinfo`    | alninfo TSV (`-i`, `.gz`/`-` ok)   | per-read summary TSV (`-o`, 33 cols)    |
| `compare`     | two readinfo TSVs (`-a`, `-b`)     | per-read comparison TSV (`-o`, 86 cols, +6 with `--compare-genomic-junctions`) |
| `compare-junctions` | two readinfo TSVs (`-a`, `-b`) | streamlined splice-focused comparison TSV (`-o`, 45 cols) |
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

# 0.5. Pre-sort each PAF by Query_Name (byte-lex). Skip if already byte-lex sorted.
#      Required so downstream readinfo groups correctly and compare's merge-join
#      sees byte-lex-sorted input. Unix `sort` uses external-sort → bounded memory.
LC_ALL=C sort -t$'\t' -k1,1 <(gzcat refA.paf.gz) | gzip > refA.sorted.paf.gz
LC_ALL=C sort -t$'\t' -k1,1 <(gzcat refB.paf.gz) | gzip > refB.sorted.paf.gz

# 1. PAF → alninfo  (one row per alignment, 35 cols; pure streaming, constant memory)
$BIN paf2alninfo -i refA.sorted.paf.gz -o refA.alninfo.tsv.gz
$BIN paf2alninfo -i refB.sorted.paf.gz -o refB.alninfo.tsv.gz

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
(`seqid`, `Query_Aln_Len`, `Query_Aln_Cov`).

**Pure streaming, constant memory** (since v0.2.7). Each PAF line is parsed and written
independently in input order — no internal collect-then-sort. For the standard pipeline
(`paf2alninfo` → `readinfo` → `compare`), pre-sort the PAF by `Query_Name` once upstream:

```bash
LC_ALL=C sort -t$'\t' -k1,1 in.paf > sorted.paf
```

Unix `sort` does external-sort with bounded memory and handles files larger than RAM.
The pre-sort satisfies both `readinfo`'s contiguity requirement and `compare`'s byte-lex
sort requirement in one pass.

**Unaligned reads are kept.** A PAF record with `Target_Name == "*"` produces a full row
with zeroed alignment statistics, allowing unaligned reads to flow through the entire
pipeline.

### `readinfo`

Groups alninfo rows by `Query_Name` (contiguous in sorted input) and collapses each group
to one summary row:

- **Best alignment** is the row with the highest `ms`, ties broken by highest `AS`. Full
  `(ms, AS)` ties fall through to alninfo input order (stable sort within each `Query_Name`
  run) — typically the aligner's emission order for that read.
- **Aggregates over all alignments of the read:** `AS_Max`, `ms_Max`, `Query_Aln_Cov_Max`,
  `Query_Aln_Len_Max`, `seqid_Max`.
- `Num_Aln` counts only aligned rows (`Target_Name != "*"`), so a read that is present but
  entirely unaligned gets a row with `Num_Aln = 0` and zeroed stats.
- `Num_Aln_MaxScore` counts alignments tied at the chosen-best sort key for this read —
  i.e., tied at **both** `ms_Max` **and** the highest `AS` among ms-tied rows. This matches
  the full `(ms desc, AS desc)` selection rule used to pick the best alignment.
  `Num_Aln_MaxScore = 1` ⇒ a single unambiguous winner under the selection rule;
  `> 1` ⇒ alignments remain indistinguishable on both `ms` and `AS`, and file order
  broke the tie. Practical note: STAR-style aligners write `ms=0` for every alignment, so
  `ms_Max = 0` and `AS` does the actual selection work — counting at `(ms, AS)` keeps
  `Num_Aln_MaxScore` informative in that case (otherwise it would equal `Num_Aln`).
- `Query_Start` / `Query_End` and `Target_Start` / `Target_End` carry the best alignment's
  query-coordinate span on the read and target-coordinate span on the reference (both
  0-based half-open, same convention as PAF / BED). Combined with `TargetChr` and `Strand`,
  this gives each read a complete BED-style alignment interval — useful for downstream
  genomic-region analysis (e.g., `bedtools merge` on filtered subsets of the compare output
  to identify regions where SetA and SetB differ).

### `compare`

A two-pointer **merge-join** over two sorted readinfo files, matching on
**(Read_Name, Read_Len)**. For each matched read it emits the 30 data columns from each
side (suffixed with `--label-a` / `--label-b`) plus 24 comparison/object columns
(`AS_Diff`, `ms_Ratio`, `seqid_Diff`, `Junction_Distance`, `N_Matched_Junctions`, `Junctions_OnlyA`, …).

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

**Genomic-junction comparison (`--compare-genomic-junctions`).** When both alignments are to
the *same* reference, you can compare junctions in **reference coordinates** instead of (or
in addition to) the query-coordinate metrics above. Pass `--compare-genomic-junctions` to
`compare` and four additional columns appear at the end of the row:

| Column                          | Meaning |
|---------------------------------|---------|
| `Genomic_N_Matched_Junctions`   | overlap on `(chrom, start, end)` sets |
| `Genomic_N_Unmatched_Junctions` | set symmetric difference (= `OnlyA + OnlyB`) |
| `Genomic_N_Junctions_OnlyA`     | only in A |
| `Genomic_N_Junctions_OnlyB`     | only in B |

The `genomic_junctions` column (always emitted by `paf2alninfo` and `readinfo`) uses
0-based half-open BED coordinates in **`((start, end), ...)`** Python-tuple-of-tuples
form, parseable with `ast.literal_eval`. The chromosome is **not** in each tuple — it's
in the sibling `TargetChr` (alninfo) / `TargetChr_<label>` (compare) column. Cross-
chromosome safety in the set comparison is still preserved: `compare` and
`compare-junctions` reconstruct full `(chrom, start, end)` keys internally by combining
each row's parsed pairs with its per-side `TargetChr`, so junctions on different contigs
cannot accidentally match.

> **Format change (v0.2.3).** The genomic-junction tuples used to include the chrom as
> the first element (e.g. `(('chr22', 100, 250), …)`). That was redundant with the
> `TargetChr` column, so it was dropped. Pre-v0.2.3 TSVs need to be regenerated from PAF
> to be readable by `compare` / `compare-junctions`.

**Strand tracking and renames (v0.2.1+).** Each side now carries a `Strand_<label>` data
column (the best alignment's strand), and the comparison block starts with a `Strand_Match`
(true/false) metric that flags strand-flips between A and B. The legacy column name
`TargetRef_1st` has been renamed to `TargetChr` (suffixed in compare output as
`TargetChr_<label>`).

**Non-overlap junction objects (v0.2.2+).** In addition to the *counts* of non-overlapping
junctions (`N_Junctions_OnlyA/B`, `Genomic_N_Junctions_OnlyA/B`), the comparison outputs
now append the actual junction **objects** that failed to overlap at the very end of each
row: `Junctions_OnlyA`, `Junctions_OnlyB` (query-coord tuples, always emitted) and
`Genomic_Junctions_OnlyA`, `Genomic_Junctions_OnlyB` (genome-coord tuples, emitted with
`--compare-genomic-junctions` in `compare`; always emitted in `compare-junctions`). These
use the same Python tuple format as the per-side `junctions` / `genomic_junctions` data
columns — parse with `ast.literal_eval` in Python.

> **Migration note (column positions).** Schema growth has shifted column positions twice
> in pre-1.0 development: first when `genomic_junctions` was added, then again when
> `Strand` and `Strand_Match` were added (v0.2.1). Scripts that filter by column *number*
> (`awk '$73 > 0'`) need updating each time; prefer column-*name* lookup using the header,
> which is robust to future schema growth:
> ```bash
> awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $c["N_Unmatched_Junctions"]>0'
> ```
> Also: if you have older readinfo TSVs with `TargetRef_1st`, rerun `readinfo` to get the
> renamed column (or rename in your scripts).

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

**Sorting requirement.** The merge-join requires both readinfo files to be in the
**same** `(Read_Name, Read_Len)` order — byte-lex order is the most convenient guarantee.
As of v0.2.7, `paf2alninfo` no longer sorts internally (it's pure streaming, order-
preserving). The cleanest fix is a one-time pre-sort upstream on the PAF, which carries
through the entire chain:

```bash
# Pre-sort PAF (plain or gzipped) by Query_Name.
LC_ALL=C sort -t$'\t' -k1,1 in.paf > sorted.paf
zcat in.paf.gz | LC_ALL=C sort -t$'\t' -k1,1 | gzip > sorted.paf.gz
```

If you already have an unsorted **alninfo** TSV, sort it on `Query_Name` (col 1) with the
header kept separately:

```bash
# alninfo: plain
(head -1 in.alninfo.tsv;
 tail -n +2 in.alninfo.tsv | LC_ALL=C sort -t$'\t' -k1,1) \
  > sorted.alninfo.tsv

# alninfo: gzipped
( zcat in.alninfo.tsv.gz | head -1;
  zcat in.alninfo.tsv.gz | tail -n +2 | LC_ALL=C sort -t$'\t' -k1,1
) | gzip > sorted.alninfo.tsv.gz
```

If you already have an unsorted **readinfo** TSV, sort it on `(Read_Name, Read_Len)` —
col 1 (string) then col 2 (numeric):

```bash
# readinfo: plain
(head -1 in.readinfo.tsv;
 tail -n +2 in.readinfo.tsv | LC_ALL=C sort -t$'\t' -k1,1 -k2,2n) \
  > sorted.readinfo.tsv

# readinfo: gzipped
( zcat in.readinfo.tsv.gz | head -1;
  zcat in.readinfo.tsv.gz | tail -n +2 | LC_ALL=C sort -t$'\t' -k1,1 -k2,2n
) | gzip > sorted.readinfo.tsv.gz
```

Notes on the sort flags:
- `LC_ALL=C` forces byte-lex order (locale-independent and deterministic).
- `-t$'\t'` sets the field separator to TAB.
- `-k1,1` sorts on column 1 as a string (`Query_Name` for alninfo, `Read_Name` for readinfo).
- `-k2,2n` (readinfo only) breaks ties by `Read_Len` numerically.
- `sort` uses external-sort under the hood, so memory stays bounded even on files larger
  than RAM. Override its scratch directory and memory cap with `-T` and `-S` if needed
  (e.g. `-T /scratch -S 8G`).

`readinfo` emits a one-time WARNING on stderr if its input alninfo is not byte-lex
sorted, flagging the most common foot-gun (a name-sorted-but-not-byte-lex aligner
output like STAR's, or a shuffled multi-threaded aligner output).

### `compare-junctions`

A **streamlined, splice-focused** variant of `compare`. Same streaming merge-join over two
readinfo files, but emits only **45 columns** instead of 86–92 — useful when the question is
"how do the splice junctions for each read differ between two alignments?" rather than full
score/indel/coverage diffs.

| Cols | Content |
|------|---------|
| 1–2   | `Read_Name`, `Read_Len` (join keys) |
| 3–30  | 14 per-side data columns × 2 sides: `TargetChr, Strand, Num_Aln, Num_Aln_MaxScore, JuncCount, seqid_Max, Query_Aln_Cov_Max, junctions, genomic_junctions, cs, Query_Start, Query_End, Target_Start, Target_End` |
| 31–41 | 11 comparison metrics: `Strand_Match` + `seqid_Diff` + `QueryAlnCov_Diff` + 4 query-junction set metrics (matched / unmatched / OnlyA / OnlyB) + 4 parallel `Genomic_*` set metrics |
| 42–45 | 4 object lists at the end: `Junctions_OnlyA`, `Junctions_OnlyB`, `Genomic_Junctions_OnlyA`, `Genomic_Junctions_OnlyB` — the actual tuples of junctions that failed to overlap (Python tuple format, parseable with `ast.literal_eval`) |

Genomic-junction metrics are always emitted (no flag) — `chrom` is embedded in each
genomic-junction tuple, so cross-chromosome compares correctly produce zero overlap.

```bash
$BIN compare-junctions \
  -a refA.readinfo.tsv.gz --label-a RefA \
  -b refB.readinfo.tsv.gz --label-b RefB \
  -o RefA_vs_RefB.junctions.tsv.gz
```

The metric values that overlap with the full `compare` output (the 8 set columns) are
identical row-for-row — `compare-junctions` is purely a column-selection variant, not a
different algorithm.

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

`test_data/` contains two Chr22-scale PAFs (gzipped, ~0.5 MB each) for end-to-end testing.

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

# Inspect the compare output header (column number → column name)
zcat < /tmp/Splice_vs_SpliceHQ.compare.tsv.gz | head -1 | tr '\t' '\n' | nl



# Streamlined splice-focused comparison (45 cols: per-side junction info + alignment span + set-overlap metrics)
time $BIN compare-junctions \
  -a /tmp/Splice.readinfo.tsv.gz   --label-a Splice \
  -b /tmp/SpliceHQ.readinfo.tsv.gz --label-b SpliceHQ \
  -o /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz

# Inspect the compare-junctions output header
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | head -1 | tr '\t' '\n' | nl

# Check all unique values in columns 25 and 29 (Checking number of unmatched junctions from query and genome perspective)
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | cut -f 25 | sort | uniq -c 
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | cut -f 29 | sort | uniq -c 


# Look at number of reads with non-concordant junction positions (query)
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $25 > 0' | wc -l 
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $25 > 0' | cut -f 1,9,18,24,25,26,27,32,33 | less -S


zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $25 > 0' | cut -f 1,3,4,12,13,9,18,24,25,26,27,32,33 | less -S


zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $25 > 0' | cut -f 1,3,4,12,13,9,18,24,25,26,27,32,33,34,35 | less -S


zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | cut -f 25 | sort | uniq -c 

# Look at number of reads with non-concordant junction positions (GENOMIC)
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $29 > 0' | wc -l 
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $29 > 0' | cut -f 1,9,18,28,29,30,31,34,35 | less -S

zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $29 > 0' | cut -f 1,3,4,12,13,29,30,31,34,35 | less -S



# Verify column counts (expect 35, 32, 86, 45)
zcat < /tmp/Splice.alninfo.tsv.gz                       | awk -F'\t' '{print NF}' | sort | uniq -c
zcat < /tmp/Splice.readinfo.tsv.gz                      | awk -F'\t' '{print NF}' | sort | uniq -c
zcat < /tmp/Splice_vs_SpliceHQ.compare.tsv.gz           | awk -F'\t' '{print NF}' | sort | uniq -c
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' '{print NF}' | sort | uniq -c
```

---

## Troubleshooting

### "`compare` matched far fewer reads than I expected" — sort-order diagnostic

`compare` and `compare-junctions` use a streaming merge-join keyed on
`(Read_Name, Read_Len)`. The algorithm runs in O(1) memory and O(N+M) time,
but it **assumes both readinfo files are sorted in the same byte-lexicographic
order**. If they aren't, matches are silently missed and the `matched` count
in the end-of-run stderr summary comes out lower than it should.

Since v0.2.7 `paf2alninfo` is pure streaming and **preserves input order** —
sortedness must be supplied by you upstream of the pipeline (or recovered
afterward, see Sorting requirement above). `readinfo` emits a one-time
WARNING on stderr if it detects a byte-lex decrease in its input alninfo's
`Query_Name` column, flagging the common foot-gun (a name-sorted-but-not-byte-lex
aligner output like STAR's, or a shuffled multi-threaded aligner output).

This diagnostic still bites you if you re-sorted a file
externally (e.g. with a non-`LC_ALL=C` locale) or assembled the inputs from
multiple sources.

To check, compute the expected match count via a set intersection of the
`(Read_Name, Read_Len)` keys. If the merge-join's `matched` count matches
this number, you're fine — any low overlap is a property of the input data.
If it's lower, the inputs aren't sorted consistently.

**Inline one-liner** (works on plain or `.gz` readinfo TSVs):

```bash
# Expected match count = full-key set intersection
LC_ALL=C comm -12 \
  <(zcat -f a.readinfo.tsv.gz | tail -n +2 | awk -F'\t' '{print $1"\t"$2}' | LC_ALL=C sort -u) \
  <(zcat -f b.readinfo.tsv.gz | tail -n +2 | awk -F'\t' '{print $1"\t"$2}' | LC_ALL=C sort -u) \
  | wc -l
```

**Or use the bundled diagnostic script** at
[`scripts/check-readinfo-overlap.sh`](scripts/check-readinfo-overlap.sh) — same
calculation plus a Name-only intersection (to spot reads that share names but
differ in `Read_Len`, e.g. due to soft-clip differences), and a printable
report:

```bash
./scripts/check-readinfo-overlap.sh a.readinfo.tsv.gz b.readinfo.tsv.gz
```

Example output:

```
Sort/overlap diagnostic
  A: a.readinfo.tsv.gz
  B: b.readinfo.tsv.gz

  rows in A (unique full-key):     466392
  rows in B (unique full-key):     496164
  intersection by (Name, Len):     1023
  intersection by Name only:       1023

  ⇒ `compare`'s 'matched' count should equal 1023. If maligno's reported
    matched count is lower than that, the two readinfo files are not sorted
    in the same byte-lexicographic order. Re-sort each with:

        LC_ALL=C sort -t$'\t' -k1,1 -k2,2n input.readinfo.tsv > sorted.tsv
        # (keep the header separately)

    or rerun the maligno pipeline from PAF level — paf2alninfo + readinfo
    produce consistently-sorted output by construction.
```

If `intersection by Name only` is larger than `intersection by (Name, Len)`,
some reads share names but have different `Read_Len` between the two files
(usually an upstream soft-clip / qlen difference). Those rows can't match
in `compare` regardless of sort order.

---

## Source layout

```
src/
├── main.rs                 — CLI dispatcher (clap subcommands)
├── paf2alninfo.rs          — PAF → per-alignment info TSV
├── readinfo.rs             — alninfo → per-read summary TSV
├── compare_streaming.rs    — streaming merge-join comparison (full)
├── compare_junctions.rs    — streamlined splice-focused comparison (45 cols)
├── record.rs               — AlnInfo struct + TSV serialisation
├── paf.rs                  — PAF record parser
├── cs_parser.rs            — cs-tag parser (PAF → stats + genomic junctions)
├── cigar_junctions.rs      — CIGAR-based intron extractor (utility, not yet wired in)
├── io_utils.rs             — open_input / open_output (gzip transparent)
├── junction.rs             — junction parsers + set-overlap stats
└── sam2paf/
    ├── mod.rs              — sam2paf CLI args + run()
    ├── convert.rs          — SAM → PAF conversion logic
    ├── cigar.rs            — CIGAR string parser
    ├── md.rs               — MD-tag iterator
    └── cs_generator.rs     — cs-tag generator (MD + CIGAR → cs string)

scripts/
└── check-readinfo-overlap.sh   — sort-order / overlap diagnostic for compare inputs
```
