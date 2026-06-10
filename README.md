<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/maligno-dark.png">
  <source media="(prefers-color-scheme: light)" srcset="docs/maligno-light.png">
  <img alt="maligno" src="docs/maligno-light.png" width="300">
</picture>

# maligno

**maligno** a unified toolkit for alignment processing and cross-file comparison. Subcommands cover the full pipeline from raw alignments (BAM/PAF) to per-read comparison statistics.

---

## Two ways to compare two PAFs

Both start from PAFs for sample A and B and produce the **identical** comparison
table. Pick based on whether you also want the intermediate per-file tables.

**1. One pass (primary) — `compare`:**

```
  A.paf ─┐
         ├─ compare ─▶ compare.tsv      maligno compare -a A.paf -b B.paf -o compare.tsv.gz
  B.paf ─┘
```

**2. Explicit intermediates — `paf2tables` then `compare-readinfo`:**

```
  A.paf ──paf2tables──▶ A.readinfo.tsv ─┐
                                         ├─ compare-readinfo ─▶ compare.tsv
  B.paf ──paf2tables──▶ B.readinfo.tsv ─┘
```
```bash
maligno paf2tables -i A.paf --readinfo A.readinfo.tsv.gz
maligno paf2tables -i B.paf --readinfo B.readinfo.tsv.gz
maligno compare-readinfo -a A.readinfo.tsv.gz -b B.readinfo.tsv.gz -o compare.tsv.gz
```

| Subcommand    | Input                              | Output                                  |
|---------------|------------------------------------|-----------------------------------------|
| **`compare`** | two PAFs (`-a`, `-b`)              | **per-read comparison TSV** (`-o`), one pass — **the primary entry point**; `--mode full` (default, 94 cols) or `--mode junctions` (47-col splice view) |
| `paf2tables`  | PAF (`-i`, `.gz`/`-` ok)           | **alninfo TSV** (`--alninfo`, 35 cols) and/or **readinfo TSV** (`--readinfo`, 33 cols), in one pass |
| `compare-readinfo` | two readinfo TSVs (`-a`, `-b`) | per-read comparison TSV (`-o`); same `--mode full`/`junctions` as `compare` |
| `sam2paf`     | SAM file or stdin (`-`)            | PAF written to stdout                   |
| `utils-readinfo` | alninfo TSV (`-i`, `.gz`/`-` ok) | per-read summary TSV (`-o`, 33 cols) — low-level utility; most users want `paf2tables --readinfo` |

All inputs/outputs transparently support gzip (`.gz` suffix) and stdin/stdout (`-`).

### Working with PAF files: `paf2tables` (start here)

`paf2tables` is the one command that turns a PAF into the downstream tables. Give
it the output path(s) you want — it writes the **alninfo** table, the **readinfo**
table, or **both in a single pass** (reading the PAF only once):

```bash
# Both tables in one pass (reads the PAF once):
maligno paf2tables -i in.paf --alninfo alninfo.tsv.gz --readinfo readinfo.tsv.gz

# Just the per-alignment table (no grouping required; works on unsorted PAF):
maligno paf2tables -i in.paf --alninfo alninfo.tsv.gz

# Just the per-read summary (best alignment per read):
maligno paf2tables -i in.paf --readinfo readinfo.tsv.gz
```

At least one of `--alninfo` / `--readinfo` must be given. The output is
**byte-identical** to the older two-command flow because every record is routed
through the exact same conversion + serialization code
(`AlnInfo::from_paf` → `write_row`) before the unchanged collapse logic runs;
when both outputs are requested the alninfo rows are simply *tee'd* off the same
stream that feeds readinfo grouping.

**Grouping requirement.** The `--readinfo` output requires the PAF be **grouped**
by `Query_Name` (each read's alignments contiguous); `--alninfo` never cares
about order. By default `paf2tables` groups contiguous runs and prints a one-shot
warning if the input isn't byte-lex sorted. Add **`--strict-grouping`** to turn a
silently-wrong non-grouped input into a hard error — it tracks completed read
names (memory ∝ number of distinct reads) and aborts if a `Query_Name` reappears
non-contiguously:

```bash
maligno paf2tables -i in.paf --readinfo readinfo.tsv.gz --strict-grouping
```

The standard way to guarantee grouping (and satisfy the downstream
`compare-readinfo` byte-lex requirement at the same time) is a single up-front sort:

```bash
LC_ALL=C sort -t$'\t' -k1,1 in.paf | maligno paf2tables -i - --readinfo readinfo.tsv.gz
```

#### Handling non-grouped input — current behavior and roadmap

readinfo correctness depends on every read's alignments being **contiguous**;
alninfo never cares. Current behavior + the roadmap for richer handling
(cheapest → most robust):

1. **Contiguous streaming + one-shot lex-decrease warning** *(current default,
   O(1) memory).* Fast, but silently wrong if a read's alignments are scattered
   and no lex-decrease trips the warning.
2. **`--strict-grouping` seen-set guard** *(now; O(#distinct reads) memory).*
   Hard-errors on non-contiguous reappearance — turns "silently wrong" into
   "loudly wrong" without buffering alignments.
3. **`--unsorted` full in-memory grouping** *(future; O(file) memory).* Buffer
   all alignments into a `HashMap<Query_Name, …>`, then collapse — order-
   independent; fine for small/medium PAFs.
4. **Internal external-sort fallback** *(future; bounded memory, any size).*
   Spill to temp files and merge-sort by `Query_Name` before grouping.
5. **Two-pass offset index on seekable input** *(future; bounded memory).* Pass 1
   indexes `Query_Name → byte offsets`; pass 2 seeks per read (regular files only).
6. **Actionable auto-detect error** *(future).* On detecting unsorted input, print
   the exact tailored `LC_ALL=C sort … | maligno paf2tables …` command to run.
7. **`--assume-grouped` fast-path** *(future).* Skip all checks for maximum
   throughput when the caller guarantees grouping.

### Primary: one-pass `compare`

`compare` collapses the whole pipeline into one streaming pass over the two PAFs:

```bash
maligno compare -a A.paf -b B.paf --label-a A --label-b B -o compare.tsv.gz
maligno compare -a A.paf -b B.paf --mode junctions -o compare_junctions.tsv.gz   # 47-col view
```

It is a strict lock-step **zip**: it requires both PAFs to list the **same
`Query_Names` in the same order**. Unlike `compare-readinfo` it does *not* require
byte-lex sorting — any shared ordering works. On the first read-name mismatch (or
if one side has extra reads) it prints an ERROR and stops; pass
`--ignore-row-mismatch` to skip reads present in only one PAF instead (requires
lex-sorted PAFs). Ideal for comparing two parameter sets / references run on the
**same** read or transcript set.

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

### Primary: one-pass `compare` from BAM/PAF

```bash
BIN=./target/release/maligno

# 0. BAM → PAF  (requires samtools; -h preserves @SQ header for contig lengths)
samtools view -h refA.bam | $BIN sam2paf -U - | gzip > refA.paf.gz
samtools view -h refB.bam | $BIN sam2paf -U - | gzip > refB.paf.gz

# 0.5. Pre-sort each PAF by Query_Name (byte-lex), so both list reads in the same order.
LC_ALL=C sort -t$'\t' -k1,1 <(gzcat refA.paf.gz) | gzip > refA.sorted.paf.gz
LC_ALL=C sort -t$'\t' -k1,1 <(gzcat refB.paf.gz) | gzip > refB.sorted.paf.gz

# 1. Compare the two PAFs in one pass (no intermediate files).
$BIN compare -a refA.sorted.paf.gz -b refB.sorted.paf.gz \
  --label-a RefA --label-b RefB -o RefA_vs_RefB.compare.tsv.gz
```

### Explicit intermediates: `paf2tables` then `compare-readinfo`

Same result, but materializes the per-file `readinfo` (and optionally `alninfo`)
tables for other analyses.

```bash
BIN=./target/release/maligno
# (steps 0 and 0.5 as above)

# 1. PAF → readinfo  (one row per read; best alignment chosen by ms, then AS, then MQ).
#    Add --alninfo <path> to also emit the 35-col per-alignment table in the same pass.
$BIN paf2tables -i refA.sorted.paf.gz --readinfo refA.readinfo.tsv.gz
$BIN paf2tables -i refB.sorted.paf.gz --readinfo refB.readinfo.tsv.gz

# 2. Compare the two readinfo files (streaming, constant memory; strict order by default).
$BIN compare-readinfo \
  -a refA.readinfo.tsv.gz --label-a RefA \
  -b refB.readinfo.tsv.gz --label-b RefB \
  -o RefA_vs_RefB.compare.tsv.gz
```

---

## How each step works

### `paf2tables` (alninfo conversion)

Parses each PAF record (12 mandatory fields + `ms:i`, `AS:i`, `cs:Z` tags), walks the
`cs` tag to accumulate match/substitution/insertion/deletion/splice statistics, computes
soft-clip lengths, junction coordinates (strand-aware), and derived scalars
(`seqid`, `Query_Aln_Len`, `Query_Aln_Cov`). This is the `--alninfo` output of
`paf2tables`.

**Pure streaming, constant memory** for the alninfo output. Each PAF line is parsed and
written independently in input order — no internal collect-then-sort. For the explicit
pipeline (`paf2tables` → `compare-readinfo`), pre-sort the PAF by `Query_Name` once upstream:

```bash
LC_ALL=C sort -t$'\t' -k1,1 in.paf > sorted.paf
```

Unix `sort` does external-sort with bounded memory and handles files larger than RAM.
The pre-sort satisfies both the `--readinfo` contiguity requirement and
`compare-readinfo`'s byte-lex sort requirement in one pass.

**Unaligned reads are kept.** A PAF record with `Target_Name == "*"` produces a full row
with zeroed alignment statistics, allowing unaligned reads to flow through the entire
pipeline.

### `utils-readinfo`

The per-read collapse step. **Most users don't call this directly** — `paf2tables --readinfo`
produces the identical readinfo table straight from a PAF. `utils-readinfo` is the low-level
utility for the case where you already have an `alninfo.tsv` (and not the PAF) and want to
collapse it to readinfo.

Groups alninfo rows by `Query_Name` (contiguous in sorted input) and collapses each group
to one summary row:

- **Best alignment** is the row with the highest `ms`, ties broken by highest `AS`, then by
  highest `MQ`. Full `(ms, AS, MQ)` ties fall through to alninfo input order (stable sort
  within each `Query_Name` run) — typically the aligner's emission order for that read.
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
- `MQ_Best` carries the mapping-quality (PAF col 12) of the best-scoring alignment — the
  same alignment from which `TargetChr`, `Strand`, `cs`, `junctions`, etc. are taken. For
  STAR-aligned data the common values are 255 (uniquely mapped), 3 (NH=2), 1 (NH=3), 0
  (NH>3). A difference in `MQ_Best` between two readinfo files surfaces reads where the two
  aligners (or parameter sets) disagree on mapping uniqueness.
- `Query_Start` / `Query_End` and `Target_Start` / `Target_End` carry the best alignment's
  query-coordinate span on the read and target-coordinate span on the reference (both
  0-based half-open, same convention as PAF / BED). Combined with `TargetChr` and `Strand`,
  this gives each read a complete BED-style alignment interval — useful for downstream
  genomic-region analysis (e.g., `bedtools merge` on filtered subsets of the compare output
  to identify regions where SetA and SetB differ).

### `compare-readinfo` (and the comparison core)

A two-pointer **merge-join** over two sorted readinfo files, matching on
**(Read_Name, Read_Len)**. This is the engine behind both `compare-readinfo` (readinfo
TSVs in) and the primary `compare` (which feeds it collapsed rows straight from PAFs).
For each matched read it emits the 31 data columns from each
side (suffixed with `--label-a` / `--label-b`) plus 30 comparison/object columns
(`AS_Diff`, `ms_Ratio`, `seqid_Diff`, `Junction_Distance`, `N_Matched_Junctions`, `Genomic_N_Matched_Junctions`, `Junctions_OnlyA`, …).

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

**Genomic-junction comparison.** When both alignments are to the *same* reference, junctions
are also compared in **reference coordinates** in addition to the query-coordinate metrics
above. These are **always emitted** (both `--mode full` and `--mode junctions`); four count
columns appear in the comparison block:

| Column                          | Meaning |
|---------------------------------|---------|
| `Genomic_N_Matched_Junctions`   | overlap on `(chrom, start, end)` sets |
| `Genomic_N_Unmatched_Junctions` | set symmetric difference (= `OnlyA + OnlyB`) |
| `Genomic_N_Junctions_OnlyA`     | only in A |
| `Genomic_N_Junctions_OnlyB`     | only in B |

The `genomic_junctions` column (always emitted in the alninfo and readinfo tables) uses
0-based half-open BED coordinates in **`((start, end), ...)`** Python-tuple-of-tuples
form, parseable with `ast.literal_eval`. The chromosome is **not** in each tuple — it's
in the sibling `TargetChr` (alninfo) / `TargetChr_<label>` (compare) column. Cross-
chromosome safety in the set comparison is still preserved: the comparison commands
reconstruct full `(chrom, start, end)` keys internally by combining
each row's parsed pairs with its per-side `TargetChr`, so junctions on different contigs
cannot accidentally match.

> **Format change (v0.2.3).** The genomic-junction tuples used to include the chrom as
> the first element (e.g. `(('chr22', 100, 250), …)`). That was redundant with the
> `TargetChr` column, so it was dropped. Pre-v0.2.3 TSVs need to be regenerated from PAF
> to be readable by `compare` / `compare-readinfo`.

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
always emitted in both `--mode full` and `--mode junctions`). These
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
> Also: if you have older readinfo TSVs with `TargetRef_1st`, rerun `paf2tables --readinfo`
> to get the renamed column (or rename in your scripts).

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
The `paf2tables` alninfo output is pure streaming and order-preserving (it does not sort
internally). The cleanest fix is a one-time pre-sort upstream on the PAF, which carries
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

The `--readinfo` collapse emits a one-time WARNING on stderr if its input is not byte-lex
sorted, flagging the most common foot-gun (a name-sorted-but-not-byte-lex aligner
output like STAR's, or a shuffled multi-threaded aligner output).

### `--mode junctions`

A **streamlined, splice-focused** view of the comparison (selected with `--mode junctions`
on either `compare` or `compare-readinfo`). Same metrics, but emits only **47 columns**
instead of 94 — useful when the question is "how do the splice junctions for each read
differ between two alignments?" rather than full score/indel/coverage diffs.

| Cols | Content |
|------|---------|
| 1–2   | `Read_Name`, `Read_Len` (join keys) |
| 3–32  | 15 per-side data columns × 2 sides: `TargetChr, Strand, MQ_Best, Num_Aln, Num_Aln_MaxScore, JuncCount, seqid_Max, Query_Aln_Cov_Max, junctions, genomic_junctions, cs, Query_Start, Query_End, Target_Start, Target_End` |
| 33–43 | 11 comparison metrics: `Strand_Match` + `seqid_Diff` + `QueryAlnCov_Diff` + 4 query-junction set metrics (matched / unmatched / OnlyA / OnlyB) + 4 parallel `Genomic_*` set metrics |
| 44–47 | 4 object lists at the end: `Junctions_OnlyA`, `Junctions_OnlyB`, `Genomic_Junctions_OnlyA`, `Genomic_Junctions_OnlyB` — the actual tuples of junctions that failed to overlap (Python tuple format, parseable with `ast.literal_eval`) |

Genomic-junction metrics are always emitted in this mode (no flag) — `chrom` is embedded in
each genomic-junction tuple, so cross-chromosome compares correctly produce zero overlap.

```bash
# From readinfo tables:
$BIN compare-readinfo --mode junctions \
  -a refA.readinfo.tsv.gz --label-a RefA \
  -b refB.readinfo.tsv.gz --label-b RefB \
  -o RefA_vs_RefB.junctions.tsv.gz

# …or directly from PAFs in one pass:
$BIN compare --mode junctions -a refA.sorted.paf.gz -b refB.sorted.paf.gz \
  --label-a RefA --label-b RefB -o RefA_vs_RefB.junctions.tsv.gz
```

The metric values that overlap with the `--mode full` output (the shared set columns) are
identical row-for-row — `--mode junctions` is purely a column-selection view of the same
computation, not a different algorithm.

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

# PAF → readinfo in one pass (add --alninfo <path> to also keep the per-alignment table).
time $BIN paf2tables -i test_data/Splice.AlnToHG38.PriAln.paf.gz   --readinfo /tmp/Splice.readinfo.tsv.gz
time $BIN paf2tables -i test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz --readinfo /tmp/SpliceHQ.readinfo.tsv.gz

# Sort each readinfo on (Read_Name, Read_Len) so the two files share byte-lex order
# (compare is strict by default — it errors on a read-name mismatch unless inputs match).
for S in Splice SpliceHQ; do
  ( zcat < /tmp/$S.readinfo.tsv.gz | { IFS= read -r h; printf '%s\n' "$h"; \
      LC_ALL=C sort -t$'\t' -k1,1 -k2,2n; } ) | gzip > /tmp/$S.readinfo.sorted.tsv.gz
done

time $BIN compare-readinfo \
  -a /tmp/Splice.readinfo.sorted.tsv.gz   --label-a Splice \
  -b /tmp/SpliceHQ.readinfo.sorted.tsv.gz --label-b SpliceHQ \
  -o /tmp/Splice_vs_SpliceHQ.compare.tsv.gz

# Inspect the compare output header (column number → column name)
zcat < /tmp/Splice_vs_SpliceHQ.compare.tsv.gz | head -1 | tr '\t' '\n' | nl



# Streamlined splice-focused comparison (47 cols: per-side junction info + alignment span + set-overlap metrics)
time $BIN compare-readinfo --mode junctions \
  -a /tmp/Splice.readinfo.sorted.tsv.gz   --label-a Splice \
  -b /tmp/SpliceHQ.readinfo.sorted.tsv.gz --label-b SpliceHQ \
  -o /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz

# Inspect the junction-view output header
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | head -1 | tr '\t' '\n' | nl

# Tip: to skip the intermediate files entirely, the primary `compare` runs the whole
# comparison in one pass (both PAFs must list reads in the same order):
#   $BIN compare -a sorted_Splice.paf.gz -b sorted_SpliceHQ.paf.gz \
#     --label-a Splice --label-b SpliceHQ -o /tmp/Splice_vs_SpliceHQ.compare.tsv.gz

# Check all unique values in columns 25 and 29 (Checking number of unmatched junctions from query and genome perspective)
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | cut -f 35 | sort | uniq -c 
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | cut -f 39 | sort | uniq -c 


# Look at number of reads with non-concordant junction positions (query)
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $35 > 0' | wc -l 
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $35 > 0' | cut -f 1,9,18,24,25,26,27,32,33 | less -S


zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $35 > 0' | cut -f 1,3,4,12,13,9,18,24,25,26,27,32,33 | less -S


zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $35 > 0' | cut -f 1,3,10,13,14,34,35,36,42,43,44,45 | column -t -s $'\t' | less -S






zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | cut -f 35 | sort | uniq -c 

# Look at number of reads with non-concordant junction positions (GENOMIC)
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $39 > 0' | wc -l 
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $39 > 0' | cut -f 1,10,18,28,29,30,31,34,35 | column -t | less -S

zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' 'NR==1 || $39 > 0' | cut -f 1,3,4,12,13,29,30,31,34,35 | column -t | less -S



# Verify column counts (expect 35, 32, 86, 45)
zcat < /tmp/Splice.alninfo.tsv.gz                       | awk -F'\t' '{print NF}' | sort | uniq -c
zcat < /tmp/Splice.readinfo.tsv.gz                      | awk -F'\t' '{print NF}' | sort | uniq -c
zcat < /tmp/Splice_vs_SpliceHQ.compare.tsv.gz           | awk -F'\t' '{print NF}' | sort | uniq -c
zcat < /tmp/Splice_vs_SpliceHQ.compare_junctions.tsv.gz | awk -F'\t' '{print NF}' | sort | uniq -c
```

---

## Troubleshooting

### "`compare-readinfo` matched far fewer reads than I expected" — sort-order diagnostic

`compare-readinfo` (both `--mode` views) uses a streaming merge-join keyed on
`(Read_Name, Read_Len)`. The algorithm runs in O(1) memory and O(N+M) time,
but it **assumes both readinfo files are sorted in the same byte-lexicographic
order**. (The primary `compare` reads PAFs directly and instead requires both
PAFs to list reads in the *same order* — see its `--ignore-row-mismatch`.) By default a read-name mismatch is a hard error (non-zero exit) so you
notice immediately; passing `--ignore-row-mismatch` reverts to skip-and-count,
where unmatched reads are dropped and tallied in the end-of-run stderr summary.

The `paf2tables` alninfo output is pure streaming and **preserves input order** —
sortedness must be supplied by you upstream of the pipeline (or recovered
afterward, see Sorting requirement above). The `--readinfo` collapse emits a
one-time WARNING on stderr if it detects a byte-lex decrease in its input's
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

  ⇒ `compare-readinfo`'s 'matched' count should equal 1023. If maligno's reported
    matched count is lower than that, the two readinfo files are not sorted
    in the same byte-lexicographic order. Re-sort each with:

        LC_ALL=C sort -t$'\t' -k1,1 -k2,2n input.readinfo.tsv > sorted.tsv
        # (keep the header separately)

    or pre-sort the PAF once (`LC_ALL=C sort -t$'\t' -k1,1`) and rerun
    `paf2tables --readinfo` — both readinfo files then share byte-lex order.
```

If `intersection by Name only` is larger than `intersection by (Name, Len)`,
some reads share names but have different `Read_Len` between the two files
(usually an upstream soft-clip / qlen difference). Those rows can't match
in the comparison regardless of sort order.

---

## Source layout

```
src/
├── main.rs                 — CLI dispatcher (clap subcommands)
├── pafcompare.rs           — PRIMARY `compare` command: paired-PAF → comparison (one pass; lock-step zip)
├── paf2tables.rs           — PAF → alninfo and/or readinfo (one pass)
├── compare_streaming.rs    — `compare-readinfo` command + shared comparison core (emit/header/ReadKey/CompareMode)
├── compare_junctions.rs    — junction-view (47-col) header/row emitters (library; used by --mode junctions)
├── readinfo.rs             — collapse library (collapse_group/ReadInfoRow/AlnRow) + utils-readinfo entry point
├── paf_groups.rs           — shared PAF → per-read group reader, with optional alninfo tee
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
