# maligno — Reference

Full reference for **maligno**, a unified toolkit for alignment processing and
cross-file comparison. Subcommands cover the full pipeline from raw alignments
(BAM/PAF) to per-read comparison statistics.

> New here? Start with the [README](../README.md) quick start. This document is
> the complete reference: every subcommand, the full column dictionary, schema
> migration notes, troubleshooting, and the source layout.

---

## Two ways to compare two PAFs

Both start from PAFs for sample A and B and produce the **identical** comparison
table. Use `compare` for the safe, do-it-all path; drop to the subcommands when
you want manual control.

**1. On-rails (primary) — `compare`:** one command that sorts both PAFs (so
order is guaranteed), verifies they carry the same read-ID
set, and writes a results directory with the per-set alninfo + readinfo tables
and the comparison table.

```
  A.paf ─┐                                            results/AvsB.{A,B}.alninfo.tsv.gz
         ├─ compare (sort → verify → tables) ─▶       results/AvsB.{A,B}.readinfo.tsv.gz
  B.paf ─┘                                            results/AvsB.compare.tsv.gz
```
```bash
maligno compare -a A.paf -b B.paf --label-a A --label-b B --outdir results/ --prefix AvsB
```

**2. Manual building blocks — `compare-pipeline paf2tables` then
`compare-pipeline merge-readinfo`:** full control (you sort and pick exactly
what to produce).

```
  A.paf ──paf2tables──▶ A.readinfo.tsv ─┐
                                         ├─ merge-readinfo ─▶ compare.tsv
  B.paf ──paf2tables──▶ B.readinfo.tsv ─┘
```
```bash
LC_ALL=C sort -t$'\t' -k1,1 <(gzcat A.paf.gz) | gzip > A.sorted.paf.gz   # (and B)
maligno compare-pipeline paf2tables -i A.sorted.paf.gz --readinfo A.readinfo.tsv.gz
maligno compare-pipeline paf2tables -i B.sorted.paf.gz --readinfo B.readinfo.tsv.gz
maligno compare-pipeline merge-readinfo -a A.readinfo.tsv.gz -b B.readinfo.tsv.gz -o compare.tsv.gz
```

| Subcommand    | Input                              | Output                                  |
|---------------|------------------------------------|-----------------------------------------|
| **`compare`** | two PAFs (`-a`, `-b`; file paths only, no stdin) | **a results directory** (`--outdir`/`--prefix`): per-set alninfo + readinfo for A and B, plus the comparison TSV — **the primary entry point**. Sorts inputs and verifies read-ID sets match. Emits the single 96-column comparison table as TSV, Parquet or both (`--format`) |
| `compare-pipeline paf2tables`  | PAF (`-i`, `.gz`/`-` ok)           | **alninfo TSV** (`--alninfo`, 35 cols) and/or **readinfo TSV** (`--readinfo`, 33 cols), in one pass |
| `compare-pipeline merge-readinfo` | two readinfo TSVs (`-a`, `-b`) | per-read comparison table (`-o`, 96 cols); the same table `compare` writes. Writes Parquet when `-o` ends in `.parquet`, TSV otherwise |
| `sam2paf`     | SAM file or stdin (`-`)            | PAF written to stdout                   |

`compare-pipeline` groups together the lower-level building blocks used
internally by `compare` (`paf2tables`, `merge-readinfo`, and `summary` — see
below) — most users only need `compare` itself, `find-aln-diff`, and
`compare-pipeline summary`.

All inputs/outputs transparently support gzip (`.gz` suffix) and stdin/stdout (`-`) —
except `compare`'s `-a`/`-b`, which require real file paths (no stdin), since
`compare` always needs two independent inputs.

### Working with PAF files: `compare-pipeline paf2tables` (start here)

`compare-pipeline paf2tables` is the one command that turns a PAF into the
downstream tables. Give it the output path(s) you want — it writes the
**alninfo** table, the **readinfo** table, or **both in a single pass**
(reading the PAF only once):

```bash
# Both tables in one pass (reads the PAF once):
maligno compare-pipeline paf2tables -i in.paf --alninfo alninfo.tsv.gz --readinfo readinfo.tsv.gz

# Just the per-alignment table (no grouping required; works on unsorted PAF):
maligno compare-pipeline paf2tables -i in.paf --alninfo alninfo.tsv.gz

# Just the per-read summary (best alignment per read):
maligno compare-pipeline paf2tables -i in.paf --readinfo readinfo.tsv.gz
```

At least one of `--alninfo` / `--readinfo` must be given. 

**Grouping requirement.** The `--readinfo` output requires the PAF be **grouped**
by `Query_Name` (each read's alignments contiguous); `--alninfo` never cares
about order. By default `paf2tables` groups contiguous runs and prints a one-shot
warning if the input isn't byte-lex sorted. Add **`--strict-grouping`** to turn a
silently-wrong non-grouped input into a hard error — it tracks completed read
names (memory ∝ number of distinct reads) and aborts if a `Query_Name` reappears
non-contiguously:

```bash
maligno compare-pipeline paf2tables -i in.paf --readinfo readinfo.tsv.gz --strict-grouping
```

The standard way to guarantee grouping (and satisfy the downstream
`merge-readinfo` byte-lex requirement at the same time) is a single up-front sort:

```bash
LC_ALL=C sort -t$'\t' -k1,1 in.paf | maligno compare-pipeline paf2tables -i - --readinfo readinfo.tsv.gz
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
   the exact tailored `LC_ALL=C sort … | maligno compare-pipeline paf2tables …` command to run.
7. **`--assume-grouped` fast-path** *(future).* Skip all checks for maximum
   throughput when the caller guarantees grouping.

### Primary: on-rails `compare`

`compare` runs the whole pipeline for you and **owns its preconditions** — you
don't have to pre-sort or worry about ordering:

```bash
maligno compare -a A.paf -b B.paf --label-a A --label-b B --outdir results/ --prefix AvsB
```

What it does, in order:
1. **Sorts** both PAFs by `Query_Name` (byte-lex) with an in-process external sort
   (`ext-sort`: buffers up to `--sort-mem`, default 1G, spilling to temp files under
   `--sort-tmp-dir`, default `--outdir`). This guarantees grouping and a consistent
   matching order — it can't silently mis-compare unsorted input.
2. **Verifies** both PAFs carry the **same `Query_Name` set** (O(1) memory). By
   default it **errors** if they differ, reporting how many IDs are shared / only
   in A / only in B (with examples). Pass `--allow-id-mismatch` to compare the
   shared intersection instead.
3. In a **single in-memory pass**, collapses both sorted PAFs in lock-step and
   feeds the merge-join directly (no readinfo written-then-reread), teeing out the
   per-set `alninfo` + `readinfo` tables and writing the comparison table:
   ```
   {prefix}.{label_a}.alninfo.tsv.gz    {prefix}.{label_b}.alninfo.tsv.gz
   {prefix}.{label_a}.readinfo.tsv.gz   {prefix}.{label_b}.readinfo.tsv.gz
   {prefix}.compare.tsv.gz
   ```
   The sorted PAFs are scratch (removed unless `--keep-sorted-paf`).

Pass **`--no-alninfo`** and/or **`--no-readinfo`** to skip writing those per-set
tables entirely (no file is created — the bytes are never serialized/compressed;
`--no-alninfo` is the biggest time/disk saver since alninfo is the largest output).
The comparison itself is unaffected.

**`--presorted`** skips the internal sort (Step 1) when your inputs are already
prepared. It only requires that both PAFs contain the **same reads in the same
relative order**, grouped by `Query_Name` — *any* consistent ordering works (e.g.
`samtools sort -n` output; byte-lex is **not** required). Instead of the upfront
set-check, the single compare pass verifies the two files line up read-for-read
and **errors on the first divergence** (leaving no partial output). Because nothing
is sorted, `--presorted` cannot be combined with `--allow-id-mismatch` (computing a
shared intersection needs a known sort order) or `--keep-sorted-paf` (no temp files are
created). Use it to avoid the sort cost when you trust your inputs are aligned.

Ideal for comparing two parameter sets / references run on the **same** read or
transcript set. **Precondition:** a `Query_Name` identifies one read/sequence
(maligno sorts by name only). The comparison table is identical to the manual
`compare-pipeline paf2tables` → `compare-pipeline merge-readinfo` path on the
same inputs.

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

### Primary: on-rails `compare` from BAM/PAF

```bash
BIN=./target/release/maligno

# 0. BAM → PAF  (requires samtools; -h preserves @SQ header for contig lengths)
samtools view -h refA.bam | $BIN sam2paf -U - | gzip > refA.paf.gz
samtools view -h refB.bam | $BIN sam2paf -U - | gzip > refB.paf.gz

# 1. Compare — sorts both PAFs, checks read-ID sets match, writes the results dir.
#    No pre-sorting needed; compare does it (--sort-mem caps the in-RAM sort buffer).
$BIN compare -a refA.paf.gz -b refB.paf.gz \
  --label-a RefA --label-b RefB --outdir results/ --prefix RefA_vs_RefB --sort-mem 2G
# → results/RefA_vs_RefB.{RefA,RefB}.alninfo.tsv.gz
#   results/RefA_vs_RefB.{RefA,RefB}.readinfo.tsv.gz
#   results/RefA_vs_RefB.compare.tsv.gz
```

### Explicit intermediates: `compare-pipeline paf2tables` then `compare-pipeline merge-readinfo`

Same result, but materializes the per-file `readinfo` (and optionally `alninfo`)
tables for other analyses.

```bash
BIN=./target/release/maligno
# (steps 0 and 0.5 as above)

# 1. PAF → readinfo  (one row per read; best alignment chosen by ms, then AS, then MQ).
#    Add --alninfo <path> to also emit the 35-col per-alignment table in the same pass.
$BIN compare-pipeline paf2tables -i refA.sorted.paf.gz --readinfo refA.readinfo.tsv.gz
$BIN compare-pipeline paf2tables -i refB.sorted.paf.gz --readinfo refB.readinfo.tsv.gz

# 2. Compare the two readinfo files (streaming, constant memory; strict order by default).
$BIN compare-pipeline merge-readinfo \
  -a refA.readinfo.tsv.gz --label-a RefA \
  -b refB.readinfo.tsv.gz --label-b RefB \
  -o RefA_vs_RefB.compare.tsv.gz
```

---

## How each step works

### `compare-pipeline paf2tables` (alninfo conversion)

Parses each PAF record (12 mandatory fields + `ms:i`, `AS:i`, `cs:Z` tags), walks the
`cs` tag to accumulate match/substitution/insertion/deletion/splice statistics, computes
soft-clip lengths, junction coordinates (strand-aware), and derived scalars
(`seqid`, `Query_Aln_Len`, `Query_Aln_Cov`). This is the `--alninfo` output of
`paf2tables`.

**Pure streaming, constant memory** for the alninfo output. Each PAF line is parsed and
written independently in input order — no internal collect-then-sort. For the explicit
pipeline (`paf2tables` → `merge-readinfo`), pre-sort the PAF by `Query_Name` once upstream:

```bash
LC_ALL=C sort -t$'\t' -k1,1 in.paf > sorted.paf
```

Unix `sort` does external-sort with bounded memory and handles files larger than RAM.
The pre-sort satisfies both the `--readinfo` contiguity requirement and
`merge-readinfo`'s byte-lex sort requirement in one pass.

**Unaligned reads are kept.** A PAF record with `Target_Name == "*"` produces a full row
with zeroed alignment statistics, allowing unaligned reads to flow through the entire
pipeline.

### Readinfo collapse (used by `compare-pipeline paf2tables --readinfo` and `compare`)

The per-read collapse step. Groups alninfo rows by `Query_Name` (contiguous in sorted input) and collapses each group
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

### `compare-pipeline merge-readinfo` (and the comparison core)

A two-pointer **merge-join** over two sorted readinfo files, matching on
**(Read_Name, Read_Len)**. This is the engine behind both `merge-readinfo` (readinfo
TSVs in) and the primary `compare` (which feeds it collapsed rows straight from PAFs).
For each matched read it emits the 31 data columns from each
side (suffixed `_A` / `_B`) plus 30 comparison/object columns
(`AS_Diff`, `ms_Ratio`, `seqid_Diff`, `Junction_Distance`, `N_Matched_Junctions`, `Genomic_N_Matched_Junctions`, `Junctions_OnlyA`, …) — 96 columns in total, including the four leading key/label columns.

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
above. These are **always emitted**; four count columns appear in the comparison block:

| Column                          | Meaning |
|---------------------------------|---------|
| `Genomic_N_Matched_Junctions`   | overlap on `(chrom, start, end)` sets |
| `Genomic_N_Unmatched_Junctions` | set symmetric difference (= `OnlyA + OnlyB`) |
| `Genomic_N_Junctions_OnlyA`     | only in A |
| `Genomic_N_Junctions_OnlyB`     | only in B |

The `genomic_junctions` column (always emitted in the alninfo and readinfo tables) uses
0-based half-open BED coordinates in **`((start, end), ...)`** Python-tuple-of-tuples
form, parseable with `ast.literal_eval`. The chromosome is **not** in each tuple — it's
in the sibling `TargetChr` (alninfo) / `TargetChr_A` & `TargetChr_B` (compare) column. Cross-
chromosome safety in the set comparison is still preserved: the comparison commands
reconstruct full `(chrom, start, end)` keys internally by combining
each row's parsed pairs with its per-side `TargetChr`, so junctions on different contigs
cannot accidentally match.

> **Format change (v0.2.3).** The genomic-junction tuples used to include the chrom as
> the first element (e.g. `(('chr22', 100, 250), …)`). That was redundant with the
> `TargetChr` column, so it was dropped. Pre-v0.2.3 TSVs need to be regenerated from PAF
> to be readable by `compare` / `compare-readinfo`.

> **Breaking schema change (v0.13.0) — fixed `_A`/`_B` side suffixes.** Per-side
> comparison columns used to be suffixed with the *dataset label*
> (`TargetChr_Splice`, `cs_SpliceHQ`, …). They now always use the fixed suffixes
> `_A` and `_B`, and two new columns — **`Label_A`** and **`Label_B`**, at
> positions 3–4 — record which dataset each side is, repeated on every row so any
> row subset stays self-describing.
>
> Why: label-suffixed names were dataset-specific (every downstream script had to
> interpolate the label) and ambiguous whenever a label itself contained an
> underscore, since `Target_Start_my_run` cannot be decomposed reliably. Column
> counts grew by two: **94 → 96** for the full table and **47 → 49** for the
> then-available `--mode junctions` view (removed in v0.14.0, see below).
>
> The `…summary.tsv` categories changed to match: `aligned_only_A` /
> `aligned_only_B` and `present_only_in_A_by_id` / `present_only_in_B_by_id`
> (previously label-interpolated), with the labels emitted as `label_A` /
> `label_B` provenance rows. Summary keys are therefore now stable across
> datasets.
>
> **Pre-v0.13.0 comparison tables are not readable** by `compare-summary` or
> `find-query-diff` in v0.13+; they exit with an error telling you to regenerate.
> Regenerate from PAF with `compare` (or from readinfo with `compare-readinfo`).
> `--label-a` / `--label-b` are now also validated: they must be non-empty,
> distinct, and free of tabs, newlines, and path separators. Underscores are fine.

> **Breaking schema change (v0.14.0) — one comparison table, columns regrouped.**
> `--mode full | junctions` has been **removed** from `compare` and
> `compare-readinfo`. There is now exactly one comparison table (96 columns), and
> the per-side and comparison blocks are **grouped by topic** rather than
> following the readinfo header order. **No column was added or removed, and no
> value changed — only positions moved.**
>
> Why: the 49-column `junctions` view contained no column that was not already in
> the 96-column view, under the same name and computed identically — it was purely
> a column selection, and it bought nothing but a smaller file (and little of
> that: `cs`, `junctions` and `genomic_junctions`, which it kept, are the bulk of
> the table's bytes). Subsetting columns is better done downstream, where it is
> not limited to one hard-coded choice of subset. Removing the flag also deletes a
> duplicated row emitter that every future schema change would have had to update
> twice.
>
> The per-side block (each column suffixed `_A` then `_B`) is now ordered: locus
> and span → alignment selection and score → identity and coverage → junction
> counts → cs-derived event counts → the three long strings (`junctions`,
> `genomic_junctions`, `cs`) last. The comparison block is ordered:
> orientation/identity → score → event diffs → query-space junction metrics →
> genomic-space junction metrics → the four object lists last. Layout:
> keys 1–4, `_A` 5–35, `_B` 36–66, comparison metrics 67–92, object lists 93–96.
>
> This decouples the comparison table's per-side order from `readinfo.rs`'s
> `READINFO_HEADER` (which is **unchanged** — the readinfo and alninfo formats are
> untouched). Columns are read by name, so the divergence is deliberate.
>
> **Impact:** scripts that select columns by *name* need no change; scripts that
> use hard-coded column *numbers* must be updated (see the migration note on
> column positions below, and the `tsvcut` / `tsvwhere` helpers in the test-data
> section). Existing tables remain readable — `compare-summary` and
> `find-query-diff` resolve every column by name.
>
> **Not affected:** `find-query-diff`'s `--compare-by all | junctions` is a
> *different* flag and is unchanged. It selects what counts as a difference (whole
> cs tag vs. query-space junction set only), not which columns are written, and it
> still suffixes its own outputs with `.junctions`. Every junction column stays in
> the table: `N_Matched_Junctions`, `N_Unmatched_Junctions`,
> `N_Junctions_OnlyA/B`, the four `Genomic_N_*`, `Junction_Distance`,
> `Junc_Dist_V2`, and the four object lists. `compare` no longer writes
> `{prefix}.compare.junctions.tsv.gz` / `.compare.junctions.summary.tsv`.

> **Bug fix changing emitted values (v0.15.0) — unmapped reads no longer report
> soft-clipping.** For an unmapped record, `N_SoftClipped_Bases_Start` was reported
> as the full read length and `N_SoftClipped_Events` as `1`. Both are now `0`.
> Affects the `alninfo`, `readinfo` and comparison tables, and in the comparison
> table also `N_SoftClipped_Bases_Start_Diff`. `N_SoftClipped_Bases_End` was already
> `0` and is unchanged.
>
> Why: soft-clip length is computed from the alignment geometry
> (`query_start`/`query_end`/`strand`), not from the cs tag. An unmapped PAF record
> carries a placeholder interval of `(0, 0)` and strand `*`, so the minus/unknown
> branch of the formula returned `query_len - 0` — mechanically concluding that the
> whole read was soft-clipped. There is no alignment, so there are no unaligned
> *ends*; `0` is now reported, consistent with every other alignment-derived field,
> all of which already came out `0`/`NaN` for unmapped rows.
>
> **Classification is unaffected.** The identity classifier reads `TargetChr`,
> `Strand`, `cs`, `Query_Start`/`Query_End` and `Target_Start`/`Target_End` — never
> soft-clip — so `query_identical`, `reference_identical`, every `…summary.tsv`
> count and all `find-query-diff` output are byte-identical across this change.
>
> **Pre-v0.15.0 tables carry the wrong values** for unmapped reads. They remain
> readable and every other column is unaffected, so regenerate only if soft-clip
> statistics on unmapped reads matter to your analysis. Reads that aligned in both
> sets were never affected.

> **New output format (v0.16.0) — Parquet.** `compare --format tsv|parquet|both`
> (default `both`) and `compare-readinfo -o …parquet` write the comparison table as
> Parquet in addition to, or instead of, the gzipped TSV. Same 96 columns, same
> names, same order — `pd.read_parquet` is a drop-in for `pd.read_csv`.
>
> Why: the table is written once and read many times, and column pruning means a
> reader touching a few of the 96 columns skips the rest of the file. Measured on
> the 507,365-row Splice-vs-SpliceHQ comparison (DuckDB 1.5.5, default CSV
> sampling, best of three) — the 16 columns `find-query-diff` needs take 2.01 s
> from `tsv.gz` and 0.55 s from Parquet (3.6×); a two-column aggregate drops from
> 0.85 s to 0.02 s (48×); a full 96-column scan gains least, 3.97 s to 2.31 s
> (1.7×). Writing is *faster* too — 4.5 s vs 10.8 s, measured with maligno itself —
> because zstd beats gzip here and the Parquet path neither formats numbers to text
> nor escapes 66 fields per row. Costs: ~31% more disk (85 vs 65 MB — six
> long-string columns dominate this table) and ~306 MB peak RSS while writing
> versus 35 MB.
>
> Benchmark caveat worth recording: an earlier version of this note quoted 15× and
> 470×, measured with DuckDB's `sample_size=-1`. That forces a full-file
> type-inference scan before any row is read — a cost neither maligno nor a normal
> reader pays — and inflated the TSV side 2–4×. Note also that decompression is not
> the bottleneck it might appear: gunzipping the whole 408 MB table takes 0.22 s.
> The cost being avoided is splitting every row into 96 fields.
>
> **Nulls mean "undefined", nothing else.** A null appears where the TSV carries
> `NaN` in a float column, or an empty `cs`. Values that mean something are kept as
> values: `TargetChr` and `Strand` stay `*` for an unmapped side (they are the
> mapping indicator), an empty junction set stays `"()"`, and a real `0` stays `0`.
> A per-side column that is absent or unparseable becomes a null rather than a
> fabricated `0`.
>
> The two serializations are related by a documented inverse, so a Parquet file can
> be turned back into the exact TSV: null → `NaN` for float columns and → `""`
> otherwise, then `escape_tsv_field` once over the per-side and object-list string
> columns. Note the asymmetry there — per-side values reach the writer already
> escaped once, and the TSV escapes them again, so Parquet holds the *less* escaped
> form.
>
> At the time, `compare-summary` and `find-query-diff` read TSV only, so
> `--format parquet` produced a table they could not consume; `both` was the
> default for that reason. (Before v0.17.0, `--format parquet` was additionally
> rejected unless `--skip-find-query-diff` was passed, because `compare` ran
> `find-query-diff` itself. That coupling is gone.) Since v0.18.1, `find-aln-diff`
> and `compare-pipeline summary` (renamed from `find-query-diff` /
> `compare-summary`) can read Parquet directly via `--input-format`; `both`
> remains the default anyway, for backward compatibility. The exact-pinned
> `arrow-array` / `arrow-schema` / `parquet` dependencies must be bumped together.

> **Breaking change (v0.17.0) — `compare` no longer runs `find-query-diff`.** It
> now produces the per-set alninfo + readinfo tables, the comparison table
> (`--format`), and `{prefix}.compare.summary.tsv`. The four query-diff outputs —
> `{prefix}.query_diff_reads.tsv[.gz]`, `{prefix}.query_diff_regions.{A,B}.bed[.gz]`
> and `{prefix}.query_diff_summary.tsv` — are no longer produced by `compare`.
>
> Why: `compare` wrote the comparison table and then **re-read the whole thing** to
> produce those four files — a second full pass, for outputs the caller may not
> want. Splitting the commands removes that pass, lets `find-query-diff` be re-run
> with different options without redoing the comparison, and makes `--format`
> orthogonal (it no longer has to guarantee a TSV for an internal consumer).
>
> **Migration** — add one command after `compare` (later renamed to
> `find-aln-diff`, and gzip later became the default output, so the flag below
> is no longer needed to match the settings the fused step used):
> ```bash
> maligno find-aln-diff -i results/AvsB.compare.tsv.gz \
>   --outdir results/ --prefix AvsB --compare-by all
> ```
> Those reproduce the four files **byte-for-byte** — verified against the
> v0.16.0 gate baseline for all four comparison scenarios. Any other `--compare-by`
> choice is now equally available (as is `--no-gzip`, for plain-text output).
>
> `--skip-find-query-diff` is removed; passing it is an unknown-argument error. The
> summary table is unaffected — it is accumulated *during* the merge pass, not by
> re-reading the table.

**Strand tracking and renames (v0.2.1+).** Each side now carries a `Strand_A` / `Strand_B` data
column (the best alignment's strand), and the comparison block starts with a `Strand_Match`
(true/false) metric that flags strand-flips between A and B. The legacy column name
`TargetRef_1st` has been renamed to `TargetChr` (suffixed in compare output as
`TargetChr_A` / `TargetChr_B`).

**Non-overlap junction objects (v0.2.2+).** In addition to the *counts* of non-overlapping
junctions (`N_Junctions_OnlyA/B`, `Genomic_N_Junctions_OnlyA/B`), the comparison outputs
now append the actual junction **objects** that failed to overlap at the very end of each
row: `Junctions_OnlyA`, `Junctions_OnlyB` (query-coord tuples, always emitted) and
`Genomic_Junctions_OnlyA`, `Genomic_Junctions_OnlyB` (genome-coord tuples, also always
emitted). These
use the same Python tuple format as the per-side `junctions` / `genomic_junctions` data
columns — parse with `ast.literal_eval` in Python.

> **Migration note (column positions).** Schema changes have shifted column positions
> several times in pre-1.0 development: when `genomic_junctions` was added, again when
> `Strand` and `Strand_Match` were added (v0.2.1), again when `Label_A`/`Label_B` were
> inserted at columns 3–4 (v0.13.0), and again when the per-side and comparison blocks
> were regrouped (v0.14.0). Scripts that filter by column *number*
> (`awk '$73 > 0'`) need updating each time; prefer column-*name* lookup using the header,
> which is robust to future schema growth:
> ```bash
> awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $c["N_Unmatched_Junctions"]>0'
> ```
> Also: if you have older readinfo TSVs with `TargetRef_1st`, rerun
> `compare-pipeline paf2tables --readinfo` to get the renamed column (or rename
> in your scripts).

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
The `compare-pipeline paf2tables` alninfo output is pure streaming and order-preserving
(it does not sort internally). The cleanest fix is a one-time pre-sort upstream on the PAF, which carries
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
time $BIN compare-pipeline paf2tables -i test_data/Splice.AlnToHG38.PriAln.paf.gz   --readinfo /tmp/Splice.readinfo.tsv.gz
time $BIN compare-pipeline paf2tables -i test_data/SpliceHQ.AlnToHG38.PriAln.paf.gz --readinfo /tmp/SpliceHQ.readinfo.tsv.gz

# Sort each readinfo on (Read_Name, Read_Len) so the two files share byte-lex order
# (compare is strict by default — it errors on a read-name mismatch unless inputs match).
for S in Splice SpliceHQ; do
  ( zcat < /tmp/$S.readinfo.tsv.gz | { IFS= read -r h; printf '%s\n' "$h"; \
      LC_ALL=C sort -t$'\t' -k1,1 -k2,2n; } ) | gzip > /tmp/$S.readinfo.sorted.tsv.gz
done

time $BIN compare-pipeline merge-readinfo \
  -a /tmp/Splice.readinfo.sorted.tsv.gz   --label-a Splice \
  -b /tmp/SpliceHQ.readinfo.sorted.tsv.gz --label-b SpliceHQ \
  -o /tmp/Splice_vs_SpliceHQ.compare.tsv.gz

# Inspect the compare output header (column number → column name)
zcat < /tmp/Splice_vs_SpliceHQ.compare.tsv.gz | head -1 | tr '\t' '\n' | nl


# Tip: the primary on-rails `compare` does all of the above (sort + per-set tables +
# comparison) in one command, writing everything to a results directory:
#   $BIN compare -a Splice.paf.gz -b SpliceHQ.paf.gz \
#     --label-a Splice --label-b SpliceHQ --outdir results/ --prefix Splice_vs_SpliceHQ

# ── Select and filter columns BY NAME ───────────────────────────────────────
# Hardcoded `cut -f N` breaks whenever the schema changes — it did in v0.13.0
# (which inserted Label_A/Label_B at columns 3-4) and again in v0.14.0 (which
# regrouped the per-side and comparison blocks). These two helpers resolve columns
# from the header instead, so they keep working across versions.

CMP=/tmp/Splice_vs_SpliceHQ.compare.tsv.gz

# tsvcut <file.gz> <Name1,Name2,...>  — print just those columns, in that order.
# (`gzip -dc` rather than `zcat`: macOS zcat rejects a plain `.gz` name.)
tsvcut () {
  gzip -dc < "$1" | awk -F'\t' -v want="$2" '
    NR==1 { n=split(want,w,","); for (i=1;i<=NF;i++) h[$i]=i
            for (j=1;j<=n;j++) {
              if (!(w[j] in h)) { print "no such column: " w[j] > "/dev/stderr"; exit 1 }
              c[j]=h[w[j]]
            } }
    { line=$c[1]; for (j=2;j<=n;j++) line=line "\t" $c[j]; print line }'
}

# tsvwhere <file.gz> <Name> <awk-test>  — keep the header + rows passing the test.
tsvwhere () {
  gzip -dc < "$1" | awk -F'\t' -v col="$2" -v test="$3" '
    NR==1 { for (i=1;i<=NF;i++) if ($i==col) k=i
            if (!k) { print "no such column: " col > "/dev/stderr"; exit 1 }
            print; next }
    { v=$k+0 } test=="pos" ? v>0 : v==0'
}

# Distribution of unmatched junctions, query space then genomic space.
tsvcut "$CMP" N_Unmatched_Junctions         | tail -n +2 | sort | uniq -c
tsvcut "$CMP" Genomic_N_Unmatched_Junctions | tail -n +2 | sort | uniq -c

# Reads whose QUERY-space junction sets disagree: count, then inspect.
tsvwhere "$CMP" N_Unmatched_Junctions pos | wc -l
tsvwhere "$CMP" N_Unmatched_Junctions pos \
  | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "Read_Name\tJuncCount_A\tJuncCount_B\tN_Matched_Junctions\tN_Junctions_OnlyA\tN_Junctions_OnlyB"; next}
                {print $h["Read_Name"]"\t"$h["JuncCount_A"]"\t"$h["JuncCount_B"]"\t"$h["N_Matched_Junctions"]"\t"$h["N_Junctions_OnlyA"]"\t"$h["N_Junctions_OnlyB"]}' \
  | column -t -s $'\t' | less -S

# Same, but the actual non-overlapping junction tuples rather than counts.
tsvwhere "$CMP" N_Unmatched_Junctions pos \
  | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "Read_Name\tJunctions_OnlyA\tJunctions_OnlyB"; next}
                {print $h["Read_Name"]"\t"$h["Junctions_OnlyA"]"\t"$h["Junctions_OnlyB"]}' \
  | column -t -s $'\t' | less -S

# Reads whose GENOMIC-space junction sets disagree.
tsvwhere "$CMP" Genomic_N_Unmatched_Junctions pos | wc -l
tsvwhere "$CMP" Genomic_N_Unmatched_Junctions pos \
  | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; print "Read_Name\tTargetChr_A\tTargetChr_B\tGenomic_N_Matched_Junctions\tGenomic_N_Junctions_OnlyA\tGenomic_N_Junctions_OnlyB"; next}
                {print $h["Read_Name"]"\t"$h["TargetChr_A"]"\t"$h["TargetChr_B"]"\t"$h["Genomic_N_Matched_Junctions"]"\t"$h["Genomic_N_Junctions_OnlyA"]"\t"$h["Genomic_N_Junctions_OnlyB"]}' \
  | column -t -s $'\t' | less -S

# Verify column counts (expect 35, 33, 96)
gzip -dc < /tmp/Splice.alninfo.tsv.gz             | awk -F'\t' '{print NF}' | sort | uniq -c
gzip -dc < /tmp/Splice.readinfo.tsv.gz            | awk -F'\t' '{print NF}' | sort | uniq -c
gzip -dc < /tmp/Splice_vs_SpliceHQ.compare.tsv.gz | awk -F'\t' '{print NF}' | sort | uniq -c
```

---

## Summary statistics (`compare-pipeline summary`)

`compare` tallies aggregate summary statistics over the per-read comparison as the
rows stream out (constant memory) and writes them to
`{prefix}.compare[.junctions].summary.tsv` (2 columns: `Category`, `Count`) plus a
stderr block. The **same** statistics can be computed from an existing comparison
table with the standalone command:

```bash
maligno compare-pipeline summary -i AvsB.compare.tsv.gz -o AvsB.compare.summary.tsv
# -i: a compare / compare-pipeline merge-readinfo table (.gz, .parquet, or - ok); -o: optional (else stderr only)
```

`compare-pipeline summary` requires the fixed `TargetChr_A` / `TargetChr_B` columns, and
reads the set names from `Label_A` / `Label_B`. Every column is resolved **by name**, so
column order does not matter. It sees only matched rows, so the `present_only_in_*_by_id`
counts are always 0 there; the built-in `compare` tally fills those from the read-ID merge.

### Classification (per matched read)

Computed from each side's **representative (best) alignment** — the row already
selected for the readinfo/compare output — using these columns: `TargetChr`,
`Strand`, `cs`, `Query_Start`, `Query_End`, `Target_Start`, `Target_End`.

- **Mapping status** (`TargetChr == "*"` or empty ⇒ unmapped): `aligned_both`,
  `aligned_only_A`, `aligned_only_B`, `aligned_neither`.
- **Query-coordinate identical** — both sides mapped, same query span
  (`Query_Start`/`Query_End`, which PAF reports in forward-read coordinates), and
  the same alignment relative to the read via **either**:
  - *same strand*: `Strand_A == Strand_B` and `cs_A == cs_B`, or
  - *reverse-complement*: `Strand_A != Strand_B` and `cs_A == cs_revcomp(cs_B)` —
    an inverted/opposite-strand alignment of the same query-to-reference
    correspondence (e.g. a locus inverted between two assemblies). `cs_revcomp`
    reverses the cs op order and complements each op length-preservingly
    (`:N`→`:N`; `=ACGT`→`=`+revcomp; `*xy`→complemented bases; `+`/`-`→revcomp of
    the sequence; intron `~gt…ag`→`~ct…ac`).

  Both `cs_A == cs_B` and `cs_A == cs_revcomp(cs_B)` are evaluated **motif-blind**:
  each cs string is first passed through `cs_strip_splice_motifs`, which replaces
  every intron's 2-letter donor/acceptor motif with a fixed `nn`/`nn` placeholder
  while keeping the intron length and every other op (matches, substitutions,
  indels) unchanged. This absorbs a real-world aligner-output limitation — e.g.
  STAR reports `~nn<len>nn` where minimap2 reports the true motif (`~ct<len>ac`)
  for the identical intron — so it no longer, by itself, makes two otherwise
  identical alignments count as different.

  Reported as `query_identical`, split into `query_identical_same_strand` and
  `query_identical_revcomp`; `query_not_identical = aligned_both - query_identical`.
- **Reference-space classification** — independent of query span, and a strict,
  literal comparison: no reverse-complement accommodation, so a real strand
  difference always counts as a different position. Two axes, both motif-blind
  on `cs`:
  - **position**: identical `TargetChr` + `Strand` + `Target_Start` (a matching
    start plus matching `cs` implies a matching `Target_End` too, since the `cs`
    tag's own operations determine the alignment's length — no need to check it
    separately).
  - **alignment**: identical `cs`.

  The four combinations are reported as `ref_same_position_same_aln` (reference-
  identical), `ref_same_position_diff_aln`, `ref_diff_position_same_aln`
  (relocated), and `ref_diff_position_diff_aln`.

### Summary TSV schema

| Category | Meaning |
|----------|---------|
| `label_A` / `label_B` | which dataset each side is (the `--label-a` / `--label-b` values). String values, not counts — provenance rows so the summary is self-describing |
| `reads_compared` | matched reads written to the comparison table |
| `aligned_both` | representative alignment mapped in both sets |
| `aligned_only_A` / `aligned_only_B` | mapped in one set, `"*"` in the other |
| `aligned_neither` | unmapped (`"*"`) in both |
| `query_identical` | query-coordinate identical (see above) |
| `query_identical_same_strand` | …via the same-strand branch |
| `query_identical_revcomp` | …via the reverse-complement branch |
| `query_not_identical` | both mapped but not query-identical |
| `ref_same_position_same_aln` | same `TargetChr`/`Strand`/`Target_Start` **and** same `cs` (motif-blind) — reference-identical |
| `ref_same_position_diff_aln` | same position, different `cs` |
| `ref_diff_position_same_aln` | same `cs`, different position — relocated |
| `ref_diff_position_diff_aln` | both position and `cs` differ |
| `present_only_in_A_by_id` / `present_only_in_B_by_id` | read present in only one set's PAF (built-in `compare` only; 0 unless `--allow-id-mismatch`) |

---

## Differing reads & regions (`find-aln-diff`)

`find-aln-diff` (renamed from `find-query-diff`) reads a `compare` /
`compare-pipeline merge-readinfo` comparison table and reports every read whose
alignment differs between A and B, in query space or reference space
(`--space`), plus merged genomic regions showing where those differing reads
cluster.

**A separate command since v0.17.0.** `compare` used to run this itself as a final
step, which meant re-reading the comparison table it had just written — a second
full pass for outputs the caller may not want. `compare` now prints the exact
command to run instead. Running it by hand with `--space query` (the default)
reproduces the previous fused-step outputs byte-for-byte, given the defaults the
fused step used (`--compare-by all`; gzip is on by default).

**Usage:**

```bash
maligno find-aln-diff -i AvsB.compare.tsv.gz --outdir results/ --prefix AvsB \
  [--space query|reference] [--no-gzip] [--compare-by all|junctions] [--emit-identical-reads]
# -i: a compare / compare-pipeline merge-readinfo table (.gz, .parquet, or - ok)
```

The fixed `TargetChr_A` / `TargetChr_B` columns are required (same as
`compare-pipeline summary`), and set names come from the `Label_A` / `Label_B`
columns.

### `--space` — which coordinate space defines a difference

| Value | "Identical" means |
|-------|-------------------|
| `query` (default) | `query_identical` (the exact same `classify()` used by `compare-pipeline summary` — a reverse-complement match still counts as identical) |
| `reference` | `ref_same_position_same_aln` (the reference-space classification above): same `TargetChr`/`Strand`/`Target_Start`, and same alignment content per `--compare-by` below. Strict and literal — no reverse-complement accommodation. |

`--space` picks the output filename stem (`query_diff`/`query_identical` vs.
`reference_diff`/`reference_identical`), so a `query`-space and
`reference`-space run at the same `--outdir`/`--prefix` never clobber each
other.

### `--compare-by` — what defines a difference, within the active `--space`

For reads mapped in **both** sets, `--compare-by` chooses the identity test; the
map-status handling (only-A / only-B / neither) is unchanged either way, and its
category names (`diff_aln_only_A` / `diff_aln_only_B`) are the same in both spaces.

| Value | Query space (`--space query`) | Reference space (`--space reference`) |
|-------|-------------------------------|----------------------------------------|
| `all` (default) | the full `cs` tag matches, motif-blind; reverse-complement counts as identical | same `TargetChr`/`Strand`/`Target_Start` **and** same `cs`, motif-blind (`ref_class.same_position() && ref_class.same_aln()`) |
| `junctions` | the **query-space** splice-junction set matches (`junction_set_stats` on `junctions_A`/`junctions_B`) | same `TargetChr`/`Strand`/`Target_Start` **and** same **genomic-coordinate** junction set (`genomic_junction_set_stats` on `genomic_junctions_A`/`genomic_junctions_B`, chrom reattached) |

Under `--compare-by junctions` (either space), any mismatch/indel/soft-clip
difference that doesn't move a splice junction no longer counts — a
differently-reported intron motif at the same position/length (e.g. STAR's `nn`
placeholder vs. minimap2's true motif) doesn't count as a difference under
`all` either, motif-blind in both spaces.

In `junctions` mode every output filename gains a `.junctions` segment (e.g.
`{prefix}.query_diff_reads.junctions.tsv.gz`) so it never clobbers an `all` run at the
same `--outdir`/`--prefix`. Both modes write leading `space<TAB><query|reference>`
and `compare_by<TAB><all|junctions>` rows in the summary TSV so the file is
self-describing.

### Categories

Derived from the shared `CompareSummary` counters — no re-derivation of the
classification logic. Row/category names follow `--space`:

| Category (`--space query`) | Category (`--space reference`) | Definition |
|----------|----------|------------|
| `diff_aln_to_both` | `reference_diff` | mapped in both sets, not identical under the active space/compare-by (`= aligned_both - <space>_identical`) |
| `diff_aln_only_A` / `diff_aln_only_B` | *(same names)* | mapped in one set only (`= aligned_only_a` / `aligned_only_b`) — identical in both spaces |
| `query_different_total` | `reference_diff_total` | sum of the three categories above |
| `query_identical_total` | `reference_identical_total` | excluded from all outputs |
| `aligned_neither` | *(same name)* | excluded (unmapped in both — no difference to report) |

Reconciliation: `reads_compared == <total> + <identical_total> + aligned_neither`.

### Outputs

Filenames use the `query_diff`/`query_identical` stem under `--space query`
(default) or `reference_diff`/`reference_identical` under `--space reference`;
the table below uses the query-space names. `[.gz]` is present by default —
`--no-gzip` omits it (covers both per-read tables and both region tables).

| File | Contents |
|------|----------|
| `{prefix}.query_diff_reads.tsv[.gz]` | one row per differing read: `Read_Name`, `outcome` (the canonical category name above), plus 8 classification booleans (see below) |
| `{prefix}.query_diff_regions.A.bed[.gz]` | merged loci over reads with an A placement (both-mapped-and-differing + `diff_aln_only_A`) |
| `{prefix}.query_diff_regions.B.bed[.gz]` | merged loci over reads with a B placement (both-mapped-and-differing + `diff_aln_only_B`) |
| `{prefix}.query_diff_summary.tsv` | the category tally above (+ stderr); never gzipped |
| `{prefix}.query_identical_reads.tsv[.gz]` | *(opt-in, `--emit-identical-reads`)* one row per identical read: `Read_Name`, `category` (`query_identical_same_strand` / `query_identical_revcomp` / `query_identical_junctions` under `--space query`; `reference_identical` under `--space reference`), plus the same 8 classification booleans — the complement of the diff-reads file. Off by default; **not** produced by `compare`'s built-in invocation. |

Under `--compare-by junctions` every filename above gains a `.junctions` segment
(e.g. `{prefix}.query_diff_reads.junctions.tsv.gz`).

#### Classification booleans

Both per-read tables above carry the same 8 boolean columns (`1`/`0`) after
`outcome`/`category`, computed the same way **regardless of the active
`--space`/`--compare-by`** — so one run shows, e.g., a read that's
query-different but reference-identical, without a second run in the other
`--space`. `false`/`0` for a read where the axis doesn't apply (e.g. all 8 are
`0` for a read mapped on only one side).

| Column | `1` when |
|---|---|
| `query_identical_same_strand` | `query_identical` and reached via the same-strand branch |
| `query_identical_revcomp` | `query_identical` and reached via the reverse-complement branch |
| `query_junctions_identical` | both mapped and the **query-space** junction sets match (`junction_set_stats`) — computed unconditionally, not just under `--compare-by junctions` |
| `ref_same_position_same_aln` | `RefClass::SamePositionSameAln` (reference-identical) |
| `ref_same_position_diff_aln` | `RefClass::SamePositionDiffAln` |
| `ref_diff_position_same_aln` | `RefClass::DiffPositionSameAln` (relocated) |
| `ref_diff_position_diff_aln` | `RefClass::DiffPositionDiffAln` |
| `ref_same_position_same_junctions` | both mapped, same position, and the **genomic-coordinate** junction sets match (`genomic_junction_set_stats`) — computed unconditionally |

For a both-mapped read, exactly one of the four `ref_*` columns is `1` (they
partition `RefClass`); all four are `0` for a read mapped on only one side.

Region-table columns: `#chrom  start  end  n_reads  n_both  n_only_<A\|B>  n_plus  n_minus`
— `n_reads = n_both + n_only_*`; `n_plus + n_minus <= n_reads`. Loci are formed by
a generic sort + single-sweep merge (`src/interval_merge.rs`), equivalent to
`bedtools merge -c -o count`, verified against a real `bedtools` oracle at both
small (~11.6K reads) and genome scale (~986K reads, 31.5K differing) — exact match
on `(chrom, start, end, n_reads)` in both cases.

**A read may appear on only one side.** A `diff_aln_only_B` read has no A
coordinate and is absent from the A region table (but still counted and listed in
the read TSV); symmetric for `diff_aln_only_A` and the B table. Bad
intervals (unparseable or `end <= start`) are skipped and counted internally
rather than aborting the run.

---

## Troubleshooting

### "`merge-readinfo` matched far fewer reads than I expected" — sort-order diagnostic

`compare-pipeline merge-readinfo` uses a streaming merge-join keyed on
`(Read_Name, Read_Len)`. The algorithm runs in O(1) memory and O(N+M) time,
but it **assumes both readinfo files are sorted in the same byte-lexicographic
order**. (The primary `compare` reads PAFs directly and instead requires both
PAFs to list reads in the *same order* — see its `--ignore-row-mismatch`.) By default a read-name mismatch is a hard error (non-zero exit) so you
notice immediately; passing `--ignore-row-mismatch` reverts to skip-and-count,
where unmatched reads are dropped and tallied in the end-of-run stderr summary.

The `compare-pipeline paf2tables` alninfo output is pure streaming and **preserves input order** —
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
[`scripts/check-readinfo-overlap.sh`](../scripts/check-readinfo-overlap.sh) — same
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

  ⇒ `merge-readinfo`'s 'matched' count should equal 1023. If maligno's reported
    matched count is lower than that, the two readinfo files are not sorted
    in the same byte-lexicographic order. Re-sort each with:

        LC_ALL=C sort -t$'\t' -k1,1 -k2,2n input.readinfo.tsv > sorted.tsv
        # (keep the header separately)

    or pre-sort the PAF once (`LC_ALL=C sort -t$'\t' -k1,1`) and rerun
    `compare-pipeline paf2tables --readinfo` — both readinfo files then share byte-lex order.
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
├── compare.rs              — PRIMARY `compare` command (on-rails: sort → verify read-IDs → tables → compare)
├── external_sort.rs        — in-process PAF external sort (ext-sort) + O(1) read-ID set check
├── paf2tables.rs           — PAF → alninfo and/or readinfo (one pass)
├── comparison_row.rs       — comparison-table schema: column lists, ComparisonRow/AlignmentRow/AlignmentDiff, TSV writers
├── parquet_out.rs          — Parquet writer + OutputFormat; Arrow schema derived from comparison_row's column lists
├── compare_streaming.rs    — `compare-pipeline merge-readinfo` command + the merge-join machinery (ReadKey/ReadInfoReader)
├── readinfo.rs             — collapse library (collapse_group/ReadInfoRow/AlnRow); utils-readinfo CLI unregistered but code kept
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
