# maligno — Reference

Complete reference for **maligno**: every command, what it reads and writes,
and how it decides what to report. For installation and a quick start, see the
[README](../README.md). For column-by-column definitions of every output file,
see [`docs/output-tables/`](output-tables/README.md).

**Contents**

- [Overview](#overview)
- Commands
  - [`compare`](#compare): end-to-end comparison of two alignment files (primary entry point)
  - [`sam2paf`](#sam2paf): SAM/BAM → PAF converter
  - [`toolkit paf2tables`](#toolkit-paf2tables): PAF → alninfo and/or readinfo tables
  - [`toolkit merge-readinfo`](#toolkit-merge-readinfo): two readinfo tables → comparison table
  - [`toolkit summary`](#toolkit-summary): comparison table → summary statistics
  - [`toolkit find-aln-diff`](#toolkit-find-aln-diff): comparison table → differing reads and regions
  - [`toolkit query-junction-diff`](#toolkit-query-junction-diff): comparison table → per-read splice-junction diff
- Shared definitions
  - [Representative alignment](#representative-alignment)
  - [Classification](#classification)
  - [Junction comparison](#junction-comparison)
  - [TSV and Parquet](#tsv-and-parquet)
- [Working with the output tables](#working-with-the-output-tables)
- [Troubleshooting](#troubleshooting)
- [Build](#build)

---

## Overview

maligno compares two sets of alignments of the **same reads** (for example,
one read set aligned with two parameter sets, two aligners, or to two
references), read by read.

There are two ways to produce the comparison table, and both give the
**identical** table:

**1. `compare` (recommended):** one command that sorts both inputs, checks
they carry the same reads, and writes a results directory.

```
  A.paf|sam|bam ─┐                                  {prefix}.compare.tsv.gz / .parquet
                 ├─ compare (sort → verify → compare) ─▶  {prefix}.compare.summary.tsv
  B.paf|sam|bam ─┘                                  {prefix}.query_diff_reads.tsv.gz
                                                    {prefix}.query_diff_regions.{A,B}.bed.gz
```

```bash
maligno compare -a A.bam -b B.bam --label-a A --label-b B --outdir results/ --prefix AvsB
```

**2. Building blocks:** `toolkit paf2tables`, then `toolkit merge-readinfo`.
You sort the inputs yourself and keep the intermediate per-read tables.

```
  A.paf ──paf2tables──▶ A.readinfo.tsv ─┐
                                         ├─ merge-readinfo ─▶ compare.tsv
  B.paf ──paf2tables──▶ B.readinfo.tsv ─┘
```

Either way, the comparison table can then be analysed further with
`toolkit summary`, `toolkit find-aln-diff` and `toolkit query-junction-diff`.

**Input and output conventions.** Inputs and outputs ending in `.gz` are
read and written as gzip. Most commands accept `-` for stdin/stdout. The
exceptions are `compare`'s `-a`/`-b` (real file paths only), Parquet inputs
(real file paths only), and `toolkit query-junction-diff` (Parquet only).

**Read identity.** maligno matches alignments across files by read name
(`Query_Name` in PAF, `QNAME` in SAM/BAM), so each name must identify one
read. For paired-end SAM/BAM, `/1` or `/2` is appended to each mate's name, so
mates are compared separately.

---

## `compare`

```bash
maligno compare -a <a.paf|a.sam|a.bam> -b <b.paf|b.sam|b.bam> \
  --label-a <A> --label-b <B> --outdir <DIR> --prefix <NAME> [options]
```

Runs the full comparison and handles its own preconditions: inputs don't
need to be pre-sorted.

### Inputs

`-a`/`--aln-a` and `-b`/`--aln-b` are each **PAF** (`.paf[.gz]`), **SAM**
(plain or gzip/BGZF) or **BAM**. The format is detected per file, so the two
can differ (e.g. BAM vs PAF). CRAM is not supported. Any BAM sort order works.

SAM/BAM inputs are converted to PAF on the fly with the same code as
`maligno sam2paf -U` (no intermediate file is written), so results are identical
to converting first. Unmapped reads are always kept. `--sam-records` picks
which records are converted:

| `--sam-records` | Records kept | Equivalent |
|---|---|---|
| `primary-supp` (default) | primary + supplementary; skip secondary | `sam2paf -p -U` |
| `primary` | primary only | `sam2paf -P -U` |
| `all` | every record | `sam2paf -U` |

This also decides which alignments the per-read `Num_Aln_tpP` / `Num_Aln_tpS`
counts can see. Under the default, secondaries are dropped, so `Num_Aln_tpS`
is `0`.

Every kept, mapped SAM/BAM record needs a `cs` tag, or an `MD` tag plus `SEQ`
to build one from. A record with neither stops the run (see
[`sam2paf`](#sam2paf), and the README's section on preprocessing STAR BAMs).

### What it does

1. **Sort.** Both inputs are sorted by read name (byte order) with an
   in-process external sort. It holds up to `--sort-mem` (default `1G`) per
   file in memory, uses `--sort-threads` threads (default `1`), and spills to
   temporary files under `--sort-tmp-dir` (default: `--outdir`).
2. **Verify.** Checks that both inputs carry the **same set of read names**.
   If they differ, `compare` stops with an error that reports how many names
   are shared, only in A, and only in B (with examples). With
   `--allow-id-mismatch` it instead compares the shared reads only; reads
   found in only one input get no row and are counted in the summary as
   `present_only_in_A_by_id` / `present_only_in_B_by_id`.
3. **Compare.** In one pass over both sorted inputs, each read's alignments
   are collapsed to a [representative alignment](#representative-alignment)
   on each side, the two are [classified](#classification), and one row per
   read is written to the comparison table. The summary statistics and the
   differing-reads and region tables are produced in the same pass.

### Outputs

Written to `--outdir`, named by `--prefix`:

| File | Written | Spec |
|---|---|---|
| `{prefix}.compare.tsv.gz` | by default (`--format tsv` or `both`) | [compare.md](output-tables/compare.md) |
| `{prefix}.compare.parquet` | by default (`--format parquet` or `both`) | [compare.md](output-tables/compare.md) |
| `{prefix}.compare.summary.tsv` | always | [compare-summary.md](output-tables/compare-summary.md) |
| `{prefix}.query_diff_reads.tsv.gz` | by default | [query-diff-reads.md](output-tables/query-diff-reads.md) |
| `{prefix}.query_diff_regions.A.bed.gz`, `.B.bed.gz` | by default | [query-diff-regions.md](output-tables/query-diff-regions.md) |
| `{prefix}.{label_a}.alninfo.tsv.gz`, `{prefix}.{label_b}.alninfo.tsv.gz` | with `--emit-alninfo` | [alninfo.md](output-tables/alninfo.md) |
| `{prefix}.{label_a}.readinfo.tsv.gz`, `{prefix}.{label_b}.readinfo.tsv.gz` | with `--emit-readinfo` | [readinfo.md](output-tables/readinfo.md) |

- `--format tsv|parquet|both` (default `both`) chooses the comparison table's
  serialization.
- The differing-reads and region tables are `toolkit find-aln-diff`'s output
  at its default settings (`--space query --compare-by all`). They are
  byte-identical to running that command on the comparison table.
  `--skip-find-aln-diff` turns them off. For the other `--space` /
  `--compare-by` modes, run [`toolkit find-aln-diff`](#toolkit-find-aln-diff)
  on the comparison table.
- The per-set alninfo and readinfo tables are off by default. They are the
  largest outputs, and no other maligno command reads them.
- An abbreviated summary is also printed to stderr.
- The sorted intermediate PAFs are deleted at the end unless you pass
  `--keep-sorted-paf`. For SAM/BAM inputs, these are the converted PAFs.

### Labels

`--label-a` / `--label-b` (defaults `SetA` / `SetB`) name the two sets. They
are written to the comparison table's `Label_A` / `Label_B` columns and to the
summary, and are used in the per-set output filenames. They must be non-empty,
must differ from each other, and must not contain tabs, newlines, `/` or `\`.
Underscores are fine. Columns are always suffixed `_A` / `_B`, whatever the
labels are.

### `--presorted`

Skips the sort (step 1) when the inputs are already prepared. Both inputs must
contain the **same reads in the same relative order**, with each read's
records contiguous. Any consistent order works; byte order is not required.
Instead of the upfront check (step 2), the comparison pass checks that the
files line up read for read and stops at the first divergence, leaving no
partial output. `--presorted` cannot be combined with `--allow-id-mismatch` or
`--keep-sorted-paf`.

For SAM/BAM inputs, name-sort both with `samtools sort -n`. Its tie-break
keeps each mate's records contiguous, which raw paired-end aligner output
generally does not. A SAM/BAM whose header says `@HD SO:coordinate` is
rejected under `--presorted`.

Output rows follow the input order. If your sort differs from maligno's,
ties for the [representative alignment](#representative-alignment) can
occasionally resolve differently (a few reads per million in testing).

### Options

| Flag | Purpose |
|---|---|
| `-a`/`--aln-a`, `-b`/`--aln-b` | input alignments for set A / B (PAF, SAM or BAM) |
| `--sam-records` | SAM/BAM only: `primary-supp` (default), `primary`, or `all` |
| `--label-a`, `--label-b` | set names (default `SetA` / `SetB`) |
| `-o`/`--outdir`, `-p`/`--prefix` | output directory and filename prefix |
| `--format` | comparison table as `tsv`, `parquet`, or `both` (default) |
| `--allow-id-mismatch` | compare the shared reads instead of erroring when the read-name sets differ |
| `--presorted` | skip the sort; inputs must already be in the same read order |
| `--emit-alninfo`, `--emit-readinfo` | also write the per-set alninfo / readinfo tables |
| `--sort-mem` | in-memory sort buffer per file (default `1G`; `K`/`M`/`G` suffix) |
| `--sort-tmp-dir` | directory for sort spill files (default: `--outdir`) |
| `--sort-threads` | sort threads (default `1`) |
| `--keep-sorted-paf` | keep the sorted intermediate PAFs |
| `--skip-find-aln-diff` | don't write the differing-reads and region tables |

---

## `sam2paf`

```bash
maligno sam2paf [options] <in.sam|in.bam|->  > out.paf
```

Converts SAM or BAM to PAF on stdout. Its output is byte-for-byte compatible
with the `sam2paf` command of minimap2's `paftools.js`. The input format is
detected automatically; `-` reads SAM text from stdin. `maligno sam2paf in.bam`
gives the same bytes as `samtools view -h in.bam | maligno sam2paf -`.

| Flag | Meaning |
|---|---|
| `-p` | primary + supplementary alignments only (skip secondary, FLAG 0x100) |
| `-P` | primary alignments only (skip secondary and supplementary, FLAG 0x800); implies `-p` |
| `-L` | write the `cs` tag in long form (`=ACGT`) instead of the default short form (`:N`) |
| `-U` | write a placeholder PAF record for each unmapped read (otherwise unmapped reads are dropped) |

Pass `-U` when the PAF will be compared: unmapped reads then appear as
`Num_Aln = 0` rows instead of disappearing.

**The `cs` requirement.** Each kept, mapped record's `cs` comes from its
`cs:Z:` tag when present, or is built from `CIGAR` + `MD` + `SEQ`. If neither
is possible (no `cs` and no `MD`, or `MD` with `SEQ` = `*`, as on many
secondary records), the run **stops with an error** naming the read. Fixes:

- align with `cs` output enabled,
- add `MD` with `samtools calmd`,
- or skip secondaries with `-p` / `-P`.

Records dropped by `-p` / `-P`, and unmapped records, are not checked. Other
per-record problems (e.g. a contig missing from the `@SQ` header lines) print
a warning and skip the record.

The `cs` strings maligno compares are text, so both sides of a comparison
should use the same form. SAM/BAM conversion produces the short form, which is
minimap2's default. Don't compare it against a PAF made with
`minimap2 --cs=long`.

---

## `toolkit paf2tables`

```bash
maligno toolkit paf2tables -i <in.paf[.gz]|-> [--alninfo <alninfo.tsv[.gz]>] [--readinfo <readinfo.tsv[.gz]>] [--strict-grouping]
```

Turns a PAF into maligno's tables, writing either or both in a single pass
over the PAF. At least one output must be given.

- **`--alninfo`** writes one row per PAF record (36 columns; see
  [alninfo.md](output-tables/alninfo.md)). Each record's `cs` tag is walked to
  count matches, substitutions, insertions, deletions and splice junctions;
  soft-clip lengths, junction coordinates and identity/coverage are derived
  from it. Rows are written in input order with constant memory, and input
  order doesn't matter. Unmapped records (`Target_Name == "*"`) get a row with
  zeroed statistics.
- **`--readinfo`** writes one row per read (36 columns; see
  [readinfo.md](output-tables/readinfo.md)): the read's
  [representative alignment](#representative-alignment) plus per-read
  aggregates.

**Grouping.** `--readinfo` requires each read's records to be **contiguous**
in the input. By default contiguous runs are grouped, and a one-time warning
is printed if the read names are not in byte order. That catches the common
cases (e.g. a name-sorted but not byte-sorted aligner output, or shuffled
multi-threaded output), but a scattered read can still slip through.
**`--strict-grouping`** turns any non-contiguous read into an error. It keeps
the set of completed read names in memory, so memory grows with the number of
reads.

The simplest way to guarantee grouping, and the byte order
`toolkit merge-readinfo` needs, is to sort the PAF once up front:

```bash
LC_ALL=C sort -t$'\t' -k1,1 in.paf | maligno toolkit paf2tables -i - --readinfo readinfo.tsv.gz
gzip -dc in.paf.gz | LC_ALL=C sort -t$'\t' -k1,1 | maligno toolkit paf2tables -i - --readinfo readinfo.tsv.gz
```

---

## `toolkit merge-readinfo`

```bash
maligno toolkit merge-readinfo -a <readinfo_a.tsv[.gz]> -b <readinfo_b.tsv[.gz]> \
  [--label-a A] [--label-b B] -o <compare.parquet|compare.tsv[.gz]> [--allow-id-mismatch]
```

Joins two readinfo tables into the alignment comparison table. It runs the same
pairing, [classification](#classification) and row writing as `compare`, so the
table is identical to the one `compare` writes from the same reads (see
[compare.md](output-tables/compare.md)). Output is Parquet when `-o` ends in
`.parquet`, otherwise TSV (`.tsv[.gz]`).

**Matching.** Reads are matched on **`Read_Name`** with a streaming merge:
constant memory, one pass over each file. If a read's `Read_Len` differs between
the two files, the output carries side A's. Both inputs must be **sorted by
`Read_Name` in byte order**. Sort the PAF before `paf2tables` (see above), or
sort existing readinfo files with the header kept on top:

```bash
( gzip -dc in.readinfo.tsv.gz | head -1
  gzip -dc in.readinfo.tsv.gz | tail -n +2 | LC_ALL=C sort -t$'\t' -k1,1
) | gzip > sorted.readinfo.tsv.gz
```

- `LC_ALL=C` forces byte order, regardless of locale.
- `-k1,1` sorts on `Read_Name`.
- `sort` works on files larger than memory. `-T <dir>` and `-S <size>` set its
  scratch directory and memory cap.

**Unmatched reads.** By default, a read present in only one file is an error.
With `--allow-id-mismatch`, only the reads present in both files are compared,
and the others are counted in the summary as `present_only_in_A_by_id` /
`present_only_in_B_by_id`. At the end, the same summary `compare` prints is
written to stderr:

```
  shared: 14   only in Splice: 1   only in SpliceHQ: 1
Comparison summary:
  label_A                            Splice
  label_B                            SpliceHQ
  reads_compared                     14
  aligned_both                       11
  aligned_only_A                     0
  aligned_only_B                     1
  aligned_neither                    2
  query_identical                    7
  query_not_identical                6
```

If `reads_compared` is lower than you expect, see
[Troubleshooting](#troubleshooting).

---

## `toolkit summary`

```bash
maligno toolkit summary -i <compare.parquet|compare.tsv[.gz]|-> [-o summary.tsv[.gz]] [--input-format auto|tsv|parquet]
```

Computes the summary statistics from an existing comparison table. The schema
is the same as `{prefix}.compare.summary.tsv`; see
[compare-summary.md](output-tables/compare-summary.md). The result is printed
to stderr, and also written to `-o` when given (`.tsv[.gz]`, or `-` for stdout).

- **Input format:** `--input-format auto` (the default) reads Parquet for a
  `.parquet` path and TSV otherwise. Parquet must be a real file, not stdin.
- **Columns:** every column is looked up by name, so column order doesn't
  matter. The set names come from `Label_A` / `Label_B`.
- **Reads in only one input:** the table only contains reads present in both
  inputs, so `present_only_in_A_by_id` / `present_only_in_B_by_id` are always
  `0` here. Only `compare` and `toolkit merge-readinfo` can fill them in.

The categories are defined in [Classification](#classification).

---

## `toolkit find-aln-diff`

```bash
maligno toolkit find-aln-diff -i <compare.parquet|compare.tsv[.gz]|-> --outdir <DIR> --prefix <STR> \
  [--space query|reference] [--compare-by all|junctions] [--emit-identical-reads] [--no-gzip] \
  [--input-format auto|tsv|parquet]
```

Reads a comparison table and reports every read whose alignment **differs**
between A and B. It also writes merged genomic regions showing where those
reads cluster. `compare` already writes this command's default-mode output;
run it standalone for the other modes, or to regenerate the output from an
existing table.

### `--space`

Which coordinate space defines a difference:

| Value | A read is identical when |
|---|---|
| `query` (default) | it is `query_identical`: same query span and same alignment relative to the read; a reverse-complement match counts as identical (see [Classification](#classification)) |
| `reference` | it is reference-identical: same `TargetChr` / `Strand` / `Target_Start` and same alignment content (per `--compare-by`). Strict and literal: no reverse-complement accommodation. |

### `--compare-by`

What counts as "same alignment" for reads mapped in **both** sets. Reads
mapped on only one side are always differences (`diff_aln_only_A` /
`diff_aln_only_B`), whatever the mode.

| Value | `--space query` | `--space reference` |
|---|---|---|
| `all` (default) | the full `cs` tag matches (motif-blind) | same position **and** same `cs` (motif-blind) |
| `junctions` | the query-space splice-junction set matches | same position **and** same genomic-coordinate junction set |

Under `junctions`, mismatches, indels and soft-clips that don't move a splice
junction don't count. Under both modes, comparison is motif-blind: a
differently reported intron motif at the same position and length is not a
difference (see [Classification](#classification)).

### Outputs

Filenames use the `query_diff` / `query_identical` stem under `--space query`
and `reference_diff` / `reference_identical` under `--space reference`. Under
`--compare-by junctions`, every filename gains a `.junctions` segment. So runs
in different modes at the same `--outdir` / `--prefix` never overwrite each
other. Read and region tables are gzipped unless `--no-gzip` is given.

| File (`--space query --compare-by all`) | Contents |
|---|---|
| `{prefix}.query_diff_reads.tsv.gz` | one row per differing read, with its category and 8 classification flags; see [query-diff-reads.md](output-tables/query-diff-reads.md) |
| `{prefix}.query_diff_regions.A.bed.gz`, `.B.bed.gz` | merged loci of the differing reads, per side; see [query-diff-regions.md](output-tables/query-diff-regions.md) |
| `{prefix}.query_diff_summary.tsv` | the [summary schema](output-tables/compare-summary.md), with `space` and `compare_by` rows after `label_A` / `label_B`; never gzipped |
| `{prefix}.query_identical_reads.tsv.gz` | with `--emit-identical-reads` only: the complementary set of identical reads, same layout as the diff-reads table except column 2 is named `category` |

Examples of other modes' names: `{prefix}.reference_diff_reads.tsv.gz`,
`{prefix}.query_diff_reads.junctions.tsv.gz`,
`{prefix}.reference_diff_regions.A.junctions.bed.gz`.

The summary counts don't depend on the mode. The number of rows in the
diff-reads file is printed on stderr. To get it from the summary, add the
mode's "not identical" row to `aligned_only_A + aligned_only_B`, e.g.
`query_not_identical + aligned_only_A + aligned_only_B` for the default mode.

### Categories

| `--space query` | `--space reference` | Meaning |
|---|---|---|
| `diff_aln_to_both` | `reference_diff` | mapped in both, not identical under the active mode |
| `diff_aln_only_A` / `diff_aln_only_B` | (same) | mapped in one set only |
| `query_identical_same_strand` / `query_identical_revcomp` / `query_identical_junctions` | `reference_identical` | identical under the active mode (identical-reads file only) |
| `neither_mapped` | (none) | mapped in neither set (identical-reads file only, and only under `--space query --compare-by all`) |

A read mapped in neither set counts as identical under
`--space query --compare-by all`: both aligners agree it doesn't map. In the
other modes it is in neither file, because those identity tests apply only to
reads mapped on both sides.

In the region tables, a read mapped on only one side contributes only to that
side's table. Intervals that can't be parsed, or where `end <= start`, are
skipped.

---

## `toolkit query-junction-diff`

```bash
maligno toolkit query-junction-diff -i <compare.parquet> --outdir <DIR> --prefix <STR>
```

For each read whose splice junctions differ between A and B, this
reconstructs every junction on each side and pairs its position on the read
(query space) with its position on the reference (genomic space). It then
reports the junctions the other side doesn't support, and how often each one
occurs across the whole table.

- **Input:** **Parquet only**, from `compare` or
  `toolkit merge-readinfo -o x.parquet`. The command reads the file twice, so
  stdin and TSV are not accepted.
- **Not run by `compare`:** this is a standalone command.
- **Coordinates:** junction positions are re-derived from each side's `cs`
  tag, so the query and genomic coordinates are correctly paired on both
  strands.

### Read selection

A read is "differing", and gets junction reconstruction, when:

```
mapped_a = TargetChr_A not empty or "*"      (mapped_b likewise)
differing =
    mapped_a && mapped_b    => N_Junctions_OnlyA > 0 || N_Junctions_OnlyB > 0
    mapped_a && !mapped_b   => JuncCount_A > 0
    !mapped_a && mapped_b   => JuncCount_B > 0
    !mapped_a && !mapped_b  => false
```

That covers two groups: reads mapped on both sides whose query-space junction
sets differ, and reads mapped on one side only that have at least one junction
there. A junction with nothing to compare against is unsupported by
definition. Every other read is still counted in the summary.

### Outputs

| File | Contents |
|---|---|
| `{prefix}.query_junction_diff.summary.tsv` | counts of reads at each selection step; see [query-junction-diff-summary.md](output-tables/query-junction-diff-summary.md) |
| `{prefix}.query_junction_diff.per_read_per_junc_info.tsv.gz` | one row per reconstructed junction, per side, per differing read; see [per-read-query-junction-diff.md](output-tables/per-read-query-junction-diff.md) |
| `{prefix}.query_junction_diff.unmatched_junctions.A.tsv.gz` | distinct junctions found in A but unsupported in B, with the number of differing reads carrying each one and the number of reads in the whole table carrying it; see [query-junction-diff-unmatched.md](output-tables/query-junction-diff-unmatched.md) |
| `{prefix}.query_junction_diff.unmatched_junctions.B.tsv.gz` | the same, for B unsupported in A |

The per-read and unmatched-junction tables are always gzipped; the summary
never is. In the per-read table, `junction_index` is 1-based and counted
separately on each side, in the read's 5'→3' order. It is not a matching key
across sides.

---

## Representative alignment

A read can have several alignments. Before comparing, each side's alignments
for a read are collapsed to one **representative (best) alignment**. This is
used by `compare` and `toolkit paf2tables --readinfo`.

- **Selection:** highest `ms`, then highest `AS`, then highest mapping quality
  (`MQ`). Full ties keep the first in input order, usually the aligner's own
  order.
- **The best alignment supplies** `TargetChr`, `Strand`, `MQ_Best`, `cs`,
  `junctions`, `genomic_junctions`, the event counts, and the spans
  `Query_Start` / `Query_End` and `Target_Start` / `Target_End`. Spans are
  0-based half-open (as in PAF and BED). Together with `TargetChr` and
  `Strand`, they give each read a BED-style interval.
- **Maxima over all the read's alignments:** `AS_Max`, `ms_Max`,
  `Query_Aln_Cov_Max`, `Query_Aln_Len_Max`, `seqid_Max`.
- **`Num_Aln`** counts the read's mapped alignments. A read that is present but
  unmapped has `Num_Aln = 0`, `TargetChr = "*"`, and zeroed or NaN statistics.
- **`Num_Aln_MaxScore`** counts the alignments tied with the best on both `ms`
  and `AS`. `1` means a clear winner; `> 1` means mapping quality or input
  order broke the tie. Aligners such as STAR write `ms = 0` for every
  alignment, so `AS` does the actual selection.
- **`Num_Aln_tpP` / `Num_Aln_tpS`** count mapped alignments by their `tp:A`
  tag: `P` or `I` (primary, including supplementary) and `S` or `i`
  (secondary). Alignments without a `tp` tag count toward neither.
- **`MQ_Best`** is the best alignment's mapping quality. For STAR, the common
  values are 255 (unique), 3 (2 loci), 1 (3 loci) and 0 (more than 3). A
  difference between sides points to reads where the two runs disagree on how
  unique the mapping is.

---

## Classification

Each read in the comparison table is classified from the two sides'
representative alignments, using `TargetChr`, `Strand`, `cs`,
`Query_Start` / `Query_End`, `Target_Start`, `junctions` and
`genomic_junctions`. These classes feed the summary statistics and
`toolkit find-aln-diff`. The summary rows are listed in
[compare-summary.md](output-tables/compare-summary.md).

**Mapping status.** A side with `TargetChr` equal to `*` or empty is
unmapped. Each read is `aligned_both`, `aligned_only_A`, `aligned_only_B` or
`aligned_neither`.

**Motif-blind `cs` comparison.** Wherever two `cs` strings are compared,
each intron's 2-letter donor/acceptor motif letters are ignored; the intron's
position and length still count. Aligners report motifs differently (STAR
writes `~nn<len>nn` where minimap2 writes the true motif, e.g. `~ct<len>ac`),
and that alone should not make identical alignments differ.

**Query-space identity (`query_identical`).** Both sides are mapped and cover
the same query span (`Query_Start` / `Query_End`, in forward-read
coordinates), and either:

- **same strand:** `Strand_A == Strand_B` and `cs_A == cs_B`
  (→ `query_identical_same_strand`), or
- **reverse complement:** `Strand_A != Strand_B` and `cs_A` equals the reverse
  complement of `cs_B` (→ `query_identical_revcomp`). This is the same
  read-to-reference correspondence on the opposite strand, e.g. a locus
  inverted between two assemblies.

Reads mapped in neither set also count as `query_identical`, since both sides
agree. So
`query_identical = query_identical_same_strand + query_identical_revcomp + aligned_neither`,
and `query_not_identical` means "mapped in both, not identical".

**Reference-space classes.** These apply to reads mapped in both sets. They
are strict and literal: there is no reverse-complement accommodation, so a
real strand difference always counts as a different position. There are two
axes:

- **position:** same `TargetChr`, `Strand` and `Target_Start`. With matching
  `cs`, the end must match too.
- **alignment:** same `cs`.

The four combinations are `ref_same_position_same_aln` (reference-identical),
`ref_same_position_diff_aln`, `ref_diff_position_same_aln` (relocated) and
`ref_diff_position_diff_aln`.

**Junction-set identity.** This is a looser test, for reads mapped in both
sets, that compares only the deduplicated splice-junction sets and ignores
mismatches, indels and soft-clips:

- **query space:** the `junctions` sets match → `query_junctions_identical`
  (otherwise `query_junctions_not_identical`).
- **genomic space:** among same-position reads, the `genomic_junctions` sets
  also match → `ref_same_position_same_junctions` (otherwise
  `ref_same_position_diff_junctions`).

---

## Junction comparison

Splice junctions are compared as **sets** (deduplicated on each side), in two
coordinate spaces. These columns are always present in the comparison table;
see [compare.md](output-tables/compare.md).

| Query space (`junctions`) | Genomic space (`genomic_junctions`) | Meaning |
|---|---|---|
| `N_Matched_Junctions` | `Genomic_N_Matched_Junctions` | junctions in both, \|A ∩ B\| |
| `N_Junctions_OnlyA` | `Genomic_N_Junctions_OnlyA` | only in A, \|A \ B\| |
| `N_Junctions_OnlyB` | `Genomic_N_Junctions_OnlyB` | only in B, \|B \ A\| |
| `N_Unmatched_Junctions` | `Genomic_N_Unmatched_Junctions` | in exactly one side (`OnlyA + OnlyB`) |
| `Junctions_OnlyA` / `Junctions_OnlyB` | `Genomic_Junctions_OnlyA` / `Genomic_Junctions_OnlyB` | the unmatched junctions themselves |

`N_Matched_Junctions + N_Junctions_OnlyA` equals A's junction count, and
likewise for B. The same holds for the genomic columns.

**Format.** Junction lists are written as Python tuples-of-tuples with
0-based half-open coordinates, e.g. `((1200, 1201), (3400, 3401))`, and parse
with `ast.literal_eval`. An empty set is `()`. `genomic_junctions` tuples are
`(start, end)` on the reference; the chromosome is in the row's `TargetChr`.
When genomic junctions are compared, each side's tuples are paired with that
side's own `TargetChr`, so junctions on different contigs never match.

---

## TSV and Parquet

The comparison table can be written as gzipped TSV, Parquet, or both, with the
same 100 columns, names and order (`pd.read_parquet` is a drop-in for
`pd.read_csv`). Parquet is much faster to read when you only need a few
columns. `toolkit summary` and `toolkit find-aln-diff` read either form;
`toolkit query-junction-diff` reads Parquet only.

**Nulls in Parquet mean "undefined", nothing else.** A Parquet null appears
where the TSV has `NaN` in a float column, or an empty `cs`. Meaningful values
stay values:

- `TargetChr` and `Strand` stay `*` for an unmapped side; they are the mapping
  indicator.
- An empty junction set stays `()`.
- A real `0` stays `0`.

A per-side value that is missing or can't be parsed becomes null, never a
made-up `0`.

**Escaping.** In the TSV, text fields escape backslash, tab, newline and
carriage return as `\\`, `\t`, `\n`, `\r`. Parquet stores the text unescaped.
This only matters for read names or labels containing those characters.

---

## Working with the output tables

**Select columns by name, not number.** Column positions can change when
columns are added, so look them up from the header. Two small shell helpers:

```bash
# tsvcut <file.gz> <Name1,Name2,...>: print just those columns, in that order.
tsvcut () {
  gzip -dc < "$1" | awk -F'\t' -v want="$2" '
    NR==1 { n=split(want,w,","); for (i=1;i<=NF;i++) h[$i]=i
            for (j=1;j<=n;j++) {
              if (!(w[j] in h)) { print "no such column: " w[j] > "/dev/stderr"; exit 1 }
              c[j]=h[w[j]]
            } }
    { line=$c[1]; for (j=2;j<=n;j++) line=line "\t" $c[j]; print line }'
}

# tsvwhere <file.gz> <Name> <pos|zero>: keep the header plus rows where that column is >0 (pos) or ==0 (zero).
tsvwhere () {
  gzip -dc < "$1" | awk -F'\t' -v col="$2" -v test="$3" '
    NR==1 { for (i=1;i<=NF;i++) if ($i==col) k=i
            if (!k) { print "no such column: " col > "/dev/stderr"; exit 1 }
            print; next }
    { v=$k+0 } test=="pos" ? v>0 : v==0'
}
```

(`gzip -dc` is used rather than `zcat` because macOS `zcat` rejects plain
`.gz` names.) Examples on a `compare` output:

```bash
CMP=results/AvsB.compare.tsv.gz

# Column number → name
gzip -dc < "$CMP" | head -1 | tr '\t' '\n' | nl

# Distribution of unmatched junctions, query space then genomic space
tsvcut "$CMP" N_Unmatched_Junctions         | tail -n +2 | sort | uniq -c
tsvcut "$CMP" Genomic_N_Unmatched_Junctions | tail -n +2 | sort | uniq -c

# Reads whose query-space junction sets disagree: count, then the junctions themselves
tsvwhere "$CMP" N_Unmatched_Junctions pos | tail -n +2 | wc -l
tsvwhere "$CMP" N_Unmatched_Junctions pos \
  | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i}
                {print $h["Read_Name"]"\t"$h["JuncCount_A"]"\t"$h["JuncCount_B"]"\t"$h["Junctions_OnlyA"]"\t"$h["Junctions_OnlyB"]}' \
  | column -t -s $'\t' | less -S

# Check column counts (alninfo 36, readinfo 36, comparison table 100)
gzip -dc < "$CMP" | awk -F'\t' '{print NF}' | sort | uniq -c
```

---

## Troubleshooting

### `merge-readinfo` compared fewer reads than expected

`toolkit merge-readinfo` assumes both readinfo files are sorted by `Read_Name`
in the **same byte order**. By default an out-of-order read is an error. With
`--allow-id-mismatch`, though, unsorted input silently lowers `reads_compared`.
Common causes are sorting without `LC_ALL=C`, or combining files from several
sources. (`compare` sorts its inputs itself and isn't affected.)

To check, count the read names the two files share. If it equals
`reads_compared`, the sort is fine and the low overlap is real.

```bash
LC_ALL=C comm -12 \
  <(gzip -dcf a.readinfo.tsv.gz | tail -n +2 | cut -f1 | LC_ALL=C sort -u) \
  <(gzip -dcf b.readinfo.tsv.gz | tail -n +2 | cut -f1 | LC_ALL=C sort -u) \
  | wc -l
```

[`scripts/check-readinfo-overlap.sh`](../scripts/check-readinfo-overlap.sh)
runs the same check and prints a report:

```bash
./scripts/check-readinfo-overlap.sh a.readinfo.tsv.gz b.readinfo.tsv.gz
```

### `compare` stops with "read-ID sets differ"

The two inputs don't contain the same reads. For example, an aligner run
without unmapped-read output (STAR without `--outSAMunmapped Within`) omits
reads that didn't align, and those differ between references. Either
regenerate the alignments with unmapped reads included, or pass
`--allow-id-mismatch` to compare the shared reads. With
`--allow-id-mismatch`, reads present on one side only are counted in the
summary but get no row.

### `sam2paf` / `compare` stops with a missing-`cs` error

A mapped record has neither a `cs` tag nor a usable `MD` tag. See
[`sam2paf`](#sam2paf) and the README's section on preprocessing STAR BAMs.

---

## Build

```bash
# Native (macOS / Linux)
cargo build --release
# → target/release/maligno

# Static Linux binary for HPC (no runtime dependencies).
# Requires the musl cross toolchain: brew install FiloSottile/musl-cross/musl-cross
# and: rustup target add x86_64-unknown-linux-musl
CC_x86_64_unknown_linux_musl=x86_64-linux-musl-gcc \
AR_x86_64_unknown_linux_musl=x86_64-linux-musl-ar \
cargo build --release --target x86_64-unknown-linux-musl
# → target/x86_64-unknown-linux-musl/release/maligno
```

The linker for the musl targets is set in `.cargo/config.toml`. The `CC_*` /
`AR_*` variables are also needed because maligno compiles a bundled htslib (C)
for SAM/BAM input.
