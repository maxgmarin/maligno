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
  - [`toolkit compare-readinfo`](#toolkit-compare-readinfo): two readinfo tables → comparison table
  - [`toolkit compare-summary`](#toolkit-compare-summary): comparison table → summary statistics
  - [`toolkit find-aln-diff`](#toolkit-find-aln-diff): comparison table → differing reads and regions
  - [`toolkit query-junction-compare`](#toolkit-query-junction-compare): comparison table → per-read splice junctions of the reads whose query-space junctions differ
  - [`toolkit genomic-junction-compare`](#toolkit-genomic-junction-compare): comparison table → A-vs-B junctions in reference coordinates, with strand and motifs
- Shared definitions
  - [Representative alignment](#representative-alignment)
  - [Transcript strand](#transcript-strand)
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
  A.paf|sam|bam ─┐                                       {prefix}.compare.parquet (and/or .tsv.gz)
                 ├─ compare (sort → verify → compare) ─▶  {prefix}.compare.summary.tsv
  B.paf|sam|bam ─┘
```

```bash
maligno compare -a A.bam -b B.bam --label-a A --label-b B --outdir results/ --prefix AvsB
```

**2. Building blocks:** `toolkit paf2tables`, then `toolkit compare-readinfo`.
You keep the intermediate per-read tables. `paf2tables` reads PAF, so SAM/BAM
is converted first with `sam2paf -p -U`, the same records `compare` keeps by
default.

```
  A.sam|bam ──sam2paf -p -U──▶ A.paf ──paf2tables──▶ A.readinfo.tsv ─┐
                                                                     ├─ compare-readinfo ─▶ compare.tsv
  B.sam|bam ──sam2paf -p -U──▶ B.paf ──paf2tables──▶ B.readinfo.tsv ─┘
```

```bash
maligno sam2paf -p -U A.bam | maligno toolkit paf2tables -i - --readinfo A.readinfo.tsv.gz
maligno sam2paf -p -U B.bam | maligno toolkit paf2tables -i - --readinfo B.readinfo.tsv.gz
maligno toolkit compare-readinfo -a A.readinfo.tsv.gz -b B.readinfo.tsv.gz --label-a A --label-b B -o AvsB.compare.parquet
```

Either way, the comparison table can then be analysed further with
`toolkit compare-summary`, `toolkit find-aln-diff`, `toolkit query-junction-compare` and
`toolkit genomic-junction-compare`.

**Input and output conventions.** Inputs and outputs ending in `.gz` are
read and written as gzip; `.gz` inputs may be single-member gzip or
multi-member (bgzip, pigz). Most commands accept `-` for stdin/stdout. The
exceptions are `compare`'s `-a`/`-b` (real file paths only), Parquet inputs
(real file paths only), and `toolkit find-aln-diff` and
`toolkit query-junction-compare` (Parquet only).

**Read identity.** maligno matches alignments across files by read name
(`query_name` in PAF, `QNAME` in SAM/BAM), so each name must identify one
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
to converting first. Unmapped reads are always kept. `--sam-records-to-convert`
picks which records are converted:

| `--sam-records-to-convert` | Records kept | Equivalent |
|---|---|---|
| `primary-supp` (default) | primary + supplementary; skip secondary | `sam2paf -p -U` |
| `primary` | primary only | `sam2paf -P -U` |
| `all` | every record | `sam2paf -U` |

This also decides which alignments the per-read `n_alns_primsupp` /
`n_alns_secondary` counts can see. Under the default, secondaries are dropped,
so `n_alns_secondary` is `0`.

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
   read is written to the comparison table. The summary statistics are
   produced in the same pass.

### Outputs

Written to `--outdir`, named by `--prefix`:

| File | Written | Spec |
|---|---|---|
| `{prefix}.compare.parquet` | by default (`--format parquet` or `both`) | [compare.md](output-tables/compare.md) |
| `{prefix}.compare.tsv.gz` | with `--format tsv` or `both` | [compare.md](output-tables/compare.md) |
| `{prefix}.compare.summary.tsv` | always | [compare-summary.md](output-tables/compare-summary.md) |
| `{prefix}.{label_a}.alninfo.tsv.gz`, `{prefix}.{label_b}.alninfo.tsv.gz` | with `--emit-alninfo` | [alninfo.md](output-tables/alninfo.md) |
| `{prefix}.{label_a}.readinfo.tsv.gz`, `{prefix}.{label_b}.readinfo.tsv.gz` | with `--emit-readinfo` | [readinfo.md](output-tables/readinfo.md) |

- `--format parquet|tsv|both` (default `parquet`) chooses the comparison table's
  serialization. Pass `--format tsv` or `both` for the gzipped TSV, e.g. for
  the shell recipes in [Working with the output tables](#working-with-the-output-tables).
- To list the reads whose alignments differ, and the regions where they
  cluster, run [`toolkit find-aln-diff`](#toolkit-find-aln-diff) on the
  Parquet table.
- The per-set alninfo and readinfo tables are off by default. They are the
  largest outputs, and no other maligno command reads them.
- An abbreviated summary is also printed to stderr.
- The sorted intermediate PAFs are deleted at the end unless you pass
  `--keep-sorted-paf`. For SAM/BAM inputs, these are the converted PAFs.

### Labels

`--label-a` / `--label-b` (defaults `SetA` / `SetB`) name the two sets. They
are written to the comparison table's `label_A` / `label_B` columns and to the
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
| `--sam-records-to-convert` | SAM/BAM only: `primary-supp` (default), `primary`, or `all` |
| `--label-a`, `--label-b` | set names (default `SetA` / `SetB`) |
| `-o`/`--outdir`, `-p`/`--prefix` | output directory and filename prefix |
| `--format` | comparison table as `parquet` (default), `tsv`, or `both` |
| `--allow-id-mismatch` | compare the shared reads instead of erroring when the read-name sets differ |
| `--presorted` | skip the sort; inputs must already be in the same read order |
| `--emit-alninfo`, `--emit-readinfo` | also write the per-set alninfo / readinfo tables |
| `--sort-mem` | in-memory sort buffer per file (default `1G`; `K`/`M`/`G` suffix) |
| `--sort-tmp-dir` | directory for sort spill files (default: `--outdir`) |
| `--sort-threads` | sort threads (default `1`) |
| `--keep-sorted-paf` | keep the sorted intermediate PAFs |

---

## `sam2paf`

```bash
maligno sam2paf [options] <in.sam|in.bam|->  > out.paf
```

Converts SAM or BAM to PAF on stdout. Its output is byte-for-byte compatible
with the `sam2paf` command of minimap2's `paftools.js`, with one addition: the
transcript-strand tags `ts:A` (minimap2) and `XS:A` (STAR, HISAT2) are copied
to the PAF record when present, so the [strand columns](#transcript-strand)
survive conversion. The input format is
detected automatically; `-` reads SAM text from stdin. `maligno sam2paf in.bam`
gives the same bytes as `samtools view -h in.bam | maligno sam2paf -`.

| Flag | Meaning |
|---|---|
| `-p`, `--primary-or-supp` | primary + supplementary alignments only (skip secondary, FLAG 0x100) |
| `-P`, `--primary-only` | primary alignments only (skip secondary and supplementary, FLAG 0x800); implies `-p` |
| `-L`, `--long-cs` | write the `cs` tag in long form (`=ACGT`) instead of the default short form (`:N`) |
| `-U`, `--unaligned` | write a placeholder PAF record for each unmapped read (otherwise unmapped reads are dropped) |

Pass `-U` when the PAF will be compared: unmapped reads then appear as
`n_alns = 0` rows instead of disappearing.

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
maligno toolkit paf2tables -i <in.paf[.gz]|-> [--alninfo <alninfo.tsv[.gz]>] [-r/--readinfo <readinfo.tsv[.gz]>] \
  [--presorted] [--sort-mem 1G] [--sort-tmp-dir DIR] [--sort-threads 1] [--keep-sorted-paf]
```

Turns a PAF into maligno's tables, writing either or both in a single pass.
At least one output must be given.

- **`--alninfo`** writes one row per PAF record (36 columns; see
  [alninfo.md](output-tables/alninfo.md)). Each record's `cs` tag is walked to
  count matches, substitutions, insertions, deletions and splice junctions;
  soft-clip lengths, junction coordinates and identity/coverage are derived
  from it. Unmapped records (`target_name == "*"`) get a row with zeroed
  statistics.
- **`--readinfo`** writes one row per read (37 columns; see
  [readinfo.md](output-tables/readinfo.md)): the read's
  [representative alignment](#representative-alignment) plus per-read
  aggregates.

**Sorting.** The PAF is first sorted by `query_name` in byte order, with the
same sort `compare` uses, so both tables are sorted by read and the readinfo
is ready for `toolkit compare-readinfo`. The tables are identical to
`compare --emit-alninfo` / `--emit-readinfo` for the same input.

- **`--sort-mem`** caps the memory the sort uses. Larger inputs are sorted in
  chunks that are written to **`--sort-tmp-dir`** and merged. The default
  directory is the folder of `--readinfo`, or of `--alninfo` when `--readinfo`
  is not given.
- **`--keep-sorted-paf`** keeps the sorted PAF as `<input stem>.sorted.paf.gz`
  in that same folder, where the stem is the input file name without `.gz`
  and `.paf`. It needs a file `--input`, not stdin.
- **`--presorted`** skips the sort for a PAF that is already sorted by
  `query_name` in byte order (as `LC_ALL=C sort -t$'\t' -k1,1` sorts it). The
  order is checked as the PAF is read: the command stops with an error at the
  first read name that is lower than the one before it, and removes the
  partial outputs.

---

## `toolkit compare-readinfo`

```bash
maligno toolkit compare-readinfo -a/--readinfo-a <readinfo_a.tsv[.gz]> -b/--readinfo-b <readinfo_b.tsv[.gz]> \
  [--label-a A] [--label-b B] -o <compare.parquet|compare.tsv[.gz]> [--allow-id-mismatch]
```

Joins two readinfo tables into the alignment comparison table. It runs the same
pairing, [classification](#classification) and row writing as `compare`, so the
table is identical to the one `compare` writes from the same reads (see
[compare.md](output-tables/compare.md)). Output is Parquet when `-o` ends in
`.parquet`, otherwise TSV (`.tsv[.gz]`).

**Required columns.** Each input must carry `query_name`, `query_len` and every
readinfo column that goes into the comparison table (all of
[readinfo.md](output-tables/readinfo.md)), in any order. A
missing column is an error naming the file and the column, as `compare`'s own
readinfo rows always carry them all.

**Matching.** Reads are matched on **`query_name`** with a streaming merge:
constant memory, one pass over each file. If a read's `query_len` differs between
the two files, the output carries side A's. Both inputs must be **sorted by
`query_name` in byte order**. `toolkit paf2tables` writes them that way. Sort
readinfo files from elsewhere with the header kept on top:

```bash
( gzip -dc in.readinfo.tsv.gz | head -1
  gzip -dc in.readinfo.tsv.gz | tail -n +2 | LC_ALL=C sort -t$'\t' -k1,1
) | gzip > sorted.readinfo.tsv.gz
```

- `LC_ALL=C` forces byte order, regardless of locale.
- `-k1,1` sorts on `query_name`.
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

## `toolkit compare-summary`

```bash
maligno toolkit compare-summary -i <compare.parquet|compare.tsv[.gz]|-> [-o summary.tsv[.gz]] [--input-format auto|tsv|parquet]
```

Computes the summary statistics from an existing comparison table. The schema
is the same as `{prefix}.compare.summary.tsv`; see
[compare-summary.md](output-tables/compare-summary.md). The result is printed
to stderr, and also written to `-o` when given (`.tsv[.gz]`, or `-` for stdout).

- **Input format:** `--input-format auto` (the default) reads Parquet for a
  `.parquet` path and TSV otherwise. Parquet must be a real file, not stdin.
- **Columns:** every column is looked up by name, so column order doesn't
  matter. The set names come from `label_A` / `label_B`.
- **Reads in only one input:** the table only contains reads present in both
  inputs, so `present_only_in_A_by_id` / `present_only_in_B_by_id` are always
  `0` here. Only `compare` and `toolkit compare-readinfo` can fill them in.

The categories are defined in [Classification](#classification).

---

## `toolkit find-aln-diff`

```bash
maligno toolkit find-aln-diff -i <compare.parquet> -o/--outdir <DIR> -p/--prefix <STR> \
  [--mode any-diff-query|any-diff-reference] [--aligned-both-only] [--emit-identical-reads] [--no-gzip]
```

Reads a comparison table (**Parquet only**, from `compare` or
`toolkit compare-readinfo -o x.parquet`) and reports every read whose core
alignment **differs** between A and B. It also writes merged genomic regions
showing where those reads cluster. MAPQ and score differences don't count.

### `--mode`

When a read mapped in both sets counts as **identical**:

| Mode | Identical when |
|---|---|
| `any-diff-query` (default) | same query span (`query_start`, `query_end`) and same `cs`, motif-blind. On opposite strands, B's `cs` must be the exact reverse complement of A's. |
| `any-diff-reference` | the same as `any-diff-query` on the **same strand**, and also the same `target_name` and `target_start`. Use it when A and B were aligned to the same reference. |

In both modes, a read mapped in only one set is a difference
(`diff_aln_only_A` / `diff_aln_only_B`), and a read mapped in neither set is
identical (`neither_mapped`): both aligners agree it doesn't map. The `cs`
comparison is motif-blind (see [Classification](#classification)).

**`--aligned-both-only`** compares only the reads mapped in both sets. The
others are left out of every output and counted in the summary as
`reads_not_aligned_both_excluded`.

### Outputs

Filenames use the `query_aln_diff` stem with `--mode any-diff-query` and
`reference_aln_diff` with `--mode any-diff-reference`, so runs in the two modes
at the same `--outdir` / `--prefix` never overwrite each other. Read and region
tables are gzipped unless `--no-gzip` is given.

| File (`--mode any-diff-query`) | Contents |
|---|---|
| `{prefix}.query_aln_diff_reads.tsv.gz` | one row per differing read, with its `outcome` and 7 classification columns; see [aln-diff-reads.md](output-tables/aln-diff-reads.md) |
| `{prefix}.query_aln_diff_regions.A.bed.gz`, `.B.bed.gz` | merged loci of the differing reads, per side; see [aln-diff-regions.md](output-tables/aln-diff-regions.md) |
| `{prefix}.query_aln_diff_summary.tsv` | the mode and parameters, then the counts of identical and differing reads; see [aln-diff-summary.md](output-tables/aln-diff-summary.md) |
| `{prefix}.query_aln_identical_reads.tsv.gz` | with `--emit-identical-reads` only: the identical reads, in the same layout |

With `--mode any-diff-reference` the names are
`{prefix}.reference_aln_diff_reads.tsv.gz`,
`{prefix}.reference_aln_diff_regions.{A,B}.bed.gz`,
`{prefix}.reference_aln_diff_summary.tsv` and
`{prefix}.reference_aln_identical_reads.tsv.gz`.

### Outcomes

The `outcome` column is the same in both modes:

| `outcome` | Meaning | File |
|---|---|---|
| `diff_aln_both` | mapped in both sets, not identical under the mode | differing reads |
| `diff_aln_only_A` / `diff_aln_only_B` | mapped in one set only | differing reads |
| `identical_same_strand` | identical, on the same strand | identical reads |
| `identical_revcomp` | identical as an exact reverse complement on opposite strands (`any-diff-query` only) | identical reads |
| `neither_mapped` | mapped in neither set | identical reads |

In the region tables, a read mapped on only one side contributes only to that
side's table. Intervals that can't be parsed, or where `end <= start`, are
skipped.

---

## `toolkit query-junction-compare`

```bash
maligno toolkit query-junction-compare -i <compare.parquet> -o/--outdir <DIR> -p/--prefix <STR>
```

For each read whose splice junctions differ between A and B, this
reconstructs every junction on each side and pairs its position on the read
(query space) with its position on the reference (genomic space). It then
reports the junctions the other side doesn't support, and how often each one
occurs across the whole table.

- **Input:** **Parquet only**, from `compare` or
  `toolkit compare-readinfo -o x.parquet`. The command reads the file twice, so
  stdin and TSV are not accepted.
- **Not run by `compare`:** this is a standalone command.
- **Coordinates:** junction positions are re-derived from each side's `cs`
  tag, so the query and genomic coordinates are correctly paired on both
  strands.

### Read selection

A read is "differing", and gets junction reconstruction, when:

```
mapped_a = target_name_A not empty or "*"      (mapped_b likewise)
differing =
    mapped_a && mapped_b    => n_query_junctions_only_A > 0 || n_query_junctions_only_B > 0
    mapped_a && !mapped_b   => n_intron_events_A > 0
    !mapped_a && mapped_b   => n_intron_events_B > 0
    !mapped_a && !mapped_b  => false
```

That covers two groups: reads mapped on both sides whose query-space junction
sets differ, and reads mapped on one side only that have at least one junction
there. A junction with nothing to compare against is unsupported by
definition. Every other read is still counted in the summary.

### Outputs

| File | Contents |
|---|---|
| `{prefix}.query_junction_compare.summary.tsv` | counts of reads at each selection step, then distinct junctions per side (unique, shared, only in A or B, unmatched); see [query-junction-compare-summary.md](output-tables/query-junction-compare-summary.md) |
| `{prefix}.query_junction_compare.per_read_per_junc_info.tsv.gz` | one row per reconstructed junction, per side, per differing read; see [per-read-query-junction-compare.md](output-tables/per-read-query-junction-compare.md) |
| `{prefix}.query_junction_compare.unmatched_junctions.A.tsv.gz` | distinct junctions found in A but unsupported in B, with the number of differing reads carrying each one and the number of reads in the whole table carrying it; see [query-junction-compare-unmatched.md](output-tables/query-junction-compare-unmatched.md) |
| `{prefix}.query_junction_compare.unmatched_junctions.B.tsv.gz` | the same, for B unsupported in A |

The per-read and unmatched-junction tables are always gzipped; the summary
never is. In the per-read table, `junction_index` is 1-based and counted
separately on each side, in the read's 5'→3' order. It is not a matching key
across sides.

---

## `toolkit genomic-junction-compare`

```bash
maligno toolkit genomic-junction-compare -i <compare.parquet|compare.tsv[.gz]> \
  [--input-format auto|tsv|parquet] -o/--outdir <DIR> -p/--prefix <STR>
```

Compares A's and B's splice junctions in **genomic (reference)
coordinates**, read by read and junction by junction, with each junction's
strand and splice-site motif. The comparison table itself compares junctions
only in query space (positions in the read); this command is where the
reference-coordinate comparison lives.

**Only meaningful when A and B were aligned to the same reference.** With
different references (another assembly, or a reverse-complemented one), the
coordinates aren't comparable; don't run it.

**What a genomic junction is.** `(target_name, intron_start, intron_end,
intron_strand)`, 0-based half-open:
- `intron_strand` is the alignment's `genomic_intron_strand`: the strand of
  its transcript on the reference, `*` when unknown (see
  [Transcript strand](#transcript-strand)). Every junction of an alignment
  gets that alignment's intron strand.
- The same coordinates on `+` and `-` are **different junctions**, and `*`
  matches only `*`. Run aligners so they write the strand tags (STAR needs
  `--outSAMstrandField intronMotif`), or every junction is `*`.

**How junctions are found.** Each side's junctions are rebuilt from its
representative alignment's `cs` tag (they equal the table's
`genomic_junctions` columns). The same walk gives each intron's **splice-site
motif**: the two letters at each end of the intron (`~gt…ag` in `cs`),
written donor-acceptor in the junction's orientation, so a canonical intron
is `gt-ag` on either strand. Aligners that don't report the bases give
`nn-nn`.

- **Input:** the comparison table from `compare` or `toolkit compare-readinfo`,
  Parquet or TSV. It reads only the columns it needs.
- **Not run by `compare`:** this is a standalone command.

### Outputs

| File | Contents |
|---|---|
| `{prefix}.genomic_junction_compare.per_read.tsv.gz` | one row per read with a junction on either side: each side's locus, strands, junctions as `(start, end)` pairs, motifs and `cs`, plus shared / only-A / only-B counts and lists; see [genomic-junction-compare-per-read.md](output-tables/genomic-junction-compare-per-read.md) |
| `{prefix}.genomic_junction_compare.per_junction.tsv.gz` | one row per distinct junction: reads supporting it on each side, reads whose A and B alignments agree on it, and each side's motif; see [genomic-junction-compare-per-junction.md](output-tables/genomic-junction-compare-per-junction.md) |
| `{prefix}.genomic_junction_compare.summary.tsv` | reads by which sides have junctions, identical vs differing junction sets, and distinct-junction totals; see [genomic-junction-compare-summary.md](output-tables/genomic-junction-compare-summary.md) |

The per-read and per-junction tables are always gzipped; the summary never is.

---

## Representative alignment

A read can have several alignments. Before comparing, each side's alignments
for a read are collapsed to one **representative (best) alignment**. This is
used by `compare` and `toolkit paf2tables --readinfo`.

- **Selection:** highest `ms_score` (the `ms:i` tag), then highest `as_score`
  (`AS:i`), then highest `mapq`, then the alignment type from `tp:A`: primary
  (`P`/`I`) before untagged before secondary (`S`/`i`), and any mapped row
  before an unmapped one. Remaining ties keep the first in input order.
  - The type key only breaks **full** score ties. Among multimappers that
    score the same (STAR gives all of them the same `AS` and MAPQ), the
    aligner's own primary represents the read; a higher-scoring secondary
    still wins.
  - minimap2 writes `tp:A` on every PAF line. For SAM/BAM input,
    [`sam2paf`](#sam2paf) derives it from the FLAG: `S` when `0x100`
    (secondary) is set, otherwise `P`, so supplementary segments are `P` too.
  - A PAF without `tp:A` tags leaves every row "untagged": the type key then
    never decides, and input order breaks the tie.
  - The representative's type is the `tp_tag` column (readinfo, and
    `tp_tag_A` / `tp_tag_B` in the comparison table).
- **Every per-alignment value comes from the representative alone:**
  `target_name`, `strand`, `query_transcript_strand`, `genomic_intron_strand`, `mapq`, `as_score`, `ms_score`, `aln_identity`,
  `query_aln_cov`, `query_aln_len`, `cs`, `query_junctions`,
  `genomic_junctions`, the event counts, and the spans `query_start` /
  `query_end` and `target_start` / `target_end`. So every per-side value in
  the comparison table, and every A-vs-B difference, describes the same two
  alignments. Spans are 0-based half-open (as in PAF and BED). Together with
  `target_name` and `strand`, they give each read a BED-style interval.
- **`n_alns`** counts the read's mapped alignments, including supplementary
  segments. A read that is present but unmapped has `n_alns = 0`,
  `target_name = "*"`, and zeroed or NaN statistics.
- **`n_alns_tied_top_score`** counts the mapped alignments tied with the
  representative on the three score keys: `ms_score`, `as_score` and `mapq`.
  The `tp:A` type is not one of them, so a primary and a secondary with equal
  scores count as tied. `1` means the representative won outright; `> 1`
  means that many alignments scored the same and the type (then input order)
  picked one. Aligners such
  as STAR write `ms = 0` for every alignment, so `as_score` does the actual
  selection.
- **`n_alns_primsupp` / `n_alns_secondary`** count mapped alignments by their `tp:A`
  tag: `P` or `I` (primary, including supplementary) and `S` or `i`
  (secondary). Alignments without a `tp` tag count toward neither.
- **`mapq`** is the representative alignment's mapping quality. For STAR, the common
  values are 255 (unique), 3 (2 loci), 1 (3 loci) and 0 (more than 3). A
  difference between sides points to reads where the two runs disagree on how
  unique the mapping is.

---

## Transcript strand

Two columns record the strand of the transcript a read came from, as
inferred by the aligner from the splice signals. Each is a column of alninfo,
readinfo (the representative alignment's value) and both per-side blocks of
the comparison table (`_A` / `_B`). Their names say which frame they are in:

| Column | Relative to | `+` | `-` | `*` |
|---|---|---|---|---|
| `query_transcript_strand` | the **read** | the read is in the transcript's sense orientation | the read is antisense to the transcript | unknown: no strand tag, or unmapped |
| `genomic_intron_strand` | the **reference** | the transcript, and so every intron of the alignment, is on the reference `+` strand | on the `-` strand | unknown |

**It is relative to the read, not the reference.** The same read keeps the
same value whichever reference, or reference orientation, it is aligned to.
For example, aligned to HG38 a read maps to `strand` `+`; aligned to the
reverse-complemented HG38 it maps to `strand` `-`, and its
`query_transcript_strand` is the same in both. This keeps it comparable between
two alignments to different references, like the query-space junctions.

**Where it comes from.** maligno reads one of two SAM/PAF tags:

| Aligner | Tag | Tag's convention | `query_transcript_strand` |
|---|---|---|---|
| minimap2, splice presets (`-x splice`, `splice:hq`, `splice:sr`) | `ts:A` | relative to the read | the tag's value |
| STAR with `--outSAMstrandField intronMotif` | `XS:A` | relative to the reference | `+` if `XS` equals the alignment `strand`, otherwise `-` |
| other aligners writing Cufflinks-style `XS:A` (e.g. HISAT2) | `XS:A` | relative to the reference | as for STAR; expected to work, not tested |

- If a record carries both tags, `ts:A` is used.
- Any other tag value, or no tag, gives `*`.
- For SAM/BAM input, [`sam2paf`](#sam2paf) copies both tags into the PAF.
- Aligners only infer a strand for **spliced** alignments, so unspliced reads
  are usually `*` (most of a short-read dataset):
  - minimap2 tags only the spliced alignment;
  - STAR tags both mates of a pair when either mate is spliced.
- STAR writes `XS:A` only when run with `--outSAMstrandField intronMotif`.
  Without it every value is `*`.

**Genomic intron strand.** `genomic_intron_strand` is `strand` when
`query_transcript_strand` is `+`, the opposite of `strand` when it is `-`, and
`*` when either is unknown. For example, `strand` `-` with
`query_transcript_strand` `-` is a transcript on the reference's `+` strand.
It is relative to the reference, so use it with genomic coordinates
(annotation strand, BED strand,
[`toolkit genomic-junction-compare`](#toolkit-genomic-junction-compare)).

For a mapped alignment, any two of `strand`, `query_transcript_strand` and
`genomic_intron_strand` give the third. Both strand columns are written so
each use reads one column directly.

---

## Classification

Each read in the comparison table is classified from the two sides'
representative alignments, using `target_name`, `strand`, `cs`,
`query_start` / `query_end`, `target_start` / `target_end` and
`query_junctions`. These classes feed the summary statistics and
`toolkit find-aln-diff`'s modes. The summary rows are listed in
[compare-summary.md](output-tables/compare-summary.md).

**Mapping status.** A side with `target_name` equal to `*` or empty is
unmapped. Each read is `aligned_both`, `aligned_only_A`, `aligned_only_B` or
`aligned_neither`.

**Motif-blind `cs` comparison.** Wherever two `cs` strings are compared,
each intron's 2-letter donor/acceptor motif letters are ignored; the intron's
position and length still count. Aligners report motifs differently (STAR
writes `~nn<len>nn` where minimap2 writes the true motif, e.g. `~ct<len>ac`),
and that alone should not make identical alignments differ.

**Query-space identity (`query_identical`).** Both sides are mapped and cover
the same query span (`query_start` / `query_end`, in forward-read
coordinates), and either:

- **same strand:** `strand_A == strand_B` and `cs_A == cs_B`
  (→ `query_identical_same_strand`), or
- **reverse complement:** `strand_A != strand_B` and `cs_A` equals the reverse
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

- **position:** same `target_name`, `strand` and `target_start`. With matching
  `cs`, the end must match too.
- **alignment:** same `cs`.

The four combinations are `ref_same_position_same_aln` (reference-identical),
`ref_same_position_diff_aln`, `ref_diff_position_same_aln` (relocated) and
`ref_diff_position_diff_aln`.

These compare positions on each side's own reference, so when A and B were
aligned to different references (another assembly, or a reverse-complemented
one) every read is a different position.

**Junction-set identity.** This is a looser test, for reads mapped in both
sets: the deduplicated `query_junctions` sets match → `query_junctions_identical`
(otherwise `query_junctions_not_identical`). Mismatches, indels and soft-clips
are ignored. Genomic junction sets are not part of the classification; they
are compared by [`toolkit genomic-junction-compare`](#toolkit-genomic-junction-compare).

---

## Junction comparison

The comparison table compares splice junctions as **sets** (deduplicated on
each side) in **query space**: each junction is its position in the read, so
the comparison doesn't depend on the reference. See
[compare.md](output-tables/compare.md).

| Column | Meaning |
|---|---|
| `n_query_junctions_shared` | junctions in both, \|A ∩ B\| |
| `n_query_junctions_only_A` | only in A, \|A \ B\| |
| `n_query_junctions_only_B` | only in B, \|B \ A\| |
| `n_query_junctions_not_shared` | in exactly one side (`only_A + only_B`) |
| `query_junctions_only_A` / `query_junctions_only_B` | the unmatched junctions themselves |

`n_query_junctions_shared + n_query_junctions_only_A` equals A's number of
distinct junctions, and likewise for B.

Junctions in **genomic** (reference) coordinates are kept per side in
`genomic_junctions_A` / `genomic_junctions_B`; comparing them between A and B
is done by [`toolkit genomic-junction-compare`](#toolkit-genomic-junction-compare).

**Format.** Junction lists are written as Python tuples-of-tuples with
0-based half-open coordinates, e.g. `((1200, 1201), (3400, 3401))`, and parse
with `ast.literal_eval`. An empty set is `()`. `genomic_junctions` tuples are
`(start, end)` on the reference; the chromosome is in the row's `target_name`.

---

## TSV and Parquet

The comparison table can be written as Parquet (`compare`'s default), gzipped
TSV, or both, with the same 93 columns, names and order (`pd.read_parquet` is
a drop-in for `pd.read_csv`). Parquet is much faster to read when you only need
a few columns. `toolkit compare-summary` and
`toolkit genomic-junction-compare` read either form;
`toolkit find-aln-diff` and `toolkit query-junction-compare` read Parquet only.

**Nulls in Parquet mean "undefined", nothing else.** A Parquet null appears
where the TSV has `NaN` in a float column, or an empty `cs`. Meaningful values
stay values:

- `target_name` and `strand` stay `*` for an unmapped side; they are the mapping
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

These recipes read the gzipped TSV comparison table, so run `compare` with
`--format tsv` or `--format both`.

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

# Distribution of unmatched query-space junctions per read
tsvcut "$CMP" n_query_junctions_not_shared | tail -n +2 | sort | uniq -c

# Reads whose query-space junction sets disagree: count, then the junctions themselves
tsvwhere "$CMP" n_query_junctions_not_shared pos | tail -n +2 | wc -l
tsvwhere "$CMP" n_query_junctions_not_shared pos \
  | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i}
                {print $h["query_name"]"\t"$h["n_intron_events_A"]"\t"$h["n_intron_events_B"]"\t"$h["query_junctions_only_A"]"\t"$h["query_junctions_only_B"]}' \
  | column -t -s $'\t' | less -S

# Check the comparison table's column count (93; alninfo has 36, readinfo 37)
gzip -dc < "$CMP" | awk -F'\t' '{print NF}' | sort | uniq -c
```

---

## Troubleshooting

### `compare-readinfo` compared fewer reads than expected

`toolkit compare-readinfo` assumes both readinfo files are sorted by `query_name`
in the **same byte order**. By default an out-of-order read is an error. With
`--allow-id-mismatch`, though, unsorted input silently lowers `reads_compared`.
Common causes are readinfo files sorted outside maligno without `LC_ALL=C`, or
combining files from several sources. (`compare` and `toolkit paf2tables` sort
by `query_name` themselves, so their output isn't affected.)

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
