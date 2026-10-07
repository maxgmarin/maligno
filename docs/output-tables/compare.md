# `compare` table

`{prefix}.compare.tsv.gz` and/or `{prefix}.compare.parquet` — one row per
read, 100 columns: both sides' representative-alignment stats plus a computed
A-vs-B comparison block. This is the **primary output** of `compare` (and of
`toolkit merge-readinfo`, which produces the identical table from
readinfo TSVs directly).

Columns are organized in groups, left to right:

| Group | Cols | What it holds |
|-------|:----:|---------------|
| **Join keys** | 1–2 | `Read_Name`, `Read_Len` |
| **Set labels** | 3–4 | `Label_A`, `Label_B`: the `--label-a` / `--label-b` values, repeated on every row |
| **Per-side data, A** | 5–37 | the representative alignment's stats for set A, each column suffixed `_A` |
| **Per-side data, B** | 38–70 | the same 33 columns for set B, suffixed `_B` |
| **Comparison metrics** | 71–96 | A-vs-B differences and ratios: `Strand_Match`, `seqid_Diff`, coverage/length diffs, score diffs (`AS_Diff`, `ms_Diff`, …), indel/soft-clip diffs, and junction-set counts in both query and genomic space |
| **Non-overlap objects** | 97–100 | the junctions that failed to overlap: `Junctions_OnlyA/B` and `Genomic_Junctions_OnlyA/B` |

Within each per-side block the columns are grouped by topic: locus and span,
alignment selection and score (including the per-type alignment counts
`Num_Aln_tpP` / `Num_Aln_tpS`), identity and coverage, junction counts,
cs-derived event counts, then the three long strings (`junctions`,
`genomic_junctions`, `cs`) last. Print the exact layout of any table with
`gzip -dc … | head -1 | tr '\t' '\n' | nl`.

Per-side columns always use the **fixed** `_A` / `_B` suffixes, never the
dataset label, so column names are identical for every comparison and stay
unambiguous even when a label contains an underscore. Which dataset each side
*is* comes from the `Label_A` / `Label_B` columns, carried on every row so any
subset of rows stays self-describing.

## Join keys & labels (columns 1–4)

| # | Column | Description |
|---|---|---|
| 1 | `Read_Name` | Read ID — join key |
| 2 | `Read_Len` | Query length — join key |
| 3 | `Label_A` | The `--label-a` value, repeated on every row |
| 4 | `Label_B` | The `--label-b` value, repeated on every row |

## Per-side data (columns 5–37 suffixed `_A`, 38–70 suffixed `_B`)

Same 33 columns, once per side — read out of each side's `readinfo` row by
name (see [readinfo.md](readinfo.md) for how each is computed; `tp_tag` and
`Num_Residue_Matches`-style raw-PAF columns are readinfo/alninfo-only and do
**not** appear here). An unmapped side has `TargetChr_* == "*"`, an empty
`cs_*`, and `NaN` for `seqid_Max_*`/`Query_Aln_Cov_Max_*`.

| # (`_A` / `_B`) | Column | Description |
|---|---|---|
| 5 / 38 | `TargetChr` | Representative alignment's reference/target name; `*` if unmapped |
| 6 / 39 | `Strand` | Representative alignment's strand; `*` if unmapped |
| 7 / 40 | `Target_Start` | Target-coordinate span start (0-based half-open) |
| 8 / 41 | `Target_End` | Target-coordinate span end |
| 9 / 42 | `Query_Start` | Query-coordinate span start (0-based half-open) |
| 10 / 43 | `Query_End` | Query-coordinate span end |
| 11 / 44 | `MQ_Best` | Representative alignment's mapping quality |
| 12 / 45 | `Num_Aln` | Total alignment count for this read on this side |
| 13 / 46 | `Num_Aln_tpP` | Mapped alignments for this read on this side whose PAF `tp:A` is `P` or `I` (SAM primary **+ supplementary**); `tpP − 1` = supplementary segments. See [readinfo.md](readinfo.md) |
| 14 / 47 | `Num_Aln_tpS` | Mapped alignments for this read on this side whose PAF `tp:A` is `S` or `i` (secondary) |
| 15 / 48 | `Num_Aln_MaxScore` | Alignments tied at the winning `(ms, AS)` score |
| 16 / 49 | `AS_Max` | Highest `AS` across this side's alignments |
| 17 / 50 | `ms_Max` | Highest `ms` across this side's alignments |
| 18 / 51 | `seqid_Max` | Highest sequence identity across this side's alignments |
| 19 / 52 | `Query_Aln_Cov_Max` | Highest query-alignment coverage |
| 20 / 53 | `Query_Aln_Len_Max` | Highest aligned query length |
| 21 / 54 | `JuncCount` | Number of query-space splice junctions in the representative alignment |
| 22 / 55 | `N_Splice_Junction_Events` | Representative alignment's intron/splice-op count |
| 23 / 56 | `N_Splice_Junction_Bases` | Representative alignment's total intron bases |
| 24 / 57 | `N_Match_Events` | Representative alignment's match-run count |
| 25 / 58 | `N_Match_Bases` | Representative alignment's matched bases |
| 26 / 59 | `N_Substitution_Events` | Representative alignment's substitution-event count |
| 27 / 60 | `N_Substitution_Bases` | Representative alignment's substituted bases |
| 28 / 61 | `N_Insertion_Events` | Representative alignment's insertion-event count |
| 29 / 62 | `N_Insertion_Bases` | Representative alignment's inserted bases |
| 30 / 63 | `N_Deletion_Events` | Representative alignment's deletion-event count |
| 31 / 64 | `N_Deletion_Bases` | Representative alignment's deleted bases |
| 32 / 65 | `N_SoftClipped_Events` | Representative alignment's soft-clip event count (`0`/`1`/`2`) |
| 33 / 66 | `N_SoftClipped_Bases_Start` | Representative alignment's start soft-clip length |
| 34 / 67 | `N_SoftClipped_Bases_End` | Representative alignment's end soft-clip length |
| 35 / 68 | `junctions` | Representative alignment's query-coordinate splice junctions (Python tuple string) |
| 36 / 69 | `genomic_junctions` | Representative alignment's splice junctions in absolute reference coordinates (Python tuple-of-tuples string; chrom is the sibling `TargetChr` column) |
| 37 / 70 | `cs` | Representative alignment's compact alignment string |

## Comparison block (columns 71–100)

A-vs-B metrics computed from the two per-side blocks above. **Sign
convention**: every `*_Diff` is **B − A**, every `*_Ratio` is **B / A**
(`NaN` when the A-side denominator is `0`) — except `Junc_Dist_V2`, whose
inner subtraction is A − B but is then made symmetric by `unsigned_abs()`.

| # | Column | Description |
|---|---|---|
| 71 | `Strand_Match` | `Strand_A == Strand_B` |
| 72 | `seqid_Diff` | `seqid_Max_B - seqid_Max_A` |
| 73 | `QueryAlnCov_Diff` | `Query_Aln_Cov_Max_B - Query_Aln_Cov_Max_A` |
| 74 | `QueryAlnLen_Diff` | `Query_Aln_Len_Max_B - Query_Aln_Len_Max_A` |
| 75 | `AS_Diff` | `AS_Max_B - AS_Max_A` |
| 76 | `ms_Diff` | `ms_Max_B - ms_Max_A` |
| 77 | `AS_Ratio` | `AS_Max_B / AS_Max_A` |
| 78 | `ms_Ratio` | `ms_Max_B / ms_Max_A` |
| 79 | `N_Substitution_Bases_Diff` | `N_Substitution_Bases_B - N_Substitution_Bases_A` |
| 80 | `N_Substitution_Bases_Ratio` | `N_Substitution_Bases_B / N_Substitution_Bases_A` |
| 81 | `N_Insertion_Bases_Diff` | `N_Insertion_Bases_B - N_Insertion_Bases_A` |
| 82 | `N_Insertion_Bases_Ratio` | `N_Insertion_Bases_B / N_Insertion_Bases_A` |
| 83 | `N_Deletion_Bases_Diff` | `N_Deletion_Bases_B - N_Deletion_Bases_A` |
| 84 | `N_Deletion_Bases_Ratio` | `N_Deletion_Bases_B / N_Deletion_Bases_A` |
| 85 | `N_SoftClipped_Bases_Start_Diff` | `N_SoftClipped_Bases_Start_B - N_SoftClipped_Bases_Start_A` |
| 86 | `N_SoftClipped_Bases_End_Diff` | `N_SoftClipped_Bases_End_B - N_SoftClipped_Bases_End_A` |
| 87 | `N_Matched_Junctions` | Query-space junction sets: size of the overlap, `\|A ∩ B\|` |
| 88 | `N_Unmatched_Junctions` | Query-space: symmetric difference, `N_Junctions_OnlyA + N_Junctions_OnlyB` |
| 89 | `N_Junctions_OnlyA` | Query-space junctions found only in A, `\|A \ B\|` |
| 90 | `N_Junctions_OnlyB` | Query-space junctions found only in B, `\|B \ A\|` |
| 91 | `Junction_Distance` | Sum of absolute differences between A's and B's query-space junction positions, paired in order (stops at the shorter list) |
| 92 | `Junc_Dist_V2` | `50 * \|JuncCount_A - JuncCount_B\|` |
| 93 | `Genomic_N_Matched_Junctions` | Genomic-coordinate junction sets: size of the overlap (always emitted; meaningful when both sides map to the same reference) |
| 94 | `Genomic_N_Unmatched_Junctions` | Genomic-coordinate: symmetric difference |
| 95 | `Genomic_N_Junctions_OnlyA` | Genomic-coordinate junctions found only in A |
| 96 | `Genomic_N_Junctions_OnlyB` | Genomic-coordinate junctions found only in B |
| 97 | `Junctions_OnlyA` | The actual query-space junction coordinates present only in A (Python tuple string) |
| 98 | `Junctions_OnlyB` | The actual query-space junction coordinates present only in B |
| 99 | `Genomic_Junctions_OnlyA` | The actual genomic-space junction coordinates present only in A (Python tuple-of-tuples string) |
| 100 | `Genomic_Junctions_OnlyB` | The actual genomic-space junction coordinates present only in B |

**Notes**
- **Junctions are compared as sets**, deduplicated per side. Query-space
  comparison uses the `junctions` columns; genomic-space comparison
  reattaches each side's `TargetChr` to its `genomic_junctions` pairs first,
  so junctions on different contigs can never falsely match.
- The `junctions` / `genomic_junctions` columns and the trailing `*_OnlyA/B`
  object lists use a Python tuple-of-tuples format; parse them in Python with
  `ast.literal_eval`. See [REFERENCE.md](../REFERENCE.md#junction-comparison)
  for the format.
- **Nulls mean undefined, not zero.** An unmapped side has a null `cs` and a
  null `seqid_Max`; a ratio over a zero denominator is null. A real `0` stays
  `0`, and `TargetChr` stays `*` for unmapped (it's the mapping indicator).
- Parquet uses the same 100 columns, same names, same order as the TSV — see
  [REFERENCE.md](../REFERENCE.md#tsv-and-parquet) for the TSV↔Parquet null/escaping mapping.
- Select columns by header name, not by position: positions can change when
  columns are added.
- Full classification logic (`query_identical`, `ref_same_position_same_aln`,
  etc.) derived from this table's columns is documented in
  [compare-summary.md](compare-summary.md) and
  [query-diff-reads.md](query-diff-reads.md).
