# `compare` table

`{prefix}.compare.tsv.gz` and/or `{prefix}.compare.parquet` — one row per
read, 96 columns: both sides' representative-alignment stats plus a computed
A-vs-B comparison block. This is the **primary output** of `compare` (and of
`compare-toolkit merge-readinfo`, which produces the identical table from
readinfo TSVs directly).

Layout: 4 join/label columns, then the 31 per-side columns once suffixed
`_A` and once suffixed `_B`, then 30 comparison columns. Per-side suffixes
are always the fixed `_A`/`_B` — never the dataset label — so column names
are identical across every comparison run; which dataset each side *is*
comes from the `Label_A`/`Label_B` data columns instead.

## Join keys & labels (columns 1–4)

| # | Column | Description |
|---|---|---|
| 1 | `Read_Name` | Read ID — join key |
| 2 | `Read_Len` | Query length — join key |
| 3 | `Label_A` | The `--label-a` value, repeated on every row |
| 4 | `Label_B` | The `--label-b` value, repeated on every row |

## Per-side data (columns 5–35 suffixed `_A`, 36–66 suffixed `_B`)

Same 31 columns, once per side — read out of each side's `readinfo` row by
name (see [readinfo.md](readinfo.md) for how each is computed; `tp_tag` and
`Num_Residue_Matches`-style raw-PAF columns are readinfo/alninfo-only and do
**not** appear here). An unmapped side has `TargetChr_* == "*"`, an empty
`cs_*`, and `NaN` for `seqid_Max_*`/`Query_Aln_Cov_Max_*`.

| # (`_A` / `_B`) | Column | Description |
|---|---|---|
| 5 / 36 | `TargetChr` | Representative alignment's reference/target name; `*` if unmapped |
| 6 / 37 | `Strand` | Representative alignment's strand; `*` if unmapped |
| 7 / 38 | `Target_Start` | Target-coordinate span start (0-based half-open) |
| 8 / 39 | `Target_End` | Target-coordinate span end |
| 9 / 40 | `Query_Start` | Query-coordinate span start (0-based half-open) |
| 10 / 41 | `Query_End` | Query-coordinate span end |
| 11 / 42 | `MQ_Best` | Representative alignment's mapping quality |
| 12 / 43 | `Num_Aln` | Total alignment count for this read on this side |
| 13 / 44 | `Num_Aln_MaxScore` | Alignments tied at the winning `(ms, AS)` score |
| 14 / 45 | `AS_Max` | Highest `AS` across this side's alignments |
| 15 / 46 | `ms_Max` | Highest `ms` across this side's alignments |
| 16 / 47 | `seqid_Max` | Highest sequence identity across this side's alignments |
| 17 / 48 | `Query_Aln_Cov_Max` | Highest query-alignment coverage |
| 18 / 49 | `Query_Aln_Len_Max` | Highest aligned query length |
| 19 / 50 | `JuncCount` | Number of query-space splice junctions in the representative alignment |
| 20 / 51 | `N_Splice_Junction_Events` | Representative alignment's intron/splice-op count |
| 21 / 52 | `N_Splice_Junction_Bases` | Representative alignment's total intron bases |
| 22 / 53 | `N_Match_Events` | Representative alignment's match-run count |
| 23 / 54 | `N_Match_Bases` | Representative alignment's matched bases |
| 24 / 55 | `N_Substitution_Events` | Representative alignment's substitution-event count |
| 25 / 56 | `N_Substitution_Bases` | Representative alignment's substituted bases |
| 26 / 57 | `N_Insertion_Events` | Representative alignment's insertion-event count |
| 27 / 58 | `N_Insertion_Bases` | Representative alignment's inserted bases |
| 28 / 59 | `N_Deletion_Events` | Representative alignment's deletion-event count |
| 29 / 60 | `N_Deletion_Bases` | Representative alignment's deleted bases |
| 30 / 61 | `N_SoftClipped_Events` | Representative alignment's soft-clip event count (`0`/`1`/`2`) |
| 31 / 62 | `N_SoftClipped_Bases_Start` | Representative alignment's start soft-clip length |
| 32 / 63 | `N_SoftClipped_Bases_End` | Representative alignment's end soft-clip length |
| 33 / 64 | `junctions` | Representative alignment's query-coordinate splice junctions (Python tuple string) |
| 34 / 65 | `genomic_junctions` | Representative alignment's splice junctions in absolute reference coordinates (Python tuple-of-tuples string; chrom is the sibling `TargetChr` column) |
| 35 / 66 | `cs` | Representative alignment's compact alignment string |

## Comparison block (columns 67–96)

A-vs-B metrics computed from the two per-side blocks above. **Sign
convention**: every `*_Diff` is **B − A**, every `*_Ratio` is **B / A**
(`NaN` when the A-side denominator is `0`) — except `Junc_Dist_V2`, whose
inner subtraction is A − B but is then made symmetric by `unsigned_abs()`.

| # | Column | Description |
|---|---|---|
| 67 | `Strand_Match` | `Strand_A == Strand_B` |
| 68 | `seqid_Diff` | `seqid_Max_B - seqid_Max_A` |
| 69 | `QueryAlnCov_Diff` | `Query_Aln_Cov_Max_B - Query_Aln_Cov_Max_A` |
| 70 | `QueryAlnLen_Diff` | `Query_Aln_Len_Max_B - Query_Aln_Len_Max_A` |
| 71 | `AS_Diff` | `AS_Max_B - AS_Max_A` |
| 72 | `ms_Diff` | `ms_Max_B - ms_Max_A` |
| 73 | `AS_Ratio` | `AS_Max_B / AS_Max_A` |
| 74 | `ms_Ratio` | `ms_Max_B / ms_Max_A` |
| 75 | `N_Substitution_Bases_Diff` | `N_Substitution_Bases_B - N_Substitution_Bases_A` |
| 76 | `N_Substitution_Bases_Ratio` | `N_Substitution_Bases_B / N_Substitution_Bases_A` |
| 77 | `N_Insertion_Bases_Diff` | `N_Insertion_Bases_B - N_Insertion_Bases_A` |
| 78 | `N_Insertion_Bases_Ratio` | `N_Insertion_Bases_B / N_Insertion_Bases_A` |
| 79 | `N_Deletion_Bases_Diff` | `N_Deletion_Bases_B - N_Deletion_Bases_A` |
| 80 | `N_Deletion_Bases_Ratio` | `N_Deletion_Bases_B / N_Deletion_Bases_A` |
| 81 | `N_SoftClipped_Bases_Start_Diff` | `N_SoftClipped_Bases_Start_B - N_SoftClipped_Bases_Start_A` |
| 82 | `N_SoftClipped_Bases_End_Diff` | `N_SoftClipped_Bases_End_B - N_SoftClipped_Bases_End_A` |
| 83 | `N_Matched_Junctions` | Query-space junction sets: size of the overlap, `\|A ∩ B\|` |
| 84 | `N_Unmatched_Junctions` | Query-space: symmetric difference, `N_Junctions_OnlyA + N_Junctions_OnlyB` |
| 85 | `N_Junctions_OnlyA` | Query-space junctions found only in A, `\|A \ B\|` |
| 86 | `N_Junctions_OnlyB` | Query-space junctions found only in B, `\|B \ A\|` |
| 87 | `Junction_Distance` | Legacy positional junction-distance metric |
| 88 | `Junc_Dist_V2` | Legacy metric: `50 * \|JuncCount_A - JuncCount_B\|` |
| 89 | `Genomic_N_Matched_Junctions` | Genomic-coordinate junction sets: size of the overlap (always emitted; meaningful when both sides map to the same reference) |
| 90 | `Genomic_N_Unmatched_Junctions` | Genomic-coordinate: symmetric difference |
| 91 | `Genomic_N_Junctions_OnlyA` | Genomic-coordinate junctions found only in A |
| 92 | `Genomic_N_Junctions_OnlyB` | Genomic-coordinate junctions found only in B |
| 93 | `Junctions_OnlyA` | The actual query-space junction coordinates present only in A (Python tuple string) |
| 94 | `Junctions_OnlyB` | The actual query-space junction coordinates present only in B |
| 95 | `Genomic_Junctions_OnlyA` | The actual genomic-space junction coordinates present only in A (Python tuple-of-tuples string) |
| 96 | `Genomic_Junctions_OnlyB` | The actual genomic-space junction coordinates present only in B |

**Notes**
- **Junctions are compared as sets**, deduplicated per side. Query-space
  comparison uses the `junctions` columns; genomic-space comparison
  reattaches each side's `TargetChr` to its `genomic_junctions` pairs first,
  so junctions on different contigs can never falsely match.
- **Nulls mean undefined, not zero.** An unmapped side has a null `cs` and a
  null `seqid_Max`; a ratio over a zero denominator is null. A real `0` stays
  `0`, and `TargetChr` stays `*` for unmapped (it's the mapping indicator).
- Parquet uses the same 96 columns, same names, same order as the TSV — see
  [REFERENCE.md](../REFERENCE.md) for the TSV↔Parquet null/escaping mapping.
- Full classification logic (`query_identical`, `ref_same_position_same_aln`,
  etc.) derived from this table's columns is documented in
  [compare-summary.md](compare-summary.md) and
  [query-diff-reads.md](query-diff-reads.md).
