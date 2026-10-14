# `compare` table

`{prefix}.compare.parquet` (default) and/or `{prefix}.compare.tsv.gz` — one row per
read, 99 columns: both sides' representative-alignment stats plus a computed
A-vs-B comparison block. This is the **primary output** of `compare` (and of
`toolkit merge-readinfo`, which produces the identical table from
readinfo TSVs directly).

Columns are organized in groups, left to right:

| Group | Cols | What it holds |
|-------|:----:|---------------|
| **Join keys** | 1–2 | `query_name`, `query_len` |
| **Set labels** | 3–4 | `label_A`, `label_B`: the `--label-a` / `--label-b` values, repeated on every row |
| **Per-side data, A** | 5–37 | the representative alignment's stats for set A, each column suffixed `_A` |
| **Per-side data, B** | 38–70 | the same 33 columns for set B, suffixed `_B` |
| **Comparison metrics** | 71–95 | A-vs-B differences and ratios: `strand_match`, `aln_identity_diff`, coverage/length diffs, score diffs (`as_score_diff`, `ms_score_diff`, …), indel/soft-clip diffs, and junction-set counts in both query and genomic space |
| **Non-overlap objects** | 96–99 | the junctions that failed to overlap: `query_junctions_only_A/B` and `genomic_junctions_only_A/B` |

Within each per-side block the columns are grouped by topic: locus and span,
alignment selection and score (including the per-type alignment counts
`n_alns_tp_p` / `n_alns_tp_s`), identity and coverage, junction counts,
cs-derived event counts, then the three long strings (`query_junctions`,
`genomic_junctions`, `cs`) last. Print the exact layout of any table with
`gzip -dc … | head -1 | tr '\t' '\n' | nl`.

**Naming.** Every column is lowercase snake_case, and counts start with `n_`.
A per-side column has the same name as the readinfo column it comes from, plus
the side suffix. A comparison column is named after the per-side column it is
computed from: `{col}_diff` is `{col}_B − {col}_A` and `{col}_ratio` is
`{col}_B / {col}_A`.

Per-side columns always use the **fixed** `_A` / `_B` suffixes, never the
dataset label, so column names are identical for every comparison and stay
unambiguous even when a label contains an underscore. Which dataset each side
*is* comes from the `label_A` / `label_B` columns, carried on every row so any
subset of rows stays self-describing.

## Join keys & labels (columns 1–4)

| # | Column | Description |
|---|---|---|
| 1 | `query_name` | Read ID — join key |
| 2 | `query_len` | Query length (side A's) |
| 3 | `label_A` | The `--label-a` value, repeated on every row |
| 4 | `label_B` | The `--label-b` value, repeated on every row |

## Per-side data (columns 5–37 suffixed `_A`, 38–70 suffixed `_B`)

Same 33 columns, once per side, read out of each side's `readinfo` row by
name. Every value describes that side's **representative alignment** (see
[readinfo.md](readinfo.md) for how it is chosen), except the four `n_alns*`
counts, which cover all of the read's alignments on that side. `tp_tag` and
the raw-PAF columns such as `n_residue_matches` are readinfo/alninfo-only and
do **not** appear here. An unmapped side has `target_name_* == "*"`, an empty
`cs_*`, a `NaN` `aln_identity_*` and a `0.0` `query_aln_cov_*`.

| # (`_A` / `_B`) | Column | Description |
|---|---|---|
| 5 / 38 | `target_name` | Reference/target name; `*` if unmapped |
| 6 / 39 | `strand` | Strand; `*` if unmapped |
| 7 / 40 | `target_start` | Target-coordinate span start (0-based half-open) |
| 8 / 41 | `target_end` | Target-coordinate span end |
| 9 / 42 | `query_start` | Query-coordinate span start (0-based half-open) |
| 10 / 43 | `query_end` | Query-coordinate span end |
| 11 / 44 | `mapq` | Mapping quality |
| 12 / 45 | `n_alns` | Number of aligned rows for this read on this side |
| 13 / 46 | `n_alns_tp_p` | Aligned rows whose PAF `tp:A` is `P` or `I` (SAM primary **+ supplementary**); `n_alns_tp_p − 1` = supplementary segments. See [readinfo.md](readinfo.md) |
| 14 / 47 | `n_alns_tp_s` | Aligned rows whose PAF `tp:A` is `S` or `i` (secondary) |
| 15 / 48 | `n_alns_tied_top_score` | Aligned rows tied with the representative on `ms_score`, `as_score` **and** `mapq`; `> 1` means input order picked the representative |
| 16 / 49 | `as_score` | Alignment score (`AS:i`) |
| 17 / 50 | `ms_score` | `ms:i` score (DP score of the max-scoring segment) |
| 18 / 51 | `aln_identity` | Alignment identity, `n_residue_matches / aln_block_len` |
| 19 / 52 | `query_aln_cov` | Fraction of the read aligned, `query_aln_len / query_len` |
| 20 / 53 | `query_aln_len` | Aligned query length, `query_end − query_start` |
| 21 / 54 | `n_query_junctions` | Number of query-space splice junctions |
| 22 / 55 | `n_intron_events` | Intron (`~`) op count |
| 23 / 56 | `n_intron_bases` | Total intron bases |
| 24 / 57 | `n_match_events` | Match-run count |
| 25 / 58 | `n_match_bases` | Matched bases |
| 26 / 59 | `n_substitution_events` | Substitution-event count |
| 27 / 60 | `n_substitution_bases` | Substituted bases |
| 28 / 61 | `n_insertion_events` | Insertion-event count |
| 29 / 62 | `n_insertion_bases` | Inserted bases |
| 30 / 63 | `n_deletion_events` | Deletion-event count |
| 31 / 64 | `n_deletion_bases` | Deleted bases |
| 32 / 65 | `n_softclip_events` | Soft-clip event count (`0`/`1`/`2`) |
| 33 / 66 | `n_softclip_bases_start` | Start soft-clip length |
| 34 / 67 | `n_softclip_bases_end` | End soft-clip length |
| 35 / 68 | `query_junctions` | Query-coordinate splice junctions (Python tuple string) |
| 36 / 69 | `genomic_junctions` | Splice junctions in absolute reference coordinates (Python tuple-of-tuples string; chrom is the sibling `target_name` column) |
| 37 / 70 | `cs` | Compact alignment string |

## Comparison block (columns 71–99)

A-vs-B metrics computed from the two per-side blocks above. **Sign
convention**: every `*_diff` is **B − A**, every `*_ratio` is **B / A**
(`NaN` when the A-side denominator is `0`).

| # | Column | Description |
|---|---|---|
| 71 | `strand_match` | `strand_A == strand_B` |
| 72 | `aln_identity_diff` | `aln_identity_B - aln_identity_A` |
| 73 | `query_aln_cov_diff` | `query_aln_cov_B - query_aln_cov_A` |
| 74 | `query_aln_len_diff` | `query_aln_len_B - query_aln_len_A` |
| 75 | `as_score_diff` | `as_score_B - as_score_A` |
| 76 | `ms_score_diff` | `ms_score_B - ms_score_A` |
| 77 | `as_score_ratio` | `as_score_B / as_score_A` |
| 78 | `ms_score_ratio` | `ms_score_B / ms_score_A` |
| 79 | `n_substitution_bases_diff` | `n_substitution_bases_B - n_substitution_bases_A` |
| 80 | `n_substitution_bases_ratio` | `n_substitution_bases_B / n_substitution_bases_A` |
| 81 | `n_insertion_bases_diff` | `n_insertion_bases_B - n_insertion_bases_A` |
| 82 | `n_insertion_bases_ratio` | `n_insertion_bases_B / n_insertion_bases_A` |
| 83 | `n_deletion_bases_diff` | `n_deletion_bases_B - n_deletion_bases_A` |
| 84 | `n_deletion_bases_ratio` | `n_deletion_bases_B / n_deletion_bases_A` |
| 85 | `n_softclip_bases_start_diff` | `n_softclip_bases_start_B - n_softclip_bases_start_A` |
| 86 | `n_softclip_bases_end_diff` | `n_softclip_bases_end_B - n_softclip_bases_end_A` |
| 87 | `n_query_junctions_shared` | Query-space junction sets: size of the overlap, `\|A ∩ B\|` |
| 88 | `n_query_junctions_not_shared` | Query-space: symmetric difference, `n_query_junctions_only_A + n_query_junctions_only_B` |
| 89 | `n_query_junctions_only_A` | Query-space junctions found only in A, `\|A \ B\|` |
| 90 | `n_query_junctions_only_B` | Query-space junctions found only in B, `\|B \ A\|` |
| 91 | `query_junction_pos_distance` | Sum of absolute differences between A's and B's query-space junction positions, paired in order (stops at the shorter list) |
| 92 | `n_genomic_junctions_shared` | Genomic-coordinate junction sets: size of the overlap (always emitted; meaningful when both sides map to the same reference) |
| 93 | `n_genomic_junctions_not_shared` | Genomic-coordinate: symmetric difference |
| 94 | `n_genomic_junctions_only_A` | Genomic-coordinate junctions found only in A |
| 95 | `n_genomic_junctions_only_B` | Genomic-coordinate junctions found only in B |
| 96 | `query_junctions_only_A` | The actual query-space junction coordinates present only in A (Python tuple string) |
| 97 | `query_junctions_only_B` | The actual query-space junction coordinates present only in B |
| 98 | `genomic_junctions_only_A` | The actual genomic-space junction coordinates present only in A (Python tuple-of-tuples string) |
| 99 | `genomic_junctions_only_B` | The actual genomic-space junction coordinates present only in B |

**Notes**
- **Junctions are compared as sets**, deduplicated per side. Query-space
  comparison uses the `query_junctions` columns; genomic-space comparison
  reattaches each side's `target_name` to its `genomic_junctions` pairs first,
  so junctions on different contigs can never falsely match.
- The `query_junctions` / `genomic_junctions` columns and the trailing
  `*_only_A/B` object lists use a Python tuple-of-tuples format; parse them in
  Python with `ast.literal_eval`. See
  [REFERENCE.md](../REFERENCE.md#junction-comparison) for the format.
- **Nulls mean undefined, not zero.** An unmapped side has a null `cs` and a
  null `aln_identity`; a ratio over a zero denominator is null. A real `0` stays
  `0`, and `target_name` stays `*` for unmapped (it's the mapping indicator).
- Parquet uses the same 99 columns, same names, same order as the TSV — see
  [REFERENCE.md](../REFERENCE.md#tsv-and-parquet) for the TSV↔Parquet null/escaping mapping.
- Select columns by header name, not by position: positions can change when
  columns are added.
- Full classification logic (`query_identical`, `ref_same_position_same_aln`,
  etc.) derived from this table's columns is documented in
  [compare-summary.md](compare-summary.md) and
  [query-diff-reads.md](query-diff-reads.md).
