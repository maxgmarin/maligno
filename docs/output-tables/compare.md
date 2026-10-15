# `compare` table

`{prefix}.compare.parquet` (default) and/or `{prefix}.compare.tsv.gz` — one row per
read, 95 columns: both sides' representative-alignment stats plus a computed
A-vs-B comparison block. This is the **primary output** of `compare` (and of
`toolkit merge-readinfo`, which produces the identical table from
readinfo TSVs directly).

Columns are organized in groups, left to right:

| Group | Cols | What it holds |
|-------|:----:|---------------|
| **Join keys** | 1–2 | `query_name`, `query_len` |
| **Set labels** | 3–4 | `label_A`, `label_B`: the `--label-a` / `--label-b` values, repeated on every row |
| **Per-side data, A** | 5–38 | the representative alignment's stats for set A, each column suffixed `_A` |
| **Per-side data, B** | 39–72 | the same 34 columns for set B, suffixed `_B` |
| **Comparison metrics** | 73–93 | A-vs-B differences and ratios: `strand_match`, `aln_identity_diff`, coverage/length diffs, score diffs (`as_score_diff`, `ms_score_diff`, …), indel/soft-clip diffs, and query-space junction-set counts |
| **Non-overlap objects** | 94–95 | the query-space junctions that failed to overlap: `query_junctions_only_A/B` |

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

## Per-side data (columns 5–38 suffixed `_A`, 39–72 suffixed `_B`)

Same 34 columns, once per side, read out of each side's `readinfo` row by
name. Every value describes that side's **representative alignment** (see
[readinfo.md](readinfo.md) for how it is chosen), except the four `n_alns*`
counts, which cover all of the read's alignments on that side. `tp_tag` and
the raw-PAF columns such as `n_residue_matches` are readinfo/alninfo-only and
do **not** appear here. An unmapped side has `target_name_* == "*"`, an empty
`cs_*`, a `NaN` `aln_identity_*` and a `0.0` `query_aln_cov_*`.

| # (`_A` / `_B`) | Column | Description |
|---|---|---|
| 5 / 39 | `target_name` | Reference/target name; `*` if unmapped |
| 6 / 40 | `strand` | Strand; `*` if unmapped |
| 7 / 41 | `transcript_strand` | Transcript strand relative to the read: `+` sense, `-` antisense, `*` unknown (no `ts:A` / `XS:A` tag, or unmapped). See [REFERENCE.md](../REFERENCE.md#transcript-strand) |
| 8 / 42 | `target_start` | Target-coordinate span start (0-based half-open) |
| 9 / 43 | `target_end` | Target-coordinate span end |
| 10 / 44 | `query_start` | Query-coordinate span start (0-based half-open) |
| 11 / 45 | `query_end` | Query-coordinate span end |
| 12 / 46 | `mapq` | Mapping quality |
| 13 / 47 | `n_alns` | Number of aligned rows for this read on this side |
| 14 / 48 | `n_alns_tp_p` | Aligned rows whose PAF `tp:A` is `P` or `I` (SAM primary **+ supplementary**); `n_alns_tp_p − 1` = supplementary segments. See [readinfo.md](readinfo.md) |
| 15 / 49 | `n_alns_tp_s` | Aligned rows whose PAF `tp:A` is `S` or `i` (secondary) |
| 16 / 50 | `n_alns_tied_top_score` | Aligned rows tied with the representative on `ms_score`, `as_score` **and** `mapq`; `> 1` means input order picked the representative |
| 17 / 51 | `as_score` | Alignment score (`AS:i`) |
| 18 / 52 | `ms_score` | `ms:i` score (DP score of the max-scoring segment) |
| 19 / 53 | `aln_identity` | Alignment identity, `n_residue_matches / aln_block_len` |
| 20 / 54 | `query_aln_cov` | Fraction of the read aligned, `query_aln_len / query_len` |
| 21 / 55 | `query_aln_len` | Aligned query length, `query_end − query_start` |
| 22 / 56 | `n_query_junctions` | Number of query-space splice junctions |
| 23 / 57 | `n_intron_events` | Intron (`~`) op count |
| 24 / 58 | `n_intron_bases` | Total intron bases |
| 25 / 59 | `n_match_events` | Match-run count |
| 26 / 60 | `n_match_bases` | Matched bases |
| 27 / 61 | `n_substitution_events` | Substitution-event count |
| 28 / 62 | `n_substitution_bases` | Substituted bases |
| 29 / 63 | `n_insertion_events` | Insertion-event count |
| 30 / 64 | `n_insertion_bases` | Inserted bases |
| 31 / 65 | `n_deletion_events` | Deletion-event count |
| 32 / 66 | `n_deletion_bases` | Deleted bases |
| 33 / 67 | `n_softclip_events` | Soft-clip event count (`0`/`1`/`2`) |
| 34 / 68 | `n_softclip_bases_start` | Start soft-clip length |
| 35 / 69 | `n_softclip_bases_end` | End soft-clip length |
| 36 / 70 | `query_junctions` | Query-coordinate splice junctions (Python tuple string) |
| 37 / 71 | `genomic_junctions` | Splice junctions in absolute reference coordinates (Python tuple-of-tuples string; chrom is the sibling `target_name` column) |
| 38 / 72 | `cs` | Compact alignment string |

## Comparison block (columns 73–95)

A-vs-B metrics computed from the two per-side blocks above. **Sign
convention**: every `*_diff` is **B − A**, every `*_ratio` is **B / A**
(`NaN` when the A-side denominator is `0`).

Junctions are compared here in **query space** only, which doesn't depend on
the reference. Comparing genomic (reference-coordinate) junction sets between
A and B is done by
[`toolkit genomic-junction-compare`](../REFERENCE.md#toolkit-genomic-junction-compare).

| # | Column | Description |
|---|---|---|
| 73 | `strand_match` | `strand_A == strand_B` |
| 74 | `aln_identity_diff` | `aln_identity_B - aln_identity_A` |
| 75 | `query_aln_cov_diff` | `query_aln_cov_B - query_aln_cov_A` |
| 76 | `query_aln_len_diff` | `query_aln_len_B - query_aln_len_A` |
| 77 | `as_score_diff` | `as_score_B - as_score_A` |
| 78 | `ms_score_diff` | `ms_score_B - ms_score_A` |
| 79 | `as_score_ratio` | `as_score_B / as_score_A` |
| 80 | `ms_score_ratio` | `ms_score_B / ms_score_A` |
| 81 | `n_substitution_bases_diff` | `n_substitution_bases_B - n_substitution_bases_A` |
| 82 | `n_substitution_bases_ratio` | `n_substitution_bases_B / n_substitution_bases_A` |
| 83 | `n_insertion_bases_diff` | `n_insertion_bases_B - n_insertion_bases_A` |
| 84 | `n_insertion_bases_ratio` | `n_insertion_bases_B / n_insertion_bases_A` |
| 85 | `n_deletion_bases_diff` | `n_deletion_bases_B - n_deletion_bases_A` |
| 86 | `n_deletion_bases_ratio` | `n_deletion_bases_B / n_deletion_bases_A` |
| 87 | `n_softclip_bases_start_diff` | `n_softclip_bases_start_B - n_softclip_bases_start_A` |
| 88 | `n_softclip_bases_end_diff` | `n_softclip_bases_end_B - n_softclip_bases_end_A` |
| 89 | `n_query_junctions_shared` | Query-space junction sets: size of the overlap, `\|A ∩ B\|` |
| 90 | `n_query_junctions_not_shared` | Query-space: symmetric difference, `n_query_junctions_only_A + n_query_junctions_only_B` |
| 91 | `n_query_junctions_only_A` | Query-space junctions found only in A, `\|A \ B\|` |
| 92 | `n_query_junctions_only_B` | Query-space junctions found only in B, `\|B \ A\|` |
| 93 | `query_junction_pos_distance` | Sum of absolute differences between A's and B's query-space junction positions, paired in order (stops at the shorter list) |
| 94 | `query_junctions_only_A` | The actual query-space junction coordinates present only in A (Python tuple string) |
| 95 | `query_junctions_only_B` | The actual query-space junction coordinates present only in B |

**Notes**
- **Junctions are compared as sets** of query-space positions (the
  `query_junctions` columns), deduplicated per side.
- The `query_junctions` / `genomic_junctions` columns and the trailing
  `query_junctions_only_A/B` object lists use a Python tuple (or tuple-of-tuples)
  format; parse them in
  Python with `ast.literal_eval`. See
  [REFERENCE.md](../REFERENCE.md#junction-comparison) for the format.
- **Nulls mean undefined, not zero.** An unmapped side has a null `cs` and a
  null `aln_identity`; a ratio over a zero denominator is null. A real `0` stays
  `0`, and `target_name` stays `*` for unmapped (it's the mapping indicator).
- Parquet uses the same 95 columns, same names, same order as the TSV — see
  [REFERENCE.md](../REFERENCE.md#tsv-and-parquet) for the TSV↔Parquet null/escaping mapping.
- Select columns by header name, not by position: positions can change when
  columns are added.
- Full classification logic (`query_identical`, `ref_same_position_same_aln`,
  etc.) derived from this table's columns is documented in
  [compare-summary.md](compare-summary.md) and
  [query-diff-reads.md](query-diff-reads.md).
