# `compare` table

`{prefix}.compare.parquet` (default) and/or `{prefix}.compare.tsv.gz` — one row per
read, 93 columns: both sides' representative-alignment stats plus a computed
A-vs-B comparison block. This is the **primary output** of `compare` (and of
`toolkit compare-readinfo`, which produces the identical table from
readinfo TSVs directly).

Columns are organized in groups, left to right:

| Group | Cols | What it holds |
|-------|:----:|---------------|
| **Join keys** | 1–2 | `query_name`, `query_len` |
| **Set labels** | 3–4 | `label_A`, `label_B`: the `--label-a` / `--label-b` values, repeated on every row |
| **Per-side data, A** | 5–39 | the representative alignment's stats for set A, each column suffixed `_A` |
| **Per-side data, B** | 40–74 | the same 35 columns for set B, suffixed `_B` |
| **Comparison metrics** | 75–91 | A-vs-B differences and ratios: `strand_match`, `aln_identity_diff`, coverage/length diffs, score diffs and ratios (`as_score_diff`, `ms_score_ratio`, …), indel/soft-clip diffs, and query-space junction-set counts |
| **Non-overlap objects** | 92–93 | the query-space junctions that failed to overlap: `query_junctions_only_A/B` |

Within each per-side block the columns are grouped by topic: locus and span,
alignment selection and score (the representative's `tp_tag`, then the
per-type alignment counts `n_alns_primsupp` / `n_alns_secondary`), identity and coverage, `cs` and the
two strand columns, cs-derived event counts, then introns and junctions
(`query_junctions`, `genomic_junctions` last). The
[readinfo table](readinfo.md) is `query_name`, `query_len` and exactly these 35
columns, in the same order. Print the exact layout of any table with
`gzip -dc … | head -1 | tr '\t' '\n' | nl`.

**Naming.** Every column is lowercase snake_case, and counts start with `n_`.
A per-side column has the same name as the readinfo column it comes from, plus
the side suffix. A comparison column is named after the per-side column it is
computed from: `{col}_diff` is `{col}_B − {col}_A` and `{col}_ratio` is
`{col}_B / {col}_A` (ratios only for the two alignment scores).

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

## Per-side data (columns 5–39 suffixed `_A`, 40–74 suffixed `_B`)

Same 35 columns, once per side, read out of each side's `readinfo` row by
name. Every value describes that side's **representative alignment** (see
[readinfo.md](readinfo.md) for how it is chosen), except the four `n_alns*`
counts, which cover all of the read's alignments on that side. The raw-PAF
columns such as `target_len` and `n_residue_matches` are alninfo-only and do
**not** appear here. An unmapped side has `target_name_* == "*"`, an empty
`cs_*`, a `NaN` `aln_identity_*` and a `0.0` `query_aln_cov_*`.

| # (`_A` / `_B`) | Column | Description |
|---|---|---|
| 5 / 40 | `target_name` | Reference/target name; `*` if unmapped |
| 6 / 41 | `strand` | Strand; `*` if unmapped |
| 7 / 42 | `target_start` | Target-coordinate span start (0-based half-open) |
| 8 / 43 | `target_end` | Target-coordinate span end |
| 9 / 44 | `query_start` | Query-coordinate span start (0-based half-open) |
| 10 / 45 | `query_end` | Query-coordinate span end |
| 11 / 46 | `mapq` | Mapping quality |
| 12 / 47 | `tp_tag` | The representative's PAF `tp:A` type (`P`, `S`, …; `*` if unmapped). `S` means a secondary outscored every primary (or the input has no primary for the read) |
| 13 / 48 | `n_alns` | Number of aligned rows for this read on this side |
| 14 / 49 | `n_alns_primsupp` | Aligned rows whose PAF `tp:A` is `P` or `I` (SAM primary **+ supplementary**); `n_alns_primsupp − 1` = supplementary segments. See [readinfo.md](readinfo.md) |
| 15 / 50 | `n_alns_secondary` | Aligned rows whose PAF `tp:A` is `S` or `i` (secondary) |
| 16 / 51 | `n_alns_tied_top_score` | Aligned rows tied with the representative on `ms_score`, `as_score` **and** `mapq` (the `tp_tag` tie-break doesn't count); `> 1` means the representative was picked among equals, primary first |
| 17 / 52 | `as_score` | Alignment score (`AS:i`) |
| 18 / 53 | `ms_score` | `ms:i` score (DP score of the max-scoring segment) |
| 19 / 54 | `aln_identity` | Alignment identity, `n_residue_matches / aln_block_len` |
| 20 / 55 | `query_aln_cov` | Fraction of the read aligned, `query_aln_len / query_len` |
| 21 / 56 | `query_aln_len` | Aligned query length, `query_end − query_start` |
| 22 / 57 | `cs` | Compact alignment string |
| 23 / 58 | `query_transcript_strand` | Transcript strand **relative to the read**: `+` sense, `-` antisense, `*` unknown (no `ts:A` / `XS:A` tag, or unmapped). See [REFERENCE.md](../REFERENCE.md#transcript-strand) |
| 24 / 59 | `genomic_intron_strand` | Transcript (and intron) strand **on the reference**: `strand` when `query_transcript_strand` is `+`, the opposite when `-`, `*` if unknown |
| 25 / 60 | `n_match_events` | Match-run count |
| 26 / 61 | `n_match_bases` | Matched bases |
| 27 / 62 | `n_substitution_events` | Substitution-event count |
| 28 / 63 | `n_substitution_bases` | Substituted bases |
| 29 / 64 | `n_insertion_events` | Insertion-event count |
| 30 / 65 | `n_insertion_bases` | Inserted bases |
| 31 / 66 | `n_deletion_events` | Deletion-event count |
| 32 / 67 | `n_deletion_bases` | Deleted bases |
| 33 / 68 | `n_softclip_events` | Soft-clip event count (`0`/`1`/`2`) |
| 34 / 69 | `n_softclip_bases_start` | Start soft-clip length |
| 35 / 70 | `n_softclip_bases_end` | End soft-clip length |
| 36 / 71 | `n_intron_events` | Intron (`~`) op count, i.e. the number of junctions |
| 37 / 72 | `n_intron_bases` | Total intron bases |
| 38 / 73 | `query_junctions` | Query-coordinate splice junctions (Python tuple string) |
| 39 / 74 | `genomic_junctions` | Splice junctions in absolute reference coordinates (Python tuple-of-tuples string; chrom is the sibling `target_name` column) |

## Comparison block (columns 75–93)

A-vs-B metrics computed from the two per-side blocks above. **Sign
convention**: every `*_diff` is **B − A**, every `*_ratio` is **B / A**
(`NaN` when the A-side denominator is `0`).

Junctions are compared here in **query space** only, which doesn't depend on
the reference. Comparing genomic (reference-coordinate) junction sets between
A and B is done by
[`toolkit genomic-junction-compare`](../REFERENCE.md#toolkit-genomic-junction-compare).

| # | Column | Description |
|---|---|---|
| 75 | `strand_match` | `strand_A == strand_B` |
| 76 | `aln_identity_diff` | `aln_identity_B - aln_identity_A` |
| 77 | `query_aln_cov_diff` | `query_aln_cov_B - query_aln_cov_A` |
| 78 | `query_aln_len_diff` | `query_aln_len_B - query_aln_len_A` |
| 79 | `as_score_diff` | `as_score_B - as_score_A` |
| 80 | `ms_score_diff` | `ms_score_B - ms_score_A` |
| 81 | `as_score_ratio` | `as_score_B / as_score_A` |
| 82 | `ms_score_ratio` | `ms_score_B / ms_score_A` |
| 83 | `n_substitution_bases_diff` | `n_substitution_bases_B - n_substitution_bases_A` |
| 84 | `n_insertion_bases_diff` | `n_insertion_bases_B - n_insertion_bases_A` |
| 85 | `n_deletion_bases_diff` | `n_deletion_bases_B - n_deletion_bases_A` |
| 86 | `n_softclip_bases_start_diff` | `n_softclip_bases_start_B - n_softclip_bases_start_A` |
| 87 | `n_softclip_bases_end_diff` | `n_softclip_bases_end_B - n_softclip_bases_end_A` |
| 88 | `n_query_junctions_shared` | Query-space junction sets: size of the overlap, `\|A ∩ B\|` |
| 89 | `n_query_junctions_not_shared` | Query-space: symmetric difference, `n_query_junctions_only_A + n_query_junctions_only_B` |
| 90 | `n_query_junctions_only_A` | Query-space junctions found only in A, `\|A \ B\|` |
| 91 | `n_query_junctions_only_B` | Query-space junctions found only in B, `\|B \ A\|` |
| 92 | `query_junctions_only_A` | The actual query-space junction coordinates present only in A (Python tuple string) |
| 93 | `query_junctions_only_B` | The actual query-space junction coordinates present only in B |

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
- Parquet uses the same 93 columns, same names, same order as the TSV — see
  [REFERENCE.md](../REFERENCE.md#tsv-and-parquet) for the TSV↔Parquet null/escaping mapping.
- Select columns by header name, not by position: positions can change when
  columns are added.
- Full classification logic (`query_identical`, `ref_same_position_same_aln`,
  etc.) derived from this table's columns is documented in
  [compare-summary.md](compare-summary.md) and
  [aln-diff-reads.md](aln-diff-reads.md).
