# `alninfo` table

`{prefix}.{label}.alninfo.tsv.gz` — **per-alignment** table: one row per PAF
alignment record (every alignment for every read, not just the best one).

Produced by `toolkit paf2tables --alninfo` and, opt-in, by
`compare --emit-alninfo` (once per side, i.e. `{prefix}.{label_a}.alninfo.tsv.gz`
and `{prefix}.{label_b}.alninfo.tsv.gz`). Off by default in `compare` — see
[REFERENCE.md](../REFERENCE.md#outputs) for why. Pure streaming output: rows are
written in PAF input order, one PAF line in, one row out (malformed lines are
skipped with a stderr warning). An unmapped PAF record (`target_name == "*"`)
still produces a full row, with alignment-derived fields zeroed/`NaN`.

The per-alignment columns use the same names as in [readinfo.md](readinfo.md),
which describes each read by its representative alignment.

37 columns.

| # | Column | Description |
|---|---|---|
| 1 | `query_name` | Read/query sequence ID (PAF col 1) |
| 2 | `query_len` | Query sequence length (PAF col 2) |
| 3 | `query_start` | Alignment's query start, 0-based half-open (PAF col 3) |
| 4 | `query_end` | Alignment's query end, 0-based half-open (PAF col 4) |
| 5 | `strand` | `+`/`-` relative to the reference; `*` if unmapped (PAF col 5) |
| 6 | `transcript_strand` | Transcript strand relative to the read: `+` sense, `-` antisense, `*` unknown or unmapped. From minimap2's `ts:A` or STAR's `XS:A`; see [REFERENCE.md](../REFERENCE.md#transcript-strand) |
| 7 | `target_name` | Reference/target sequence name (chromosome); `*` if unmapped (PAF col 6) |
| 8 | `target_len` | Target sequence length (PAF col 7) |
| 9 | `target_start` | Alignment's target start, 0-based half-open (PAF col 8) |
| 10 | `target_end` | Alignment's target end, 0-based half-open (PAF col 9) |
| 11 | `n_residue_matches` | Number of matching bases in the alignment (PAF col 10) |
| 12 | `aln_block_len` | Alignment block length (PAF col 11) |
| 13 | `mapq` | Mapping quality (PAF col 12) |
| 14 | `ms_score` | minimap2's `ms:i` tag: DP score of the max-scoring segment in the alignment; `0` if absent |
| 15 | `as_score` | Alignment score (`AS:i` tag); `0` if absent |
| 16 | `cs` | Compact alignment string (`cs:Z` tag); empty for unmapped |
| 17 | `n_match_events` | Number of contiguous match runs walked out of `cs` |
| 18 | `n_match_bases` | Total matched bases |
| 19 | `n_substitution_events` | Number of substitution events |
| 20 | `n_substitution_bases` | Total substituted bases |
| 21 | `n_insertion_events` | Number of insertion events |
| 22 | `n_insertion_bases` | Total inserted bases |
| 23 | `n_deletion_events` | Number of deletion events |
| 24 | `n_deletion_bases` | Total deleted bases |
| 25 | `n_intron_events` | Number of intron/splice (`~`) ops in `cs` |
| 26 | `n_intron_bases` | Total intron bases |
| 27 | `n_softclip_bases_start` | Unaligned bases clipped at the read's start (strand-aware); `0` for unmapped |
| 28 | `n_softclip_bases_end` | Unaligned bases clipped at the read's end; `0` for unmapped |
| 29 | `n_softclip_events` | `0`, `1`, or `2` — how many ends are clipped |
| 30 | `query_junctions` | Query-coordinate splice-junction positions, plus-strand-adjusted, sorted ascending. Python tuple string, e.g. `(108, 359)`; `()` if none |
| 31 | `n_query_junctions` | Number of entries in `query_junctions` |
| 32 | `target_start_1based` | `target_start + 1` — 1-based convenience column |
| 33 | `aln_identity` | Alignment identity, `n_residue_matches / aln_block_len`; `NaN` if unmapped |
| 34 | `query_aln_len` | `query_end - query_start` |
| 35 | `query_aln_cov` | `query_aln_len / query_len`; `NaN` if `query_len` is `0` |
| 36 | `genomic_junctions` | Splice junctions in absolute reference coordinates, 0-based half-open. Python tuple-of-tuples string, e.g. `((100, 250), (400, 800))`; `()` if none. Chromosome is **not** embedded — it's the sibling `target_name` column |
| 37 | `tp_tag` | PAF `tp:A` alignment-type tag: `P` primary, `S` secondary, `I` inversion of the primary. `*` when absent or unmapped |

**Notes**
- `query_junctions` / `genomic_junctions` use a Python tuple-of-tuples format —
  parse with `ast.literal_eval` in Python.
- `aln_identity` is `NaN`, not `0`, when undefined (unmapped);
  `query_aln_cov` is `NaN` only when `query_len` is `0`.
