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

37 columns: the 12 standard PAF fields in PAF order, then the remaining
columns in the same order as the per-side block of the
[comparison table](compare.md).

| # | Column | Description |
|---|---|---|
| 1 | `query_name` | Read/query sequence ID (PAF col 1) |
| 2 | `query_len` | Query sequence length (PAF col 2) |
| 3 | `query_start` | Alignment's query start, 0-based half-open (PAF col 3) |
| 4 | `query_end` | Alignment's query end, 0-based half-open (PAF col 4) |
| 5 | `strand` | `+`/`-` relative to the reference; `*` if unmapped (PAF col 5) |
| 6 | `target_name` | Reference/target sequence name (chromosome); `*` if unmapped (PAF col 6) |
| 7 | `target_len` | Target sequence length (PAF col 7) |
| 8 | `target_start` | Alignment's target start, 0-based half-open (PAF col 8) |
| 9 | `target_end` | Alignment's target end, 0-based half-open (PAF col 9) |
| 10 | `n_residue_matches` | Number of matching bases in the alignment (PAF col 10) |
| 11 | `aln_block_len` | Alignment block length (PAF col 11) |
| 12 | `mapq` | Mapping quality (PAF col 12) |
| 13 | `tp_tag` | PAF `tp:A` alignment-type tag: `P` primary, `S` secondary, `I` inversion of the primary. `*` when absent or unmapped |
| 14 | `as_score` | Alignment score (`AS:i` tag); `0` if absent |
| 15 | `ms_score` | minimap2's `ms:i` tag: DP score of the max-scoring segment in the alignment; `0` if absent |
| 16 | `aln_identity` | Alignment identity, `n_residue_matches / aln_block_len`; `NaN` if unmapped |
| 17 | `query_aln_cov` | `query_aln_len / query_len`; `NaN` if `query_len` is `0` |
| 18 | `query_aln_len` | `query_end - query_start` |
| 19 | `cs` | Compact alignment string (`cs:Z` tag); empty for unmapped |
| 20 | `query_transcript_strand` | Strand of the transcript **relative to the read**: `+` sense, `-` antisense, `*` unknown or unmapped. From minimap2's `ts:A` or STAR's `XS:A`; see [REFERENCE.md](../REFERENCE.md#transcript-strand) |
| 21 | `genomic_intron_strand` | Strand of the transcript, and so of every intron, **on the reference**: `strand` when `query_transcript_strand` is `+`, the opposite when `-`; `*` if unknown or unmapped |
| 22 | `n_match_events` | Number of contiguous match runs walked out of `cs` |
| 23 | `n_match_bases` | Total matched bases |
| 24 | `n_substitution_events` | Number of substitution events |
| 25 | `n_substitution_bases` | Total substituted bases |
| 26 | `n_insertion_events` | Number of insertion events |
| 27 | `n_insertion_bases` | Total inserted bases |
| 28 | `n_deletion_events` | Number of deletion events |
| 29 | `n_deletion_bases` | Total deleted bases |
| 30 | `n_softclip_events` | `0`, `1`, or `2` — how many ends are clipped |
| 31 | `n_softclip_bases_start` | Unaligned bases clipped at the read's start (strand-aware); `0` for unmapped |
| 32 | `n_softclip_bases_end` | Unaligned bases clipped at the read's end; `0` for unmapped |
| 33 | `n_intron_events` | Number of intron/splice (`~`) ops in `cs`, i.e. the number of junctions |
| 34 | `n_intron_bases` | Total intron bases |
| 35 | `query_junctions` | Query-coordinate splice-junction positions, plus-strand-adjusted, sorted ascending. Python tuple string, e.g. `(108, 359)`; `()` if none |
| 36 | `genomic_junctions` | Splice junctions in absolute reference coordinates, 0-based half-open. Python tuple-of-tuples string, e.g. `((100, 250), (400, 800))`; `()` if none. Chromosome is **not** embedded — it's the sibling `target_name` column |
| 37 | `target_start_1based` | `target_start + 1` — 1-based convenience column |

**Notes**
- `query_junctions` / `genomic_junctions` use a Python tuple-of-tuples format —
  parse with `ast.literal_eval` in Python.
- `aln_identity` is `NaN`, not `0`, when undefined (unmapped);
  `query_aln_cov` is `NaN` only when `query_len` is `0`.
- Any two of `strand`, `query_transcript_strand` and `genomic_intron_strand`
  give the third (for a mapped row); both strand columns are kept so each use
  reads one column directly. See [REFERENCE.md](../REFERENCE.md#transcript-strand).
