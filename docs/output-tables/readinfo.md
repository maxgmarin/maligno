# `readinfo` table

`{prefix}.{label}.readinfo.tsv.gz` — **per-read** table: one row per read,
describing its representative alignment, plus counts over all of that read's
alignments.

Produced by `toolkit paf2tables --readinfo` and, opt-in, by
`compare --emit-readinfo` (once per side). Off by default in `compare`.
Requires the PAF be grouped by `query_name` (contiguous runs); built by the
**readinfo collapse** step shared by both commands.

**Representative-alignment selection**: within a read's group of alignments,
the representative is the one with the highest `ms_score`, ties broken by
highest `as_score`, then by highest `mapq`; remaining ties fall through to
input order (stable sort).

**Every per-alignment column comes from the representative alone** (locus,
strand, scores, identity, coverage, junctions, `cs`, coordinates, event counts,
`tp_tag`), and has the same name as in [alninfo.md](alninfo.md). Only the four
`n_alns*` columns are computed over the read's whole group of alignments.

36 columns.

| # | Column | Description |
|---|---|---|
| 1 | `query_name` | Read/query ID — join key |
| 2 | `query_len` | Query sequence length |
| 3 | `target_name` | Representative alignment's reference/target name; `*` if unaligned |
| 4 | `strand` | Representative alignment's strand; `*` if unaligned |
| 5 | `mapq` | Representative alignment's mapping quality |
| 6 | `as_score` | Representative alignment's alignment score (`AS:i`) |
| 7 | `ms_score` | Representative alignment's `ms:i` score (DP score of the max-scoring segment) |
| 8 | `query_aln_cov` | Representative alignment's `query_aln_len / query_len`; `NaN` if `query_len` is `0` |
| 9 | `query_aln_len` | Representative alignment's aligned query length, `query_end − query_start` |
| 10 | `aln_identity` | Representative alignment's identity, `n_residue_matches / aln_block_len`; `NaN` if unaligned |
| 11 | `query_junctions` | Representative alignment's query-coordinate splice junctions (Python tuple string) |
| 12 | `n_alns` | Number of aligned rows (`target_name != "*"`) for this read; `0` if the read is present but entirely unaligned |
| 13 | `n_alns_tied_top_score` | Number of aligned rows tied with the representative on **all** selection keys (`ms_score`, `as_score` **and** `mapq`). `1` ⇒ the representative won outright; `> 1` ⇒ input order picked it among equals; `0` for an unaligned read |
| 14 | `n_query_junctions` | Number of entries in `query_junctions` |
| 15 | `n_match_events` | Representative alignment's match-run count |
| 16 | `n_match_bases` | Representative alignment's matched bases |
| 17 | `n_substitution_events` | Representative alignment's substitution-event count |
| 18 | `n_substitution_bases` | Representative alignment's substituted bases |
| 19 | `n_insertion_events` | Representative alignment's insertion-event count |
| 20 | `n_insertion_bases` | Representative alignment's inserted bases |
| 21 | `n_deletion_events` | Representative alignment's deletion-event count |
| 22 | `n_deletion_bases` | Representative alignment's deleted bases |
| 23 | `n_intron_events` | Representative alignment's intron (`~`) op count |
| 24 | `n_intron_bases` | Representative alignment's total intron bases |
| 25 | `n_softclip_bases_start` | Representative alignment's start soft-clip length |
| 26 | `n_softclip_bases_end` | Representative alignment's end soft-clip length |
| 27 | `n_softclip_events` | Representative alignment's soft-clip event count (`0`/`1`/`2`) |
| 28 | `cs` | Representative alignment's compact alignment string |
| 29 | `genomic_junctions` | Representative alignment's splice junctions in absolute reference coordinates (Python tuple-of-tuples string; chrom is the sibling `target_name` column) |
| 30 | `query_start` | Representative alignment's query-coordinate span start (0-based half-open) |
| 31 | `query_end` | Representative alignment's query-coordinate span end |
| 32 | `target_start` | Representative alignment's target-coordinate span start (0-based half-open) |
| 33 | `target_end` | Representative alignment's target-coordinate span end |
| 34 | `tp_tag` | Representative alignment's PAF `tp:A` tag (`P`/`S`/`I`); `*` if absent or unaligned |
| 35 | `n_alns_tp_p` | Number of the read's **aligned** rows whose `tp:A` is `P` or `I` (SAM primary **+ supplementary**); `0` for an unaligned read |
| 36 | `n_alns_tp_s` | Number of the read's **aligned** rows whose `tp:A` is `S` or `i` (secondary); `0` for an unaligned read |

**Notes**
- `query_start`/`query_end` + `target_start`/`target_end` + `target_name` +
  `strand` together give each read a complete BED-style alignment interval —
  useful for `bedtools`-style downstream region analysis.
- A `mapq` difference between two `readinfo` files (same reads, different
  aligner/parameters) flags reads where the two runs disagree on mapping
  uniqueness.
- **`n_alns` counts every aligned row the input contains for the read**,
  including supplementary segments, which cover other parts of the read rather
  than competing placements. `n_alns_tied_top_score` is counted over the same
  rows.
- **`n_alns_tp_p` / `n_alns_tp_s`** split `n_alns` by PAF `tp:A` type:
  - `n_alns_tp_p` counts `P` and `I` rows. In PAF (and in maligno's `sam2paf`
    output), **supplementary alignments are `tp:A:P` too**, so it is SAM
    primary + supplementary, and `n_alns_tp_p > 1` means the read is split
    into `n_alns_tp_p − 1` extra segments (chimeric / split alignment).
  - `n_alns_tp_s` counts `S` and `i` rows: secondary (alternative) placements.
  - A mapped row with no `tp` tag counts toward neither, so
    `n_alns_tp_p + n_alns_tp_s <= n_alns`; the remainder is untagged rows.
  - The representative alignment, `n_alns` and `n_alns_tied_top_score` don't
    depend on these counts.
  - **Counts reflect the input as filtered.** On SAM/BAM input, `compare`'s
    default `--sam-records-to-convert primary-supp` (`sam2paf -p`) drops
    secondary records, so `n_alns_tp_s` is always `0` unless you pass
    `--sam-records-to-convert all`; `--sam-records-to-convert primary` (`-P`)
    also drops supplementary records, so `n_alns_tp_p <= 1`.
    For minimap2 PAF, secondaries appear only if minimap2 emitted them.
