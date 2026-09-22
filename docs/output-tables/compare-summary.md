# `compare.summary.tsv`

`{prefix}.compare.summary.tsv` — aggregate summary statistics over the
comparison table, tallied as `compare` streams its rows (constant memory).
2 columns: `Category`, `Count`. The same schema is written by
`compare-toolkit summary` (computed from an existing `compare.tsv`/`.parquet`)
and, with `space`/`compare_by` provenance rows prepended, by
`compare-toolkit find-aln-diff`.

Every `Category` below is one **row** in the file, in this order:

| Category | Meaning |
|---|---|
| `label_A` | Which dataset side A is (the `--label-a` value). String, not a count — provenance row |
| `label_B` | Which dataset side B is (the `--label-b` value) |
| `reads_compared` | Matched reads written to the comparison table |
| `aligned_both` | Representative alignment mapped in both sets |
| `aligned_only_A` | Mapped in A, `"*"` in B |
| `aligned_only_B` | Mapped in B, `"*"` in A |
| `aligned_neither` | Unmapped (`"*"`) in both |
| `query_identical` | Both mapped with the same query span and same alignment relative to the read (same-strand `cs` match, or reverse-complement `cs` match), motif-blind — **or** neither side mapped at all (`aligned_neither`; both aligners agreeing a read doesn't map is agreement, not disagreement) |
| `query_identical_same_strand` | …of the both-mapped identical reads, reached via the same-strand branch (excludes `aligned_neither`) |
| `query_identical_revcomp` | …of the both-mapped identical reads, reached via the reverse-complement branch (excludes `aligned_neither`) |
| `query_not_identical` | Both mapped but not `query_identical` — i.e. `aligned_both - query_identical_same_strand - query_identical_revcomp` (not gated on `aligned_neither`, which is never part of this count) |
| `query_junctions_identical` | Both mapped, **query-space** splice-junction sets match — a looser criterion than `query_identical` (ignores mismatches/indels/soft-clips) |
| `query_junctions_not_identical` | Both mapped but query-space junction sets differ |
| `ref_same_position_same_aln` | Same `TargetChr`/`Strand`/`Target_Start` **and** same `cs` (motif-blind) — reference-identical |
| `ref_same_position_diff_aln` | Same position, different `cs` |
| `ref_diff_position_same_aln` | Same `cs`, different position — relocated |
| `ref_diff_position_diff_aln` | Both position and `cs` differ |
| `ref_same_position_same_junctions` | Among same-position reads, **genomic-coordinate** junction sets also match |
| `ref_same_position_diff_junctions` | Among same-position reads, genomic-coordinate junction sets differ |
| `present_only_in_A_by_id` | Read present in only A's PAF (built-in `compare` tally only; `0` unless `--allow-id-mismatch`) |
| `present_only_in_B_by_id` | Read present in only B's PAF |

**Notes**
- Classification is computed from each side's **representative (best)
  alignment** using `TargetChr`, `Strand`, `cs`, `Query_Start`/`Query_End`,
  `Target_Start`/`Target_End`, `junctions`, `genomic_junctions` — the same
  columns from the [compare table](compare.md).
- Reference-space classification is strict and literal: no
  reverse-complement accommodation, so a real strand difference always
  counts as a different position.
- `query_junctions_identical`/`ref_same_position_same_junctions` are looser,
  complementary axes to the `query_identical`/`ref_*_aln` categories —
  `None`/not-counted unless both sides are mapped.
- `query_identical = query_identical_same_strand + query_identical_revcomp +
  aligned_neither` always holds — `aligned_neither`'s full count is folded
  into `query_identical` without a dedicated sub-category row of its own
  (its value is already exactly the `aligned_neither` row above).
- `compare` prints an abbreviated version of this file to stderr
  (`label_A`/`label_B`, `reads_compared`, `aligned_*`, `query_identical`,
  `query_not_identical`) — the on-disk TSV always carries the full row set
  above.
