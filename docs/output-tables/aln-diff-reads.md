# `{query|reference}_aln_diff_reads.tsv`

`{prefix}.query_aln_diff_reads.tsv.gz` (`--mode any-diff-query`) /
`{prefix}.reference_aln_diff_reads.tsv.gz` (`--mode any-diff-reference`) — one
row per read whose core alignment **differs** between set A and set B, plus 7
classification columns. Written by `toolkit find-aln-diff`; see
[REFERENCE.md](../REFERENCE.md#toolkit-find-aln-diff) for the two modes.

With `--emit-identical-reads`, `{prefix}.query_aln_identical_reads.tsv.gz` /
`{prefix}.reference_aln_identical_reads.tsv.gz` list the reads that don't
differ, in the same layout.

9 columns.

| # | Column | Description |
|---|---|---|
| 1 | `query_name` | Read ID |
| 2 | `outcome` | Differing reads: `diff_aln_both` (mapped in both sets, not identical under the mode), `diff_aln_only_A` (mapped in A only) or `diff_aln_only_B` (mapped in B only). Identical reads: `identical_same_strand`, `identical_revcomp` (`any-diff-query` only) or `neither_mapped` |
| 3 | `query_identical_same_strand` | `1` if both mapped with the same query span and `cs` on the same strand |
| 4 | `query_identical_revcomp` | `1` if both mapped with the same query span, on opposite strands, with B's `cs` the exact reverse complement of A's |
| 5 | `query_junctions_identical` | `1` if both mapped and the **query-space** splice-junction sets match |
| 6 | `ref_same_position_same_aln` | `1` if same `target_name`/`strand`/`target_start` and same `cs` |
| 7 | `ref_same_position_diff_aln` | `1` if same position, different `cs` |
| 8 | `ref_diff_position_same_aln` | `1` if same `cs`, different position (relocated) |
| 9 | `ref_diff_position_diff_aln` | `1` if both position and `cs` differ |

**Notes**
- Columns 3–9 are computed the same way in both modes, so one run can show,
  e.g., a read that differs in reference space but is identical in query
  space. All 7 are `0` for a read mapped on one side or neither.
- For a both-mapped read, exactly one of columns 6–9 is `1`.
- `cs` comparisons are motif-blind: intron donor/acceptor letters are ignored.
- With `--aligned-both-only`, only reads mapped in both sets appear in either
  table.
- The genomic regions where these reads cluster are in
  [aln-diff-regions.md](aln-diff-regions.md).
