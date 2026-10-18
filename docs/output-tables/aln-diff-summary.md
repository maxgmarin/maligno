# `{query|reference}_aln_diff_summary.tsv`

`{prefix}.query_aln_diff_summary.tsv` (`--mode any-diff-query`) /
`{prefix}.reference_aln_diff_summary.tsv` (`--mode any-diff-reference`) — the
mode and parameters of a `toolkit find-aln-diff` run, then how many reads it
found identical and differing. Always uncompressed. The same block is printed
to stderr. 2 columns: `category`, `count`.

For the overall classification of the table (mapping status, query and
reference identity, junction identity), use `toolkit compare-summary`; see
[compare-summary.md](compare-summary.md).

| Category | Meaning |
|---|---|
| `label_A` / `label_B` | the set names, from the table's `label_A` / `label_B` columns |
| `mode` | `any-diff-query` or `any-diff-reference` |
| `aligned_both_only` | `true` with `--aligned-both-only`, otherwise `false` |
| `reads_in_table` | rows in the comparison table |
| `reads_not_aligned_both_excluded` | reads left out by `--aligned-both-only` (mapped on one side or neither); `0` without it |
| `reads_compared` | `reads_in_table - reads_not_aligned_both_excluded` |
| `reads_identical` | reads whose alignments are the same under the mode |
| `identical_same_strand` | …identical on the same strand |
| `identical_revcomp` | …identical as an exact reverse complement on opposite strands (`any-diff-query` only; `0` in `any-diff-reference`) |
| `neither_mapped` | …mapped in neither set |
| `reads_differ` | reads whose alignments differ: the rows of the differing-reads table |
| `diff_aln_both` | …mapped in both sets |
| `diff_aln_only_A` / `diff_aln_only_B` | …mapped in one set only |
| `n_regions_A` / `n_regions_B` | rows in each region table |

**Invariants**
- `reads_compared = reads_identical + reads_differ`
- `reads_identical = identical_same_strand + identical_revcomp + neither_mapped`
- `reads_differ = diff_aln_both + diff_aln_only_A + diff_aln_only_B`
- Without `--aligned-both-only`, `diff_aln_only_A` / `diff_aln_only_B` /
  `neither_mapped` equal `compare-summary`'s `aligned_only_A` /
  `aligned_only_B` / `aligned_neither`. In `any-diff-query` mode,
  `diff_aln_both` equals its `query_not_identical`.
