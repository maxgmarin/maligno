# `{query|reference}_aln_diff_regions.{A,B}.bed`

`{prefix}.query_aln_diff_regions.A.bed.gz` / `.B.bed.gz` (`--mode any-diff-query`),
or `{prefix}.reference_aln_diff_regions.A.bed.gz` / `.B.bed.gz`
(`--mode any-diff-reference`) — merged genomic loci where the differing reads
from [`aln_diff_reads.tsv`](aln-diff-reads.md) cluster, one file per side (A's
own placement / B's own placement). Written by `toolkit find-aln-diff`.

Each differing read contributes its representative alignment's span on that
side (`target_start` to `target_end`, introns included). Spans are merged with
a sort and single sweep, like `bedtools merge -c -o count`: overlapping or
touching spans join one locus. A read mapped on both sides contributes to
**both** files; a `diff_aln_only_A` read only to the A file, and a
`diff_aln_only_B` read only to the B file. Identical reads never contribute.

8 columns.

| # | Column | Description |
|---|---|---|
| 1 | `#chrom` | Chromosome/contig of the merged locus |
| 2 | `start` | Locus start, 0-based half-open |
| 3 | `end` | Locus end |
| 4 | `n_diff_aln_total` | Differing reads whose span falls in this locus |
| 5 | `n_diff_aln_both` | Of those, how many are mapped in **both** A and B |
| 6 | `n_diff_aln_only_<A\|B>` | Of those, how many are mapped on **only this side** (`n_diff_aln_only_A` in the A file, `n_diff_aln_only_B` in the B file) |
| 7 | `n_diff_aln_plus_strand` | Of `n_diff_aln_total`, how many are on the `+` strand on this side |
| 8 | `n_diff_aln_minus_strand` | Of `n_diff_aln_total`, how many are on the `-` strand on this side |

**Invariants**
- `n_diff_aln_total = n_diff_aln_both + n_diff_aln_only_<A|B>`
- `n_diff_aln_plus_strand + n_diff_aln_minus_strand <= n_diff_aln_total`

**Notes**
- Rows are sorted by `#chrom` (byte order), then `start`.
- Spans that can't be parsed, or where `end <= start`, are skipped and their
  number is printed on stderr.
- When most reads differ, the loci merge into roughly one per contig. The
  tables are most informative when the differing set is small.
