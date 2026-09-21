# `query_diff_regions.{A,B}.bed`

`{prefix}.query_diff_regions.A.bed.gz` / `{prefix}.query_diff_regions.B.bed.gz`
— merged genomic loci where the differing reads from
[`query_diff_reads.tsv`](query-diff-reads.md) cluster, one file per
coordinate side (A's own placement / B's own placement). Written by
`compare` by default (same fused pass as the reads table) and by standalone
`compare-toolkit find-aln-diff`.

Loci are formed by a generic sort + single-sweep interval merge
(`bedtools merge -c -o count` equivalent), over the genomic interval of
every *differing* read that has a placement on this side — a both-mapped
read contributes to **both** the A and B region tables; a `diff_aln_only_A`
read contributes only to the A table (symmetric for `diff_aln_only_B`).
Identical reads never contribute an interval and are entirely absent from
these tables.

10 columns.

| # | Column | Description |
|---|---|---|
| 1 | `#chrom` | Chromosome/contig of the merged locus |
| 2 | `start` | Locus start, 0-based half-open |
| 3 | `end` | Locus end |
| 4 | `n_diff_aln_total` | Total differing reads whose interval falls in this locus |
| 5 | `n_diff_aln_to_both` | Of those, how many mapped on **both** A and B (as opposed to only this side) |
| 6 | `n_diff_aln_only_<A\|B>` | Of those, how many mapped on **only this side** (column name is `n_diff_aln_only_A` in the A-side file, `n_diff_aln_only_B` in the B-side file) |
| 7 | `n_diff_aln_both_junctions_differ` | Of `n_diff_aln_to_both`, how many have differing **query-space** splice junctions |
| 8 | `n_diff_aln_both_junctions_same` | Of `n_diff_aln_to_both`, how many have matching query-space splice junctions (so the difference is elsewhere — mismatch/indel/soft-clip) |
| 9 | `n_diff_aln_plus_strand` | Of `n_diff_aln_total`, how many are `+` strand on this side |
| 10 | `n_diff_aln_minus_strand` | Of `n_diff_aln_total`, how many are `-` strand on this side |

**Invariants**
- `n_diff_aln_total = n_diff_aln_to_both + n_diff_aln_only_<A|B>`
- `n_diff_aln_to_both = n_diff_aln_both_junctions_differ + n_diff_aln_both_junctions_same`
- `n_diff_aln_plus_strand + n_diff_aln_minus_strand <= n_diff_aln_total`

**Notes**
- `n_diff_aln_both_junctions_differ`/`n_diff_aln_both_junctions_same` use the
  **query-space** junction comparison regardless of the active
  `--space`/`--compare-by` that selected these reads as differing — so under
  `--space reference` a read can be junction-identical in query space yet
  still be counted here (it was flagged as a *reference*-space difference).
- Bad intervals (unparseable coordinates, or `end <= start`) are skipped and
  counted internally rather than aborting the run.
- Verified against a real `bedtools merge` oracle at both small (~11.6K
  reads) and genome scale (~986K reads, 31.5K differing) — exact match on
  `(chrom, start, end, n_diff_aln_total)` in both cases.
