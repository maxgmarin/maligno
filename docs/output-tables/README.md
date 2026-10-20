# Output-table specs

Column-by-column specifications for every file `maligno compare` (and its
`toolkit` building blocks) can write. Each doc gives, for its file:
what produces it, when it's written, and a table with one row per column.

| File | Spec |
|---|---|
| `{prefix}.{label}.alninfo.tsv.gz` | [alninfo.md](alninfo.md) |
| `{prefix}.{label}.readinfo.tsv.gz` | [readinfo.md](readinfo.md) |
| `{prefix}.compare.parquet` / `.compare.tsv.gz` | [compare.md](compare.md) |
| `{prefix}.compare.summary.tsv` | [compare-summary.md](compare-summary.md) |
| `{prefix}.{query\|reference}_aln_diff_reads.tsv.gz` (and `_aln_identical_reads`) | [aln-diff-reads.md](aln-diff-reads.md) |
| `{prefix}.{query\|reference}_aln_diff_regions.{A,B}.bed.gz` | [aln-diff-regions.md](aln-diff-regions.md) |
| `{prefix}.{query\|reference}_aln_diff_summary.tsv` | [aln-diff-summary.md](aln-diff-summary.md) |
| `{prefix}.query_junction_compare.summary.tsv` | [query-junction-compare-summary.md](query-junction-compare-summary.md) |
| `{prefix}.query_junction_compare.per_read_per_junc_info.tsv.gz` | [per-read-query-junction-compare.md](per-read-query-junction-compare.md) |
| `{prefix}.query_junction_compare.unmatched_junctions.{A,B}.tsv.gz` | [query-junction-compare-unmatched.md](query-junction-compare-unmatched.md) |
| `{prefix}.query_junction_compare.all_junctions.{A,B}.tsv.gz` | [query-junction-compare-all-junctions.md](query-junction-compare-all-junctions.md) |
| `{prefix}.genomic_junction_compare.per_read.tsv.gz` | [genomic-junction-compare-per-read.md](genomic-junction-compare-per-read.md) |
| `{prefix}.genomic_junction_compare.per_junction.tsv.gz` | [genomic-junction-compare-per-junction.md](genomic-junction-compare-per-junction.md) |
| `{prefix}.genomic_junction_compare.summary.tsv` | [genomic-junction-compare-summary.md](genomic-junction-compare-summary.md) |

These are per-column references. For how each command works and the
classification logic behind the derived columns, see the main
[REFERENCE.md](../REFERENCE.md).
