# `genomic_junction_compare.per_read.tsv.gz`

`{prefix}.genomic_junction_compare.per_read.tsv.gz` — written by
[`toolkit genomic-junction-compare`](../REFERENCE.md#toolkit-genomic-junction-compare).
One row per read that has **at least one genomic junction on either side**
(reads with none on both sides are only counted in the
[summary](genomic-junction-compare-summary.md)). Join it to the comparison
table on `query_name`.

A genomic junction is `(target_name, intron_start, intron_end, genome_strand)`.
Each side's junctions are rebuilt from its representative alignment's `cs`
tag, so their coordinates equal that side's `genomic_junctions` column in the
[comparison table](compare.md). Unlike that column, the lists here write each
junction as a complete tuple, `('chr11', 61300932, 61302256, '-')`: contig,
start, end (0-based half-open), and genome strand.

21 columns.

| # | Column | Description |
|---|---|---|
| 1 | `query_name` | Read ID |
| 2 | `target_name_A` | A's reference/target name; `*` if unmapped |
| 3 | `target_name_B` | B's reference/target name |
| 4 | `strand_A` | A's alignment strand (read relative to the reference) |
| 5 | `strand_B` | B's alignment strand |
| 6 | `genome_strand_A` | Strand of A's transcript on the reference, so of every A junction: `strand_A` when `transcript_strand_A` is `+`, the opposite when it is `-`; `*` if unknown or unmapped. See [REFERENCE.md](../REFERENCE.md#transcript-strand) |
| 7 | `genome_strand_B` | The same for B |
| 8 | `n_genomic_junctions_A` | Number of junctions in A's alignment |
| 9 | `n_genomic_junctions_B` | Number of junctions in B's alignment |
| 10 | `n_genomic_junctions_shared` | Junctions in both sides, \|A ∩ B\| |
| 11 | `n_genomic_junctions_not_shared` | Junctions in exactly one side, `n_genomic_junctions_only_A + n_genomic_junctions_only_B` |
| 12 | `n_genomic_junctions_only_A` | Junctions only in A, \|A \ B\| |
| 13 | `n_genomic_junctions_only_B` | Junctions only in B, \|B \ A\| |
| 14 | `genomic_junctions_A` | A's junctions in reference order, each `(contig, start, end, genome_strand)`, e.g. `(('chr11', 61300932, 61302256, '-'),)`; `()` if none |
| 15 | `genomic_junctions_B` | B's junctions |
| 16 | `junction_motifs_A` | A's splice-site motifs, one per junction in the same order, e.g. `('gt-ag', 'gc-ag')` |
| 17 | `junction_motifs_B` | B's splice-site motifs |
| 18 | `genomic_junctions_only_A` | The junctions only in A, as `(contig, start, end, genome_strand)` tuples, sorted |
| 19 | `genomic_junctions_only_B` | The junctions only in B |
| 20 | `cs_A` | A's compact alignment string |
| 21 | `cs_B` | B's compact alignment string |

**Notes**
- **The strand is part of a junction.** The same coordinates with different
  genome strands are different junctions, and `*` matches only `*`. So when
  `genome_strand_A` and `genome_strand_B` differ, every junction is "only A"
  or "only B", and the same `(start, end)` can appear in both
  `genomic_junctions_only_*` lists, with different strands in their tuples.
- **Motifs** are the 2 + 2 intron-edge letters from `cs`, written
  donor-acceptor in the junction's orientation: on a `-` strand junction the
  reference letters are reverse-complemented, so a canonical intron reads
  `gt-ag` on either strand. With `genome_strand` `*` they are written
  reference-forward. `nn-nn` means the aligner didn't report them (STAR, or
  `cs` that `sam2paf` built from `MD`).
- Lists use Python tuple syntax; parse them with `ast.literal_eval`.
