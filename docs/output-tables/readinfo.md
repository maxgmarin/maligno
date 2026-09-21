# `readinfo` table

`{prefix}.{label}.readinfo.tsv.gz` — **per-read** table: one row per read,
carrying its best (representative) alignment plus a few aggregates over all
of that read's alignments.

Produced by `compare-toolkit paf2tables --readinfo` and, opt-in, by
`compare --emit-readinfo` (once per side). Off by default in `compare`.
Requires the PAF be grouped by `Query_Name` (contiguous runs); built by the
**readinfo collapse** step shared by both commands.

**Best-alignment selection**: within a read's group of alignments, the
representative row is the one with the highest `ms`, ties broken by highest
`AS`, then by highest `MQ`; remaining ties fall through to input order
(stable sort). All per-alignment columns below (locus, strand, `cs`,
junctions, coordinates, event counts, `tp_tag`) come from that single
winning row — only `AS_Max`/`ms_Max`/`Query_Aln_Cov_Max`/`Query_Aln_Len_Max`/
`seqid_Max`/`Num_Aln`/`Num_Aln_MaxScore` are aggregated over the *whole*
group.

34 columns.

| # | Column | Description |
|---|---|---|
| 1 | `Read_Name` | Read/query ID — join key |
| 2 | `Read_Len` | Query sequence length — join key (matches with `Read_Name`) |
| 3 | `TargetChr` | Best alignment's reference/target name; `*` if unaligned |
| 4 | `Strand` | Best alignment's strand; `*` if unaligned |
| 5 | `MQ_Best` | Best alignment's mapping quality |
| 6 | `AS_Max` | Highest `AS` across **all** alignments of this read |
| 7 | `ms_Max` | Highest `ms` across all alignments of this read |
| 8 | `Query_Aln_Cov_Max` | Highest `Query_Aln_Cov` across all alignments |
| 9 | `Query_Aln_Len_Max` | Highest `Query_Aln_Len` across all alignments |
| 10 | `seqid_Max` | Highest `seqid` across all alignments |
| 11 | `junctions` | Best alignment's query-coordinate splice junctions (Python tuple string) |
| 12 | `Num_Aln` | Count of aligned rows (`Target_Name != "*"`) for this read; `0` if the read is present but entirely unaligned |
| 13 | `Num_Aln_MaxScore` | Count of alignments tied at the chosen-best sort key (`ms` **and** `AS`). `1` ⇒ unambiguous winner; `>1` ⇒ tie broken by `MQ`/file order |
| 14 | `JuncCount` | Number of entries in `junctions` |
| 15 | `N_Match_Events` | Best alignment's match-run count |
| 16 | `N_Match_Bases` | Best alignment's matched bases |
| 17 | `N_Substitution_Events` | Best alignment's substitution-event count |
| 18 | `N_Substitution_Bases` | Best alignment's substituted bases |
| 19 | `N_Insertion_Events` | Best alignment's insertion-event count |
| 20 | `N_Insertion_Bases` | Best alignment's inserted bases |
| 21 | `N_Deletion_Events` | Best alignment's deletion-event count |
| 22 | `N_Deletion_Bases` | Best alignment's deleted bases |
| 23 | `N_Splice_Junction_Events` | Best alignment's intron/splice-op count |
| 24 | `N_Splice_Junction_Bases` | Best alignment's total intron bases |
| 25 | `N_SoftClipped_Bases_Start` | Best alignment's start soft-clip length |
| 26 | `N_SoftClipped_Bases_End` | Best alignment's end soft-clip length |
| 27 | `N_SoftClipped_Events` | Best alignment's soft-clip event count (`0`/`1`/`2`) |
| 28 | `cs` | Best alignment's compact alignment string |
| 29 | `genomic_junctions` | Best alignment's splice junctions in absolute reference coordinates (Python tuple-of-tuples string; chrom is the sibling `TargetChr` column) |
| 30 | `Query_Start` | Best alignment's query-coordinate span start (0-based half-open) |
| 31 | `Query_End` | Best alignment's query-coordinate span end |
| 32 | `Target_Start` | Best alignment's target-coordinate span start (0-based half-open) |
| 33 | `Target_End` | Best alignment's target-coordinate span end |
| 34 | `tp_tag` | Best alignment's PAF `tp:A` tag (`P`/`S`/`I`); `*` if absent or unaligned |

**Notes**
- `Query_Start`/`Query_End` + `Target_Start`/`Target_End` + `TargetChr` +
  `Strand` together give each read a complete BED-style alignment interval —
  useful for `bedtools`-style downstream region analysis.
- A `MQ_Best` difference between two `readinfo` files (same reads, different
  aligner/parameters) flags reads where the two runs disagree on mapping
  uniqueness.
- Column count bumped 33 → 34 when `tp_tag` was added (v0.25.0).
