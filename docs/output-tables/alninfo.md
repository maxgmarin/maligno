# `alninfo` table

`{prefix}.{label}.alninfo.tsv.gz` — **per-alignment** table: one row per PAF
alignment record (every alignment for every read, not just the best one).

Produced by `compare-toolkit paf2tables --alninfo` and, opt-in, by
`compare --emit-alninfo` (once per side, i.e. `{prefix}.{label_a}.alninfo.tsv.gz`
and `{prefix}.{label_b}.alninfo.tsv.gz`). Off by default in `compare` — see
[REFERENCE.md](../REFERENCE.md) for why. Pure streaming output: rows are
written in PAF input order, one PAF line in, one row out (malformed lines are
skipped with a stderr warning). An unmapped PAF record (`Target_Name == "*"`)
still produces a full row, with alignment-derived fields zeroed/`NaN`.

36 columns.

| # | Column | Description |
|---|---|---|
| 1 | `Query_Name` | Read/query sequence ID (PAF col 1) |
| 2 | `Query_Len` | Query sequence length (PAF col 2) |
| 3 | `Query_Start` | Alignment's query start, 0-based half-open (PAF col 3) |
| 4 | `Query_End` | Alignment's query end, 0-based half-open (PAF col 4) |
| 5 | `Strand` | `+`/`-` relative to the reference; `*` if unmapped (PAF col 5) |
| 6 | `Target_Name` | Reference/target sequence name (chromosome); `*` if unmapped (PAF col 6) |
| 7 | `Target_Len` | Target sequence length (PAF col 7) |
| 8 | `Target_Start` | Alignment's target start, 0-based half-open (PAF col 8) |
| 9 | `Target_End` | Alignment's target end, 0-based half-open (PAF col 9) |
| 10 | `Num_Residue_Matches` | Number of matching bases in the alignment (PAF col 10) |
| 11 | `Aln_Block_Len` | Alignment block length (PAF col 11) |
| 12 | `MQ` | Mapping quality (PAF col 12) |
| 13 | `ms` | Chaining/minimizer score (`ms:i` tag); `0` if absent |
| 14 | `AS` | Alignment score (`AS:i` tag); `0` if absent |
| 15 | `cs` | Compact alignment string (`cs:Z` tag); empty for unmapped |
| 16 | `N_Match_Events` | Number of contiguous match runs walked out of `cs` |
| 17 | `N_Match_Bases` | Total matched bases |
| 18 | `N_Substitution_Events` | Number of substitution events |
| 19 | `N_Substitution_Bases` | Total substituted bases |
| 20 | `N_Insertion_Events` | Number of insertion events |
| 21 | `N_Insertion_Bases` | Total inserted bases |
| 22 | `N_Deletion_Events` | Number of deletion events |
| 23 | `N_Deletion_Bases` | Total deleted bases |
| 24 | `N_Splice_Junction_Events` | Number of intron/splice (`~`) ops in `cs` |
| 25 | `N_Splice_Junction_Bases` | Total intron bases |
| 26 | `N_SoftClipped_Bases_Start` | Unaligned bases clipped at the read's start (strand-aware); `0` for unmapped |
| 27 | `N_SoftClipped_Bases_End` | Unaligned bases clipped at the read's end; `0` for unmapped |
| 28 | `N_SoftClipped_Events` | `0`, `1`, or `2` — how many ends are clipped |
| 29 | `junctions` | Query-coordinate splice-junction positions, plus-strand-adjusted, sorted ascending. Python tuple string, e.g. `(108, 359)`; `()` if none |
| 30 | `splice_junction_count` | Number of entries in `junctions` |
| 31 | `Target_Start_1based` | `Target_Start + 1` — 1-based convenience column |
| 32 | `seqid` | Sequence identity, `Num_Residue_Matches / Aln_Block_Len`; `NaN` if unmapped |
| 33 | `Query_Aln_Len` | `Query_End - Query_Start` |
| 34 | `Query_Aln_Cov` | `Query_Aln_Len / Query_Len`; `NaN` if `Query_Len` is `0` |
| 35 | `genomic_junctions` | Splice junctions in absolute reference coordinates, 0-based half-open. Python tuple-of-tuples string, e.g. `((100, 250), (400, 800))`; `()` if none. Chromosome is **not** embedded — it's the sibling `Target_Name` column |
| 36 | `tp_tag` | PAF `tp:A` alignment-type tag: `P` primary, `S` secondary, `I` inversion of the primary. `*` when absent or unmapped |

**Notes**
- `junctions` / `genomic_junctions` use a Python tuple-of-tuples format —
  parse with `ast.literal_eval` in Python.
- `seqid` and `Query_Aln_Cov` are `NaN`, not `0`, when undefined (unmapped, or
  a zero-length denominator).
- Column count bumped 35 → 36 when `tp_tag` was added (v0.25.0).
