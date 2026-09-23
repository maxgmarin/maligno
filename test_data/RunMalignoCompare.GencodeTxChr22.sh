#!/usr/bin/env bash
set -euo pipefail

maligno compare \
  -a Splice.AlnToHG38.PriAln.paf.gz \
  -b SpliceHQ.AlnToHG38.PriAln.paf.gz \
  --label-a Splice --label-b SpliceHQ \
  --outdir test_results --prefix Splice_vs_SpliceHQ \
  --emit-alninfo --emit-readinfo

maligno compare-toolkit query-junction-diff \
  -i test_results/Splice_vs_SpliceHQ.compare.parquet \
  --outdir test_results --prefix Splice_vs_SpliceHQ
