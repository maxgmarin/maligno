#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Check how many read names two readinfo tables share.
#
# Use it when `toolkit merge-readinfo` reports fewer `reads_compared` than
# expected: if `reads_compared` is lower than the shared count printed here, the
# two files are not sorted by Read_Name in the same byte order.
#
# Usage: check-readinfo-overlap.sh a.readinfo.tsv[.gz] b.readinfo.tsv[.gz]
# -----------------------------------------------------------------------------
set -euo pipefail

if [ $# -ne 2 ]; then
  echo "usage: $0 a.readinfo.tsv[.gz] b.readinfo.tsv[.gz]" >&2
  exit 1
fi

a="$1"
b="$2"

decompress() {
  case "$1" in
    *.gz) zcat < "$1" ;;
    *)    cat    "$1" ;;
  esac
}

# Unique Read_Name values (column 1) from each file.
a_names=$(mktemp); b_names=$(mktemp)
trap 'rm -f "$a_names" "$b_names"' EXIT

decompress "$a" | tail -n +2 | cut -f1 | LC_ALL=C sort -u > "$a_names"
decompress "$b" | tail -n +2 | cut -f1 | LC_ALL=C sort -u > "$b_names"

n_a=$(wc -l < "$a_names" | tr -d ' ')
n_b=$(wc -l < "$b_names" | tr -d ' ')
n_shared=$(LC_ALL=C comm -12 "$a_names" "$b_names" | wc -l | tr -d ' ')

cat <<EOF
Read-name overlap
  A: $a
  B: $b

  read names in A:      $n_a
  read names in B:      $n_b
  shared read names:    $n_shared

  ⇒ \`toolkit merge-readinfo\`'s reads_compared should equal $n_shared. If it is lower,
    the two files are not sorted by Read_Name in the same byte order. Re-sort each with:

        LC_ALL=C sort -t\$'\\t' -k1,1 input.readinfo.tsv > sorted.tsv
        # (keep the header separately)

    or regenerate both from PAFs sorted with \`LC_ALL=C sort -t\$'\\t' -k1,1\`
    using \`toolkit paf2tables --readinfo\`.
EOF
