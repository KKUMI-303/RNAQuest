sed -i '1c\>chr1' tests/data/reference/mini.fa#!/usr/bin/env bash
set -euo pipefail

REF="data/reference"
OUT="tests/data/reference"

CHROM="1"
GTF_CHROM="chr1"
START=74626055
END=75626054

mkdir -p "$OUT"

# Index the original genome if needed
samtools faidx "$REF/mm39.fa"

# Extract the genomic interval
samtools faidx "$REF/mm39.fa" \
    "${CHROM}:${START}-${END}" \
    > "$OUT/mini.fa"

# Rename the FASTA header to match the GTF chromosome name
sed -i '1c\>chr1' tests/data/reference/mini.fa

# Extract features fully contained within the selected interval.
# Convert their coordinates to the mini-reference coordinate system.
awk -F '\t' -v OFS='\t' \
    -v chrom="$GTF_CHROM" \
    -v start="$START" -v end="$END" '
    /^#/ { next }
    $1 == chrom && $4 >= start && $5 <= end {
        $4 = $4 - start + 1
        $5 = $5 - start + 1
        print
    }
' "$REF/annotation.gtf" > "$OUT/mini.gtf"

echo "Mini-reference created:"
grep '^>' "$OUT/mini.fa"
echo "FASTA size:"
du -h "$OUT/mini.fa"
echo "GTF feature counts:"
awk -F '\t' '{n[$3]++} END {
    for (type in n) print type, n[type]
}' "$OUT/mini.gtf"
