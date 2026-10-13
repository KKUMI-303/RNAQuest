#!/usr/bin/env bash
set -euo pipefail

N=20000
OUT="tests/data/raw/mock"

mkdir -p "$OUT"

make_mock() {
    sample="$1"
    input="$2"
    seed="$3"

    seqtk sample -s"$seed" \
        "tests/data/raw/${input}_1.fastq" \
        "$N" > "$OUT/${sample}_R1.fastq"

    seqtk sample -s"$seed" \
        "tests/data/raw/${input}_2.fastq" \
        "$N" > "$OUT/${sample}_R2.fastq"

    gzip "$OUT/${sample}_R1.fastq"
    gzip "$OUT/${sample}_R2.fastq"
}

make_mock "ConditonA1" "SRR24142039" 101
make_mock "ConditonA2" "SRR24142040" 102
make_mock "ConditonA3" "SRR24142041" 103

make_mock "ConditionB1" "SRR24142033" 104
make_mock "ConditionB2" "SRR24142034" 105
make_mock "ConditionB3" "SRR24142035" 106
