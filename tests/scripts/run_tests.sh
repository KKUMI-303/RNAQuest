#!/usr/bin/env bash

set -e

echo "========================================"
echo "Running RNA-seq pipeline tests"
echo "========================================"

echo ""
echo "1. Running Snakemake..."
snakemake \
    --configfile tests/config/config.yaml \
    --cores 2

echo ""
echo "2. Validating DESeq2 outputs..."
python tests/scripts/validate.py

echo ""
echo "========================================"
echo "All tests passed!"
echo "========================================"
