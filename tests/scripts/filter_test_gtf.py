#!/usr/bin/env python3

import csv
import re
import statistics

counts_file = "tests/results/counts/gene_counts.txt"
gtf_file = "tests/data/reference/mini.gtf"
output_gtf = "tests/data/reference/mini.expressed.gtf"

# Read featureCounts and identify genes with variable counts.
variable_genes = set()

with open(counts_file) as handle:
    lines = (line for line in handle if not line.startswith("#"))
    reader = csv.reader(lines, delimiter="\t")
    header = next(reader)

    # featureCounts columns:
    # Geneid, Chr, Start, End, Strand, Length, then sample counts.
    sample_columns = header[6:]

    if len(sample_columns) != 6:
        raise ValueError(
            f"Expected 6 sample columns, found {len(sample_columns)}"
        )

    for row in reader:
        if len(row) < 7:
            continue

        gene_id = row[0]
        values = [float(x) for x in row[6:]]

        if len(values) == 6 and statistics.pvariance(values) > 0:
            variable_genes.add(gene_id)

if len(variable_genes) < 3:
    raise RuntimeError(
        f"Only {len(variable_genes)} variable genes found. "
        "Inspect the counts before filtering the annotation."
    )

# Retain every GTF feature belonging to a selected gene.
written = 0

with open(gtf_file) as source, open(output_gtf, "w") as target:
    for line in source:
        if line.startswith("#"):
            target.write(line)
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            continue

        match = re.search(r'gene_id "([^"]+)"', fields[8])

        if match and match.group(1) in variable_genes:
            target.write(line)
            written += 1

print(f"Variable genes retained: {len(variable_genes)}")
print(f"GTF features retained: {written}")
print(f"Output: {output_gtf}")
