import csv
import re

gtf = "tests/data/reference/mini.gtf"
output = "tests/data/reference/tx2gene.csv"

with open(gtf) as source, open(output, "w", newline="") as target:
    writer = csv.writer(target)
    writer.writerow(["transcript_id", "gene_id"])

    seen = set()
    for line in source:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or fields[2] != "transcript":
            continue

        tx = re.search(r'transcript_id "([^"]+)"', fields[8])
        gene = re.search(r'gene_id "([^"]+)"', fields[8])

        if tx and gene and tx.group(1) not in seen:
            writer.writerow([tx.group(1), gene.group(1)])
            seen.add(tx.group(1))

print(f"Wrote {len(seen)} transcript-to-gene mappings to {output}")
