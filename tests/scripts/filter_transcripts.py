import re

fasta_path = "tests/data/reference/transcripts.fa"
gtf_path = "tests/data/reference/mini.gtf"
output_path = "tests/data/reference/transcripts.matched.fa"

# Collect transcript IDs from GTF transcript records.
gtf_ids = set()

with open(gtf_path) as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9 or fields[2] != "transcript":
            continue

        match = re.search(r'transcript_id "([^"]+)"', fields[8])

        if match:
            gtf_ids.add(match.group(1))

# Filter FASTA records without requiring external libraries.
kept = 0
excluded = 0
write_record = False

with open(fasta_path) as source, open(output_path, "w") as target:
    for line in source:
        if line.startswith(">"):
            transcript_id = line[1:].strip().split()[0]
            write_record = transcript_id in gtf_ids

            if write_record:
                kept += 1
                target.write(line)
            else:
                excluded += 1
        elif write_record:
            target.write(line)

print(f"Retained transcripts: {kept}")
print(f"Excluded transcripts: {excluded}")
print(f"Output: {output_path}")
