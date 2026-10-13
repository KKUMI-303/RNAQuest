#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(polyester)
    library(Biostrings)
})

set.seed(20261012)

# -------------------------
# Paths and parameters
# -------------------------
REF_DIR <- "tests/data/reference"
OUT_DIR <- "tests/data/raw/mock"
SIM_DIR <- tempfile("polyester_")

TRANSCRIPTS <- file.path(REF_DIR, "transcripts.matched.fa")
GTF <- file.path(REF_DIR, "mini.gtf")
METADATA <- "tests/data/metadata.csv"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SIM_DIR, recursive = TRUE, showWarnings = FALSE)

on.exit(unlink(SIM_DIR, recursive = TRUE), add = TRUE)

READ_LENGTH <- 100
FRAGMENT_LENGTH <- 200
FRAGMENT_SD <- 20
SEED <- 20261012

# -------------------------
# Load transcript sequences
# -------------------------
tx <- readDNAStringSet(TRANSCRIPTS)
names(tx) <- sub("\\s.*$", "", names(tx))

# Exclude transcripts too short for the intended fragments.
keep <- width(tx) >= FRAGMENT_LENGTH
tx <- tx[keep]

if (length(tx) == 0) {
    stop("No transcripts long enough for simulation.")
}

cat("Transcripts available:", length(tx), "\n")

# -------------------------
# Map transcripts to genes
# -------------------------
gtf <- read.delim(
    GTF,
    header = FALSE,
    sep = "\t",
    comment.char = "#",
    quote = "",
    stringsAsFactors = FALSE
)

# Transcript rows contain transcript_id and gene_id.
gtf <- gtf[gtf[[3]] == "transcript", ]

get_attribute <- function(x, key) {
    pattern <- paste0(key, ' "([^"]+)"')
    matches <- regmatches(x, regexec(pattern, x))

    vapply(
        matches,
        function(m) {
            if (length(m) >= 2) m[2] else NA_character_
        },
        character(1)
    )
}

tx_map <- data.frame(
    transcript_id = get_attribute(gtf[[9]], "transcript_id"),
    gene_id = get_attribute(gtf[[9]], "gene_id"),
    gene_type = get_attribute(gtf[[9]], "gene_type"),
    stringsAsFactors = FALSE
)

tx_map <- unique(tx_map)
tx_map <- tx_map[!is.na(tx_map$transcript_id) &
                 !is.na(tx_map$gene_id), ]

idx <- match(names(tx), tx_map$transcript_id)

if (anyNA(idx)) {
    stop(
        "Some transcript FASTA IDs could not be mapped to the GTF. ",
        "Check the FASTA headers and GTF transcript_id attributes."
    )
}

gene_ids <- tx_map$gene_id[idx]
gene_types <- tx_map$gene_type[idx]

# Select genes with protein-coding transcripts.
coding_genes <- sort(unique(
    gene_ids[!is.na(gene_types) &
             gene_types == "protein_coding"]
))

if (length(coding_genes) < 12) {
    stop("Fewer than 12 protein-coding genes are available.")
}

# Six genes will increase in Condition B.
# Six different genes will increase in Condition A.
up_in_B <- coding_genes[1:6]
up_in_A <- coding_genes[7:12]

cat("Genes upregulated in Condition B:", length(up_in_B), "\n")
cat("Genes upregulated in Condition A:", length(up_in_A), "\n")

# -------------------------
# Define expression changes
# -------------------------
# Rows: transcripts
# Columns: conditions
fold_changes <- matrix(
    1,
    nrow = length(tx),
    ncol = 2,
    dimnames = list(names(tx), c("ConditionA", "ConditionB"))
)

fold_changes[gene_ids %in% up_in_B, "ConditionB"] <- 4
fold_changes[gene_ids %in% up_in_A, "ConditionA"] <- 4

# More reads for longer transcripts.
reads_per_transcript <- pmax(
    30,
    round(width(tx) / READ_LENGTH * 20)
)

# Write the filtered transcript FASTA used by Polyester.
SIM_TX <- file.path(SIM_DIR, "transcripts.fa")
writeXStringSet(tx, SIM_TX, format = "fasta")

# -------------------------
# Simulate six paired-end samples
# -------------------------
simulate_experiment(
    fasta = SIM_TX,
    reads_per_transcript = reads_per_transcript,
    num_reps = c(3, 3),
    fold_changes = fold_changes,
    outdir = SIM_DIR,
    paired = TRUE,
    readlen = READ_LENGTH,
    fraglen = FRAGMENT_LENGTH,
    fragsd = FRAGMENT_SD,
    error_rate = 0.001,
    seed = SEED
)

# -------------------------
# Convert FASTA reads to FASTQ
# -------------------------
# Polyester writes FASTA reads. Use a fixed high-quality score
# for each base because FASTA output has no quality information.
fasta_to_fastq_gz <- function(input_file, output_file) {
    reads <- readDNAStringSet(input_file)

    if (length(reads) == 0) {
        stop("No reads generated in ", input_file)
    }

    seqs <- as.character(reads)
    ids <- names(reads)

    if (is.null(ids)) {
        ids <- paste0("read_", seq_along(reads))
    }

    quality <- strrep("I", nchar(seqs))

    lines <- as.vector(rbind(
        paste0("@", ids),
        seqs,
        "+",
        quality
    ))

    con <- gzfile(output_file, open = "wt")
    on.exit(close(con), add = TRUE)
    writeLines(lines, con)

    invisible(length(reads))
}

samples <- c(
    "ConditionA1", "ConditionA2", "ConditionA3",
    "ConditionB1", "ConditionB2", "ConditionB3"
)

conditions <- c(
    rep("ConditionA", 3),
    rep("ConditionB", 3)
)

# Remove prior generated fixtures to prevent stale files.
unlink(list.files(
    OUT_DIR,
    pattern = "^Condition[AB][123]_[12]\\.fastq\\.gz$",
    full.names = TRUE
))

for (i in seq_along(samples)) {
    prefix <- sprintf("sample_%02d", i)

    for (mate in 1:2) {
        fasta_file <- file.path(
            SIM_DIR,
            paste0(prefix, "_", mate, ".fasta")
        )

        if (!file.exists(fasta_file)) {
            stop("Expected Polyester output not found: ", fasta_file)
        }

        fastq_file <- file.path(
            OUT_DIR,
            paste0(samples[i], "_", mate, ".fastq.gz")
        )

        n <- fasta_to_fastq_gz(fasta_file, fastq_file)

        cat(samples[i], "mate", mate, ":", n, "reads\n")
    }
}

# -------------------------
# Write metadata
# -------------------------
metadata <- data.frame(
    sample = samples,
    condition = conditions
)

write.csv(metadata, METADATA, row.names = FALSE)

cat("\nSynthetic dataset created successfully.\n")
cat("FASTQ directory:", OUT_DIR, "\n")
cat("Metadata:", METADATA, "\n")
