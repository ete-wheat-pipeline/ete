#!/usr/bin/env Rscript

# _8_repetitiveness_score_.R
#
# Repetitiveness validation of the putative ETEs. Genomic DNA reads mapped
# (and filtered) against (A) all genes and (B) the TE library give a mean
# depth per sequence; the TE depth quartiles set the thresholds:
#   depth <= Q1(TE)        -> repetitivity_score = +1
#   Q1(TE) < depth <= Q3   -> repetitivity_score =  0
#   depth >  Q3(TE)        -> repetitivity_score = -1
#
# Usage:
#   Rscript _8_repetitiveness_score_.R <bam_TE> <bam_genes> <ETE_list>
#
# Example:
#   Rscript _8_repetitiveness_score_.R mapped_TEonly.filtered.bam \
#       mapped_genes.filtered.bam all_putative_ETEs_ids.txt
#
# The result (repetitivity_result.tsv) is saved in the directory from
# which Rscript is run.

# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste(
      "Usage:",
      "Rscript _8_repetitiveness_score_.R <bam_TE> <bam_genes> <ETE_list>",
      sep = "\n"
    ),
    call. = FALSE
  )
}

bam_te_file <- args[1]
bam_genes_file <- args[2]
ete_list_file <- args[3]

# Current directory: the output will be saved here
output_dir <- getwd()
output_file <- file.path(output_dir, "repetitivity_result.tsv")

# ---------------------------------------------------------------------------
# Check that the input files exist
# ---------------------------------------------------------------------------
input_files <- c(bam_te_file, bam_genes_file, ete_list_file)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    paste(
      "The following input files were not found:",
      paste(missing_files, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------------------
library(tidyverse)
library(GenomicAlignments)

# ---------------------------------------------------------------------------
# Reference B: TEs only
# ---------------------------------------------------------------------------
reads_te <- readGAlignments(bam_te_file)

# ---------------------------------------------------------------------------
# Reference A: genes (including ETEs)
# ---------------------------------------------------------------------------
reads_genes <- readGAlignments(bam_genes_file)

ete_tsv <- read_tsv(
  ete_list_file,
  col_names = "query_name"
)

genes_of_interest <- unique(ete_tsv$query_name)

reads_ete <- reads_genes[
  as.character(seqnames(reads_genes)) %in% genes_of_interest
]

# ---------------------------------------------------------------------------
# Coverage metrics
# ---------------------------------------------------------------------------
compute_coverage_metrics <- function(reads_subset) {
  cov_obj <- coverage(reads_subset)

  tibble::tibble(
    seqname       = names(cov_obj),
    length        = sapply(cov_obj, length),
    covered_bases = sapply(cov_obj, function(x) sum(x > 0)),
    total_bases   = sapply(cov_obj, sum),
    coverage      = covered_bases / length,
    depth         = total_bases / length
  )
}

te_cov <- compute_coverage_metrics(reads_te) %>%
  filter(grepl("^(LTR|TIR)_", seqname))

gene_cov <- compute_coverage_metrics(reads_genes)

# The ETEs are taken from gene_cov, so ETEs with depth = 0
# stay in the table.
ete_cov <- gene_cov %>%
  filter(seqname %in% genes_of_interest)

# ---------------------------------------------------------------------------
# Thresholds from the TE quartiles
# ---------------------------------------------------------------------------
q1_te <- quantile(te_cov$depth, 0.25, na.rm = TRUE)
q3_te <- quantile(te_cov$depth, 0.75, na.rm = TRUE)

cat("Q1 TE =", q1_te, "  Q3 TE =", q3_te, "\n")
summary(te_cov$depth)

# ---------------------------------------------------------------------------
# Score assignment
# ---------------------------------------------------------------------------
result <- ete_cov %>%
  mutate(
    repetitivity_score = case_when(
      depth <= q1_te ~  1,
      depth >  q3_te ~ -1,
      TRUE           ~  0
    )
  )

table(result$repetitivity_score)

# ---------------------------------------------------------------------------
# Save the result
# ---------------------------------------------------------------------------
write_tsv(result, output_file)

cat("\nResult saved to:\n", output_file, "\n", sep = "")
