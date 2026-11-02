#!/usr/bin/env Rscript
#
# _10_score_rpkm_table.R
#
# Consolidated version of the expression-scoring script. It takes only two
# input files:
#   1) a FASTA of the target genes (e.g. GENE_putative_ete.fasta) - the
#      headers (after ">") are used as the list of loci to validate
#   2) the expression CSV (e.g. expression_all_cds_bwa.csv), with a
#      "Feature ID" column and one "<library> RPKM" column per library
#
# For every target gene:
#   - among its isoforms (Feature IDs sharing the same locus), take the one
#     with the highest RPKM across all libraries
#   - assign the score: RPKM = 0 -> -1 | 0 < RPKM <= 1 -> 0 | RPKM > 1 -> +1
#
# The number of RPKM columns is detected automatically (not hard-coded),
# so it works with any number of libraries (9 in the published run).
#
# Usage:
#   Rscript _10_score_rpkm_table.R GENE_putative_ete.fasta expression_all_cds_bwa.csv [output.tsv]
#
# Example:
#   Rscript _10_score_rpkm_table.R GENE_putative_ete.fasta expression_all_cds_bwa.csv expression_result.tsv

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript _10_score_rpkm_table.R <target_genes.fasta> <expression.csv> [output.tsv]")
}
fasta_target <- args[1]
expression_csv <- args[2]
out_tsv <- if (length(args) >= 3) args[3] else "expression_result.tsv"

########################################
# 1) List of the target loci from the FASTA headers
########################################
fasta_lines <- readLines(fasta_target)
header_lines <- fasta_lines[startsWith(fasta_lines, ">")]
target_loci <- header_lines %>%
  str_remove("^>") %>%
  str_split(" ", simplify = TRUE) %>%   # in case the header has a description after the first space
  .[, 1] %>%
  unique()

cat("Target genes read from the FASTA:", length(target_loci), "\n")

########################################
# 2) Expression table: extract the RPKM columns (whatever the number of
#    libraries) and compute the locus by removing the isoform suffix
########################################
expression_table <- read_csv(expression_csv, show_col_types = FALSE)

list_of_columns <- colnames(expression_table)[grepl("RPKM", colnames(expression_table))]
if (length(list_of_columns) == 0) {
  stop("No RPKM column found in ", expression_csv,
       " (columns with 'RPKM' in the name are expected).")
}
cat("RPKM columns detected (", length(list_of_columns), "):\n  ",
    paste(list_of_columns, collapse = ", "), "\n", sep = "")

expression_rpkm <- expression_table %>% select(`Feature ID`, all_of(list_of_columns))

expressed_genes <- expression_rpkm %>%
  mutate(locus = sub("\\..*$", "", `Feature ID`)) %>%           # group the isoforms by locus
  mutate(max_expression = pmax(!!!syms(list_of_columns), na.rm = TRUE))  # maximum RPKM across libraries, per isoform

########################################
# 3) Keep only the target genes, with the most expressed isoform of each locus
########################################
expression_target <- expressed_genes %>% filter(locus %in% target_loci)

loci_not_found <- setdiff(target_loci, unique(expression_target$locus))
if (length(loci_not_found) > 0) {
  cat("WARNING:", length(loci_not_found),
      "target genes not found in the expression table (absent from the mapping reference).\n")
  cat("  Examples:", paste(head(loci_not_found, 5), collapse = ", "), "\n")
}

expressed_genes_result <- expression_target %>%
  arrange(desc(max_expression)) %>%
  group_by(locus) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(locus, max_expression) %>%
  mutate(
    score = case_when(
      max_expression == 0 ~ -1,
      max_expression > 0 & max_expression <= 1 ~ 0,
      max_expression > 1 ~ 1
    )
  )

expressed_genes_result

cat("\nScore distribution:\n")
print(table(expressed_genes_result$score))

write_tsv(expressed_genes_result, out_tsv)
cat("\nWritten:", out_tsv, "\n")
