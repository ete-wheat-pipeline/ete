#!/usr/bin/env Rscript
#
# _10_build_rpkm_table.R
#
# Combines the per-library counts (counts/<library>.counts.tsv: one row per
# reference sequence, "<Feature ID>\t<reads>") into a single table with one
# RPKM column per library, in the same format as the CLC CSV export used
# previously: a "Feature ID" column plus one "<library> RPKM" column per
# library. This way the scoring script (_10_score_rpkm_table.R) can read it
# directly. (An earlier version read pipeline_out/counts/*.fragment_counts.tsv
# written by bwa_expression_pipeline.sh; neither file is part of this project.
# In the published run the counts were produced by the samtools loop
# documented in the expression-validation section of 0INFO.md.)
#
# RPKM_gene = (fragments_mapped_to_gene * 1e9) / (gene_length_bp * total_fragments_mapped_in_library)
# (standard definition, Mortazavi et al. 2008, the same reported in the CLC
# RNA-Seq reports: "Reads Per Kilobase of exon Model")
#
# Usage:
#   Rscript _10_build_rpkm_table.R <counts_dir> <file.fai> <output.csv>
#
# Example:
#   Rscript _10_build_rpkm_table.R counts ref/Chr1A_cds.fasta.fai expression_all_cds_bwa.csv

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
counts_dir <- if (length(args) >= 1) args[1] else "counts"
fasta_fai  <- if (length(args) >= 2) args[2] else "reference/CDS.fasta.fai"
out_csv    <- if (length(args) >= 3) args[3] else "expression_all_cds_bwa.csv"

if (!file.exists(fasta_fai)) {
  stop(".fai file not found: ", fasta_fai,
       "\nCreate it with: samtools faidx <reference.fasta>")
}

# Lengths of the transcripts/references from the .fai (samtools faidx: name, length, offset, linebases, linewidth)
lengths <- read_tsv(fasta_fai,
                     col_names = c("Feature ID", "length", "offset", "linebases", "linewidth"),
                     col_types = cols_only(`Feature ID` = col_character(),
                                            length = col_double()))

count_files <- list.files(counts_dir, pattern = "\\.counts\\.tsv$", full.names = TRUE)
if (length(count_files) == 0) {
  stop("No *.counts.tsv file found in: ", counts_dir)
}

rpkm_list <- map(count_files, function(f) {
  lib <- str_remove(basename(f), "\\.counts\\.tsv$")
  counts <- read_tsv(f, col_names = c("Feature ID", "reads"), col_types = "cd")
  total_mapped <- sum(counts$reads)

  if (total_mapped == 0) {
    warning("Library '", lib, "': 0 total mapped reads, RPKM will be NA/Inf.")
  }

  lengths %>%
    left_join(counts, by = "Feature ID") %>%
    mutate(reads = replace_na(reads, 0)) %>%
    transmute(`Feature ID`,
              !!paste0(lib, " RPKM") := (reads * 1e9) / (length * total_mapped))
})

expression_all <- reduce(rpkm_list, full_join, by = "Feature ID")
write_csv(expression_all, out_csv)

cat("Written:", out_csv, "\n")
cat("Columns:", paste(colnames(expression_all), collapse = ", "), "\n")
cat("\nNext step: Rscript _10_score_rpkm_table.R <GENE_putative_ete.fasta> ", out_csv,
    " expression_result.tsv\n", sep = "")
