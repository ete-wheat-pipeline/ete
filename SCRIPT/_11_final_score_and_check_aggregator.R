#!/usr/bin/env Rscript
# =============================================================================
# _11_final_score_and_check_aggregator.R
#
# Merges the four validation scores of the putative ETEs (similarity to
# TEs, siRNA, repetitiveness, expression) with the check information,
# derived from the files already produced in the previous steps:
#   - HMM coverage per isoform: 6_SCANNING_RESULT/coverage_summary.tsv
#   - isoform filter: isoforms of the gene passing the coverage filter out
#     of the total isoforms (headers of fine_cleaned_pep.fasta)
#   - LTR/TIR lineages hit and ambiguity: "targets" column of
#     coverage_summary.tsv
# Replaces filtro_aa_perc.py, filtro_isoforme.sh and
# filtro_ambiguous_e_hit_check.R (earlier scripts, not part of this project).
#
# Usage (from the 11_RESULTS directory):
#   Rscript ../SCRIPT/_11_final_score_and_check_aggregator.R ..
#
#   1st argument: project directory                  [default: ..]
#   2nd argument: output directory                   [default: .]
#
# Output:
#   score_per_locus.tsv                          one row per ETE locus
#   final_validation_result_filterable.tsv       one row per isoform with hits
#   isoform_labels.tsv                           query_name / score_label /
#                                                iso_filter
#
# Stops if a module has no locus in common with the ETE list (empty file or
# file from another run). A locus missing even one of the four scores has
# SCORE_SUM = NA; n_score tells how many scores it has.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
})
# rtracklayer loads S4Vectors, IRanges and other packages that mask some
# dplyr functions (rename, count, first, slice, desc...): dplyr:: is used
# where needed

args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args) >= 1) args[1] else ".."
out  <- if (length(args) >= 2) args[2] else "."

f <- c(
  ete       = file.path(base, "6_SCANNING_RESULT", "all_putative_ETEs_ids.txt"),
  coverage  = file.path(base, "6_SCANNING_RESULT", "coverage_summary.tsv"),
  pep       = file.path(base, "0_STARTING_FILES", "fine_cleaned_pep.fasta"),
  gff       = file.path(base, "0_STARTING_FILES", "fine_cleaned_gff.gff3"),
  sim       = file.path(base, "7_TE_SIMILARITY_VALIDATION", "similarity_score_final.tsv"),
  sirna     = file.path(base, "8_SIRNA_VALIDATION", "sirna_score.tsv"),
  repet     = file.path(base, "9_REPETITIVENESS_VALIDATION", "repetitivity_result.tsv"),
  expr      = file.path(base, "10_EXPRESSION_VALIDATION", "expression_result.tsv")
)

not_found <- f[!file.exists(f)]
if (length(not_found) > 0) {
  stop("Files not found:\n  ", paste(not_found, collapse = "\n  "), call. = FALSE)
}

read_chr <- function(path) {
  read_tsv(path, col_types = cols(.default = col_character()), progress = FALSE)
}
gene_of <- function(id) sub("\\.[^.]+$", "", id)

# --------------------------------------------------------------- ETE list --
ete <- tibble(locus = str_trim(read_lines(f[["ete"]]))) %>%
  filter(locus != "") %>%
  distinct()

# ---------------------------------------------------------- four scores --
# similarity: the documented procedure adds the header with sed '1i'; the
# "locus similarity_score" row is removed if present.
sim <- read_tsv(f[["sim"]], col_names = c("locus", "similarity_score"),
                col_types = "cc", progress = FALSE) %>%
  filter(locus != "locus") %>%
  mutate(similarity_score = as.numeric(similarity_score))

sirna <- read_chr(f[["sirna"]]) %>%
  transmute(locus,
            sirna_score = as.numeric(sirna_score),
            sirna_category = category)

# _9_repetitiveness_score_.R calls the gene column "seqname"
repet <- read_chr(f[["repet"]])
if (!"locus" %in% names(repet) && "seqname" %in% names(repet)) {
  repet <- dplyr::rename(repet, locus = seqname)
}
repet <- repet %>%
  transmute(locus,
            depth = as.numeric(depth),
            repetitivity_score = as.numeric(repetitivity_score))

expr <- read_chr(f[["expr"]]) %>%
  transmute(locus,
            max_expression = as.numeric(max_expression),
            score = as.numeric(score))

modules <- list(
  similarity     = list(data = sim,   file = f[["sim"]],   col = "similarity_score"),
  siRNA          = list(data = sirna, file = f[["sirna"]], col = "sirna_score"),
  repetitiveness = list(data = repet, file = f[["repet"]], col = "repetitivity_score"),
  expression     = list(data = expr,  file = f[["expr"]],  col = "score")
)

cat("ETE loci in the list:", nrow(ete), "\n\n")
cat(sprintf("  %-15s %8s %10s %12s %10s\n",
            "module", "rows", "in common", "no score", "foreign"))
for (m in names(modules)) {
  x <- modules[[m]]$data
  dup <- sum(duplicated(x$locus))
  if (dup > 0) {
    stop("Module ", m, " has ", dup, " repeated loci: ", modules[[m]]$file,
         call. = FALSE)
  }
  common <- sum(ete$locus %in% x$locus)
  if (common == 0) {
    stop("Module ", m, " has no locus in common with the ETE list ",
         "(empty file or from another run?): ", modules[[m]]$file, call. = FALSE)
  }
  no_score <- nrow(ete) - sum(ete$locus %in% x$locus[!is.na(x[[modules[[m]]$col]])])
  cat(sprintf("  %-15s %8d %10d %12d %10d\n",
              m, nrow(x), common, no_score, sum(!x$locus %in% ete$locus)))
}

# ---------------------------------------------------- gene coordinates --
gff <- import(f[["gff"]], format = "gff3", feature.type = "gene")
genes <- tibble(locus  = as.character(mcols(gff)$ID),
                chr    = as.character(seqnames(gff)),
                start  = start(gff),
                end    = end(gff),
                strand = as.character(strand(gff))) %>%
  mutate(gene_len = end - start + 1L) %>%
  distinct(locus, .keep_all = TRUE)

n_no_coord <- sum(!ete$locus %in% genes$locus)
if (n_no_coord > 0) {
  cat("\nWARNING:", n_no_coord,
      "ETE loci not found among the GFF genes (coordinates NA).\n")
}

# ------------------------------------- coverage per isoform (step 6) --
cov <- read_chr(f[["coverage"]]) %>%
  transmute(locus = gene_of(query),
            query,
            qlen = as.integer(qlen),
            n_hit = as.integer(n_hits),
            n_blocks = as.integer(n_blocks),
            covered_ali_len = as.integer(covered_ali_len),
            percentage_ali_len = as.numeric(percentage_ali_len),
            passes_coverage = as.logical(pass),
            targets)

# ------------------------------- isoform filter: passing / total --
pep_id <- read_lines(f[["pep"]], progress = FALSE)
pep_id <- sub("^>\\s*(\\S+).*$", "\\1", pep_id[startsWith(pep_id, ">")])

iso_tot  <- tibble(locus = gene_of(unique(pep_id))) %>%
  dplyr::count(locus, name = "n_tot")
iso_pass <- cov %>%
  filter(passes_coverage) %>%
  distinct(query, .keep_all = TRUE) %>%
  dplyr::count(locus, name = "n_pass")

iso <- iso_tot %>%
  left_join(iso_pass, by = "locus") %>%
  mutate(n_pass = coalesce(n_pass, 0L),
         isoform_filter = if_else(n_pass == n_tot, "yes", "no"),
         isoforms_hit = paste0(n_pass, "/", n_tot),
         numeric_ratio = n_pass / n_tot) %>%
  dplyr::select(locus, isoform_filter, isoforms_hit, numeric_ratio)

# ------------------------------- LTR/TIR lineages and ambiguity --
LTR_lineage <- c("Ale", "Alesia", "Angela", "Bianca", "Bryco", "Gymco I",
                 "Gymco II", "Gymco III", "Gymco IV", "Ikeros", "Ivana", "Lyco",
                 "Osser", "SIRE", "TAR", "Tork", "Athila", "Chlamyvir", "CRM",
                 "Galadriel", "Ogre", "Phygy", "Reina", "Retand", "Selgy",
                 "Tat I", "Tat II", "Tat III", "Tcn1", "Tekay")
TIR_lineage <- c("CACTA", "hAT", "Kolobok", "Merlin", "Mutator", "Novosib", "P",
                 "PIF", "Harbinger", "PiggyBac", "Sola1", "Sola2", "Tc1", "Mariner")

fam <- cov %>%
  dplyr::select(query, targets) %>%
  separate_rows(targets, sep = ",") %>%
  mutate(lineage = sub("_win_.*$", "", targets),
         class   = case_when(lineage %in% LTR_lineage ~ "LTR",
                             lineage %in% TIR_lineage ~ "TIR",
                             TRUE ~ "UNKNOWN")) %>%
  group_by(query) %>%
  summarise(LTR_families     = paste(sort(unique(lineage[class == "LTR"])), collapse = ";"),
            TIR_families     = paste(sort(unique(lineage[class == "TIR"])), collapse = ";"),
            UNKNOWN_families = paste(sort(unique(lineage[class == "UNKNOWN"])), collapse = ";"),
            ambiguous        = if_else(any(class == "LTR") & any(class == "TIR"), "yes", "no"),
            .groups = "drop")

n_unknown <- sum(fam$UNKNOWN_families != "")
if (n_unknown > 0) {
  cat("\nWARNING:", n_unknown,
      "isoforms hit by lineages classified neither as LTR nor as TIR (UNKNOWN_families column).\n")
}

# ------------------------------------------------------ per-locus table --
per_locus <- ete %>%
  left_join(genes, by = "locus") %>%
  left_join(sirna, by = "locus") %>%
  left_join(repet, by = "locus") %>%
  left_join(sim,   by = "locus") %>%
  left_join(expr,  by = "locus") %>%
  mutate(n_score = (!is.na(sirna_score)) + (!is.na(repetitivity_score)) +
                   (!is.na(similarity_score)) + (!is.na(score)),
         SCORE_SUM = sirna_score + repetitivity_score + similarity_score + score) %>%
  left_join(iso, by = "locus")

# ---------------------------------------------------- per-isoform table --
final <- per_locus %>%
  left_join(cov, by = "locus") %>%
  left_join(fam, by = "query") %>%
  dplyr::select(locus, chr, start, end, strand, gene_len, qlen,
                sirna_score, repetitivity_score, similarity_score, score,
                SCORE_SUM, n_score,
                query, n_hit, n_blocks, covered_ali_len, percentage_ali_len,
                passes_coverage, isoform_filter, isoforms_hit, numeric_ratio,
                everything())

isoform_labels <- final %>%
  filter(!is.na(query)) %>%
  transmute(query_name  = query,
            score_label = if_else(is.na(SCORE_SUM), "NA", as.character(SCORE_SUM)),
            iso_filter  = case_when(isoform_filter == "yes" ~ "yes",
                                    isoform_filter == "no" ~ "no",
                                    TRUE ~ "NA")) %>%
  distinct(query_name, .keep_all = TRUE)

# ------------------------------------------------------------------ output --
dir.create(out, showWarnings = FALSE, recursive = TRUE)
write_tsv(per_locus, file.path(out, "score_per_locus.tsv"), na = "NA")
write_tsv(final, file.path(out, "final_validation_result_filterable.tsv"), na = "NA")
write_tsv(isoform_labels, file.path(out, "isoform_labels.tsv"), na = "NA")

cat("\nSCORE_SUM per locus:\n")
print(table(SCORE_SUM = per_locus$SCORE_SUM, useNA = "ifany"))

cl4 <- per_locus %>% filter(SCORE_SUM == 4)
cat("\nLoci with SCORE_SUM = 4                 :", nrow(cl4), "\n")
cat("  of which with isoform_filter = yes    :", sum(cl4$isoform_filter == "yes", na.rm = TRUE), "\n")
cat("Class-4 isoforms passing 25 aa and 10 %:",
    sum(final$SCORE_SUM == 4 & final$passes_coverage, na.rm = TRUE), "\n")
cat("Ambiguous LTR/TIR isoforms:", sum(fam$ambiguous == "yes"), "\n")

cat("\nWritten to", normalizePath(out), ":\n",
    " score_per_locus.tsv\n",
    " final_validation_result_filterable.tsv\n",
    " isoform_labels.tsv\n")
