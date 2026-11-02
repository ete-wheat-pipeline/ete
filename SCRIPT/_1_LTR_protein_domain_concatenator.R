#!/usr/bin/env Rscript
# Extracts the concatenated amino-acid sequences of the protein domains of
# each LTR element from a DANTE_LTR GFF3 file, writing one FASTA per
# Final_Classification.
#
# Element selection (aligned with 1_dante_domains2fasta.py):
#   1. drop elements with Rank exactly "D" (keep DL, DLT, DLP, DLTP, ...)
#   2. drop elements with at least one domain with Relat_Interruptions > 0
#   3. keep only elements whose Ndomains equals the maximum observed for
#      their lineage, computed on the elements surviving steps 1-2
#   4. if a domain occurs more than once in the same element (nested
#      elements / duplications), keep only the best hit (highest score)
#
# Concatenation order of the domains: the GFF "domain_order" field (the
# DANTE_LTR source of truth for the 5'->3' order of the protein). Fallback
# by genomic coordinates (increasing on +, decreasing on -) ONLY for the
# elements in which domain_order is missing, with an on-screen warning.
#
# Sequence field: Region_Seq (not Query_Seq).
#
# FASTA header:
#   >ElementID_start_end_Final_Classification
# where start/end are those of the "transposable_element" feature (the
# whole element), not those of the individual domains. The chromosome is
# not included (it is already part of the ID) and the strand is omitted.
#
# Usage: ./_1_LTR_protein_domain_concatenator.R input.gff3 [output_dir]

suppressMessages({
  library(tidyverse)
  library(rtracklayer)
})

args   <- commandArgs(TRUE)
input  <- args[1]
outdir <- ifelse(length(args) >= 2, args[2], ".")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

stopifnot(!is.na(input), file.exists(input))

ltr <- import(input)
ltr <- as_tibble(ltr)

# GFF attributes with variable cardinality (Parent, etc.) are imported as
# CharacterList. They must be converted to character IMMEDIATELY: doing it
# after match()/filter() have already been used silently gives wrong results.
list_cols <- c("ID", "Parent", "Name", "Final_Classification", "Rank",
               "Region_Seq", "Relat_Interruptions", "Ndomains", "domain_order")
for (col in intersect(list_cols, colnames(ltr))) {
  if (is(ltr[[col]], "List")) {
    ltr[[col]] <- vapply(ltr[[col]], function(x) paste(as.character(x), collapse = ","),
                         character(1))
  } else {
    ltr[[col]] <- as.character(ltr[[col]])
  }
}

# --- filter 1: drop Rank exactly "D" (keep DL, DLT, DLP, DLTP, ...) --------
ltr_TE <- ltr %>%
  filter(type == "transposable_element", Rank != "D") %>%
  mutate(Ndomains = as.numeric(Ndomains))

# Protein domains that belong ONLY to these elements
ltr_domain <- ltr %>%
  filter(type == "protein_domain", Parent %in% ltr_TE$ID) %>%
  mutate(score = suppressWarnings(as.numeric(score)),
         score = ifelse(is.na(score), 0, score),
         Relat_Interruptions = suppressWarnings(as.numeric(Relat_Interruptions)))

df <- ltr_domain %>%
  mutate(domain_order = suppressWarnings(as.numeric(domain_order))) %>%
  select(seqnames, start, end, strand, score,
         Final_Classification, Name, Parent, Region_Seq, Relat_Interruptions,
         domain_order)

# --- filter 2: drop elements with at least one interrupted domain ----------
# (equivalent to "any(Relat_Interruptions > 0)" in the Python version; with
# non-negative values it is the same as the "mean == 0" used previously)
bad_parents <- df %>%
  filter(Relat_Interruptions > 0) %>%
  distinct(Parent) %>%
  pull(Parent)

df <- df %>% filter(!Parent %in% bad_parents)
ltr_TE <- ltr_TE %>% filter(!ID %in% bad_parents)

# --- deduplication: if a domain occurs more than once in the same element, -
# keep only the best hit (highest score), as in the Python version
df <- df %>%
  group_by(Parent, Name) %>%
  slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
  ungroup()

# --- filter 3: maximum Ndomains per lineage, computed on the survivors -----
# (elements left after filters 1-2, i.e. those still present in df)
survivors <- ltr_TE %>% filter(ID %in% unique(df$Parent))

max_nd <- survivors %>%
  group_by(Final_Classification) %>%
  summarize(max_ndomains = max(Ndomains, na.rm = TRUE), .groups = "drop")

keep_ids <- survivors %>%
  left_join(max_nd, by = "Final_Classification") %>%
  filter(Ndomains == max_ndomains) %>%
  pull(ID)

df <- df %>% filter(Parent %in% keep_ids)

# --- coordinates of the WHOLE element (transposable_element feature) -------
# Used for the FASTA header: not the start/end of the individual domains but
# the interval of the element the domains belong to (Parent == ID).
te_coords <- ltr_TE %>%
  filter(ID %in% keep_ids) %>%
  transmute(Parent   = ID,
            te_start = start,
            te_end   = end)

# --- concatenation order: domain_order (DANTE source of truth) -------------
# DANTE_LTR already records the 5'->3' order of the protein in
# "domain_order": using it avoids assuming that the spatial order on the
# genome (coordinates) always matches the order in the protein, an
# assumption that can fail for duplicated/nested domains or non-canonical
# topologies.
# If domain_order is missing for AT LEAST one domain of an element, the
# coordinate fallback is used ONLY for that element (increasing on +,
# decreasing on -), and a warning is printed so that it can be checked.

parents_missing_order <- df %>%
  filter(is.na(domain_order)) %>%
  distinct(Parent) %>%
  pull(Parent)

if (length(parents_missing_order) > 0) {
  warning("domain_order missing for ", length(parents_missing_order),
          " element(s): the coordinate fallback is used for these only. ",
          "IDs: ", paste(head(parents_missing_order, 10), collapse = ", "),
          if (length(parents_missing_order) > 10) ", ..." else "")
}

# order_key: value used to sort row by row. For Parents with a complete
# domain_order, domain_order is used; for the others (fallback) start is
# used, with the sign inverted on the minus strand - so a single ascending
# arrange on order_key covers all cases.
df <- df %>%
  mutate(order_key = ifelse(
    Parent %in% parents_missing_order,
    ifelse(strand == "-", -start, start),
    domain_order
  ))

concat_domains <- function(data) {
  data %>%
    arrange(Parent, order_key) %>%
    group_by(Parent) %>%
    summarize(Final_Classification = dplyr::first(Final_Classification),
              concat = paste(Region_Seq, collapse = ""),
              .groups = "drop")
}

df_concat <- concat_domains(df) %>%
  left_join(te_coords, by = "Parent") %>%
  mutate(header = paste(Parent, te_start, te_end, Final_Classification, sep = "_"))

# one FASTA per Final_Classification, generated dynamically
# (file-name prefix "concatenated_domains_")
df_concat %>%
  group_by(Final_Classification) %>%
  group_walk(~ {
    safe_name  <- gsub("[^A-Za-z0-9]+", "_", .y$Final_Classification)
    fasta_path <- file.path(outdir, paste0("concatenated_domains_", safe_name, ".fasta"))
    write_lines(paste0(">", .x$header, "\n", .x$concat),
                fasta_path)
  })

message("Done. Files written to: ", normalizePath(outdir))
