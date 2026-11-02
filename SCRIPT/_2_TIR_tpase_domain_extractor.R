#!/usr/bin/env Rscript
# =============================================================================
# _2_TIR_tpase_domain_extractor.R
#
# Combines two outputs of the DANTE pipeline:
#   - DANTE_TIR.gff3  -> TIR elements (Class_II) with boundaries, TIRs, TSDs
#   - DANTE.gff3       -> protein domains (e.g. TPase) already translated from
#                         the genome and aligned to REXdb, with quality
#                         metrics (Identity, Relat_Length, Relat_Interruptions)
#
# The genome is not needed: the amino-acid sequence is already present in
# DANTE.gff3 (Region_Seq field), correctly translated by DANTE itself
# (with frameshift correction), with identity and the number of
# interruptions (stop codons/frameshifts) already computed -> the "few stop
# codons" filter becomes a direct filter on Relat_Interruptions / the
# count of '*'.
#
# Region_Seq vs Query_Seq: Query_Seq is the portion aligned to the REXdb
# profile (DB_Seq) used to compute identity/similarity, and it may contain
# pure alignment gaps '-' (not real amino acids) to match the positions of
# DB_Seq. Region_Seq is instead the real protein sequence, without gaps, of
# the whole genomic region annotated by the feature: it is the correct one
# to use as the final amino-acid sequence.
#
# Strategy:
#   1. Filter the "good" TIR elements (TIRs present, clean, not degenerate)
#      from DANTE_TIR.gff3 (no genome needed).
#   2. Filter the "good" domains from DANTE.gff3 (identity, coverage,
#      very few interruptions).
#   3. Associate each TIR element with the domains CONTAINED in its genomic
#      coordinates (configurable minimum overlap fraction) AND whose family
#      classification matches.
#   4. For each element with >= 1 associated domain, keep the best domain
#      (highest Identity) and clean/extract its amino-acid sequence;
#      each domain can be assigned to ONE element only.
#   5. Write one FASTA per lineage/classification (e.g. one for
#      PIF_Harbinger, one for EnSpm_CACTA, ...), with header
#      ID_chr_start_end_classification.
#
# NOTES on the fixes relative to the previous version:
#   [FIX 1] Missing classifications (NA) in the domains were propagated into
#           logical indexing, and in R subsetting a data.frame with a logical
#           vector containing NA does NOT drop the row but inserts an all-NA
#           row: sequences written as "NA" could end up in the FASTA.
#           Final_Classification is now required to be non-NA in the filter
#           and all selections use which() (NA-safe).
#   [FIX 2] The same domain could be the best hit of several overlapping/
#           nested TIR elements, duplicating the same protein under
#           different IDs (downstream pseudo-replication). Each domain is now
#           assigned to one element only, resolving conflicts by decreasing
#           overlap fraction and, on ties, by the shorter element (the one
#           that fits the domain most tightly).
#   [FIX 3] The association required an overlap of just 1 bp, so a domain
#           could be attached to the wrong element through an end-to-end
#           contact. A minimum fraction of the DOMAIN length must now fall
#           inside the element (min_domain_overlap_frac; 1.0 = strict
#           containment).
#
# Usage:
#   Rscript _2_TIR_tpase_domain_extractor.R
# (edit the parameters in CONFIG, or pass command-line arguments, see the end)
# =============================================================================

## =============================================================================
## CONFIG - edit these parameters
## =============================================================================

CONFIG <- list(
  tir_gff     = "DANTE_TIR_final.gff3",      # TIR elements (DANTE_TIR)
  domain_gff  = "DANTE.gff3",                # protein domains (DANTE)
  out_dir     = ".",                         # output directory (default: current dir); one FASTA per lineage
  out_prefix  = "TIR_domain_proteins",       # file-name prefix, e.g. TIR_domain_proteins_Class_II_..._PIF_Harbinger.fasta
  out_table   = "TIR_domain_elements.tsv",   # single summary table (all lineages together), in the current dir

  # --- Filter 1: quality of the annotated TIRs (from DANTE_TIR.gff3 only) ---
  tir_min_len    = 10,     # minimum length (bp) of tir_seq5/tir_seq3
  tir_max_n_frac = 0.10,   # maximum fraction of N tolerated inside the TIRs

  # --- Filter 2 (optional): minimum element length, from coordinates only ---
  elem_min_len = 300,

  # --- Filter 3: quality of the protein domain (from DANTE.gff3) ---
  domain_names          = c("TPase"),  # accepted domain names (TIR/Class_II -> TPase)
  min_identity           = 0.50,       # minimum identity vs REXdb (DANTE uses 0.35 by default)
  min_relat_length       = 0.90,       # minimum domain coverage (close to full length)
  max_relat_interruptions = 0,         # maximum number of interruptions (stop codons/frameshifts)
  max_stop_in_seq         = 0,         # maximum direct count of '*' in the kept sequence

  # --- How to associate domain <-> TIR element ---
  require_same_strand = FALSE,  # if TRUE, the domain must be on the same strand as the element

  # [FIX 3] Minimum fraction of the DOMAIN length that must fall inside the
  # TIR element for the association to be accepted.
  #   1.0  = strict containment (the domain must lie entirely inside)
  #   0.9  = tolerates a small overhang due to slightly different boundaries
  #          between DANTE and DANTE_TIR (they are two different tools)
  # Low values (e.g. 0.1) reproduce the old permissive behaviour.
  min_domain_overlap_frac = 0.90
)

## =============================================================================
## FUNCTIONS
## =============================================================================

# Robust parsing of column 9 (attributes) of a GFF3
parse_gff_attributes <- function(attr_string) {
  fields <- strsplit(attr_string, ";", fixed = TRUE)[[1]]
  fields <- fields[nzchar(fields)]
  keys <- sub("=.*$", "", fields)
  vals <- sub("^[^=]*=", "", fields)
  setNames(as.list(vals), keys)
}

read_gff_generic <- function(gff_file) {
  lines <- readLines(gff_file)
  lines <- lines[!startsWith(lines, "#")]
  lines <- lines[nzchar(lines)]
  split_lines <- strsplit(lines, "\t")
  n_ok <- vapply(split_lines, length, integer(1)) == 9
  if (any(!n_ok)) {
    warning(sum(!n_ok), " lines discarded in ", gff_file, " (not 9 columns).")
    split_lines <- split_lines[n_ok]
  }
  list(
    chr    = vapply(split_lines, `[`, character(1), 1),
    type   = vapply(split_lines, `[`, character(1), 3),
    start  = as.integer(vapply(split_lines, `[`, character(1), 4)),
    end    = as.integer(vapply(split_lines, `[`, character(1), 5)),
    strand = vapply(split_lines, `[`, character(1), 7),
    attrs  = lapply(vapply(split_lines, `[`, character(1), 9), parse_gff_attributes)
  )
}

get_attr_vec <- function(attrs, name) {
  vapply(attrs, function(a) { v <- a[[name]]; if (is.null(v)) NA_character_ else v }, character(1))
}

# --- read DANTE_TIR.gff3 -> data.frame of TIR elements ---
read_dante_tir <- function(gff_file) {
  g <- read_gff_generic(gff_file)
  data.frame(
    chr            = g$chr,
    type           = g$type,
    start          = g$start,
    end            = g$end,
    strand         = g$strand,
    ID             = get_attr_vec(g$attrs, "ID"),
    Classification = get_attr_vec(g$attrs, "Classification"),
    tir_seq5       = get_attr_vec(g$attrs, "tir_seq5"),
    tir_seq3       = get_attr_vec(g$attrs, "tir_seq3"),
    tsd            = get_attr_vec(g$attrs, "tsd"),
    stringsAsFactors = FALSE
  )
}

# --- read DANTE.gff3 -> data.frame of protein domains ---
read_dante_domains <- function(gff_file) {
  g <- read_gff_generic(gff_file)
  data.frame(
    chr                 = g$chr,
    type                = g$type,
    start               = g$start,
    end                 = g$end,
    strand              = g$strand,
    Name                = get_attr_vec(g$attrs, "Name"),
    Final_Classification = get_attr_vec(g$attrs, "Final_Classification"),
    Region_Seq          = get_attr_vec(g$attrs, "Region_Seq"),
    Identity            = as.numeric(get_attr_vec(g$attrs, "Identity")),
    Similarity          = as.numeric(get_attr_vec(g$attrs, "Similarity")),
    Relat_Length        = as.numeric(get_attr_vec(g$attrs, "Relat_Length")),
    Relat_Interruptions = as.numeric(get_attr_vec(g$attrs, "Relat_Interruptions")),
    stringsAsFactors = FALSE
  )
}

n_fraction <- function(seq_chr) {
  seq_chr <- toupper(seq_chr)
  n_count <- lengths(regmatches(seq_chr, gregexpr("N", seq_chr)))
  n_count / nchar(seq_chr)
}

# Extracts the final "family name" from the TIR classification (underscore-separated)
# e.g. "Class_II_Subclass_1_TIR_PIF_Harbinger" -> "PIF_Harbinger"
family_from_tir_classification <- function(x) sub(".*_TIR_", "", x)

# Extracts the final "family name" from the domain classification (pipe-separated)
# e.g. "Class_II|Subclass_1|TIR|PIF/Harbinger" -> "PIF/Harbinger"
family_from_domain_classification <- function(x) sub(".*\\|TIR\\|", "", x)

# Normalizes for comparison (removes everything non-alphanumeric, lower case)
normalize_family <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

# Cleans a DANTE Region_Seq for use as the final protein.
# Region_Seq (unlike Query_Seq) NEVER contains alignment gaps '-' (those in
# Query_Seq are padding for column correspondence with DB_Seq, not real
# amino acids): here we still remove any '-' for safety, remove the
# frameshift-join separator '/', and normalize the case (lower case in DANTE
# marks residues in the extended region around the core alignment, but they
# are still translated from the genome).
clean_region_seq <- function(x) {
  x <- gsub("-", "", x, fixed = TRUE)
  x <- gsub("/", "", x, fixed = TRUE)
  toupper(x)
}

count_stops <- function(x) lengths(regmatches(x, gregexpr("\\*", x, fixed = FALSE)))

## =============================================================================
## MAIN PIPELINE
## =============================================================================

run_pipeline <- function(cfg) {

  ## ---- 1. TIR elements: filter on present/clean TIRs + length ----
  message(">> Reading TIR elements: ", cfg$tir_gff)
  tir <- read_dante_tir(cfg$tir_gff)
  message("   Total TIR elements: ", nrow(tir))

  # [FIX 1] which() instead of logical indexing: NA-safe.
  keep_tir <- which(
    !is.na(tir$tir_seq5) & !is.na(tir$tir_seq3) &
      !is.na(tir$Classification) & !is.na(tir$ID) &
      nchar(tir$tir_seq5) >= cfg$tir_min_len &
      nchar(tir$tir_seq3) >= cfg$tir_min_len &
      n_fraction(tir$tir_seq5) <= cfg$tir_max_n_frac &
      n_fraction(tir$tir_seq3) <= cfg$tir_max_n_frac &
      (tir$end - tir$start + 1) >= cfg$elem_min_len
  )

  tir_f <- tir[keep_tir, , drop = FALSE]
  rownames(tir_f) <- NULL
  tir_f$family   <- normalize_family(family_from_tir_classification(tir_f$Classification))
  tir_f$elem_len <- tir_f$end - tir_f$start + 1
  message("   Passing the TIR + length filter: ", nrow(tir_f), " / ", nrow(tir))

  if (nrow(tir_f) == 0) stop("No TIR element left after filtering. Check the parameters.")

  ## ---- 2. Protein domains: filter on identity/coverage/interruptions ----
  message(">> Reading protein domains: ", cfg$domain_gff)
  dom <- read_dante_domains(cfg$domain_gff)
  message("   Total domains: ", nrow(dom))

  # [FIX 1] Final_Classification required to be non-NA: without this check
  # the family was NA and the comparison "== e$family" returned NA,
  # generating all-NA ghost rows in the subsetting further down.
  keep_dom <- which(
    dom$Name %in% cfg$domain_names &
      !is.na(dom$Final_Classification) &
      !is.na(dom$Region_Seq) &
      !is.na(dom$Identity)            & dom$Identity            >= cfg$min_identity &
      !is.na(dom$Relat_Length)        & dom$Relat_Length        >= cfg$min_relat_length &
      !is.na(dom$Relat_Interruptions) & dom$Relat_Interruptions <= cfg$max_relat_interruptions
  )

  dom_f <- dom[keep_dom, , drop = FALSE]
  rownames(dom_f) <- NULL
  dom_f$family <- normalize_family(family_from_domain_classification(dom_f$Final_Classification))

  # additional direct check on the number of '*' in the cleaned sequence
  dom_f$protein <- clean_region_seq(dom_f$Region_Seq)
  dom_f$n_stops <- count_stops(dom_f$protein)
  dom_f <- dom_f[which(dom_f$n_stops <= cfg$max_stop_in_seq), , drop = FALSE]
  rownames(dom_f) <- NULL

  # [FIX 2] unique domain identifier, needed to guarantee that a domain is
  # not assigned to more than one TIR element.
  dom_f$dom_uid  <- seq_len(nrow(dom_f))
  dom_f$dom_len  <- dom_f$end - dom_f$start + 1

  message("   Passing the domain filter (Name in {", paste(cfg$domain_names, collapse = ","),
          "}, Identity>=", cfg$min_identity, ", Relat_Length>=", cfg$min_relat_length,
          ", Relat_Interruptions<=", cfg$max_relat_interruptions,
          ", direct stops<=", cfg$max_stop_in_seq, "): ", nrow(dom_f), " / ", nrow(dom))

  if (nrow(dom_f) == 0) stop("No domain left after filtering. Consider relaxing the thresholds.")

  ## ---- 3. Association TIR element <-> domain ----
  ## [FIX 3] a 1-bp overlap is no longer enough: at least
  ## min_domain_overlap_frac of the domain length must fall inside the element.
  message(">> Associating TIR elements <-> domains (minimum domain overlap: ",
          cfg$min_domain_overlap_frac, ")...")

  assign_list <- vector("list", nrow(tir_f))
  for (i in seq_len(nrow(tir_f))) {
    e <- tir_f[i, ]

    # [FIX 1] which() makes the selection immune to residual NAs
    idx <- which(
      dom_f$chr == e$chr &
        !is.na(dom_f$family) & dom_f$family == e$family &
        dom_f$start <= e$end & dom_f$end >= e$start
    )
    if (length(idx) == 0) next

    if (cfg$require_same_strand) {
      idx <- idx[dom_f$strand[idx] == e$strand]
      if (length(idx) == 0) next
    }

    # [FIX 3] fraction of the domain length contained in the element
    ovl <- pmin(dom_f$end[idx], e$end) - pmax(dom_f$start[idx], e$start) + 1
    ovl_frac <- ovl / dom_f$dom_len[idx]

    ok <- which(ovl_frac >= cfg$min_domain_overlap_frac)
    if (length(ok) == 0) next
    idx      <- idx[ok]
    ovl_frac <- ovl_frac[ok]

    # keep the domain with the highest Identity (ties -> highest Relat_Length)
    ord <- order(-dom_f$Identity[idx], -dom_f$Relat_Length[idx])
    idx      <- idx[ord]
    ovl_frac <- ovl_frac[ord]

    assign_list[[i]] <- data.frame(
      elem_i   = i,
      dom_uid  = dom_f$dom_uid[idx[1]],
      ovl_frac = ovl_frac[1],
      elem_len = e$elem_len,
      stringsAsFactors = FALSE
    )
  }

  assign_df <- do.call(rbind, assign_list)
  if (is.null(assign_df) || nrow(assign_df) == 0) {
    stop("No TIR element has an associated domain. Check the classification/coordinate correspondence, or lower min_domain_overlap_frac.")
  }
  message("   TIR elements with at least one associated domain: ", nrow(assign_df), " / ", nrow(tir_f))

  ## [FIX 2] one domain -> one element only.
  ## Conflicts resolved by: decreasing overlap fraction, then shorter element
  ## (fits the domain more tightly), then index, for determinism.
  n_before <- nrow(assign_df)
  assign_df <- assign_df[order(assign_df$dom_uid,
                               -assign_df$ovl_frac,
                               assign_df$elem_len,
                               assign_df$elem_i), , drop = FALSE]
  assign_df <- assign_df[!duplicated(assign_df$dom_uid), , drop = FALSE]
  n_dropped <- n_before - nrow(assign_df)
  if (n_dropped > 0) {
    message("   [dedup] ", n_dropped, " assignments removed: the same domain was ",
            "the best hit of several overlapping/nested TIR elements.")
  }

  # restore the per-element order
  assign_df <- assign_df[order(assign_df$elem_i), , drop = FALSE]

  final  <- tir_f[assign_df$elem_i, , drop = FALSE]
  hit_df <- dom_f[match(assign_df$dom_uid, dom_f$dom_uid), , drop = FALSE]
  rownames(final) <- NULL
  rownames(hit_df) <- NULL

  final$domain_name           <- hit_df$Name
  final$domain_classification <- hit_df$Final_Classification
  final$domain_identity       <- hit_df$Identity
  final$domain_relat_length   <- hit_df$Relat_Length
  final$domain_interruptions  <- hit_df$Relat_Interruptions
  final$domain_start          <- hit_df$start
  final$domain_end            <- hit_df$end
  final$domain_overlap_frac   <- round(assign_df$ovl_frac, 4)
  final$protein               <- hit_df$protein
  final$protein_len_aa        <- nchar(hit_df$protein)

  message("   Final TIR elements (after domain dedup): ", nrow(final))

  ## ---- 4. Output: underscore between header fields; one FASTA per lineage ----
  headers <- sprintf("%s_%s_%d_%d_%s",
                     final$ID, final$chr, final$start, final$end, final$domain_classification)

  if (!dir.exists(cfg$out_dir)) dir.create(cfg$out_dir, recursive = TRUE)

  # grouping/file names use the TIR Classification (already underscore-safe)
  lineage_safe <- gsub("[^A-Za-z0-9]+", "_", final$Classification)
  lineages <- sort(unique(lineage_safe))

  message(">> Writing ", length(lineages), " FASTA files (one per lineage) to: ", cfg$out_dir)
  for (lin in lineages) {
    idx <- which(lineage_safe == lin)
    fasta_lines <- character(2 * length(idx))
    fasta_lines[c(TRUE, FALSE)] <- paste0(">", headers[idx])
    fasta_lines[c(FALSE, TRUE)] <- final$protein[idx]
    out_file <- file.path(cfg$out_dir, paste0(cfg$out_prefix, "_", lin, ".fasta"))
    writeLines(fasta_lines, out_file)
    message("   - ", out_file, " (", length(idx), " sequences)")
  }

  out_tab <- final[, c("ID", "chr", "start", "end", "strand", "Classification",
                       "domain_name", "domain_classification", "domain_start", "domain_end",
                       "domain_identity", "domain_relat_length", "domain_interruptions",
                       "domain_overlap_frac", "protein_len_aa")]
  write.table(out_tab, cfg$out_table, sep = "\t", row.names = FALSE, quote = FALSE)
  message(">> Summary table (all lineages): ", cfg$out_table)

  message("\n== DONE == ", nrow(final), " TIR elements with a high-quality transposase domain, out of ",
          nrow(tir), " TIR elements in the original GFF.")
  invisible(final)
}

## =============================================================================
## EXECUTION (with optional command-line override)
##   Rscript _2_TIR_tpase_domain_extractor.R DANTE_TIR.gff3 DANTE.gff3 out_dir out.tsv out_prefix
## =============================================================================
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) CONFIG$tir_gff    <- args[1]
  if (length(args) >= 2) CONFIG$domain_gff <- args[2]
  if (length(args) >= 3) CONFIG$out_dir    <- args[3]
  if (length(args) >= 4) CONFIG$out_table  <- args[4]
  if (length(args) >= 5) CONFIG$out_prefix <- args[5]

  run_pipeline(CONFIG)
}
