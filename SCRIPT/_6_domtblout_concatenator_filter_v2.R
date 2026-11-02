#!/usr/bin/env Rscript
# =============================================================================
# _6_domtblout_concatenator_filter_v2.R
#
# 1) recursively reads all the domtblout files (hmmscan) in a directory
# 2) merges them into a single domtblout
# 3) for each query name computes the coverage of the hits and filters:
#      - hits of the same query are merged if they overlap or are
#        contiguous (even if they come from different target names); if at
#        least one aa is left uncovered the blocks stay separate
#      - covered_ali_len = sum of the lengths of the merged blocks
#      - percentage_ali_len = covered_ali_len / qlen * 100
#      - different query names of the same gene are NEVER summed
#      - no per-row filter on E-value / i-Evalue
#      - pass: covered_ali_len >= min_len  AND  percentage_ali_len >= min_pct
# 4) filters proteome, CDS, mRNA (optional) and genes, keeping only the
#    sequences that pass
#
# A gene with at least one passing isoform is a putative ETE.
#
# Usage:
#   Rscript _6_domtblout_concatenator_filter_v2.R -d ../directory_hmm_scan \
#       -p ../starting_file/genome_pep.fasta \
#       -c ../starting_file/genome_cds.fasta \
#       -m ../starting_file/genome_mrna.fasta \
#       -g ../starting_file/genome_gene.fasta
#
# Base R only, no package to install.
# =============================================================================

usage <- function() {
  cat("
Usage: Rscript _6_domtblout_concatenator_filter_v2.R -d DIR -p PEP.fasta -c CDS.fasta -g GENE.fasta [-m MRNA.fasta] [options]

Required:
  -d, --dir      directory with the sub-directories of hmmscan results
  -p, --pep      FASTA of the proteome   (header ID.X)
  -c, --cds      FASTA of the CDS        (header ID.X)
  -g, --gene     FASTA of the genes      (header ID, without isoform)

Optional:
  -m, --mrna     FASTA of the mRNAs      (header ID.X, as for pep and CDS)
  -o, --outdir   output directory                    [default: putative_ete_out]
  -x, --pattern  regex on the names of the files to read [default: domtblout]
  --min-len      minimum covered_ali_len (aa)        [default: 25]
  --min-pct      minimum percentage_ali_len (%)      [default: 10]
  -h, --help     show this message

Output (in --outdir):
  merged.domtblout          merged domtblout
  coverage_summary.tsv      coverage of every query name (with pass flag)
  PEP_putative_ete.fasta
  CDS_putative_ete.fasta
  MRNA_putative_ete.fasta   (only if -m is given)
  GENE_putative_ete.fasta
\n")
}

# -------------------------------------------------------------- arguments --
parse_args <- function(argv) {
  opts <- list(dir = NULL, pep = NULL, cds = NULL, mrna = NULL, gene = NULL,
               outdir = "putative_ete_out", pattern = "domtblout",
               min_len = 25, min_pct = 10)
  map <- c("-d" = "dir",     "--dir" = "dir",
           "-p" = "pep",     "--pep" = "pep",
           "-c" = "cds",     "--cds" = "cds",
           "-m" = "mrna",    "--mrna" = "mrna",
           "-g" = "gene",    "--gene" = "gene",
           "-o" = "outdir",  "--outdir" = "outdir",
           "-x" = "pattern", "--pattern" = "pattern",
           "--min-len" = "min_len", "--min-pct" = "min_pct")
  if (length(argv) == 0) { usage(); quit(status = 1) }
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a %in% c("-h", "--help")) { usage(); quit(status = 0) }
    if (!a %in% names(map)) stop("Unknown option: ", a, call. = FALSE)
    if (i == length(argv)) stop("Missing value for ", a, call. = FALSE)
    opts[[map[[a]]]] <- argv[i + 1]
    i <- i + 2
  }
  opts$min_len <- as.numeric(opts$min_len)
  opts$min_pct <- as.numeric(opts$min_pct)
  if (is.na(opts$min_len) || is.na(opts$min_pct))
    stop("--min-len and --min-pct must be numeric", call. = FALSE)
  for (k in c("dir", "pep", "cds", "gene"))
    if (is.null(opts[[k]])) stop("Missing required argument --", k, call. = FALSE)
  if (!dir.exists(opts$dir)) stop("Directory not found: ", opts$dir, call. = FALSE)
  for (k in c("pep", "cds", "gene"))
    if (!file.exists(opts[[k]])) stop("File not found: ", opts[[k]], call. = FALSE)
  # optional mRNA: if given it must exist
  if (!is.null(opts$mrna) && !file.exists(opts$mrna))
    stop("File not found: ", opts$mrna, call. = FALSE)
  opts
}

# ---------------------------------------------- steps 1 + 2: read / merge --
merge_domtblout <- function(dir, pattern, outdir, merged_file) {
  files <- list.files(dir, pattern = pattern, recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  # avoid re-reading previous outputs if outdir is inside the input directory
  out_abs   <- normalizePath(outdir, mustWork = FALSE)
  files_abs <- normalizePath(files, mustWork = FALSE)
  files <- files[!startsWith(files_abs, paste0(out_abs, .Platform$file.sep))]
  if (length(files) == 0)
    stop("No '", pattern, "' file found in ", dir, call. = FALSE)
  message("[1/4] Found ", length(files), " domtblout files")

  con <- file(merged_file, open = "w")
  on.exit(close(con), add = TRUE)
  header_written <- FALSE
  n_lines <- 0L
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    is_data <- nzchar(lines) & !startsWith(lines, "#")
    if (!header_written) {
      # initial comments (the 3 header lines) of the first file
      first_data <- if (any(is_data)) which(is_data)[1] else length(lines) + 1L
      hdr <- lines[seq_len(first_data - 1L)]
      hdr <- hdr[startsWith(hdr, "#")]
      writeLines(hdr, con)
      header_written <- TRUE
    }
    if (any(is_data)) {
      writeLines(lines[is_data], con)
      n_lines <- n_lines + sum(is_data)
    } else {
      message("      (no hit in ", f, ")")
    }
  }
  writeLines(paste0("# Merged from ", length(files), " domtblout files"), con)
  message("[2/4] Merge completed: ", n_lines, " hits -> ", merged_file)
  invisible(n_lines)
}

# ------------------------------------------------------- read the domtblout --
read_domtblout <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- sub("\r+$", "", lines)
  lines <- lines[nzchar(trimws(lines)) & !startsWith(lines, "#")]
  if (length(lines) == 0) stop("The merged domtblout contains no hits", call. = FALSE)

  # --- FIX: discard truncated/malformed lines BEFORE scan() -----------------
  # a valid line has at least 22 fields separated by spaces/tabs
  ok <- grepl("^(?:\\S+[ \t]+){21}\\S+", lines, perl = TRUE)
  if (any(!ok)) {
    bad <- file.path(dirname(path), "malformed_lines.txt")
    writeLines(lines[!ok], bad)
    message("      WARNING: ", sum(!ok), " malformed lines discarded -> ", bad)
    lines <- lines[ok]
    if (length(lines) == 0) stop("No valid line in the domtblout", call. = FALSE)
  }
  # ---------------------------------------------------------------------------

  # keep only the first 22 fields (the 23rd, the description, contains spaces)
  lines <- sub("^((?:\\S+\\s+){21}\\S+).*$", "\\1", lines, perl = TRUE)
  # columns: 1 target name | 4 query name | 6 qlen | 18 ali from | 19 ali to
  what <- list(character(), NULL, NULL, character(), NULL, integer(),
               NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
               integer(), integer(), NULL, NULL, NULL)
  # FIX: multi.line = FALSE -> scan() does not complete a record with the next line
  x <- scan(text = lines, what = what, quiet = TRUE, quote = "",
            comment.char = "", multi.line = FALSE)
  d <- data.frame(target = x[[1]], query = x[[4]], qlen = x[[6]],
                  ali_from = pmin(x[[18]], x[[19]]),
                  ali_to   = pmax(x[[18]], x[[19]]),
                  stringsAsFactors = FALSE)
  # FIX: no NA may survive parsing (it would end up in pass_q)
  if (anyNA(d)) stop("Missing fields in the domtblout after parsing", call. = FALSE)
  d
}

# ------------------------------------------------ step 3: coverage / filter --
compute_coverage <- function(d, min_len, min_pct) {
  d <- d[order(d$query, d$ali_from, d$ali_to), ]
  n <- nrow(d)

  # maximum end reached up to the previous row, within the same query
  cmax      <- ave(d$ali_to, d$query, FUN = cummax)
  prev_cmax <- c(-Inf, cmax[-n])
  new_query <- c(TRUE, d$query[-1] != d$query[-n])

  # new block if the query changes or at least 1 aa is left uncovered
  # (start > previous_end + 1; start == previous_end + 1 = contiguous -> merged)
  new_block <- new_query | (d$ali_from > prev_cmax + 1)
  block_id  <- cumsum(new_block)

  blocks <- data.frame(
    query = d$query[new_block],
    start = d$ali_from[new_block],
    end   = as.integer(tapply(d$ali_to, block_id, max)),
    stringsAsFactors = FALSE)
  blocks$len <- blocks$end - blocks$start + 1L

  cov <- data.frame(query = unique(d$query), stringsAsFactors = FALSE)
  cov$qlen    <- d$qlen[match(cov$query, d$query)]
  cov$n_hits  <- as.integer(table(d$query)[cov$query])
  cov$n_blocks <- as.integer(table(blocks$query)[cov$query])
  cov$covered_ali_len <- as.integer(rowsum(blocks$len, blocks$query)[cov$query, 1])
  cov$percentage_ali_len <- cov$covered_ali_len / cov$qlen * 100
  cov$targets <- vapply(split(d$target, d$query)[cov$query],
                        function(t) paste(unique(t), collapse = ","), character(1))
  cov$pass <- cov$covered_ali_len >= min_len & cov$percentage_ali_len >= min_pct
  cov
}

# ------------------------------------------------------------------ FASTA ----
fasta_ids <- function(path) {
  lines <- readLines(path, warn = FALSE)
  sub("^>\\s*(\\S+).*$", "\\1", lines[startsWith(lines, ">")])
}

filter_fasta <- function(infile, outfile, keep_ids, label) {
  lines  <- readLines(infile, warn = FALSE)
  is_hdr <- startsWith(lines, ">")
  grp    <- cumsum(is_hdr)                       # record index for every line
  ids    <- sub("^>\\s*(\\S+).*$", "\\1", lines[is_hdr])
  keep_rec  <- ids %in% keep_ids
  keep_line <- grp > 0 & keep_rec[pmax(grp, 1L)]
  writeLines(lines[keep_line], outfile)
  n_out <- sum(keep_rec)
  missing <- setdiff(keep_ids, ids)
  message(sprintf("      %-5s %d/%d records written -> %s", label, n_out, length(ids), outfile))
  if (length(missing) > 0)
    warning(label, ": ", length(missing), " IDs not found in the FASTA (e.g. ",
            paste(head(missing, 3), collapse = ", "), ")", call. = FALSE)
  invisible(n_out)
}

# ------------------------------------------------------------------- main ----
main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(opts$outdir, showWarnings = FALSE, recursive = TRUE)
  merged_file <- file.path(opts$outdir, "merged.domtblout")

  # steps 1-2
  merge_domtblout(opts$dir, opts$pattern, opts$outdir, merged_file)

  # step 3
  d   <- read_domtblout(merged_file)
  cov <- compute_coverage(d, opts$min_len, opts$min_pct)
  cov_out <- cov[order(cov$query), ]
  cov_out$percentage_ali_len <- round(cov_out$percentage_ali_len, 2)
  write.table(cov_out, file.path(opts$outdir, "coverage_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  pass_q <- cov$query[cov$pass]
  message(sprintf("[3/4] %d query names in total, %d pass the filter (covered >= %g aa AND pct >= %g%%)",
                  nrow(cov), length(pass_q), opts$min_len, opts$min_pct))
  if (length(pass_q) == 0) stop("No query passes the filter", call. = FALSE)

  # step 4
  message("[4/4] Filtering the FASTA files")
  filter_fasta(opts$pep, file.path(opts$outdir, "PEP_putative_ete.fasta"), pass_q, "PEP")
  filter_fasta(opts$cds, file.path(opts$outdir, "CDS_putative_ete.fasta"), pass_q, "CDS")

  # mRNA (optional): same headers as pep and CDS -> same IDs (query names)
  if (!is.null(opts$mrna))
    filter_fasta(opts$mrna, file.path(opts$outdir, "MRNA_putative_ete.fasta"), pass_q, "mRNA")

  # genes: ID without isoform -> remove the last ".X" from the query name
  # (if the whole ID is already a gene ID and the shortened one is not, keep the whole one)
  gene_ids   <- fasta_ids(opts$gene)
  stripped   <- sub("\\.[^.]+$", "", pass_q)
  gene_query <- ifelse(stripped %in% gene_ids | !(pass_q %in% gene_ids), stripped, pass_q)
  gene_query <- unique(gene_query)
  filter_fasta(opts$gene, file.path(opts$outdir, "GENE_putative_ete.fasta"), gene_query, "GENE")

  message("Done. Output in: ", normalizePath(opts$outdir))
}

main()
