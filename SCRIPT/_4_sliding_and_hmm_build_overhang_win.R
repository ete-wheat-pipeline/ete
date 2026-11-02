#!/usr/bin/env Rscript
#
# _4_sliding_and_hmm_build_overhang_win.R
#
# Variant of pipeline_lineages.R for ALREADY ALIGNED input: it skips the
# MAFFT step and goes straight to sliding windows -> split per window ->
# hmmbuild + hmmpress (HMMER), for every aligned FASTA found in an input
# directory.
#
# LINEAGE NAME
# ------------
# Obtained by removing the "_ali" suffix (if present) and the ".fasta"
# extension. It follows the convention of the files already in use:
#   Ale_ali.fasta      -> Ale
#   Angela_ali.fasta   -> Angela
#   EnSpm_CACTA.TPase_ali.fasta -> EnSpm_CACTA.TPase
# Unlike the previous version, the last field after "_" is NO LONGER used:
# with "Ale_ali.fasta" that rule would have returned "ali" as the lineage
# name, which is wrong.
#
# MODEL AND WINDOW NAMES
# ----------------------
# Each window is saved ONLY inside the directory of its own lineage (no
# directory or file shared between different lineages), and the HMM model
# is named "<lineage>_win_<start>-<end>", not "win_<start>-<end>". Reason:
# different lineages generate the same window labels (e.g. all of them have
# a "1-100" window), so without the prefix the models would end up with
# duplicated names as soon as profiles of several lineages are merged into
# the same HMM database (hmmbuild does not notice this by itself: the
# conflict only appears later, when hmmscan can no longer attribute a hit
# to the correct lineage). With the prefix, model names stay unique
# however they are combined afterwards.
#
# The number of windows is computed dynamically by reading the headers
# actually produced by seqkit sliding.
#
# Output created in the working directory from which the script is run:
#   ./<lineage>/ali.fasta                    <- raw output of seqkit sliding
#   ./<lineage>/<start>-<end>.fasta          <- FASTA of the single window
#   ./<lineage>/<start>-<end>.hmm(.h3f/.h3i/.h3m/.h3p)  <- pressed HMM,
#       with NAME = <lineage>_win_<start>-<end>
#
# NB: here the alignment is given as input, so it is NOT copied/duplicated
# into the parent directory (unlike the pipeline with MAFFT, where the
# alignment is a generated output).
#
# Exit status: 0 if every step succeeded; 1 if any step failed or any input
# file or lineage was skipped (all lineages are still attempted first, and
# the failed steps are listed at the end).
#
# Usage:
#   ./_4_sliding_and_hmm_build_overhang_win.R <input_dir> [n_cores]
#
# If n_cores is not given, 14 is used by default.

suppressMessages({
  library(tidyverse)
  library(Biostrings)
})

args <- commandArgs(TRUE)

if (length(args) < 1) {
  stop("Usage: ./_4_sliding_and_hmm_build_overhang_win.R <input_dir> [n_cores]", call. = FALSE)
}

input_folder <- args[1]
threads <- if (length(args) >= 2) args[2] else "14"

# Sliding-window parameters: used both for the seqkit call and to compute
# the optional final window with the overhang (see step 5 below), so they
# are kept in one place.
window_size <- 100L
step_size <- 20L

if (!dir.exists(input_folder)) {
  stop(paste0("Error: directory '", input_folder, "' does not exist."), call. = FALSE)
}

fasta_files <- list.files(input_folder, pattern = "\\.fasta$", full.names = TRUE)

if (length(fasta_files) == 0) {
  stop(paste0("No .fasta file found in '", input_folder, "'."), call. = FALSE)
}

# Failed or skipped steps, reported at the end (non-zero exit status).
failures <- character(0)
record_failure <- function(what) failures <<- c(failures, what)

# --- helper: runs a system command, reports a failure without stopping
#     the whole pipeline (only the current lineage is skipped) ---
run_cmd <- function(cmd, cmd_args, stdout_file = NULL, step_name = "") {
  status <- system2(
    cmd, cmd_args,
    stdout = if (is.null(stdout_file)) "" else stdout_file,
    stderr = ""
  )
  if (is.null(status) || status != 0) {
    warning(paste0("!!! ", step_name, " failed (exit code ", status, ")"), call. = FALSE)
    record_failure(step_name)
    return(FALSE)
  }
  if (!is.null(stdout_file) &&
      (!file.exists(stdout_file) || file.info(stdout_file)$size == 0)) {
    warning(paste0("!!! ", step_name, " produced an empty file: ", stdout_file), call. = FALSE)
    record_failure(step_name)
    return(FALSE)
  }
  TRUE
}

# --- helper: given a set of sequences already cut for one window (same
#     width for all) and its "start-end" label, writes the window FASTA and
#     runs hmmbuild + hmmpress. Used both for the regular windows and for
#     the optional final window with the overhang (see step 5), so the
#     build logic stays identical and is not duplicated. ---
build_window <- function(seq_set, range_label, lineage_dir, lineage, threads) {
  window_fasta <- file.path(lineage_dir, paste0(range_label, ".fasta"))
  writeXStringSet(seq_set, window_fasta)

  hmm_out <- file.path(lineage_dir, paste0(range_label, ".hmm"))
  model_name <- paste0(lineage, "_win_", range_label)  # e.g. Ale_win_1-100

  ok <- run_cmd(
    "hmmbuild",
    c("-n", model_name, "--cpu", threads, hmm_out, window_fasta),
    step_name = paste0("hmmbuild ", range_label, " (", lineage, ")")
  )
  if (!ok) return(invisible(FALSE))

  run_cmd(
    "hmmpress",
    c("-f", hmm_out),
    step_name = paste0("hmmpress ", range_label, " (", lineage, ")")
  )
  invisible(TRUE)
}

base_dir <- getwd()

message("Found ", length(fasta_files), " aligned FASTA files as input. Threads: ", threads)

# --- lineage names already assigned, to catch collisions BEFORE silently
#     overwriting the directory of a previous lineage (e.g. two input files
#     that reduce to the same name) ---
lineage_seen <- character(0)

for (fasta_file in fasta_files) {

  # --- lineage name: basename without "_ali" (if present) and without
  #     the ".fasta" extension ---
  lineage <- basename(fasta_file)
  lineage <- sub("\\.fasta$", "", lineage)
  lineage <- sub("_ali$", "", lineage)

  if (lineage %in% lineage_seen) {
    warning(paste0("!!! Lineage name '", lineage, "' already used by another input ",
                    "file: '", basename(fasta_file), "' would be written to the SAME ",
                    "directory, overwriting the windows already produced. File skipped."),
            call. = FALSE)
    record_failure(paste0("lineage name collision: ", basename(fasta_file)))
    next
  }
  lineage_seen <- c(lineage_seen, lineage)

  message("\n========== Lineage: ", lineage, " (", basename(fasta_file), ") ==========")

  # quick check: is the FASTA actually aligned? (all sequences must have
  # the same length)
  seqs_check <- tryCatch(readAAStringSet(fasta_file), error = function(e) NULL)
  if (is.null(seqs_check) || length(seqs_check) == 0) {
    warning(paste0("!!! Cannot read ", fasta_file, " - lineage skipped"), call. = FALSE)
    record_failure(paste0("unreadable input: ", basename(fasta_file)))
    next
  }
  if (length(unique(width(seqs_check))) != 1) {
    warning(paste0("!!! WARNING: the sequences in ", basename(fasta_file),
                    " have different lengths: the file does not look aligned. Continuing anyway."),
            call. = FALSE)
  }
  # NB: unlike the previous version, seqs_check is NO LONGER discarded here:
  # it is needed later (step 5) to build the final window with the overhang,
  # which cannot be derived from the seqkit sliding output (that only
  # contains window_size-aa pieces, not the complete alignment).

  # 1) Output directory of the lineage: exclusive to this lineage, no file
  #    or name shared with the other lineages
  lineage_dir <- file.path(base_dir, lineage)
  dir.create(lineage_dir, showWarnings = FALSE, recursive = TRUE)

  sliding_out <- file.path(lineage_dir, "ali.fasta")

  # 2) Sliding windows with seqkit (directly on the input FASTA)
  message(">>> seqkit sliding...")
  ok <- run_cmd(
    "seqkit",
    c("sliding", "-s", as.character(step_size), "-W", as.character(window_size), fasta_file),
    stdout_file = sliding_out,
    step_name = paste0("seqkit sliding (", lineage, ")")
  )
  if (!ok) next

  # 3) Read the sliding output and discover the windows dynamically
  whole <- tryCatch(readAAStringSet(sliding_out), error = function(e) NULL)
  if (is.null(whole) || length(whole) == 0) {
    warning(paste0("!!! Cannot read the sliding windows for ", lineage,
                    " - lineage skipped"), call. = FALSE)
    record_failure(paste0("unreadable sliding windows: ", lineage))
    next
  }

  window_labels <- str_extract(names(whole), "sliding:\\d+-\\d+$")
  window_labels <- unique(na.omit(window_labels))

  if (length(window_labels) == 0) {
    warning(paste0("!!! No 'sliding:start-end' header found for ", lineage,
                    " - lineage skipped"), call. = FALSE)
    record_failure(paste0("no sliding headers: ", lineage))
    next
  }

  # sort the windows numerically by start position
  starts <- as.integer(sub("sliding:(\\d+)-\\d+$", "\\1", window_labels))
  window_labels <- window_labels[order(starts)]
  # NB: starts must be recomputed AFTER sorting window_labels, otherwise it
  # stays in the original (pre-sort) order and no longer matches the order
  # of window_labels. The (sorted) starts vector is also needed further
  # down to identify the last regular window.
  starts <- as.integer(sub("sliding:(\\d+)-\\d+$", "\\1", window_labels))

  message(">>> ", length(window_labels), " windows found for ", lineage)

  # 4) Split per window + hmmbuild + hmmpress
  for (label in window_labels) {

    range_label <- sub("^sliding:", "", label)   # e.g. "1-100"
    pattern <- paste0(label, "$")                  # e.g. "sliding:1-100$" (anchored)

    matches <- grepl(pattern, names(whole))
    if (!any(matches)) next

    # files and model stay ONLY inside lineage_dir: no window or name is
    # shared with the directory of another lineage
    build_window(whole[matches], range_label, lineage_dir, lineage, threads)
  }

  # 5) Final window with the overhang (the residues left after the last
  #    "full" window of window_size aa).
  #
  #    Example with window_size=100, step_size=20 and an alignment 145 long:
  #    the regular windows are 1-100, 21-120, 41-140 (the last full window
  #    starts at 41 and ends at 140, 5 aa before the end). Here ONE extra
  #    window is added that continues the regular series from the next
  #    start (41 + step_size = 61) but runs to the end of the alignment
  #    instead of to 61+window_size-1=160: hence 61-145 (85 aa, shorter than
  #    the others because it is what is left).
  #
  #    This window is NOT generated by seqkit sliding (which, even in
  #    -g/--greedy mode, would produce several overlapping final windows
  #    with labels not aligned to the real length), but built here directly
  #    from the complete alignment (seqs_check) with Biostrings::subseq, so
  #    the "start-end" label always reflects the interval actually extracted.
  ali_len <- min(width(seqs_check))

  last_start <- max(starts)
  last_label <- window_labels[which.max(starts)]
  last_end <- as.integer(sub("^sliding:\\d+-(\\d+)$", "\\1", last_label))

  if (last_end < ali_len) {
    extra_start <- last_start + step_size
    extra_end <- ali_len

    if (extra_start <= ali_len) {
      extra_range_label <- paste0(extra_start, "-", extra_end)
      message(">>> Final window with the overhang: ", extra_range_label)

      extra_seqs <- tryCatch(
        subseq(seqs_check, start = extra_start, end = extra_end),
        error = function(e) NULL
      )

      if (is.null(extra_seqs)) {
        warning(paste0("!!! Cannot build the final window ", extra_range_label,
                        " for ", lineage, " - window skipped"), call. = FALSE)
        record_failure(paste0("final window ", extra_range_label, " (", lineage, ")"))
      } else {
        names(extra_seqs) <- paste0(names(seqs_check), "_sliding:", extra_range_label)
        build_window(extra_seqs, extra_range_label, lineage_dir, lineage, threads)
      }
    }
  }

  message(">>> Completed: ", lineage)
}

if (length(failures) > 0) {
  message("\nPipeline finished WITH ", length(failures), " failed or skipped step(s):")
  for (f in failures) message("  - ", f)
  quit(status = 1)
}

message("\nPipeline completed for all lineages.")
