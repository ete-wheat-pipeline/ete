#!/usr/bin/env bash
# =============================================================================
# _11_te_derived_fusion.sh
#
# Splits the class-4 ETE genes (all of them, no isoform filter) into
# TE-derived_host_gene and TE-gene-fusion.
#
# Isoforms considered (from the aggregator output):
#   SCORE_SUM == 4                 class-4 ETE
#   passes_coverage == TRUE        the isoform passes the coverage filter (>= 25 aa and >= 10 %)
# Coverage = covered_ali_len / qlen * 100 (without the rounding of the
# percentage_ali_len column).
#
# A gene is a TE-derived_host_gene only if:
#   1. all its passing isoforms have coverage >= THRESHOLD, and
#   2. all its isoforms in the scanned proteome were recognized by hmmscan,
#      i.e. they appear in merged.domtblout.
# Otherwise the whole gene, with all its isoforms, is a TE-gene-fusion.
# The printed summary gives the reason:
#   coverage_below_threshold   all passing isoforms are < THRESHOLD
#   mixed                      isoforms above and below THRESHOLD
#   unrecognized_isoforms      all >= THRESHOLD, but at least one isoform of
#                              the proteome is missing from merged.domtblout
#
# Usage (from the 11_RESULTS directory, after _11_final_score_and_check_aggregator.R):
#   bash ../SCRIPT/_11_te_derived_fusion.sh [aggregator_file] [threshold] [proteome]
#   defaults: final_validation_result_filterable.tsv  70
#             ../0_STARTING_FILES/fine_cleaned_pep.fasta (the proteome given to hmmscan)
#
# Output:
#   class4_te_derived_fusion.tsv        one row per passing isoform
#   unrecognized_isoforms.tsv           isoforms missing from the domtable that
#                                       moved a gene into TE-gene-fusion
#   summary on screen
# =============================================================================

set -euo pipefail

IN="${1:-final_validation_result_filterable.tsv}"
THRESHOLD="${2:-70}"
PEP="${3:-../0_STARTING_FILES/fine_cleaned_pep.fasta}"
DOM="../6_SCANNING_RESULT/merged.domtblout"
OUT="class4_te_derived_fusion.tsv"
OUT_NR="unrecognized_isoforms.tsv"

for f in "$IN" "$PEP" "$DOM"; do
  [ -s "$f" ] || { echo "File not found or empty: $f" >&2; exit 1; }
done

REC=$(mktemp); ISO=$(mktemp)
trap 'rm -f "$REC" "$ISO"' EXIT

# isoforms recognized by hmmscan (column 4 of merged.domtblout)
grep -v '^#' "$DOM" | awk '{print $4}' | sort -u > "$REC"
# isoforms of the scanned proteome (first field of the header)
grep '^>' "$PEP" | sed 's/^>//; s/[[:space:]].*//' > "$ISO"

awk -F'\t' -v s="$THRESHOLD" -v out="$OUT" -v out_nr="$OUT_NR" '
# 1. aggregator file: passing class-4 isoforms
FILENAME == ARGV[1] {
  if (FNR == 1) {
    for (i = 1; i <= NF; i++) col[$i] = i
    n = split("locus query qlen covered_ali_len SCORE_SUM passes_coverage", needed, " ")
    for (i = 1; i <= n; i++)
      if (!(needed[i] in col)) { print "Missing column in " FILENAME ": " needed[i] > "/dev/stderr"; error = 1; exit 1 }
    next
  }
  if ($col["SCORE_SUM"] != "4" || $col["passes_coverage"] != "TRUE") next
  gene = $col["locus"]
  pct  = $col["covered_ali_len"] * 100 / $col["qlen"]
  r++
  row_gene[r] = gene
  row[r] = sprintf("%s\t%s\t%s\t%s\t%.2f", gene, $col["query"], $col["qlen"], $col["covered_ali_len"], pct)
  genes4[gene] = 1
  if (pct < s) below[gene] = 1; else above[gene] = 1
  next
}
# 2. isoforms present in merged.domtblout
FILENAME == ARGV[2] { rec[$1] = 1; next }
# 3. isoforms of the proteome: those of class-4 genes missing from the domtable
{
  iso = $1; gene = iso; sub(/\.[^.]+$/, "", gene)
  if (!(gene in genes4)) next
  in_proteome[gene] = 1
  if (!(iso in rec)) { nonrec[gene]++; nr++; nr_gene[nr] = gene; nr_iso[nr] = iso }
}
END {
  if (error) exit 1

  # class and reason for every gene
  for (g in genes4) {
    if ((g in below) && (g in above)) { class[g] = "TE-gene-fusion"; reason[g] = "mixed" }
    else if (g in below)              { class[g] = "TE-gene-fusion"; reason[g] = "coverage_below_threshold" }
    else if (g in nonrec)             { class[g] = "TE-gene-fusion"; reason[g] = "unrecognized_isoforms" }
    else                              { class[g] = "TE-derived_host_gene"; reason[g] = "-" }
    ng++; genes[class[g]]++; nreason[reason[g]]++
    if (!(g in in_proteome)) no_pep++
  }

  # per-isoform table
    print "locus\tquery\tqlen\tcovered_ali_len\tcoverage_pct\tclassification" > out
  for (k = 1; k <= r; k++) {
    g = row_gene[k]
    print row[k] "\t" class[g] > out
    iso_cl[class[g]]++; tot++
  }

  # unrecognized isoforms that moved a gene
  print "locus\tunrecognized_isoform" > out_nr
  for (k = 1; k <= nr; k++)
    if (reason[nr_gene[k]] == "unrecognized_isoforms") { print nr_gene[k] "\t" nr_iso[k] > out_nr; n_iso_moved++ }

  printf "Coverage threshold: %s%%\n\n", s
  printf "%-22s %10s %8s\n", "class", "isoforms", "genes"
  printf "%-22s %10d %8d\n", "TE-derived_host_gene", iso_cl["TE-derived_host_gene"], genes["TE-derived_host_gene"]
  printf "%-22s %10d %8d\n", "TE-gene-fusion", iso_cl["TE-gene-fusion"], genes["TE-gene-fusion"]
  printf "%-22s %10d %8d\n", "total", tot, ng
  print ""
  print "TE-gene-fusion genes by reason:"
  printf "  coverage below threshold:                     %d\n", nreason["coverage_below_threshold"]
  printf "  mixed (isoforms above and below threshold):   %d\n", nreason["mixed"]
  printf "  >= threshold but isoforms not recognized:     %d  (%d isoforms)\n", nreason["unrecognized_isoforms"], n_iso_moved
  if (no_pep) printf "\nWARNING: %d class-4 genes without isoforms in the proteome (different IDs or wrong proteome?)\n", no_pep
  if (tot == 0) print "\nWARNING: no class-4 isoform passing the coverage filter." > "/dev/stderr"
  printf "\nWritten: %s  %s\n", out, out_nr
}' "$IN" "$REC" "$ISO"
