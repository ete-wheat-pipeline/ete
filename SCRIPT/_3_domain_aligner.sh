#!/usr/bin/env bash
# =============================================================================
# _3_domain_aligner.sh
#
# Takes all the FASTA files produced by the R scripts (TIR/LTR...) in a
# directory, groups them by "lineage" based on the file name, and for each
# lineage runs:
#
#   mafft --thread N --reorder --maxiterate 1000 --globalpair <input> > <lineage>_ali.fasta
#
# If several files map to the same lineage (e.g. the same lineage extracted
# from several subgenomes) they are concatenated BEFORE aligning, so one
# alignment per lineage is produced instead of overwriting the output.
#
# Usage:
#   ./_3_domain_aligner.sh <fasta_dir> [n_threads] [output_dir]
#
#   fasta_dir   : directory containing the .fasta/.fa files to align (required)
#   n_threads   : threads for mafft (default: 4)
#   output_dir  : where to write the *_ali.fasta files (default: current
#                 directory, i.e. where the script is launched)
# =============================================================================

set -euo pipefail

usage() {
  echo "Usage: $0 <fasta_dir> [n_threads] [output_dir]"
  echo "  fasta_dir  : directory with the .fasta/.fa files to align"
  echo "  n_threads  : number of mafft threads (default: 4)"
  echo "  output_dir : destination directory for the *_ali.fasta files (default: current dir)"
  exit 1
}

[ $# -ge 1 ] || usage
IN_DIR="$1"
THREADS="${2:-4}"
OUT_DIR="${3:-.}"

[ -d "$IN_DIR" ] || { echo "Error: '$IN_DIR' is not a directory." >&2; exit 1; }
[[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "Error: n_threads must be a positive integer (got: '$THREADS')." >&2; exit 1; }
command -v mafft >/dev/null 2>&1 || { echo "Error: 'mafft' not found in PATH." >&2; exit 1; }

mkdir -p "$OUT_DIR"

# -----------------------------------------------------------------------
# Map: keyword recognized in the file name (case-insensitive, as an
# isolated token) -> lineage name used in the output (<lineage>_ali.fasta)
# Includes aliases to handle both the DANTE_TIR notation (underscore) and
# the DANTE notation (pipe/slash) used in the names/headers generated so far.
# -----------------------------------------------------------------------
declare -A LINEAGE_MAP=(
  # Class I - Copia
  [ALE]="Ale"
  [ANGELA]="Angela"
  [BIANCA]="Bianca"
  [IKEROS]="Ikeros"
  [IVANA]="Ivana"
  [SIRE]="SIRE"
  [TAR]="TAR"
  # Class I - Gypsy
  [ATHILA]="Athila"
  [CRM]="CRM"
  [OGRE]="Ogre"
  [REINA]="Reina"
  [RETAND]="Retand"
  [TEKAY]="Tekay"
  [TORK]="Tork"
  # Class II - TIR
  [CACTA]="CACTA"
  [ENSPM]="CACTA"        # EnSpm/CACTA
  [HARBINGER]="Harbinger"
  [PIF]="Harbinger"      # PIF/Harbinger
  [MUTATOR]="Mutator"
  [MUDR]="Mutator"       # MuDR/Mutator
  [HAT]="hAT"
  [MARINER]="Mariner"
  [TC1]="Mariner"        # Tc1/Mariner
)

shopt -s nullglob
FASTA_FILES=("$IN_DIR"/*.fasta "$IN_DIR"/*.fa "$IN_DIR"/*.faa)
shopt -u nullglob

if [ ${#FASTA_FILES[@]} -eq 0 ]; then
  echo "No .fasta/.fa/.faa file found in '$IN_DIR'." >&2
  exit 1
fi

echo "Found ${#FASTA_FILES[@]} FASTA files in '$IN_DIR'"
echo "Output in: '$OUT_DIR'  |  mafft threads: $THREADS"
echo

# group the files by the lineage recognized in the name (a file can end up
# in one lineage only: the loop stops at the first match)
declare -A LINEAGE_FILES
SKIPPED=()

for f in "${FASTA_FILES[@]}"; do
  base=$(basename "$f")
  name_no_ext="${base%.*}"
  upper_name=$(printf '%s' "$name_no_ext" | tr '[:lower:]' '[:upper:]')

  matched=""
  for key in "${!LINEAGE_MAP[@]}"; do
    if [[ "$upper_name" =~ (^|[^A-Z0-9])${key}([^A-Z0-9]|$) ]]; then
      matched="${LINEAGE_MAP[$key]}"
      break
    fi
  done

  if [ -z "$matched" ]; then
    SKIPPED+=("$base")
    continue
  fi

  LINEAGE_FILES["$matched"]+="${f}"$'\n'
done

if [ ${#SKIPPED[@]} -gt 0 ]; then
  echo "Files ignored (no known lineage in the name):"
  printf '  - %s\n' "${SKIPPED[@]}"
  echo
fi

if [ ${#LINEAGE_FILES[@]} -eq 0 ]; then
  echo "No file matches a known lineage. Nothing to align." >&2
  exit 1
fi

# sort the lineages alphabetically only for a more readable output
mapfile -t LINEAGES_SORTED < <(printf '%s\n' "${!LINEAGE_FILES[@]}" | sort)

for lineage in "${LINEAGES_SORTED[@]}"; do
  files_list="${LINEAGE_FILES[$lineage]}"
  n_files=$(printf '%s' "$files_list" | grep -c . || true)

  combined="$(mktemp "${OUT_DIR}/.${lineage}_combined.XXXXXX.fasta")"
  : > "$combined"
  while IFS= read -r file; do
    [ -n "$file" ] && cat "$file" >> "$combined"
  done <<< "$files_list"

  n_seqs=$(grep -c '^>' "$combined" || true)

  if [ "$n_seqs" -lt 2 ]; then
    echo "[skip] lineage=$lineage : only $n_seqs sequence(s) found, at least 2 are needed for an alignment."
    rm -f "$combined"
    continue
  fi

  out_file="$OUT_DIR/${lineage}_ali.fasta"
  log_file="$OUT_DIR/${lineage}_ali.log"

  echo "[align] lineage=$lineage  input_files=$n_files  sequences=$n_seqs  ->  $out_file"
  mafft --thread "$THREADS" --reorder --maxiterate 1000 --globalpair "$combined" > "$out_file" 2> "$log_file"

  rm -f "$combined"
done

echo
echo "== DONE =="
