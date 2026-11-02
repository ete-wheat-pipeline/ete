#!/usr/bin/env bash
# =============================================================================
# _5_massive_hmm_search.sh
#
# Scans a proteome with every sliding-window HMM built by
# _4_sliding_and_hmm_build_overhang_win.R:
#
#   for each <hmm_dir>/<family>/<window>.hmm (pressed with hmmpress):
#     hmmscan --cpu N --max -E 1e-3 --tblout --domtblout <window>.hmm <pep.fasta>
#
# Output, written in the CURRENT directory, one sub-directory per family:
#   <family>/<window>.tblout
#   <family>/<window>.domtblout
#   <family>/<window>.hmmscan.out
#
# Usage:
#   _5_massive_hmm_search.sh <pep.fasta> <hmm_dir> [n_cpu]
#
#   n_cpu defaults to 28 (the value used for the published results). The
#   number of worker threads does not change hmmscan results.
# =============================================================================
set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage:"
    echo "  $0 <pep.fasta> <hmm_dir> [n_cpu]"
    echo "Example:"
    echo "  $0 /path/to/pep.fasta /path/to/hmm_input 8"
    exit 1
fi

PEP_FASTA="$(realpath "$1")"
INPUT_DIR="$(realpath "$2")"
CPU="${3:-28}"
OUTPUT_DIR="$PWD"

if [ ! -f "$PEP_FASTA" ]; then
    echo "Error: FASTA file not found:"
    echo "$PEP_FASTA"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: input directory not found:"
    echo "$INPUT_DIR"
    exit 1
fi

[[ "$CPU" =~ ^[0-9]+$ ]] || { echo "Error: n_cpu must be a positive integer (got: '$CPU')." >&2; exit 1; }

echo "FASTA database : $PEP_FASTA"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "hmmscan --cpu  : $CPU"
echo

n_models=0
for d in "$INPUT_DIR"/*/; do
    [ -d "$d" ] || continue
    family="$(basename "$d")"
    family_output="$OUTPUT_DIR/$family"
    mkdir -p "$family_output"

    for f in "$d"/*.hmm; do
        [ -f "$f" ] || continue
        model="$(basename "${f%.hmm}")"
        echo "[$family] $model"

        hmmscan \
            --cpu "$CPU" \
            --max \
            -E 1e-3 \
            --tblout "$family_output/${model}.tblout" \
            --domtblout "$family_output/${model}.domtblout" \
            "$f" \
            "$PEP_FASTA" \
            > "$family_output/${model}.hmmscan.out"

        # the tabular files always contain a header, so an empty file means
        # hmmscan did not complete
        for out in "$family_output/${model}.tblout" "$family_output/${model}.domtblout"; do
            [ -s "$out" ] || { echo "Error: empty output $out" >&2; exit 1; }
        done
        n_models=$((n_models + 1))
    done
done

if [ "$n_models" -eq 0 ]; then
    echo "Error: no *.hmm file found under $INPUT_DIR/<family>/" >&2
    exit 1
fi

echo "Analysis completed: $n_models models scanned."
