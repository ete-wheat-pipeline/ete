#!/usr/bin/env python3
"""
_8_filter_sirna.py - Filters a FASTA file of small RNAs (e.g. from plantsmallrnagenes):
  1. removes the sequences whose headers indicate a microRNA
  2. keeps only sequences exactly 24 nt long (changeable with -l)
  3. removes redundant sequences (keeps the first occurrence)

Usage:
    python3 ./_8_filter_sirna.py -i sirna_file.fasta -o sirna_filtered.fasta
"""

import argparse
import re
import sys

# Default pattern to recognize miRNAs in the header:
# "miRNA", "microRNA", "micro-RNA", "miR166", "mir-156", "MIR159a", etc.
DEFAULT_MIRNA_PATTERN = r"micro-?rna|mirna|(?<![a-z])mir[-_]?\d+|(?<![a-z])mir(?![a-z])"


def read_fasta(handle):
    """Generator of (header, sequence) - also handles multi-line sequences."""
    header, seq = None, []
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq)
            header, seq = line[1:], []
        else:
            seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def main():
    ap = argparse.ArgumentParser(description="Filter siRNAs: no miRNA, fixed length, no duplicates.")
    ap.add_argument("-i", "--input", required=True, help="input FASTA")
    ap.add_argument("-o", "--output", required=True, help="filtered output FASTA")
    ap.add_argument("-l", "--length", type=int, default=24, help="required length in nt (default: 24)")
    ap.add_argument("-p", "--mirna-pattern", default=DEFAULT_MIRNA_PATTERN,
                    help="regex (case-insensitive) to recognize miRNAs in the header")
    args = ap.parse_args()

    mirna_re = re.compile(args.mirna_pattern, re.IGNORECASE)
    seen = set()
    n_tot = n_mirna = n_len = n_dup = n_out = 0

    with open(args.input) as fin, open(args.output, "w") as fout:
        for header, seq in read_fasta(fin):
            n_tot += 1

            # 1) miRNA filter (on the header)
            if mirna_re.search(header):
                n_mirna += 1
                continue

            # normalize: upper case and U -> T, so RNA/DNA and upper/lower case
            # do not hide duplicates
            seq = seq.upper().replace("U", "T")

            # 2) length filter
            if len(seq) != args.length:
                n_len += 1
                continue

            # 3) redundancy filter (sequence identity, keep the first one)
            if seq in seen:
                n_dup += 1
                continue
            seen.add(seq)

            fout.write(f">{header}\n{seq}\n")
            n_out += 1

    sys.stderr.write(
        f"Total sequences:                {n_tot}\n"
        f"Discarded (miRNA):              {n_mirna}\n"
        f"Discarded (length != {args.length}):     {n_len}\n"
        f"Discarded (duplicates):         {n_dup}\n"
        f"Sequences in output:            {n_out}\n"
    )


if __name__ == "__main__":
    main()
