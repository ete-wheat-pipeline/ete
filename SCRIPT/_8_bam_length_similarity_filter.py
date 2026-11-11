#!/usr/bin/env python3
"""
_8_bam_length_similarity_filter.py

Keeps only the mapped reads whose alignment covers at least 90 % of the
read length and has at least 90 % identity ((aligned_len - NM) / aligned_len).
Used for the repetitiveness (stage 8) and expression (stage 10) validations.

Usage:
    python3 _8_bam_length_similarity_filter.py <input.bam> <output.bam>
"""
import sys
import pysam

if len(sys.argv) != 3:
    sys.exit("Usage: python3 _8_bam_length_similarity_filter.py <input.bam> <output.bam>")

bam_in = sys.argv[1]
bam_out = sys.argv[2]
LENGTH_FRACTION_MIN = 0.9
SIMILARITY_FRACTION_MIN = 0.9

with pysam.AlignmentFile(bam_in, "rb") as inbam:
    with pysam.AlignmentFile(bam_out, "wb", template=inbam) as outbam:
        for read in inbam:
            if read.is_unmapped:
                continue
            read_len = read.infer_read_length()
            aligned_len = read.query_alignment_length
            if not read_len or not aligned_len:
                continue
            try:
                nm = read.get_tag("NM")
            except KeyError:
                continue

            length_fraction = aligned_len / read_len
            similarity_fraction = (aligned_len - nm) / aligned_len

            if length_fraction >= LENGTH_FRACTION_MIN and similarity_fraction >= SIMILARITY_FRACTION_MIN:
                outbam.write(read)
