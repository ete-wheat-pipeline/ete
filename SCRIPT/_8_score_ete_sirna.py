#!/usr/bin/env python3
"""
_8_score_ete_sirna.py

Assigns the "siRNA coverage" score to the putative ETEs following the
protocol of Ventimiglia et al. 2023 (Plant J 113:734-748), and computes the
background statistics that the original protocol does not report but that
are needed to know whether the score really discriminates.

SCORING SCHEME (Ventimiglia, siRNA module):
    +1  -> the ETE is NOT targeted by any siRNA            (good sign)
     0  -> targeted, but only in an intronic region        (mild penalty)
    -1  -> targeted in an EXON (UTRs included)             (severe penalty)

Note (translated from the original comment): there is NO -1 level in the
siRNA module of the paper, and there is no tblastn step of the HMM matching
regions against the siRNAs (that belongs to the "similarity with known TEs"
module, with a nucleotide subject).
[Translator's note: the first half of this statement conflicts with the
scheme above, which assigns -1 to exonic targeting. The scheme implemented
below is the one that produced the published results; the statement is
kept for the author to confirm or correct.]

CRITICAL POINT - MULTIMAPPING:
BWA samse reports ONE primary alignment in column 3, and the additional
alignments in the XA:Z: tag. With "-n 500" the XA tag contains up to 500
loci. Ignoring XA means losing the vast majority of the hits on a TE-rich
genome such as wheat. This script reads BOTH.

Output columns (read by _11_final_score_and_check_aggregator.R):
    locus  sirna_score  category
with category in {not_targeted, intronic, exonic}.

Usage:
    python3 _8_score_ete_sirna.py \\
        --ete-list putative_ete_gene_list.txt \\
        --sam-ete-genes bwa_out/sirna_vs_ETE_genes.sam \\
        --sam-ete-exons bwa_out/sirna_vs_ETE_exons.sam \\
        --sam-all-genes bwa_out/sirna_vs_all_genes.sam \\
        --all-genes-fasta genes.fa \\
        -o sirna_score_ventimiglia.tsv
"""

import argparse
import os
import sys


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------

def strip_isoform(name: str) -> str:
    """
    Normalizes a reference ID to the bare gene ID.

    Handles TWO different formats that appear in this pipeline:
      - bedtools getfasta -name : "TraesCS1A03G0000200::Chr1A:40097-70338(-)"
      - gffread -w             : "TraesCS1A03G0000200.1"
    In both cases it returns "TraesCS1A03G0000200".

    Without the cut on "::" the IDs of the genomic SAM files would never
    match the ETE list, and the scoring would fail entirely.
    """
    return name.split("::")[0].split(".")[0]


def parse_sam_references(sam_path: str, strip_iso: bool = True):
    """
    Returns the set of references (genes) that received at least one hit.

    Reads BOTH column 3 (primary alignment) AND the XA:Z: tag
    (multiple alignments produced by 'bwa samse -n N').

    XA format: XA:Z:ref,pos,CIGAR,NM;ref,pos,CIGAR,NM;...
    """
    hit_refs = set()
    n_lines = 0
    n_primary = 0
    n_xa_entries = 0

    if not os.path.exists(sam_path):
        sys.exit(f"ERROR: SAM file not found: {sam_path}")

    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            n_lines += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue

            rname = fields[2]
            if rname != "*":
                n_primary += 1
                hit_refs.add(strip_isoform(rname) if strip_iso else rname)

            # XA tag in the optional fields (from column 12 onwards)
            for tag in fields[11:]:
                if tag.startswith("XA:Z:"):
                    payload = tag[5:]
                    for entry in payload.split(";"):
                        if not entry:
                            continue
                        ref = entry.split(",")[0]
                        if ref and ref != "*":
                            n_xa_entries += 1
                            hit_refs.add(strip_isoform(ref) if strip_iso else ref)
                    break

    print(f"  {os.path.basename(sam_path)}:")
    print(f"    alignment lines       : {n_lines}")
    print(f"    primary hits (col 3)  : {n_primary}")
    print(f"    extra hits from XA tag: {n_xa_entries}")
    print(f"    distinct references   : {len(hit_refs)}")
    return hit_refs


def read_id_list(path: str, strip_iso: bool = True):
    """Reads a list of IDs (one per line). Ignores empty lines and comments."""
    ids = []
    seen = set()
    if not os.path.exists(path):
        sys.exit(f"ERROR: file not found: {path}")
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            gid = line.split()[0]
            if strip_iso:
                gid = strip_isoform(gid)
            if gid not in seen:
                seen.add(gid)
                ids.append(gid)
    return ids


def count_fasta_headers(path: str, strip_iso: bool = True):
    """Returns the set of IDs (first token of the header) of a FASTA."""
    ids = set()
    if not os.path.exists(path):
        return ids
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                gid = line[1:].split()[0]
                ids.add(strip_isoform(gid) if strip_iso else gid)
    return ids


# --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------

def fisher_and_or(a, b, c, d):
    """
    2x2 table:
        a = ETE positive      b = ETE negative
        c = non-ETE positive  d = non-ETE negative
    Returns (odds_ratio, p_value) using scipy if available.
    """
    try:
        from scipy.stats import fisher_exact
        orr, p = fisher_exact([[a, b], [c, d]])
        return orr, p
    except ImportError:
        # fallback: odds ratio only, no p-value
        if b == 0 or c == 0 or d == 0:
            return float("nan"), None
        orr = (a / b) / (c / d)
        return orr, None


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="siRNA scoring of the putative ETEs following Ventimiglia et al. 2023"
    )
    ap.add_argument("--ete-list", required=True,
                    help="List of the putative ETE loci (one ID per line)")
    ap.add_argument("--sam-ete-genes", required=True,
                    help="SAM: siRNA vs putative ETE genes (whole gene)")
    ap.add_argument("--sam-ete-exons", required=True,
                    help="SAM: siRNA vs exons (mRNA) of the putative ETEs")
    ap.add_argument("--sam-all-genes", default=None,
                    help="SAM: siRNA vs all genes (for the background, optional)")
    ap.add_argument("--all-genes-fasta", default=None,
                    help="FASTA of all genes (for the background total)")
    ap.add_argument("-o", "--output", default="sirna_score_ventimiglia.tsv",
                    help="Output TSV file")
    ap.add_argument("--no-strip-isoform", action="store_true",
                    help="Do not remove the isoform suffix (.1, .2) from the IDs")
    args = ap.parse_args()

    strip_iso = not args.no_strip_isoform

    print("=" * 68)
    print("PARSING THE SAM FILES")
    print("=" * 68)

    pos_ete_genes = parse_sam_references(args.sam_ete_genes, strip_iso)
    pos_ete_exons = parse_sam_references(args.sam_ete_exons, strip_iso)

    pos_all_genes = None
    if args.sam_all_genes:
        pos_all_genes = parse_sam_references(args.sam_all_genes, strip_iso)

    print()
    print("=" * 68)
    print("ID CONSISTENCY CHECKS")
    print("=" * 68)

    ete_list = read_id_list(args.ete_list, strip_iso)
    ete_set = set(ete_list)
    print(f"  Putative ETEs in the list: {len(ete_list)}")

    # Critical check: if the IDs do not match, everything would silently
    # become +1. Better to notice immediately.
    overlap_genes = len(pos_ete_genes & ete_set)
    print(f"  Positive ETEs (genes) matching the list: "
          f"{overlap_genes} / {len(pos_ete_genes)}")
    if pos_ete_genes and overlap_genes == 0:
        print()
        print("  *** WARNING: NO ID MATCHES. ***")
        print("  The IDs in the SAM do not correspond to those in the ETE list.")
        print("  Without this match every gene would wrongly be scored +1.")
        print(f"  Example IDs from the SAM : {list(pos_ete_genes)[:3]}")
        print(f"  Example IDs from the list: {ete_list[:3]}")
        sys.exit(1)
    elif pos_ete_genes and overlap_genes < 0.5 * len(pos_ete_genes):
        print("  *** WARNING: less than 50% of the IDs match. Check the formats. ***")

    # Logical consistency: exons are a subset of the gene, so every exonic
    # positive should also be positive on the whole gene.
    exon_not_in_gene = (pos_ete_exons & ete_set) - pos_ete_genes
    if exon_not_in_gene:
        print(f"  Note: {len(exon_not_in_gene)} loci positive on the exons but NOT "
              f"on the whole gene.")
        print("  This can happen if the two FASTA files have different origins; check.")

    print()
    print("=" * 68)
    print("SCORING (Ventimiglia scheme: +1 / 0 / -1)")
    print("=" * 68)

    rows = []
    counts = {1: 0, 0: 0, -1: 0}
    for gid in ete_list:
        if gid not in pos_ete_genes:
            score = 1          # not targeted
            cat = "not_targeted"
        elif gid in pos_ete_exons:
            score = -1         # hit in an exon (UTRs included)
            cat = "exonic"
        else:
            score = 0          # intronic hit only
            cat = "intronic"
        counts[score] += 1
        rows.append((gid, score, cat))

    with open(args.output, "w") as out:
        out.write("locus\tsirna_score\tcategory\n")
        for gid, score, cat in rows:
            out.write(f"{gid}\t{score}\t{cat}\n")

    total = len(rows)
    print(f"  +1  not targeted   : {counts[1]:>7}  ({100*counts[1]/total:.2f}%)")
    print(f"   0  intronic only  : {counts[0]:>7}  ({100*counts[0]/total:.2f}%)")
    print(f"  -1  exonic         : {counts[-1]:>7}  ({100*counts[-1]/total:.2f}%)")
    print(f"  TOTAL              : {total:>7}")
    print(f"\n  Written: {args.output}")

    # ----------------------------------------------------------------------
    # Background: are the ETEs really more targeted than ordinary genes?
    # ----------------------------------------------------------------------
    if pos_all_genes is not None and args.all_genes_fasta:
        print()
        print("=" * 68)
        print("BACKGROUND CHECK (ETE vs non-ETE genes)")
        print("=" * 68)

        all_genes = count_fasta_headers(args.all_genes_fasta, strip_iso)
        if not all_genes:
            print("  Cannot read the FASTA of all genes, skipping.")
            return

        tot_genes = len(all_genes)
        pos_genes = len(pos_all_genes & all_genes)

        tot_ete = len(ete_set & all_genes) or len(ete_set)
        pos_ete = len(pos_all_genes & ete_set)

        tot_non = tot_genes - tot_ete
        pos_non = pos_genes - pos_ete
        neg_ete = tot_ete - pos_ete
        neg_non = tot_non - pos_non

        if tot_ete <= 0 or tot_non <= 0:
            print("  Counts not valid for the test, skipping.")
            return

        print(f"  ETE      : {pos_ete:>7} / {tot_ete:>7} = {100*pos_ete/tot_ete:.2f}%")
        print(f"  non-ETE  : {pos_non:>7} / {tot_non:>7} = {100*pos_non/tot_non:.2f}%")
        print(f"  all      : {pos_genes:>7} / {tot_genes:>7} = {100*pos_genes/tot_genes:.2f}%")

        orr, p = fisher_and_or(pos_ete, neg_ete, pos_non, neg_non)
        print()
        print(f"  Odds ratio : {orr:.3f}")
        if p is not None:
            print(f"  Fisher p   : {p:.3e}")
        else:
            print("  Fisher p   : (scipy not available)")

        expected = tot_ete * pos_non / tot_non if tot_non else 0
        excess = pos_ete - expected
        print(f"  ETE positives expected at the background rate: {expected:.0f}")
        print(f"  Observed excess: {excess:.0f}")
        if pos_ete > 0:
            print(f"  Fraction of the positives truly 'in excess': "
                  f"{100*excess/pos_ete:.1f}%")

        print()
        if orr == orr and orr < 1.5:
            print("  READING: low odds ratio. Even if the p-value were")
            print("  significant (an effect of the very large n), the siRNA score")
            print("  discriminates poorly between ETEs and ordinary genes. Consider")
            print("  using the siRNA density per kb with quartile thresholds")
            print("  instead of simple presence/absence.")
        elif orr == orr:
            print("  READING: appreciable enrichment in the ETEs relative to the background.")


if __name__ == "__main__":
    main()
