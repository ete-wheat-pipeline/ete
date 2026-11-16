# PETER Protocol Overview

PETER provides a genome-wide framework to search for ETEs among annotated plant genes by combining the detection of TE-protein signatures along with structural and transcriptional evidence. The workflow begins with a curated collection of full-length TEs isolated from high-quality genome sequence. TE intact coding domains are grouped by superfamily or lineage and used to build Hidden Markov Model (HMM) profiles from overlapping alignment windows. These profiles are searched against the complete predicted proteome, and genes whose encoded proteins contain matching regions spanning at least 25 amino acids and 10% of their length are designated as putative ETEs.
Each putative ETE is then assessed using four complementary attributes: nucleotide similarity to known TEs, genomic repetitiveness, 24-nt siRNA targeting, and expression. The first two assess whether candidates still closely resemble TEs by sequence similarity and genomic abundance. Over time, ETEs are expected to diverge from related TEs and generally occur at low copy numbers. The other two attributes examine the putative ETEs activity. Mapping 24-nt siRNAs identifies potential targeting by a class of small RNAs involved in TE silencing through RNA-directed DNA methylation (Blevins et al. 2015), whereas RNA-seq provides evidence of transcription in the sampled tissues. ETEs, as conventional genes, show higher expression and are less targeted by siRNA than TEs (Cowan et al. 2005). The scores assigned to these four attributes are combined for each putative ETE: those receiving +1 in every analysis reach a cumulative score of +4 and are considered validated ETEs under the PETER criteria. Finally, the extent of the HMM matching regions within each validated protein is used to distinguish TE-derived genes from TE-gene fusion events.
The protocol can be used with a single representative protein per gene or with all annotated isoforms. Screening all isoforms is recommended at the discovery stage because exons derived from TE exaptation may be isoform specific. 
PETER is applicable to any genome assembly; nevertheless, assembly completeness and contiguity directly influence detection performance, with high-contiguity assemblies generally providing more reliable gene predictions and, particularly, more comprehensive TE catalogues.


<img width="671" height="914" alt="immagine" src="https://github.com/user-attachments/assets/c7f6b077-dd5e-4a92-ba80-eb4cec8db0a2" />


References:

Blevins T, Podicheti R, Mishra V, Marasco M, Wang J, Rusch D, Tang H, Pikaard CS. 2015. Identification of Pol IV and RDR2-dependent precursors of 24 nt siRNAs guiding de novo DNA methylation in Arabidopsis. eLife. 4. 

Cowan RK, Hoen DR, Schoen DJ, Bureau TE. 2005. MUSTANG is a novel family of domesticated transposase genes found in diverse angiosperms. Mol Biol Evol. 22:2084–2089. 

# How to run PETER pipeline

You need only this file, `SCRIPT/` and `0_STARTING_FILES/`.

Rules for every step below:
- Copy the commands **in order**, into **one** terminal (Linux, or Ubuntu under WSL2 on Windows).
- Run them **from the folder that contains this file**.
- Start from a folder that contains nothing else (no old result folders).
- On a remote server, run inside `screen` so that a dropped connection does not stop the run.

## Before you begin: rebuild the genome
GitHub does not accept files larger than 100 MB, so the genome `Chr1A.fasta` (612 MB) is stored compressed and split into four parts: `0_STARTING_FILES/Chr1A.fasta.gz.part_aa` ... `part_ad`. Rebuild it once, from the folder that contains this file (works with both `git clone` and "Download ZIP"; no extra software needed):
```bash
cat 0_STARTING_FILES/Chr1A.fasta.gz.part_* | gunzip > 0_STARTING_FILES/Chr1A.fasta
sha256sum 0_STARTING_FILES/Chr1A.fasta
```
The checksum must be `0757ee471d54348cf1c641141a1a055dd169095dbd636227ab4072257990b904`. If it differs, one of the parts is missing or incomplete: download the repository again.

## Before you begin: input data
`0_STARTING_FILES/` must contain:
```
Chr1A.fasta   (rebuilt in the previous section)
fine_cleaned_gff.gff3
fine_cleaned_genes.fasta
fine_cleaned_mRNA.fasta
fine_cleaned_cds.fasta
fine_cleaned_pep.fasta
siRNA24_filtered_unique_DNA.fasta
DNA_READS/subsampled_500k/library_R1_500k.fastq.gz
DNA_READS/subsampled_500k/library_R2_500k.fastq.gz
RNA_READS/<run>_1.500k.fastq.gz and <run>_2.500k.fastq.gz   (9 SRA runs)
```
```bash
ls -R 0_STARTING_FILES
```

## Software and execution

### Install conda (skip if `conda --version` works)
```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3
~/miniforge3/bin/conda init bash
exec bash
```

### Create the environments (once per computer, one command at a time)
```bash
conda create -y -n dante -c conda-forge -c bioconda -c petrnovak dante=0.2.12
conda create -y -n dante_ltr -c conda-forge -c bioconda -c petrnovak dante_ltr=0.6.4.0
conda create -y -n dante_tir -c conda-forge -c r -c bioconda -c petrnovak dante_tir=0.3.2
conda create -y -n ete -c conda-forge -c bioconda mafft=7.526 hmmer==3.3 bwa=0.7.19 samtools seqkit=2.1.0 blast python=3.12 pysam scipy r-base r-tidyverse bioconductor-rtracklayer bioconductor-biostrings bioconductor-genomicalignments
```

### Settings (repeat in every new terminal)
`T` is the number of CPU threads. DANTE_LTR uses up to about 1 GB of RAM per thread, so use `T=2` on a PC with little memory.
```bash
export LC_ALL=C
T=10
```

## TE reference preparation and protein screening

### Step 0. Build the full-length TE reference
```bash
conda activate dante
mkdir -p 1_DANTE_suite/DANTE
cd 1_DANTE_suite/DANTE
dante -q ../../0_STARTING_FILES/Chr1A.fasta -o DANTE_output.gff3 -c $T
dante_gff_output_filtering.py -dg DANTE_output.gff3
cp domains_filtered.gff DANTE_filtered_output.gff3
cd ../..

conda activate dante_ltr
mkdir -p 1_DANTE_suite/DANTE_LTR
cd 1_DANTE_suite/DANTE_LTR
dante_ltr -g ../DANTE/DANTE_filtered_output.gff3 -s ../../0_STARTING_FILES/Chr1A.fasta -o DANTE_LTR -M 1 -c $T
dante_ltr_to_library -g DANTE_LTR.gff3 -s ../../0_STARTING_FILES/Chr1A.fasta -o ltr_library -c $T
cd ../..

conda activate dante_tir
cd 1_DANTE_suite
dante_tir.py -g DANTE/DANTE_filtered_output.gff3 -f ../0_STARTING_FILES/Chr1A.fasta -o TIR -c $T
cd ..
```
The pipeline uses `1_DANTE_suite/DANTE_LTR/ltr_library/TE_DLplus.fasta` and `1_DANTE_suite/TIR/DANTE_TIR_final.*`.

`dante_ltr_to_library` ends with an mmseqs2 clustering that needs several GB of RAM. If only that last step fails, `TE_DLplus.fasta` is already complete.

### Step 1. Select intact TE coding regions
```bash
conda activate ete
mkdir -p 2_EXTRACTED_DOMAINS
cd 2_EXTRACTED_DOMAINS
Rscript ../SCRIPT/_1_LTR_protein_domain_concatenator.R ../1_DANTE_suite/DANTE_LTR/DANTE_LTR.gff3
Rscript ../SCRIPT/_2_TIR_tpase_domain_extractor.R ../1_DANTE_suite/TIR/DANTE_TIR_final.gff3 ../1_DANTE_suite/DANTE/DANTE_filtered_output.gff3
cd ..
```

### Step 2. Build HMMs from overlapping alignment windows
```bash
mkdir -p 3_ALIGNMENT
cd 3_ALIGNMENT
bash ../SCRIPT/_3_domain_aligner.sh ../2_EXTRACTED_DOMAINS/ $T
cd ..

conda activate ete
mkdir -p 4_ALIGNMENT_PARSING_HMM_BUILD
cd 4_ALIGNMENT_PARSING_HMM_BUILD
Rscript ../SCRIPT/_4_sliding_and_hmm_build_overhang_win.R ../3_ALIGNMENT $T
cd ..
```

### Step 3. Identify HMM-positive protein isoforms
```bash
mkdir -p 5_HMM_PROTEOME_SCANNING
cd 5_HMM_PROTEOME_SCANNING
bash ../SCRIPT/_5_massive_hmm_search.sh ../0_STARTING_FILES/fine_cleaned_pep.fasta ../4_ALIGNMENT_PARSING_HMM_BUILD $T
cd ..

Rscript SCRIPT/_6_domtblout_concatenator_filter_v2.R -d 5_HMM_PROTEOME_SCANNING -p 0_STARTING_FILES/fine_cleaned_pep.fasta -c 0_STARTING_FILES/fine_cleaned_cds.fasta -g 0_STARTING_FILES/fine_cleaned_genes.fasta -m 0_STARTING_FILES/fine_cleaned_mRNA.fasta -o 6_SCANNING_RESULT
```

### Step 4. Extract putative ETE genes
```bash
grep '>' 6_SCANNING_RESULT/GENE_putative_ete.fasta | sed 's/^>//; s/ .*//' | sort -u > 6_SCANNING_RESULT/all_putative_ETEs_ids.txt
wc -l 6_SCANNING_RESULT/all_putative_ETEs_ids.txt
```

## Validation procedure

The four validations below are independent of each other: each one gives every putative ETE a score of +1, 0 or -1, and Step 9 adds them up (`SCORE_SUM`).
- **Structural validation** (Steps 5-6): similarity to known TEs and repetitiveness in the genome, both measured against the DANTE TE library.
- **Functional validation** (Steps 7-8): targeting by 24-nt siRNAs and expression in RNA-seq data.

### Step 5. Assess similarity to known TEs
```bash
conda activate ete
mkdir -p 7_TE_SIMILARITY_VALIDATION
cd 7_TE_SIMILARITY_VALIDATION
cat ../1_DANTE_suite/DANTE_LTR/ltr_library/TE_DLplus.fasta ../1_DANTE_suite/TIR/DANTE_TIR_final.fasta > all_TE.fasta
grep -c '>' ../1_DANTE_suite/TIR/DANTE_TIR_final.fasta   # must print the number of TIR elements (176 in the published run), not an error
makeblastdb -in all_TE.fasta -dbtype nucl
blastn -task dc-megablast -query ../6_SCANNING_RESULT/GENE_putative_ete.fasta -db all_TE.fasta -evalue 1e-5 -dust no -outfmt 6 -num_threads $T -out step1_genes_vs_TE.tsv
awk -F'\t' '$3 >= 70 && $4 >= 70 {print $1}' step1_genes_vs_TE.tsv | sort -u > step1_hits.txt
comm -23 ../6_SCANNING_RESULT/all_putative_ETEs_ids.txt step1_hits.txt > score_plus1.txt
sed 's/^/^/; s/$/\\./' step1_hits.txt > step1_hits_regex.txt
seqkit grep -n -r -f step1_hits_regex.txt ../6_SCANNING_RESULT/CDS_putative_ete.fasta -o isoforms_step1_hits.fasta
blastn -task dc-megablast -query isoforms_step1_hits.fasta -db all_TE.fasta -evalue 1e-5 -dust no -outfmt 6 -num_threads $T -out step2_isoforms_vs_TE.tsv
awk -F'\t' '$3 >= 70 && $4 >= 70 {print $1}' step2_isoforms_vs_TE.tsv | sed 's/\.[0-9]*$//' | sort -u > step2_hits.txt
comm -23 step1_hits.txt step2_hits.txt > score_0.txt
cp step2_hits.txt score_m1.txt
{ awk '{print $1"\t1"}' score_plus1.txt; awk '{print $1"\t0"}' score_0.txt; awk '{print $1"\t-1"}' score_m1.txt; } | sort -k1,1 > similarity_score_final.tsv
sed -i '1i locus\tsimilarity_score' similarity_score_final.tsv
cd ..
```

### Step 6. Assess repetitiveness from genomic read depth
`bwa mem` keeps the original thread count (`-t 28`): its batch size depends on the number of threads and slightly changes the paired-end mapping.
```bash
conda activate ete
mkdir -p 8_REPETITIVENESS_VALIDATION
cd 8_REPETITIVENESS_VALIDATION
cp ../1_DANTE_suite/DANTE_LTR/ltr_library/TE_DLplus.fasta ../1_DANTE_suite/TIR/DANTE_TIR_final.fasta ../0_STARTING_FILES/fine_cleaned_genes.fasta .
awk '/^>/{print ">LTR_" ++i; next}{print}' TE_DLplus.fasta > cleaned_LTR.fasta
awk '/^>/{print ">TIR_" ++i; next}{print}' DANTE_TIR_final.fasta > cleaned_TIR.fasta
cat cleaned_LTR.fasta cleaned_TIR.fasta > TE_only.fasta
bwa index -a bwtsw TE_only.fasta
bwa index fine_cleaned_genes.fasta
R1=../0_STARTING_FILES/DNA_READS/subsampled_500k/library_R1_500k.fastq.gz
R2=../0_STARTING_FILES/DNA_READS/subsampled_500k/library_R2_500k.fastq.gz
bwa mem -t 28 -A 1 -B 1 -O 1,1 -E 1,1 -k 15 -T 20 fine_cleaned_genes.fasta $R1 $R2 > mapped_genes.sam
samtools sort -@ $T -o mapped_genes.sorted.bam mapped_genes.sam
samtools index mapped_genes.sorted.bam
bwa mem -t 28 -A 1 -B 1 -O 1,1 -E 1,1 -k 15 -T 20 TE_only.fasta $R1 $R2 > mapped_TEonly.sam
samtools sort -@ $T -o mapped_TEonly.sorted.bam mapped_TEonly.sam
samtools index mapped_TEonly.sorted.bam
python ../SCRIPT/_8_bam_length_similarity_filter.py mapped_genes.sorted.bam mapped_genes.filtered.bam
python ../SCRIPT/_8_bam_length_similarity_filter.py mapped_TEonly.sorted.bam mapped_TEonly.filtered.bam
samtools index mapped_genes.filtered.bam
samtools index mapped_TEonly.filtered.bam
Rscript ../SCRIPT/_8_repetitiveness_score_.R mapped_TEonly.filtered.bam mapped_genes.filtered.bam ../6_SCANNING_RESULT/all_putative_ETEs_ids.txt
cd ..
```

### Step 7. Test 24-nt siRNA targeting
```bash
conda activate ete
mkdir -p 9_SIRNA_VALIDATION
cd 9_SIRNA_VALIDATION
cp ../0_STARTING_FILES/siRNA24_filtered_unique_DNA.fasta .
bwa index ../6_SCANNING_RESULT/GENE_putative_ete.fasta
bwa index ../6_SCANNING_RESULT/MRNA_putative_ete.fasta
bwa aln -t $T ../6_SCANNING_RESULT/GENE_putative_ete.fasta siRNA24_filtered_unique_DNA.fasta > sirna_vs_ETE_genes.sai
bwa samse -n 500 ../6_SCANNING_RESULT/GENE_putative_ete.fasta sirna_vs_ETE_genes.sai siRNA24_filtered_unique_DNA.fasta > sirna_vs_ETE_genes.sam
bwa aln -t $T ../6_SCANNING_RESULT/MRNA_putative_ete.fasta siRNA24_filtered_unique_DNA.fasta > sirna_vs_ETE_exons.sai
bwa samse -n 500 ../6_SCANNING_RESULT/MRNA_putative_ete.fasta sirna_vs_ETE_exons.sai siRNA24_filtered_unique_DNA.fasta > sirna_vs_ETE_exons.sam
python ../SCRIPT/_9_score_ete_sirna.py --ete-list ../6_SCANNING_RESULT/all_putative_ETEs_ids.txt --sam-ete-genes sirna_vs_ETE_genes.sam --sam-ete-exons sirna_vs_ETE_exons.sam -o sirna_score.tsv
cd ..
```

### Step 8. Assess expression from RNA-seq data
`bwa mem` keeps the original thread count (`-t 22`) for the same reason as in Step 6.
```bash
conda activate ete
mkdir -p 10_EXPRESSION_VALIDATION/ref 10_EXPRESSION_VALIDATION/2_mapped 10_EXPRESSION_VALIDATION/3_filtered 10_EXPRESSION_VALIDATION/counts
cd 10_EXPRESSION_VALIDATION
cp ../0_STARTING_FILES/fine_cleaned_cds.fasta ref/Chr1A_cds.fasta
bwa index ref/Chr1A_cds.fasta
samtools faidx ref/Chr1A_cds.fasta
for r1 in ../0_STARTING_FILES/RNA_READS/*_1.500k.fastq.gz; do
  r2=${r1%_1.500k.fastq.gz}_2.500k.fastq.gz
  lib=$(basename ${r1%_1.500k.fastq.gz})
  bwa mem -t 22 ref/Chr1A_cds.fasta $r1 $r2 | samtools view -b -F 4 - | samtools sort -@ $T -o 2_mapped/$lib.sorted.bam -
  samtools index 2_mapped/$lib.sorted.bam
  python ../SCRIPT/_8_bam_length_similarity_filter.py 2_mapped/$lib.sorted.bam 3_filtered/$lib.filtered.bam
  samtools index 3_filtered/$lib.filtered.bam
  samtools view -F 0x904 -q 10 3_filtered/$lib.filtered.bam | cut -f3 | sort | uniq -c | awk 'BEGIN{OFS="\t"} $2!="*"{print $2,$1}' > counts/$lib.counts.tsv
done
Rscript ../SCRIPT/_10_build_rpkm_table.R counts ref/Chr1A_cds.fasta.fai expression_all_cds_bwa.csv
Rscript ../SCRIPT/_10_score_rpkm_table.R ../6_SCANNING_RESULT/GENE_putative_ete.fasta expression_all_cds_bwa.csv expression_result.tsv
cd ..
```

### Step 9. Validate putative ETEs using the cumulative score
```bash
conda activate ete
mkdir -p 11_RESULTS
cd 11_RESULTS
Rscript ../SCRIPT/_11_final_score_and_check_aggregator.R ..
cd ..
```

### Step 10. Characterize validated ETE
```bash
cd 11_RESULTS
bash ../SCRIPT/_11_te_derived_fusion.sh
cd ..
```



**Pipeline results:**

| ETE | Number |
|---|---:|
| Total class 4 loci | 219 |
| TE-derived host genes | 6 |
| TE–gene fusions | 213 |

**Score distribution across ETE loci:**

| Score | Number of loci |
|---:|---:|
| -4 | 0 |
| -3 | 2 |
| -2 | 13 |
| -1 | 24 |
| 0 | 79 |
| 1 | 57 |
| 2 | 231 |
| 3 | 46 |
| 4 | 219 |
| **Total** | **671** |
