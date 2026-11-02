# 0INFO: run the pipeline

You need only this file, `SCRIPT/` and `0_STARTING_FILES/`.

Rules for every step below:
- Copy the commands **in order**, into **one** terminal (Linux, or Ubuntu under WSL2 on Windows).
- Run them **from the folder that contains this file**.
- Start from a folder that contains nothing else (no old result folders).
- On a remote server, run inside `screen` so that a dropped connection does not stop the run.

The detailed guide (verification status, troubleshooting, demonstration) is `docs/GUIDE_DETAILED.md`.

## 0. Input files
`0_STARTING_FILES/` must contain:
```
Chr1A.fasta
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

## 1. Install conda (skip if `conda --version` works)
```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3
~/miniforge3/bin/conda init bash
exec bash
```

## 2. Create the environments (once per computer, one command at a time)
```bash
conda create -y -n dante -c conda-forge -c bioconda -c petrnovak dante=0.2.12
conda create -y -n dante_ltr -c conda-forge -c bioconda -c petrnovak dante_ltr=0.6.4.0
conda create -y -n dante_tir -c conda-forge -c r -c bioconda -c petrnovak dante_tir=0.3.2
conda create -y -n ete -c conda-forge -c bioconda mafft=7.526 hmmer==3.3 bwa=0.7.19 samtools seqkit=2.1.0 blast python=3.12 pysam scipy r-base r-tidyverse bioconductor-rtracklayer bioconductor-biostrings bioconductor-genomicalignments
```

## 3. Settings (repeat in every new terminal)
`T` is the number of CPU threads. DANTE_LTR uses up to about 1 GB of RAM per thread, so use `T=2` on a PC with little memory.
```bash
export LC_ALL=C
T=10
```

## 4. TE annotation: DANTE, DANTE_LTR, DANTE_TIR
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
dante_tir.py -g /home/samuel/Chr1A_github/1_DANTE_suite/DANTE/DANTE_filtered_output.gff3 -f /home/samuel/Chr1A_github/0_STARTING_FILES/Chr1A.fasta -o TIR -c $T
cd ..
```
The pipeline uses `1_DANTE_suite/DANTE_LTR/ltr_library/TE_DLplus.fasta` and `1_DANTE_suite/TIR/DANTE_TIR_final.*`.

`dante_ltr_to_library` ends with an mmseqs2 clustering that needs several GB of RAM. If only that last step fails, `TE_DLplus.fasta` is already complete.

## 5. Domain extraction and alignment
```bash
conda activate ete
mkdir -p 2_EXTRACTED_DOMAINS
cd 2_EXTRACTED_DOMAINS
Rscript ../SCRIPT/_1_LTR_protein_domain_concatenator.R ../1_DANTE_suite/DANTE_LTR/DANTE_LTR.gff3
Rscript ../SCRIPT/_2_TIR_tpase_domain_extractor.R ../1_DANTE_suite/TIR/DANTE_TIR_final.gff3 ../1_DANTE_suite/DANTE/DANTE_filtered_output.gff3
cd ..

mkdir -p 3_ALIGNMENT
cd 3_ALIGNMENT
bash ../SCRIPT/_3_domain_aligner.sh ../2_EXTRACTED_DOMAINS/ $T
cd ..
```

## 6. HMMs, proteome scan and putative ETEs
```bash
conda activate ete
mkdir -p 4_ALIGNMENT_PARSING_HMM_BUILD
cd 4_ALIGNMENT_PARSING_HMM_BUILD
Rscript ../SCRIPT/_4_sliding_and_hmm_build_overhang_win.R ../3_ALIGNMENT $T
cd ..

mkdir -p 5_HMM_PROTEOME_SCANNING
cd 5_HMM_PROTEOME_SCANNING
bash ../SCRIPT/_5_massive_hmm_search.sh ../0_STARTING_FILES/fine_cleaned_pep.fasta ../4_ALIGNMENT_PARSING_HMM_BUILD $T
cd ..

Rscript SCRIPT/_6_domtblout_concatenator_filter_v2.R -d 5_HMM_PROTEOME_SCANNING -p 0_STARTING_FILES/fine_cleaned_pep.fasta -c 0_STARTING_FILES/fine_cleaned_cds.fasta -g 0_STARTING_FILES/fine_cleaned_genes.fasta -m 0_STARTING_FILES/fine_cleaned_mRNA.fasta -o 6_SCANNING_RESULT
grep '>' 6_SCANNING_RESULT/GENE_putative_ete.fasta | sed 's/^>//; s/ .*//' | sort -u > 6_SCANNING_RESULT/all_putative_ETEs_ids.txt
wc -l 6_SCANNING_RESULT/all_putative_ETEs_ids.txt
```

## 7. Validation 1: similarity to known TEs
```bash
conda activate ete
mkdir -p 7_TE_SIMILARITY_VALIDATION
cd 7_TE_SIMILARITY_VALIDATION
cat ../1_DANTE_suite/DANTE_LTR/ltr_library/TE_DLplus.fasta ../1_DANTE_suite/DANTE_TIR/TIR/DANTE_TIR_final.fasta > all_TE.fasta
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

## 8. Validation 2: siRNA targeting
```bash
conda activate ete
mkdir -p 8_SIRNA_VALIDATION
cd 8_SIRNA_VALIDATION
cp ../0_STARTING_FILES/siRNA24_filtered_unique_DNA.fasta ../0_STARTING_FILES/fine_cleaned_genes.fasta .
bwa index fine_cleaned_genes.fasta
bwa index ../6_SCANNING_RESULT/GENE_putative_ete.fasta
bwa index ../6_SCANNING_RESULT/MRNA_putative_ete.fasta
bwa aln -t $T fine_cleaned_genes.fasta siRNA24_filtered_unique_DNA.fasta > sirna_vs_all_genes.sai
bwa samse -n 500 fine_cleaned_genes.fasta sirna_vs_all_genes.sai siRNA24_filtered_unique_DNA.fasta > sirna_vs_all_genes.sam
bwa aln -t $T ../6_SCANNING_RESULT/GENE_putative_ete.fasta siRNA24_filtered_unique_DNA.fasta > sirna_vs_ETE_genes.sai
bwa samse -n 500 ../6_SCANNING_RESULT/GENE_putative_ete.fasta sirna_vs_ETE_genes.sai siRNA24_filtered_unique_DNA.fasta > sirna_vs_ETE_genes.sam
bwa aln -t $T ../6_SCANNING_RESULT/MRNA_putative_ete.fasta siRNA24_filtered_unique_DNA.fasta > sirna_vs_ETE_exons.sai
bwa samse -n 500 ../6_SCANNING_RESULT/MRNA_putative_ete.fasta sirna_vs_ETE_exons.sai siRNA24_filtered_unique_DNA.fasta > sirna_vs_ETE_exons.sam
python ../SCRIPT/_8_score_ete_sirna.py --ete-list ../6_SCANNING_RESULT/all_putative_ETEs_ids.txt --sam-ete-genes sirna_vs_ETE_genes.sam --sam-ete-exons sirna_vs_ETE_exons.sam --sam-all-genes sirna_vs_all_genes.sam --all-genes-fasta fine_cleaned_genes.fasta -o sirna_score.tsv
cd ..
```

## 9. Validation 3: repetitiveness (genomic reads)
`bwa mem` keeps the original thread count (`-t 28`): its batch size depends on the number of threads and slightly changes the paired-end mapping.
```bash
conda activate ete
mkdir -p 9_REPETITIVENESS_VALIDATION
cd 9_REPETITIVENESS_VALIDATION
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
python ../SCRIPT/_9_bam_length_similarity_filter.py mapped_genes.sorted.bam mapped_genes.filtered.bam
python ../SCRIPT/_9_bam_length_similarity_filter.py mapped_TEonly.sorted.bam mapped_TEonly.filtered.bam
samtools index mapped_genes.filtered.bam
samtools index mapped_TEonly.filtered.bam
Rscript ../SCRIPT/_9_repetitiveness_score_.R mapped_TEonly.filtered.bam mapped_genes.filtered.bam ../6_SCANNING_RESULT/all_putative_ETEs_ids.txt
cd ..
```

## 10. Validation 4: expression (RNA-seq)
`bwa mem` keeps the original thread count (`-t 22`) for the same reason as in step 9.
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
  python ../SCRIPT/_9_bam_length_similarity_filter.py 2_mapped/$lib.sorted.bam 3_filtered/$lib.filtered.bam
  samtools index 3_filtered/$lib.filtered.bam
  samtools view -F 0x904 -q 10 3_filtered/$lib.filtered.bam | cut -f3 | sort | uniq -c | awk 'BEGIN{OFS="\t"} $2!="*"{print $2,$1}' > counts/$lib.counts.tsv
done
Rscript ../SCRIPT/_10_build_rpkm_table.R counts ref/Chr1A_cds.fasta.fai expression_all_cds_bwa.csv
Rscript ../SCRIPT/_10_score_rpkm_table.R ../6_SCANNING_RESULT/GENE_putative_ete.fasta expression_all_cds_bwa.csv expression_result.tsv
cd ..
```

## 11. Final results
```bash
conda activate ete
mkdir -p 11_RESULTS
cd 11_RESULTS
Rscript ../SCRIPT/_11_final_score_and_check_aggregator.R ..
bash ../SCRIPT/_11_te_derived_fusion.sh
cd ..
```
Output, in `11_RESULTS/`:
- `score_per_locus.tsv`: one row per ETE locus;
- `risultato_finale_validazione_filtrabile.tsv`: one row per isoform;
- `classe4_te_derived_fusion.tsv`: class-4 split.

Column meanings are in `docs/DATA_DICTIONARY.md`.

**Published run, for comparison:**

| Result | Value |
|---|---|
| DANTE filtered domains | 129,180 |
| DANTE_LTR elements | 30,439 |
| TIR elements | 176 |
| window HMMs | 522 |
| putative ETE loci (step 6) | 671 |
| `SCORE_SUM = 4` (class 4) | 219: 6 TE-derived host genes, 213 TE–gene fusions |
