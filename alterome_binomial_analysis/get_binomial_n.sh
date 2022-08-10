#!/bin/bash -l

# qsub -P cancergrp -N submit -j y get_binomial_n.sh 

module load bedtools

pdc=../../../noncoding_cancer/data/icgc_pdc/promoter_SNV
collab=../../../noncoding_cancer/data/icgc_collab/promoter_SNV
gencode=../../../noncoding_cancer/data/gencode/no_cds_promoter_v19_hs37.bed

counter=0
for file in $pdc/*.vcf;do
	fname=$(basename "$file" .vcf)
	bedtools intersect -wa -wb -a $gencode -b $file > bed_promoter_ids/all_samples_$fname.bed
	counter=$((counter+1))
	echo $counter
done

counter=0
for file in $collab/*.vcf;do
	fname=$(basename "$file" .vcf)
	bedtools intersect -wa -wb -a $gencode -b $file > bed_promoter_ids/all_samples_$fname.bed
	counter=$((counter+1))
	echo $counter
done

cat bed_promoter_ids/all_samples_* | cut -f 8 | sort | uniq -c | sort -nr > number_variants_per_promoter.txt

