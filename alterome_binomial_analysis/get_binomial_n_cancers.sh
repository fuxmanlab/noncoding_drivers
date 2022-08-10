#!/bin/bash -l

# sh get_binomial_n_cancers.sh biliary

module load bedtools

cancer=$1

pdc=../../../noncoding_cancer/data/icgc_pdc/cancers/$cancer/promoter_SNV
collab=../../../noncoding_cancer/data/icgc_collab/cancers/$cancer/promoter_SNV
gencode=../../../noncoding_cancer/data/gencode/no_cds_promoter_v19_hs37.bed

counter=0
for file in $pdc/*.vcf;do
	fname=$(basename "$file" .vcf)
	bedtools intersect -wa -wb -a $gencode -b $file > bed_promoter_ids/${cancer}_$fname.bed
	counter=$((counter+1))
	echo $counter
done

counter=0
for file in $collab/*.vcf;do
	fname=$(basename "$file" .vcf)
	bedtools intersect -wa -wb -a $gencode -b $file > bed_promoter_ids/${cancer}_$fname.bed
	counter=$((counter+1))
	echo $counter
done

cat bed_promoter_ids/${cancer}_* | cut -f 8 | sort | uniq -c | sort -nr > binomial_n_cancers/$cancer.txt

