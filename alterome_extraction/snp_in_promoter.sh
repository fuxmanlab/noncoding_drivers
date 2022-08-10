#!/bin/bash -l

# usage:
# sh scr/snp_in_promoter.sh vcf

module load bedtools

promoter=/restricted/projectnb/cancergrp/noncoding_cancer/data/gencode/no_cds_promoter_v19_hs37.bed 
vcfDir=$1

mkdir -p promoter_SNV
mkdir -p promoter_SNV_uniq

counter=0
for file in $vcfDir/*.vcf;do
	bedtools intersect -wa -wb -a $file -b $promoter > tmp.vcf
	grep "#" $file > header.txt
	name=$(basename $file )
	cat header.txt tmp.vcf > promoter_SNV/$name
	Rscript scr/get_unique_vcf.R tmp.vcf tmp_uniq
	cat header.txt tmp_uniq > promoter_SNV_uniq/$name
	counter=$((counter+1))
	echo $counter $name
done
	
rm header.txt
rm tmp.vcf
rm tmp_uniq
