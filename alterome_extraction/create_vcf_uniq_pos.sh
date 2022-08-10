#!/bin/bash -l
mkdir -p vcf_unique_pos
counter=0
for f in promoter_SNV/*;do
	counter=$((counter+1))
	echo $counter
	fname=$(basename $f)
	bedtools merge -i $f | awk 'BEGIN {OFS="\t"}; {$2 = $2 +1; print $0};' > tmp
	bedtools intersect -a vcf/$fname -b tmp > vcf_unique_pos/$fname
done

rm tmp
