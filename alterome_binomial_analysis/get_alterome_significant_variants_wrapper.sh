#!/bin/bash -l

cancer=$1

split -l 100 filtered_fdr/intersection/${cancer}_gain_intersection.txt ${cancer}_gain_split
split -l 100 filtered_fdr/intersection/${cancer}_loss_intersection.txt ${cancer}_loss_split

mv ${cancer}_gain_split* filtered_fdr/pwm_promoter/
mv ${cancer}_loss_split* filtered_fdr/pwm_promoter/

#echo "pdc alterome"
#cat /restricted/projectnb/cancergrp/noncoding_cancer/data/icgc_pdc/promoter_alterome_table/* > filtered_fdr/alterome.bed
#echo "collab alterome"
#cat /restricted/projectnb/cancergrp/noncoding_cancer/data/icgc_collab/promoter_alterome_table/* >> filtered_fdr/alterome.bed


for f in filtered_fdr/pwm_promoter/${cancer}_*;do
	bf=$(basename "$f")
	echo $bf
	qsub -P cancergrp -l h_rt=12:00:00 -j y -N $bf -v file=$bf,cancer=$cancer get_alterome_significant_variants.sh
done
