#!/bin/bash -l

module load bedtools
module load R/3.4.3

promoters=../../noncoding_cancer/data/gencode/no_cds_promoter_sequences_v19_hs37.bed


awk -v pwm=$pwm '{ if ($5 == pwm ) {print} }' p_alterome/chr*_$chr.bed > tmp_${pwm}_${cancer}
bedtools intersect -wa -wb -a tmp_${pwm}_${cancer} -b ../../noncoding_cancer/data/gencode/no_cds_promoter_v19_hs37.bed > pwm_${pwm}_${cancer}.bed
Rscript pwm_p_calculation_cancers.R pwm_${pwm}_${cancer}.bed ../tf_info_scored.csv $promoters freq_cancers/$cancer.csv probabilities_cancers/$cancer

rm tmp_${pwm}_${cancer}
rm pwm_${pwm}_${cancer}.bed
