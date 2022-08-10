#!/bin/bash -l

module load bedtools
module load R/3.4.3

promoters=../../noncoding_cancer/data/gencode/no_cds_promoter_sequences_v19_hs37.bed


awk -v pwm=$pwm '{ if ($5 == pwm ) {print} }' p_alterome/chr*_$chr.bed > tmp_$pwm
bedtools intersect -wa -wb -a tmp_$pwm -b $promoters > pwm_$pwm.bed
Rscript pwm_p_calculation.R pwm_$pwm.bed ../tf_info_scored.csv $promoters number_variants_per_nucleotide_promoter_divided_by_nucleotide_fre.csv probabilities

rm tmp_$pwm
rm pwm_$pwm.bed
