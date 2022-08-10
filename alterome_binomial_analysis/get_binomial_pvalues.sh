#!/bin/bash -l

# qsub -P cancergrp -l h_rt=12:00:00 -j y -N pvalues get_binomial_pvalues.sh

module load R/3.5.2


echo $file
echo $counter
for e in $file/*.txt;do
	outName=$(basename "$file")
	pwmN=$(basename $e .txt)
	pwm=${pwmN##*_}.csv
	for f in ../probabilities/pwm_${pwm};do
		echo "Rscript get_binomial.R $f number_variants_per_promoter.txt $e $counter pvalues/$outName"
		Rscript get_binomial.R $f number_variants_per_promoter.txt $e $counter pvalues/$outName
	done
done
