#!/bin/bash -l

# use binomial

module load R/3.5.2


echo $file
echo $counter

outName=$(basename "$file")
mkdir -p pvalues_cancers/${cancer}_$outName

for e in $file/*.txt;do
	outName=$(basename "$file")
	pwmN=$(basename $e .txt)
	pwm=${pwmN##*_}.csv
	for f in ../probabilities_cancers/$cancer/pwm_${pwm};do
		echo "Rscript get_binomial.R $f binomial_n_cancers/${cancer}.txt $e $counter pvalues/${cancer}_$outName"
		Rscript get_binomial.R $f binomial_n_cancers/${cancer}.txt $e $counter pvalues_cancers/${cancer}_$outName
	done
done

cat pvalues_cancers/${cancer}_$outName/*.out > pvalues_cancers/${cancer}_${outName}.txt
