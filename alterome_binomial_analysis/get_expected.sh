#!/bin/bash -l

module load R/3.4.3


# qsub -P cancergrp -l h_rt=12:00:00 -j y -N expected get_expected.sh

counter=0
for f in ../probabilities/*;do

	Rscript get_expected.R $f number_variants_per_promoter.txt 1 expected_mutations_by_pwm/loss_1
	Rscript get_expected.R $f number_variants_per_promoter.txt 2 expected_mutations_by_pwm/loss_2
	Rscript get_expected.R $f number_variants_per_promoter.txt 3 expected_mutations_by_pwm/loss_3
	
	Rscript get_expected.R $f number_variants_per_promoter.txt 4 expected_mutations_by_pwm/gain_1
	Rscript get_expected.R $f number_variants_per_promoter.txt 5 expected_mutations_by_pwm/gain_2
	Rscript get_expected.R $f number_variants_per_promoter.txt 6 expected_mutations_by_pwm/gain_3
	
	counter=$((counter+1))
	echo $counter
	
done

cat expected_mutations_by_pwm/loss_1/* > expected_mutations_by_pwm/loss_1.txt
cat expected_mutations_by_pwm/loss_2/* > expected_mutations_by_pwm/loss_2.txt
cat expected_mutations_by_pwm/loss_3/* > expected_mutations_by_pwm/loss_3.txt

cat expected_mutations_by_pwm/gain_1/* > expected_mutations_by_pwm/gain_1.txt
cat expected_mutations_by_pwm/gain_2/* > expected_mutations_by_pwm/gain_2.txt
cat expected_mutations_by_pwm/gain_3/* > expected_mutations_by_pwm/gain_3.txt

