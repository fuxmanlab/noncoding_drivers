#!/bin/bash -l

module load R/3.5.1

# qsub -P cancergrp -l h_rt=12:00:00 -v cancer=biliary -j y -N expected get_expected_cancers.sh

counter=0

mkdir -p expected_mutations_by_pwm_cancers/$cancer

mkdir -p expected_mutations_by_pwm_cancers/$cancer/loss_1
mkdir -p expected_mutations_by_pwm_cancers/$cancer/loss_2
mkdir -p expected_mutations_by_pwm_cancers/$cancer/loss_3

mkdir -p expected_mutations_by_pwm_cancers/$cancer/gain_1
mkdir -p expected_mutations_by_pwm_cancers/$cancer/gain_2
mkdir -p expected_mutations_by_pwm_cancers/$cancer/gain_3

n=binomial_n_cancers/${cancer}.txt

for f in ../probabilities_cancers/$cancer/*;do

	Rscript get_expected_cancers.R $f $n 1 expected_mutations_by_pwm_cancers/$cancer/loss_1
	Rscript get_expected_cancers.R $f $n 2 expected_mutations_by_pwm_cancers/$cancer/loss_2
	Rscript get_expected_cancers.R $f $n 3 expected_mutations_by_pwm_cancers/$cancer/loss_3
	
	Rscript get_expected_cancers.R $f $n 4 expected_mutations_by_pwm_cancers/$cancer/gain_1
	Rscript get_expected_cancers.R $f $n 5 expected_mutations_by_pwm_cancers/$cancer/gain_2
	Rscript get_expected_cancers.R $f $n 6 expected_mutations_by_pwm_cancers/$cancer/gain_3
	
	counter=$((counter+1))
	echo $counter
	
done

cat expected_mutations_by_pwm_cancers/$cancer/loss_1/* > expected_mutations_by_pwm_cancers/${cancer}_loss_1.txt
cat expected_mutations_by_pwm_cancers/$cancer/loss_2/* > expected_mutations_by_pwm_cancers/${cancer}_loss_2.txt
cat expected_mutations_by_pwm_cancers/$cancer/loss_3/* > expected_mutations_by_pwm_cancers/${cancer}_loss_3.txt

cat expected_mutations_by_pwm_cancers/$cancer/gain_1/* > expected_mutations_by_pwm_cancers/${cancer}_gain_1.txt
cat expected_mutations_by_pwm_cancers/$cancer/gain_2/* > expected_mutations_by_pwm_cancers/${cancer}_gain_2.txt
cat expected_mutations_by_pwm_cancers/$cancer/gain_3/* > expected_mutations_by_pwm_cancers/${cancer}_gain_3.txt

