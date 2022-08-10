#!/bin/bash -l

# qsub -P cancergrp -N submit -j y  scr/query_alterome_wrapper.sh

module load bedtools

mkdir -p promoter_alterome

alterome=/restricted/projectnb/cancergrp/tfbs_profile/promoters/p_alterome

for s in $alterome/*.bed;do
	sample=$(basename $s .vcf)
	echo "Submited job: $sample"
	qsub -P cancergrp -j y -N job_${sample} -l eth_speed=10 -v s=$s scr/query_alterome.sh
done


