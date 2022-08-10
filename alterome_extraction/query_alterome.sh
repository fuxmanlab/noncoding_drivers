#!/bin/bash -l

set -x

module load bedtools

# Definitions


time cp -v $s $TMPDIR

mkdir $TMPDIR/outputBed

sampleFile=$TMPDIR/$(basename $s)
sample=$(basename $sampleFile .bed)

# Intersection
for b in  promoter_SNV/*.vcf;do
	vcfName=$(basename $b .vcf)
	bedtools intersect -wa -wb -a $sampleFile -b $b > $TMPDIR/outputBed/${vcfName}_${sample}.tmp
	awk -v OFS="\t" '{ if ($2 == $18) { print } }'  $TMPDIR/outputBed/${vcfName}_${sample}.tmp > $TMPDIR/outputBed/$vcfName.bed
done


tar -zcvf $TMPDIR/$sample.tar $TMPDIR/outputBed/*.bed 

cp $TMPDIR/$sample.tar /restricted/projectnb/cancergrp/noncoding_cancer/data/icgc_pdc/promoter_alterome/
