#!/bin/bash -l


count=0
while read f1 f2 f3 f4 f5 f6; do
	count=$(( $count + 1 ))
	echo $count
    awk -v promoter=$f4 -v tf=$f5 '{ if (($4==tf) && ($8==promoter))  print }' filtered_fdr/alterome.bed >> filtered_fdr/variants/split_files/$file
done < filtered_fdr/pwm_promoter/$file

