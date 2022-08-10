#!/bin/bash -l



for d in pwm_promoter_matrix/*;do
	for f in $d/*;do
		awk '{sum+=$1} END {print sum}' $f >> tmp_x
		echo $(basename "$f" .txt) >> tmp_file
	done
	outName=$(basename "$d")
	paste tmp_file tmp_x > observed_mutations_by_pwm/$outName.txt
	rm -rf tmp_*
done
