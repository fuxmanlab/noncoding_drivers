#!/bin/bash -l

mkdir -p ../cancers
rm -rf ../cancers/$fname/freq.txt
for f in /restricted/projectnb/cancergrp/noncoding_cancer/data/pcawg/samples_cancer/*;do
	fname=$(basename "$f" .txt )
	echo $fname
	mkdir -p ../cancers/$fname
	mkdir -p ../cancers/$fname/promoter_alterome_table/
	mkdir -p ../cancers/$fname/promoter_SNV
	mkdir -p ../cancers/$fname/promoter_SNV_uniq
	while IFS=, read -r col1 col2
	do
	    cp ../promoter_alterome_table/$col1* ../cancers/$fname/promoter_alterome_table
		cp ../promoter_SNV/$col1* ../cancers/$fname/promoter_SNV/
		cp ../promoter_SNV_uniq/$col1* ../cancers/$fname/promoter_SNV_uniq/
	done < $f
	grep -v "#" ../cancers/$fname/promoter_SNV_uniq/* | LC_ALL=C grep -P "A\tC"  | wc -l >> ../cancers/$fname/freq.txt
	grep -v "#" ../cancers/$fname/promoter_SNV_uniq/* | LC_ALL=C grep -P "A\tG"  | wc -l >> ../cancers/$fname/freq.txt
	grep -v "#" ../cancers/$fname/promoter_SNV_uniq/* | LC_ALL=C grep -P "A\tT"  | wc -l >> ../cancers/$fname/freq.txt
	grep -v "#" ../cancers/$fname/promoter_SNV_uniq/* | LC_ALL=C grep -P "C\tA"  | wc -l >> ../cancers/$fname/freq.txt
	grep -v "#" ../cancers/$fname/promoter_SNV_uniq/* | LC_ALL=C grep -P "C\tG"  | wc -l >> ../cancers/$fname/freq.txt
	grep -v "#" ../cancers/$fname/promoter_SNV_uniq/* | LC_ALL=C grep -P "C\tT"  | wc -l >> ../cancers/$fname/freq.txt
	grep -v "#" ../cancers/$fname/promoter_SNV_uniq/* | LC_ALL=C grep -P "G\tA"  | wc -l >> ../cancers/$fname/freq.txt
	grep -v "#" ../cancers/$fname/promoter_SNV_uniq/* | LC_ALL=C grep -P "G\tC"  | wc -l >> ../cancers/$fname/freq.txt
	grep -v "#" ../cancers/$fname/promoter_SNV_uniq/* | LC_ALL=C grep -P "G\tT"  | wc -l >> ../cancers/$fname/freq.txt
	grep -v "#" ../cancers/$fname/promoter_SNV_uniq/* | LC_ALL=C grep -P "T\tA"  | wc -l >> ../cancers/$fname/freq.txt
	grep -v "#" ../cancers/$fname/promoter_SNV_uniq/* | LC_ALL=C grep -P "T\tC"  | wc -l >> ../cancers/$fname/freq.txt
	grep -v "#" ../cancers/$fname/promoter_SNV_uniq/* | LC_ALL=C grep -P "T\tG"  | wc -l >> ../cancers/$fname/freq.txt
done
