#!/bin/bash -l

mkdir -p promoter_alterome_untar
for f in promoter_alterome/*.tar;do
	fname=$(basename "$f" .tar )
	tar xvzf $f
	mv scratch promoter_alterome_untar/$fname
done

ls -1 promoter_alterome_untar/chr1_001/5130784.1.y/outputBed/ > samples_untar.txt

counter=0
while read p; do
	counter=$((counter+1))
	echo $counter
 	cat promoter_alterome_untar/chr*_*/*/*/$p > promoter_alterome_untar/$p
done < samples_untar.txt
