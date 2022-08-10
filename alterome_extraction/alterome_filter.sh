#!/bin/bash -l

# qsub -P cancergrp -l h_rt=06:00:00 -j y -N alterome -t 1-828 scr/alterome_filter.sh

module load R/3.5.2 

# Intersection
Rscript scr/alterome_filter.R promoter_alterome_untar/bed_files ../../../tfbs_profile/tf_info_scored.csv promoter_alterome_table $SGE_TASK_ID

