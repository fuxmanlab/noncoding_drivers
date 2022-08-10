#!/bin/bash -l

# qsub -P cancergrp -l h_rt=06:00:00 -j y -N pwm_promoter -t 1-1963 get_pwm_promoter_matrix.sh

pdc=../../../noncoding_cancer/data/icgc_pdc/promoter_alterome_table
collab=../../../noncoding_cancer/data/icgc_collab/promoter_alterome_table

cat $pdc/* $collab/* | awk '{if ($16==1) print}' | awk -v tf=$SGE_TASK_ID '{if ($4==tf) print}' | cut -f 8 | sort | uniq -c | sort -nr > pwm_promoter_matrix/loss_1/$SGE_TASK_ID.txt

cat $pdc/* $collab/* | awk '{if ($17==1) print}' | awk -v tf=$SGE_TASK_ID '{if ($4==tf) print}' | cut -f 8 | sort | uniq -c | sort -nr > pwm_promoter_matrix/loss_2/$SGE_TASK_ID.txt

cat $pdc/* $collab/* | awk '{if ($18==1) print}' | awk -v tf=$SGE_TASK_ID '{if ($4==tf) print}' | cut -f 8 | sort | uniq -c | sort -nr > pwm_promoter_matrix/loss_3/$SGE_TASK_ID.txt

cat $pdc/* $collab/* | awk '{if ($19==1) print}' | awk -v tf=$SGE_TASK_ID '{if ($4==tf) print}' | cut -f 8 | sort | uniq -c | sort -nr > pwm_promoter_matrix/gain_1/$SGE_TASK_ID.txt

cat $pdc/* $collab/* | awk '{if ($20==1) print}' | awk -v tf=$SGE_TASK_ID '{if ($4==tf) print}' | cut -f 8 | sort | uniq -c | sort -nr > pwm_promoter_matrix/gain_2/$SGE_TASK_ID.txt

cat $pdc/* $collab/* | awk '{if ($21==1) print}' | awk -v tf=$SGE_TASK_ID '{if ($4==tf) print}' | cut -f 8 | sort | uniq -c | sort -nr > pwm_promoter_matrix/gain_3/$SGE_TASK_ID.txt
