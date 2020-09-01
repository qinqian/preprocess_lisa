#!/bin/bash

ann=/data/home/qqin/01_Projects/Programming/dc2/scripts/samples_directory.tab


#hg38
for i in `cut -f 1 ../best_dc_tfcr_basedon_frip_peak_dhs.xls | sed 1d`
#mm10
#for i in `cut -f 1 ../best_dc_tfcr_basedon_frip_peak_dhs_mm.xls | sed 1d`
do
  #/data5/DC_results/Result_new/dataset37126/attic/37126_gene_score_5fold.txt
  summmit=$(ls $(grep "^$i"$'\t' $ann | cut -f 2)/attic/*5fold*txt)
  cp $summmit .
done
