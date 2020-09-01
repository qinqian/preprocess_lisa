#!/bin/bash

ann=/data/home/qqin/01_Projects/Programming/dc2/scripts/samples_directory.tab
#hg38
windows=/data/home/qqin/MARGE/PhaseA_motif_deltarp/hg38_window100bp.bed
#mm10
#windows=/data/home/qqin/seqpos2/mm10_window100bp.bed

#hg38
for i in `cut -f 1 ../best_dc_tfcr_basedon_frip_peak_dhs.xls | sed 1d`
#mm10
#for i in `cut -f 1 ../best_dc_tfcr_basedon_frip_peak_dhs_mm.xls | sed 1d`
do
  summmit=$(ls $(grep "^$i"$'\t' $ann | cut -f 2)/*summits*bed)
  echo $summmit
  if [ ! -s ${i}.index ];then
     echo $i
     bedtools intersect -wa -u -b $summmit -a $windows | cut -f 4 > ${i}.index
     python -c "import numpy; a=numpy.loadtxt('${i}.index', dtype='int32'); numpy.save('${i}.index.bin', a)"
  fi
done
