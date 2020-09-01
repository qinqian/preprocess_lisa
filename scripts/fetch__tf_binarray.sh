#!/bin/bash

# a peak or summit bed file
peak=$1
window=../../PhaseA_motif_deltarp/hg38_window100bp.bed

bedtools intersect -wa -u -b $peak -a $window | cut -f 4 > ${peak}.index
python -c "import numpy; a=numpy.loadtxt('${peak}.index', dtype='int32'); numpy.save('${peak}.index.bin', a)"
