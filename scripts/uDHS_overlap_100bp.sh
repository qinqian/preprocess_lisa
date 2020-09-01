#!/bin/bash

## TODO!!!: create unionDHS and TFBS union later to redo this
intersectBed -wa -u -a ../PhaseA_motif_deltarp/hg38_window100bp.bed  -b ~/12_data/human_unionDHS_fc5_75merge_split.bed | cut -f 4 > hg38_100bp_uDHS.bed.index
python -c "import numpy; a=numpy.loadtxt('hg38_100bp_uDHS.bed.index', dtype='int32'); numpy.save('hg38_100bp_uDHS.bed.index', a)"

#for i in *peaks*bed
#do 
#   echo $i
#   intersectBed -wa -u -a ../PhaseA_motif_deltarp/hg38_window100bp.bed  -b $i | cut -f 4 > ${i}.index
#   python -c "import numpy; a=numpy.loadtxt('${i}.index', dtype='int32'); numpy.save('${i}.index', a)"
#done
