import h5py
import os
import numpy as np

index = []
# human
#with open('bigWigAverageOverBed.originalorder.hg') as sortedf:
# mouse
with open('bigWigAverageOverBed.originalorder.mm') as sortedf:
     for line in sortedf:
          index.append(int(line.strip())-1)

# human
#with h5py.File(os.path.expanduser("~/12_data/MARGE/hg38_window1kb_H3K27ac_averageonallbase.h5")) as store:
#with h5py.File("/data/home/qqin/12_data/hg38_window1kb_H3K4me3.h5") as store:
#with h5py.File("/data/home/qqin/12_data/hg38_window1kb_DNase.h5") as store:
#with h5py.File("/data/home/qqin/12_data/hg38_window1kb_H3K4me1.h5") as store:
#with h5py.File('/data/home/qqin/12_data/hg38_window1kb_H3K4me2.h5') as store:
#with h5py.File('/data/home/qqin/12_data/hg38_window1kb_H3K27me3.h5') as store:
#with h5py.File('/data/home/qqin/12_data/hg38_window1kb_H3K36me3.h5') as store:
#with h5py.File('/data/home/qqin/12_data/hg38_window1kb_H3K9me3.h5') as store:
#with h5py.File('/data/home/qqin/12_data/hg38_window1kb_H3K9ac.h5') as store:
# mouse
# with h5py.File("/data/home/qqin/12_data/mm10_window1kb_H3K27ac.h5") as store:
# with h5py.File("/data/home/qqin/12_data/mm10_window1kb_DNase.h5") as store:
# with h5py.File("/data/home/qqin/12_data/mm10_window1kb_H3K4me3.h5") as store:
# with h5py.File("/data/home/qqin/12_data/mm10_window1kb_H3K4me1.h5") as store:
#with h5py.File("/data/home/qqin/12_data/mm10_window1kb_H3K4me2.h5") as store:
#with h5py.File("/data/home/qqin/12_data/mm10_window1kb_H3K27me3.h5") as store:
#with h5py.File("/data/home/qqin/12_data/mm10_window1kb_H3K9me3.h5") as store:
#with h5py.File("/data/home/qqin/12_data/mm10_window1kb_H3K9ac.h5") as store:
with h5py.File("/data/home/qqin/12_data/mm10_window1kb_H3K36me3.h5") as store:
    target = store['Count'][...]
    shape = target.shape
    #del store['Count']
    #del store['OrderCount']

target = target[np.argsort(index)]
# human
#with h5py.File(os.path.expanduser("~/12_data/MARGE/hg38_window1kb_H3K27ac_averageonallbase.h5")) as store:
#with h5py.File('/data/home/qqin/12_data/hg38_window1kb_H3K4me3.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/hg38_window1kb_DNase.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/hg38_window1kb_H3K4me1.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/hg38_window1kb_H3K4me2.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/hg38_window1kb_H3K27me3.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/hg38_window1kb_H3K36me3.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/hg38_window1kb_H3K9me3.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/hg38_window1kb_H3K9ac.h5', 'a') as store:
# mouse
#with h5py.File('/data/home/qqin/12_data/mm10_window1kb_H3K27ac.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/mm10_window1kb_DNase.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/mm10_window1kb_H3K4me3.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/mm10_window1kb_H3K4me1.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/mm10_window1kb_H3K4me2.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/mm10_window1kb_H3K27me3.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/mm10_window1kb_H3K9me3.h5', 'a') as store:
#with h5py.File('/data/home/qqin/12_data/mm10_window1kb_H3K9ac.h5', 'a') as store:
with h5py.File('/data/home/qqin/12_data/mm10_window1kb_H3K36me3.h5', 'a') as store:
    store["OrderCount"] = target
