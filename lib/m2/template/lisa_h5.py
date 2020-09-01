import h5py, _bw
import numpy as np
import os
from pkg_resources import resource_filename
import gzip

meta = resource_filename("_bw", "%s.tss" % (snakemake.params.sp))
_bw.getrp(snakemake.input.bw, meta, snakemake.output.rp, 1e4, 0, 0)

# single sample RP to HDF5
tmp = []; n = 0
with open(snakemake.output.rp) as inf:
    for line in inf:
        n += 1
        line = line.split()[4]
        tmp.append(float(line.strip()))
    tmp = np.array(tmp, dtype='float32')
    with h5py.File(snakemake.output.h5, 'a') as store:
        RP = store.create_dataset("RP", dtype=np.float32, shape=(n, 1), compression='gzip', shuffle=True, fletcher32=True)
        RP[:,0] = tmp
        store.flush()

os.system('bigWigAverageOverBed %s %s stdout > %s' % (snakemake.input.bw, snakemake.params.bin1kb, snakemake.output.ct))

tmp = []
n = 0
with open(snakemake.output.ct) as inf:
    for line in inf:
        n += 1
        line = line.split()[4] # average on all bases!
        tmp.append(float(line.strip()))
    target = np.array(tmp, dtype=np.float32)

# reorder
index = []
with gzip.open(snakemake.params.order, 'rb') as sortedf:
    for line in sortedf:
        index.append(int(line.strip())-1)

target = target[np.argsort(index)]

with h5py.File(snakemake.output.h5, "a") as store:
    ct = store.create_dataset("OrderCount", dtype=np.float32, shape=(n, 1), compression='gzip', shuffle=True, fletcher32=True)
    ct[:,0] = target
    store.flush()
