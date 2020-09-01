import h5py
from glob import glob
import numpy as np
import operator

def parse(x, return_refseq=False):
    """
    cut -f 1,2,3,4 9_beta_peakfold5.tab | grep -v "#" | sort -k1,1d -k2,2n -k3,3n -k4,4d | uniq | less
    """
    all_refseq = []
    with open(x) as fin:
        for i in fin:
            if i.startswith('#'): continue
            elems = i.split()
            elems[1] = int(elems[1])
            elems[2] = int(elems[2])
            all_refseq.append(elems)
            elems[4] = float(elems[4])
    all_refseq.sort(key=operator.itemgetter(0,1,2,3))
    if return_refseq:
        return np.array(list(map(lambda x: str.encode(":".join([x[0], str(x[1]), str(x[2]), x[3], x[6]]), 'utf-8'), all_refseq)))
    else:
        return np.array(list(map(lambda x:x[4], all_refseq)))


fs = glob("*tab")
refseqs = parse(fs[0], True)
with h5py.File("mm10_beta_peak5fold.h5", 'a') as store:
    RP = store.create_dataset("RP", dtype=np.float32, shape=(len(refseqs), len(fs)), compression='gzip', shuffle=True, fletcher32=True)
    RefSeq = store.create_dataset("RefSeq", shape=(len(refseqs), ), dtype='S200', compression='gzip', shuffle=True, fletcher32=True)
    ids = store.create_dataset("IDs", shape=(len(fs),), dtype='S10', compression='gzip', shuffle=True, fletcher32=True)
    #99_beta_peakfold5.tab
    ids[...] = list(map(lambda x: str.encode(x.split("_")[0], 'utf-8'), fs))
    RefSeq[...] = refseqs
    store.flush()
    for i,f in enumerate(fs):
        RP[:,i] = parse(f)
    store.flush()
