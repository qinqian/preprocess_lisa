#!/usr/bin/env python

import numpy as np
import os
import h5py
import pickle

import glob
fs = glob.glob('*npy')

def load_dnase_index(f):
    idx = os.path.basename(f).split('.')[0]
    Y = np.load(f)
    return idx, Y

##with h5py.File('marge2_tf_100bp.h5', 'a') as f:
##with h5py.File('marge2_tf_100bp_mm.h5', 'a') as f:
with h5py.File('mm10_lisa_tf_100bp_all_nonhm_nonca_peak5fold.h5', 'a') as f:
#with h5py.File('hg38_lisa_tf_100bp_all_nonhm_nonca_peak5fold.h5', 'a') as f:
    ids = f.create_dataset("IDs", shape=(len(fs),), dtype='S25', compression='gzip', shuffle=True, fletcher32=True)
    for index, fi in enumerate(fs):
        i, Y = load_dnase_index(fi)
        ids[index] = str.encode(i, 'utf-8')
        """ in marge2, Y need to -1"""
        f[i] = Y

