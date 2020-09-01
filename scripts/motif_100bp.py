#!/usr/bin/env python

import numpy as np
import os
import h5py
import pickle

import glob
#fs = glob.glob('cistrome/*bin.npy')
fs = glob.glob('cistrome/*100bp')
percentage_cutoff = 98
#window_map = np.load(os.path.expanduser('~/MARGE/PhaseA_motif_deltarp/hg38_100to1000window.out.npy'))

def load_motif_meta(map_motifsym):
    with open('motif_meta.txt') as inf:
        for line in inf:
            line = line.strip().split()
            sym = line[0].upper()
            map_motifsym[sym] = map_motifsym.get(sym,[]) + [line[1]]

def load_motif_sites(motif_annotation, f):
    idx = os.path.basename(f).split('.')[0]
    motif = motif_annotation[idx]
    #score = np.load(f)
    score = []
    with open(f) as inf:
        for line in inf:
            line = line.strip().split()
            score.append(float(line[-1]))
    score = np.array(score)
    #print score.shape
    #score = score.astype(np.float32)
    cutoff = np.percentile(score, percentage_cutoff)
    Y, = np.where(score > cutoff) # 0-based regions
    return (motif, idx, Y)
    #return (motif, idx, 'debug')

motif_table = {}
load_motif_meta(motif_table)

##with h5py.File('marge2_motif_100bp_%s.h5'%percentage_cutoff, 'a') as f:
with h5py.File('mm10_marge2_motif_100bp_%s.h5'%percentage_cutoff, 'a') as f:
    ids = f.create_dataset("IDs", shape=(len(fs),), dtype=np.dtype((str, 25)), compression='gzip', shuffle=True, fletcher32=True)
    TF = f.create_dataset("TFs", shape=(len(fs),), dtype=np.dtype((str, 25)), compression='gzip', shuffle=True, fletcher32=True)
    for index, fi in enumerate(fs):
        motif, i, Y = load_motif_sites(motif_table, fi)
        ids[index] = i
        TF[index] = motif
        f[i] = Y
