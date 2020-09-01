#!/usr/bin/env python

import numpy as np
import os
import h5py
import pickle

import glob
#fs = ['cistrome/MS00011.pwm.100bp.bin.npy', 'cistrome/MS00023.pwm.100bp.bin.npy']
fs = glob.glob('cistrome/*bin.npy')
percentage_cutoff = 97
window_map = np.load(os.path.expanduser('~/MARGE/PhaseA_motif_deltarp/hg38_100to1000window.out.npy'))

def load_motif_meta(map_motifsym):
    with open('/data/home/qqin/MARGE/PhaseA_motif_deltarp/motif_meta.txt') as inf:
        for line in inf:
            line = line.strip().split()
            sym = line[0].upper()
            map_motifsym[sym] = map_motifsym.get(sym,[]) + [line[1]]

def load_motif_sites(motif_annotation, f):
    idx = os.path.basename(f).split('.')[0]
    motif = motif_annotation[idx]
    score = np.load(f)
    score = score.astype(np.float32)
    cutoff = np.percentile(score, percentage_cutoff)
    # map from 100bp window to 1kb
    index, = np.where(score >= cutoff)
    # 0-based
    index = np.unique(window_map[index])
    Y = np.zeros(3209513).astype(np.int32)
    Y[index] = 1
    return (motif, idx, Y)

motif_table = {}
load_motif_meta(motif_table)

with h5py.File('marge2_motif_100to1kb_to_gene%s.h5'%percentage_cutoff, 'a') as f:
    ids = f.create_dataset("IDs", shape=(len(fs),), dtype=np.dtype((str, 25)), compression='gzip', shuffle=True, fletcher32=True)
    #ids = f['IDs']
    TF = f.create_dataset("TFs", shape=(len(fs),), dtype=np.dtype((str, 25)), compression='gzip', shuffle=True, fletcher32=True)
    #TF = f['TFs']
    # 3209513 bins for hg38 1kb window 
    MR = f.create_dataset("IsMotifRegion%s"%percentage_cutoff,  dtype=np.int32, shape=(3209513, len(fs)), compression='gzip', shuffle=True, fletcher32=True)
    for index, f in enumerate(fs):
        motif, i, Y = load_motif_sites(motif_table, f)
        ids[index] = i
        TF[index] = motif
        MR[:,index] = Y
