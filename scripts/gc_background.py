#!/usr/bin/env python

"""
test for hg38 H3K27ac GC background
"""
from marge2 import Model
import pickle
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as p
import seaborn as sns

import h5py
import os
import sys
import numpy as np
import pandas as pd

h5 = '/data/home/qqin/12_data/hg38_window1kb_H3K27ac_averageonallbase.h5'
tssbin = '/data/home/qqin/12_data/hg38_window1kb_tss.bed'
gc = '/data/home/qqin/MARGE/PhaseD_GC/hg38_gc5Base.bin.npy'

class Weight:
    def __init__(self, bl=1000):
        padding = int(1e5)            # TSS +/- 100kb
        assert bl > 0
        assert (2*padding+bl)%bl == 0
        
        self.bl = bl
        self.binN = (2*padding+bl)/bl      # bin number
        
        distances = np.array([ z + bl/2 for z in range(-padding-bl/2, padding+bl/2, bl) ], dtype=np.float32)
        self.alpha = -math.log(1.0/3.0)*10          # 1e5/1e4, 1e4: half decay
        self.weight = self.balance_weight(distances)     # weight
        
    def get_weight(self):
        return self.weight
    def get_binnum(self):
        return self.binN
    def balance_weight(self, distances):
        weight = np.exp(-np.fabs(distances) * self.alpha/1e5)
        return  2*weight/ (1+weight)
    
def gc_hist():
    global gc
    gc = np.load(gc)
    plt.figure()
    plt.hist(gc, 50)
    plt.savefig('hg38_gc5Base_1kb.png', dpi=300)

def chrom_bin_bounary(c):
    chrom_bin = {}
    out = os.path.basename(c+'.pkl')
    if os.path.exists(out):
        with open(out, 'rb') as p:
            chrom_bin = pickle.load(p)
    else:
        with open(c) as inf:
            for line in inf:
                line = line.strip().split()
                chrom_bin[line[0]] = int(line[-1])-1
        with open(out, 'wb') as p:
            pickle.dump(chrom_bin, p)
    return chrom_bin

def rp_from_1kb_count(x):
    """ GCRP = W * S * GC * M
    """
    W = Weight(bl=1000)
    tssb = pd.read_table('/data/home/qqin/12_data/hg38_window1kb_tss.bed', header=None)
    bins = (tssb.iloc[:,-1]-1).values
    chrs = (tssb.iloc[:,0]).values
    tss = (tssb.iloc[:,1]).values
    bc = ((tssb.iloc[:,6] + tssb.iloc[:,7])//2).values

    # case study for GC GR
    with open('GR.symbol.fs_H3K27ac_withpromoter') as handler:
        GR_model = pickle.load(handler)

    # print GR_model.X
    # print GR_model.symbols
    ids = GR_model.sid

    diff = np.zeros(len(chrs), dtype=np.int32)
    refseq = tssb.apply(lambda x: str(x[3]).split(':')[0], axis=1).values
    bc_2d = []; idx_2d = []
    print GR_model.refseqs[1:3]
    print refseq[1:3]
    diff[np.in1d(refseq,GR_model.refseqs)] = 1
    print np.sum(diff)
    for i in range(-100, 100):
        # bin center
        bc_2d.append(bc + i*W.bl) # 1000 bin length
        # TSS bin index array
        idx_2d.append(bins + i) 
    bc_2d = np.vstack(bc_2d).T
    idx_2d = np.vstack(idx_2d).T

    # boundary mask bins
    M = chrom_bin_bounary('/data/home/qqin/seqpos2/hg38_window1kb.bed')
    M = np.array([M[c] for c in chrs])
    print M.shape
    left_mask = idx_2d >= 0 
    right_mask = (idx_2d.T - M).T <= 0

    print left_mask.shape
    print right_mask.shape
    M = left_mask & right_mask
    W = W.balance_weight((bc_2d.T-tss).T)
    idx_2d = M * idx_2d
    GC = gc[idx_2d]         #          gene x bin
    ## output hdf5

    samples = len(ids)
    f = h5py.File('hg38_1kb_H3K27ac_RP_GR.h5', 'a')
    RP = f.create_dataset("RP", dtype=np.float32, shape=(len(tss), samples), compression='gzip', shuffle=True, fletcher32=True)
    GCRP = f.create_dataset("GCRP", dtype=np.float32, shape=(len(tss), samples), compression='gzip', shuffle=True, fletcher32=True)
    Refseqs = f.create_dataset("refseqs", dtype=np.dtype((str, 100)), shape=(len(tss),), compression='gzip', shuffle=True, fletcher32=True)
    Diff = f.create_dataset("diff", dtype=np.int32, shape=(len(tss), ), compression='gzip', shuffle=True, fletcher32=True)
    #GCbias = f.create_dataset("GCbias", dtype=np.float32, shape=(len(tss), samples), compression='gzip', shuffle=True, fletcher32=True)
    Refseqs = refseq
    Diff[:] = diff
    print np.sum(Diff[...])

    # S signal 
    with h5py.File(h5) as h:
        allids = h['IDs'][...]
        S = h['OrderCount']
        map_ids = {}
        for ii, i in enumerate(allids):
            map_ids[i] = ii
        idx = np.array([ map_ids[i] for i in ids ])
        # hdf index must be increasing order
        idx_sort = np.sort(idx)
        S = S[:, idx_sort]
        # map back to original iid order
        map_ids = {}
        for ii, i in enumerate(idx_sort):
            map_ids[i] = ii
        idx = np.array([ map_ids[i] for i in idx ])
        S = S[:, idx]
        for i, s in enumerate(S.T): # for each sample
            s = s[idx_2d]           # gene x bin
            GCRP[:,i] = np.sum(1000* GC * s * W * M, axis=1)
            RP[:,i] = np.sum(1000*s * W * M, axis=1)
            print GCRP.shape
        #GCbias = 1.0*GCRP[...]/RP[...]
    f.close()

#gc_hist()
#rp_from_1kb_count('a')

def draw_GCbias_hist():
    with h5py.File('hg38_1kb_H3K27ac_RP_GR.h5') as s:
        gcb = 1.0 * s['GCRP'][...] / s['RP'][...]
        gcb[np.isnan(gcb)] = -1.0
        diff = s['diff'][...]

    d_logic = (diff == 1)
    print np.sum(d_logic)
    from matplotlib.backends.backend_pdf import PdfPages
    #pdf_pages = PdfPages('GR_1kb_GCbias_count.pdf')
    pdf_pages = PdfPages('GR_1kb_GCbias.pdf')
    for col in gcb.T:
        fig = plt.figure(figsize=(12.27, 8.69), dpi=100)
        #plt.hist(col[~d_logic], bins = 100, label='background', facecolor='blue', alpha=0.3)
        #plt.hist(col[d_logic],  bins = 100, label='diff', facecolor='red', alpha=0.3)
        plt.hist(col[~d_logic], bins = 100, label='background', facecolor='blue', alpha=0.3, normed=1)
        plt.hist(col[d_logic],  bins = 100, label='diff', facecolor='red', alpha=0.3, normed=1)
        #print(np.max(col[d_logic]))
        #print(np.min(col[d_logic]))
        plt.legend(loc='upper center', shadow=True, fontsize='x-large')
        plt.xlim(-2, 100)
        plt.tight_layout()
        plt.close('all')
        pdf_pages.savefig(fig)
        plt.clf() # plf.cla() 
    pdf_pages.close()

draw_GCbias_hist()

