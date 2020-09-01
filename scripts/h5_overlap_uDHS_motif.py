#!/usr/bin/env python

import h5py
import os
import numpy as np
import argparse

p = argparse.ArgumentParser()
p.add_argument('--tf', required =True, help='numpy binary 100bp window occurence array')
p.add_argument('--uDHS', default =True, action='store_false', help='numpy binary 100bp window occurence array')

args = p.parse_args()

if args.uDHS:
    label = 'dhs'
    dhs = np.load('hg38_100bp_uDHS.out.npy')-1 # 1-based to 0-based index
else:
    label = ''

tf = np.load(args.tf) - 1 # 1-based to 0-based
print(tf)

#m = np.load("hg38_100to1000window.out.npy")
#row = len(m)
#

f = open("%s_overlap_motif.txt.%s"%(args.tf, label), 'a')
with h5py.File(os.path.expanduser('~/seqpos2/marge2_motif_100bp_%s.h5' % 99)) as store: # 0-based index from np.where > percentage_cutoff percentile
    TF = store['TFs'][...] 
    f.write('\t'.join(list(map(lambda x: x.decode('utf-8'), TF)))+'\n')

for percentage_cutoff in [97,98,99]:
    with h5py.File(os.path.expanduser('~/seqpos2/marge2_motif_100bp_%s.h5' % percentage_cutoff)) as store: # 0-based index from np.where > percentage_cutoff percentile
        # Python3 dtype S25
        # NOTE: http://docs.h5py.org/en/latest/strings.html
        ids = store['IDs'][...] 

        #row = mapping[-1]+1 # map to 1kb
        result = []
        for i, key in enumerate(ids):
            motif_index = store[key][...] # 0-based
            if args.uDHS:
                motif_index = np.intersect1d(motif_index, dhs)
            total = len(motif_index) # background motif sites
            result.append(len(np.intersect1d(motif_index, tf))*1.0/total)
            #print(result[i])
        f.write("%d\t%s\n" % (percentage_cutoff, '\t'.join(list(map(str, result)))))

f.close()
