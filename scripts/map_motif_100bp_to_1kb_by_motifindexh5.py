#!/usr/bin/env python
##python3

import h5py
import numpy as np


mapping = np.load('mm10_100to1000window.out.npy')

percentage_cutoff = 99
with h5py.File('mm10_marge2_motif_100bp_%s.h5' % percentage_cutoff) as store: # 0-based index from np.where > percentage_cutoff percentile
    with h5py.File('mm10_marge2_motif_100to1kb_to_gene%s.h5' % percentage_cutoff, 'a') as store_out:
        # Python3 dtype S25
        # NOTE: http://docs.h5py.org/en/latest/strings.html
        ids = store_out.create_dataset("IDs", shape=(len(store['IDs']),), dtype="S25", compression='gzip', shuffle=True, fletcher32=True)
        #ids = f['IDs']
        TF = store_out.create_dataset("TFs", shape=(len(store['IDs']),), dtype="S25", compression='gzip', shuffle=True, fletcher32=True)
        ids = store['IDs'][...] 
        store_out.flush()
        TF = store['TFs'][...] 
        store_out.flush()
        # mm10
        print(2730901) # 1kb mm10 2730901, 100bp mm10 27308748
        row = mapping[-1]+1
        print(row)
        MR = store_out.create_dataset("IsMotifRegion%s"%percentage_cutoff,  dtype=np.int32, shape=(row, len(store["IDs"])), compression='gzip', shuffle=True, fletcher32=True)
        for i, key in enumerate(store["IDs"]):
            map1kb = mapping[store[key][...]]
            tmp = np.zeros(row)
            tmp[map1kb] = 1
            MR[:,i] = tmp
            store_out.flush()
    
