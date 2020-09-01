#!/usr/bin/env python

import numpy as np
import _bw
import glob
import os
import h5py
import pickle


# get meta information
factors = {}
stat = {}
with open('margeFactor.csv') as inf:
    inf.readline()
    for line in inf:
        line = line.strip().split(',')
        factors[line[0]] = line[9]
        stat[line[9]] = stat.get(line[9], 0) + 1
print(stat)

# get shape 0 for h5py numpy array
f = os.popen('wc -l 1025_treat.bw.tab')
a = f.read()
f.close()
l = int(a.split()[0])

# get bigWig directory
with open(os.path.expanduser('~/01_Projects/Programming/dc2/scripts/samples_directory.txt'), 'rb') as inf:
    directory = pickle.load(inf)

# use _bw to get rp for each of the files
status = 0
shape2 = 0
files = []
for key in factors:
    fin =  glob.glob(os.path.join(directory[int(key)], '*treat*bw'))
    if fin:
        print(fin[0])
        # check file existence
        prefix = os.path.basename(fin[0])
        if not os.path.exists('%s.tab' % prefix):
            # check status
            status = os.system('bigWigInfo %s' % fin[0])
            if status == 0:
                if factors[key] == 'DNase':
                    os.system('bigWigAverageOverBed -stats=%s.stats.ra %s human_unionDHS_fc5_75merge_split.bed stdout | cut -f 4,5,6 > %s.tab  ' % (prefix, fin[0], prefix))
                else:
                    os.system('bigWigAverageOverBed -sampleAroundCenter=1000 -stats=%s.stats.ra %s human_unionDHS_fc5_75merge_split.bed stdout | cut -f 4,5,6 > %s.tab ' % (prefix, fin[0], prefix))
                if os.path.getsize('%s.tab' % prefix) > 0:
                    shape2 += 1
                    files.append('%s.tab' % prefix)
        else:
            # get shape 2
            if os.path.getsize('%s.tab' % prefix) > 0:
                shape2 += 1
                files.append('%s.tab' % prefix)
            else:
                status = os.system('bigWigInfo %s' % fin[0])
                if status == 0:
                    prefix = os.path.basename(fin[0])
                    if factors[key] == 'DNase':
                        os.system('bigWigAverageOverBed -stats=%s.stats.ra %s human_unionDHS_fc5_75merge_split.bed stdout | cut -f 4,5,6 > %s.tab  ' % (prefix, fin[0], prefix))
                    else:
                        os.system('bigWigAverageOverBed -sampleAroundCenter=1000 -stats=%s.stats.ra %s human_unionDHS_fc5_75merge_split.bed stdout | cut -f 4,5,6 > %s.tab ' % (prefix, fin[0], prefix))
                    if os.path.getsize('%s.tab' % prefix) > 0:
                        shape2 += 1
                        files.append('%s.tab' % prefix)

f = h5py.File('margeCount.h5', 'a')
if not ('Count' in f.keys()):
    # DNase = f.create_dataset("RPGroup/DNase", dtype=np.float32, shape=(l, stat['DNase']), compression='gzip', shuffle=True, fletcher32=True)
    RP = f.create_dataset("Count", dtype=np.float32, shape=(l, shape2), compression='gzip', shuffle=True, fletcher32=True)
    ids = f.create_dataset("IDs", shape=(shape2,), dtype=np.dtype((str, 6)), compression='gzip', shuffle=True, fletcher32=True)
else:
    #RP = f["/RPGroup/DNase"]
    RP = f["Count"]
    ids = f['IDs']

ids[...] = np.array(map(lambda x: x.replace('_treat.bw.tab', ''), files))
f.flush()
print(len(ids))
f.close()

n=0
print(shape2)
print(len(files))
with h5py.File('margeCount.h5', 'a') as store:
    for f in files:
       print(n)
       tmp = []
       with open(f) as inf:
           for line in inf:
               line = line.split()[2]
               tmp.append(float(line.strip()))
       tmp = np.array(tmp, dtype='float32')
       store["Count"][:,n] = tmp
       n += 1
