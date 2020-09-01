#!/usr/bin/env python

import numpy as np
import glob
import os, sys
import h5py
import pickle

if len(sys.argv) <= 1:
    print >>sys.stderr, "run_uDHS_count_h5_genome_window1kb.py <factor>"
    sys.exit(0)

factor = sys.argv[1]
# get meta information
factors = {}
stat = {}
with open('margeFactor.csv') as inf:
    inf.readline()
    for line in inf:
        line = line.strip().split(',')
        if line[9] == factor:
            factors[line[0]] = line[9]
            stat[line[9]] = stat.get(line[9], 0) + 1
print(stat)

# get shape 0 for h5py numpy array
f = os.popen('wc -l hg38_window1kb.bed')
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
                os.system('bigWigAverageOverBed %s hg38_window1kb.bed stdout > %s.tab  ' % (fin[0], prefix))
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
                    os.system('bigWigAverageOverBed %s hg38_window1kb.bed stdout > %s.tab  ' % (fin[0], prefix))
                    if os.path.getsize('%s.tab' % prefix) > 0:
                        shape2 += 1
                        files.append('%s.tab' % prefix)

f = h5py.File('hg38_window1kb_%s.h5' % factor, 'a')
if not ('Count' in f.keys()):
    RP = f.create_dataset("Count", dtype=np.float32, shape=(l, shape2), compression='gzip', shuffle=True, fletcher32=True)
    ids = f.create_dataset("IDs", shape=(shape2,), dtype=np.dtype((str, 10)), compression='gzip', shuffle=True, fletcher32=True)
else:
    RP = f["Count"]
    ids = f['IDs']

ids[...] = np.array(map(lambda x: x.replace('_treat.bw.tab', ''), files))
f.flush()
print(len(ids))
f.close()

n=0
print(shape2)
print(len(files))
with h5py.File('hg38_window1kb_%s.h5' % factor, 'a') as store:
    for f in files:
       print(f)
       print(n)
       tmp = []
       with open(f) as inf:
           for line in inf:
               line = line.split()[4] # average on all bases
               tmp.append(float(line.strip()))
       tmp = np.array(tmp, dtype='float32')
       store["Count"][:,n] = tmp
       n += 1
