#!/usr/bin/env python

import numpy as np
import _bw
import glob
import os, sys
import h5py
import pickle

if len(sys.argv) <= 1:
    print >>sys.stderr, "run_rp.py <factor>"
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
f = os.popen('wc -l histonerp/test/hg38.tss')
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
        if not os.path.exists(os.path.basename(fin[0])+'.txt'):
            # check status
            status = os.system('bigWigInfo %s' % fin[0])
            if status == 0:
                # without promoter +/- 1000bp
                #if factors[key] == 'H3K4me3':
                #    _bw.getrp(fin[0], 'histonerp/test/hg38.tss', os.path.basename(fin[0])+'_rdown1kb.txt', 1e3, -1000, 1000)
                #else:
                #    _bw.getrp(fin[0], 'histonerp/test/hg38.tss', os.path.basename(fin[0])+'_rdown1kb.txt', 1e4, -1000, 1000)

                # with promoter +/- 1000bp
                #if factors[key] == 'H3K4me3':
                #    _bw.getrp(fin[0], 'histonerp/test/hg38.tss', os.path.basename(fin[0])+'_rdown1kb.txt', 1e3, 0, 0)
                #else:
                #    _bw.getrp(fin[0], 'histonerp/test/hg38.tss', os.path.basename(fin[0])+'_rdown1kb.txt', 1e4, 0, 0)

                _bw.getrp(fin[0], 'histonerp/test/hg38.tss', os.path.basename(fin[0])+'.txt', 1e4, 0, 0)

                if os.path.getsize(os.path.basename(fin[0])+'.txt') > 0:
                    shape2 += 1
                    files.append(os.path.basename(fin[0])+'.txt')
        else:
            # get shape 2
            if os.path.getsize(os.path.basename(fin[0])+'.txt') > 0:
                shape2 += 1
                files.append(os.path.basename(fin[0])+'.txt')

f = h5py.File('margeRP_%s.h5' % factor, 'a')
if not ('RP' in f.keys()):
    RP = f.create_dataset("RP", dtype=np.float32, shape=(l, shape2), compression='gzip', shuffle=True, fletcher32=True)
    RefSeq = f.create_dataset("RefSeq", shape=(l, ), dtype=np.dtype((str, 200)), compression='gzip', shuffle=True, fletcher32=True)
    ids = f.create_dataset("IDs", shape=(shape2,), dtype=np.dtype((str, 10)), compression='gzip', shuffle=True, fletcher32=True)
else:
    RP = f["RP"]
    if "RefSeq" in f.keys():
        del f["RefSeq"]
    if "IDs" in f.keys():
        del f["IDs"]
    RefSeq = f.create_dataset("RefSeq", shape=(l, ), dtype=np.dtype((str, 200)), compression='gzip', shuffle=True, fletcher32=True)
    ids = f.create_dataset("IDs", shape=(shape2,), dtype=np.dtype((str, 10)), compression='gzip', shuffle=True, fletcher32=True)

ids[...] = np.array(map(lambda x: x.replace('_treat.bw.txt', ''), files))
f.flush()
print(ids)
tmp = []
refseqs = []
n = 0
with open(files[0]) as inf:
    for line in inf:
        line = line.split()
        tmp.append(float(line[4]))
        refseqs.append(':'.join(line[:4]))
tmp = np.array(tmp, dtype=np.float32)
RP[:,n] = tmp
RefSeq[...] = np.array(refseqs)
print(RefSeq)
f.flush()
f.close()

with h5py.File('margeRP_%s.h5' % factor, 'a') as store:
    for f in files[1:]:
       print(n)
       n += 1
       tmp = []
       with open(f) as inf:
           for line in inf:
               line = line.split()[4]
               tmp.append(float(line.strip()))
       tmp = np.array(tmp, dtype='float32')
       store["RP"][:,n] = tmp
