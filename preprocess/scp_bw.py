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
#with open('margeFactor.csv') as inf:
with open('../uDHSRead/margeK27ac.csv') as inf:
    inf.readline()
    for line in inf:
        line = line.strip().split(',')
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
        status = os.system('bigWigInfo %s' % fin[0])
        if status == 0:
            os.system('scp -P 22 %s qinq@compbio3.tongji.edu.cn:' % fin[0])
