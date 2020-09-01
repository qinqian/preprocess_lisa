import h5py
import sys
import os
import numpy as np

h5 = sys.argv[1]
index = []
## human
with open('bigWigAverageOverBed.originalorder.hg') as sortedf:
     for line in sortedf:
          index.append(int(line.strip())-1)

# human
with h5py.File(h5) as store:
    target = store['Count'][...]
    shape = target.shape
    #del store['Count']
    #del store['OrderCount']

print shape
target = target[np.argsort(index)]
# human
with h5py.File(h5) as store:
    store["OrderCount"] = target
