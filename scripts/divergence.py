import pandas as pd
import numpy as np
import scipy.stats as stats
import argparse

from numpy.linalg import norm

def JSD(P, Q):
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return 0.5 * (stats.entropy(_P, _M) + stats.entropy(_Q, _M))


p = argparse.ArgumentParser()

p.add_argument('-f', required=True, help='factor')
p.add_argument('--kl', default=False, action='store_true', help='factor')


args = p.parse_args()


pq = pd.read_csv(args.f, encoding='utf-8', index_col=0)
meta = pd.read_table('/data/home/qqin/seqpos2/cistrome.txt',index_col=0, encoding='utf-8')
meta = meta.ix[:, ['symbol', 'dbd', "entrez"]]
dc_meta = pd.read_table('/data/home/qqin/01_Projects/Programming/dc2/scripts/best_dc_tfcr_basedon_frip_peak_dhs.xls', index_col=0, encoding='utf-8')

dc_meta = dc_meta.ix[:, ['FactorName', 'cell_line', 'cell_type']]

dc_meta.FactorName = dc_meta.apply(lambda x: 'ChIP_seq_%s_in_%s_%s'%(x[0], x[1], x[2]), axis=1)

rows = np.concatenate([meta.index.values.astype(str), dc_meta.index.values.astype(str)], axis=0)
print(rows)
print(meta.shape)
print(dc_meta.shape)
meta = np.vstack([meta.values, dc_meta.values])

meta = pd.DataFrame(meta, index=rows, columns=['symbol', 'annotation1', 'annotation2'])

print(pq.head())
columns = [ i[2:-1] for i in pq.columns[1:-1].astype(str).tolist() ]
symbols = meta.ix[columns,'symbol'].values.tolist()
symbols.insert(0, 'chipseq')

symbols = np.array(symbols)
pq = pq/pq.sum(axis=0)

# entropy actually can calculate KL divergence
if args.kl:
    final=pq.iloc[:,:-1].apply(lambda x: stats.entropy(x, pq.iloc[:,-1],base=2), axis=0)
else: # JS divergence
    final=pq.iloc[:,:-1].apply(lambda x: JSD(x, pq.iloc[:,-1]), axis=0)

final.index=symbols
final.sort_values(inplace=True)

rank = final.rank(ascending=False)
final = pd.concat([final, rank], axis=1)
final.to_csv('%s_kl_rank.csv' % args.f)
