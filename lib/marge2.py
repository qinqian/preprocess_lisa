from sklearn.linear_model import LogisticRegression, Lasso, LinearRegression
import sys
from sklearn.cluster import KMeans,MiniBatchKMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import make_scorer, roc_auc_score, average_precision_score, roc_curve, r2_score

from marge2_conf import Config

import math
import h5py
from scipy.stats import ks_2samp
import pandas as pd
import numpy as np
import os
import theano

def annotate(x):
    result = ''
    sep= ''
    if not pd.isnull(x.cell_line):
       result += x.cell_line
       sep='-'
    else:
       sep=''
    if not pd.isnull(x.cell_type):
       result += sep+x.cell_type
       sep='-'
    else:
       sep=''
    if not pd.isnull(x.tissue):
       result += sep+x.tissue
    return result

def parse_simple_file(f):
    """ one line per file
    """
    with open(f) as fin:
         content = set([line.strip().upper() for line in fin]) # to uppercase
    return np.array(list(content))

class Weight:
    def __init__(self, bl=1000):
        padding = int(1e5)            # TSS +/- 100kb
        assert bl > 0
        assert (2*padding+bl)%bl == 0

        self.bl = bl
        self.binN = (2*padding+bl)/bl      # bin number

        distances = np.array([ z + bl/2 for z in range(int(-padding-bl/2), int(padding+bl/2), bl) ], dtype=np.float32)
        self.alpha = -math.log(1.0/3.0)*10          # 1e5/1e4, 1e4: half decay
        self.weight = self.balance_weight(distances)     # weight
    def get_weight(self):
        return self.weight
    def get_binnum(self):
        return self.binN
    def balance_weight(self, distances):
        weight = np.exp(-np.fabs(distances) * self.alpha/1e5)
        return  2*weight/ (1+weight)

class Annotation(object):
    def __init__(self, c, gene_file, bw=None, bl=1000):
        self.c = c
        self.chrom = get_chrom_len(c)
        self.genes = parse_simple_file(gene_file)
        self.annotation = c.get_annotation
        self.f = open(c.get_annotation)
        self.sample = bw
        self.weight = Weight(bl)

    def get_gene_set(self):
        return self.genes

    def _get_bw_signal(self, chrom, bwstart, bwend, bins):
        """ get bigwig signal summary from BigWig file
        """
        cmd = 'bigWigSummary {} {} {} {} {}'
        sysf = os.popen(cmd.format(self.sample, chrom, bwstart, bwend, bins))
        bwsum = sysf.read().strip().split()
        score = []
        for bs in bwsum:
            if bs == 'n/a':
                score.append(0)
            else:
                score.append(float(bs))
        bwsum = np.array(score, np.float32)
        return bwsum

    def _get_h5_signal(self, index):
        """ get bigwig signal summary from HDF5
        when self.sample is list or integer
        """
        h5 = HDF(self.c)
        return h5.count_hdf_reader(iid=self.sample, index=index)

    def get_signal(self, *args):
        if isinstance(self.sample, str):
            chrom, bwstart, bwend, bins = args
            signal = self._get_bw_signal(chrom, bwstart, bwend, bins)
            return signal
        elif isinstance(self.sample, int) or isinstance(self.sample, list):
            index, = args
            return self._get_h5_signal(index)
        else:
            raise Exception

    def get_delta_rp(self, tss, region):
        """ delete a region or regions RP
        do not consider boundary problem
        merge into get_rp later
        region.label determine how many bins to divide into
        """
        assert tss.chrom == region.chrom
        # delete regions
        if region.label > 1:
            bin_len, distances = region.get_bin_center(region.label)
            weight = self.weight.balance_weight(tss.start-distances)
            signal = self.get_signal(*region.tuple())
            return weight, np.dot(weight, bin_len * signal)
        # delete a region
        if region.label == 1:
            weight = self.weight.balance_weight([tss.start-region.get_center()])
            signal = self.get_signal(*region.tuple())
            return weight, len(region) * weight * signal

    def get_rp(self, gene, region_to_del=None):
        """
        Calculate RP or DeltaRP based on one TSS
        gene: a Region object represent one TSS
        region_to_del:
          a Region object to be deleted,
          the label property determine how many bins to divide into, e.g. 1 means 1bp, 1000 means 1000bp window width
        """
        if region_to_del == None:
            delta = 0
        else:
            _, delta = self.get_delta_rp(gene, region_to_del)
        if isinstance(self.sample, str):
            l, r, region = gene.get_context(self.chrom, self.weight.bl)
            bwsum  = self.get_signal(*region.tuple())
            bwsum = np.concatenate([np.zeros(l), bwsum, np.zeros(r)])
            weight = self.weight.get_weight()
            rp = self.weight.bl * np.dot(bwsum, weight)
        elif isinstance(self.sample, int):
            rp = HDF(self.c).rp_hdf_reader(iid=self.sample, gene=gene.label)
        return rp, delta, rp - delta, 1.0*(delta)/rp

    def __iter__(self):
        return self

    def __next__(self):
        try:
            line = next(self.f)
        except StopIteration:
            self.f.close()
            raise StopIteration()
        line = line.strip().split()
        ###rs = line[3].split(':')[0]
        rs = line[3].upper() ## support cluster score, to uppercase
        return (Region(line[0], int(line[1]), int(line[2]), rs),
                Region(line[-4], int(line[-3]), int(line[-2]), int(line[-1])))


class Region(object):
    def __init__(self, chrom, start, end, *args, **kwargs):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.label = args[0]

    def __str__(self):
        return "%s\t%s\t%s\t%s" %(self.chrom, self.start, self.end, self.label)

    def __len__(self):
        return self.end - self.start

    def get_bin_center(self, bin_num):
        """ divide the region into bins
        fetch the center of each bin
        """
        bin_len = (self.end - self.start)//bin_num
        if bin_len == 1:
            return 1, np.arange(self.start, self.end)

        elif bin_len > 1:
            half_bin = bin_len // 2
            distances = np.linspace(self.start + half_bin,
                                    self.end   - half_bin,
                                    bin_num)
            return bin_len, distances

    def get_center(self):
        """ get the center f the region """
        return (self.start+self.end)//2

    def clip_region(self, chrom, bl=1000, padding=int(1e5)):
        """ generate a weight vector and check rp score
        """
        center = self.start
        bwstart = center - padding - bl/2 if center > padding + bl/2 else 0
        bwend   = center + padding + bl/2 if center + padding + bl/2 <= chrom[self.chrom] else chrom[self.chrom]
        bins = (bwend - bwstart) / bl
        left = padding+bl/2+bwstart-center
        right = padding+bl/2-(bwend-center)
        left   = left/bl
        right  = right/bl
        return left, right, Region(self.chrom, bwstart, bwend, bins)

    def get_context(self, chrom, bl=1000):
        if isinstance(self.label, str):
            return self.clip_region(chrom, bl=bl)
        else:
            raise TypeError

    def tuple(self):
        return (self.chrom, self.start, self.end, self.label)

class HDF(object):
    def __init__(self, c):
        self.c = c
        # python3 pandas http://stackoverflow.com/questions/18171739/unicodedecodeerror-when-reading-csv-file-in-pandas-with-python
        self.meta = pd.read_csv(self.c.get_meta, encoding="ISO-8859-1")

    def filter_low_quality(self, ids, factor="H3K27ac"):
        # previous for human k27ac
        # selector = (self.meta['UniquelyMappedRatio'] > 0.5) & (self.meta['MappedReadsNumber'] > 5e6) & (self.meta['AllPeaksNumber'] > 1000) & (self.meta['PBC'] > 0.85) & (self.meta['FRiP'] > 0.1) & (self.meta['FactorName'] == factor)
        # lower the standard to preserve more samples
        # selector = (self.meta['UniquelyMappedRatio'] > 0.4) & (self.meta['MappedReadsNumber'] > 4e6) & (self.meta['AllPeaksNumber'] > 50) & (self.meta['PBC'] > 0.75) & (self.meta['FRiP'] > 0.01) & (self.meta['FactorName'] == factor)

        selector = (self.meta['UniquelyMappedRatio'] > 0.4) & (self.meta['MappedReadsNumber'] > 4e6) & (self.meta['AllPeaksNumber'] > 50) & (self.meta['PBC'] > 0.7) & (self.meta['FactorName'] == factor)

        sids = self.meta.ix[selector, 'X']
        sids = set(map(str, list(sids)))
        ## ids is from hdf5 file, in python3 
        ## important NOTE!!!
        ## http://docs.h5py.org/en/latest/strings.html
        ## Variable-length ASCII (Python 2 str, Python 3 bytes)
        ## Variable-length UTF-8 (Python 2 unicode, Python 3 str)
        ## we stored the hdf5 all based on Python2 str, i.e. Python3 bytes
        ## in Python3, we convert bytes to str http://stackoverflow.com/questions/606191/convert-bytes-to-a-python-string
        iids = set()
        for i in ids:
            try:
                iids.add(i.decode('utf-8').split('_')[0]) # for the weird bytes and ID name
            except AttributeError:
                iids.add(i.split('_')[0]) # for the weird bytes and ID name

        final = list(sids & iids)
        return final

    def get_rp_XY_withfold(self, fold, histone):
        """ input fold change of pandas DataFrame,
             sample ids: numpy array
             refseqs:    numpy array, the first refseq
             symbols:    numpy array, use the first refseq for a gene symbol
             X:          numpy array, float32
             Y:          numpy array, int32, continuous vector
        """
        if isinstance(fold, pd.DataFrame):
            common_ids = self.rp_count_common_id(histone)
            # filter by quality
            common_ids = self.filter_low_quality(common_ids, factor=histone)
            ids, ref, sym, X = self.rp_hdf_reader(factor=histone, iid=common_ids, gene=None)

            irna = fold.index
            ## TODO: add heterogenous gene name match
            if len(np.intersect1d(ref, irna))>5:   # not so few overlap
                print('input refseq ...')
                Y = fold.ix[ref,1].values
            elif len(np.intersect1d(sym, irna))>5: # not so few overlap
                print('input symbol ...')
                Y = fold.ix[sym,1].values
            else:
                raise BaseException

            # remove not matching gene/refseq
            index, = np.where(~np.isnan(Y))
            Y = Y[index]
            ref = ref[index]
            sym = sym[index]
            X = X[index]
            map_dict = {}
            for i, (r, s) in enumerate(zip(ref, sym)):
                # get the first refseq for a symbol
                map_dict.setdefault(s, (i, r))
            # variable after removing replicates
            si = []
            ss = []
            rr = []
            Ym = []
            for s, (i, r) in map_dict.items():
                si.append(i)
                Ym.append(Y[i])
                rr.append(r)
                ss.append(s)
            Y = np.array(Ym, dtype=np.float32)
            X = X[si]
            return np.array(ids), np.array(rr), np.array(ss), X, Y

    def get_rp_XY(self, gene_set, histone, other_h5, DNase=False):
        """ get X and Y for training model

        gene_set: a gene list of refseq/symbol

        return:
             sample ids: numpy array
             refseqs:    numpy array, the first refseq
             symbols:    numpy array, use the first refseq for a gene symbol
             X:          numpy array, float32
             Y:          numpy array, int32, binary vector
        """
        if not DNase:
            # common id between rp and count
            common_ids = self.rp_count_common_id(histone)
            # filter by quality
            common_ids = self.filter_low_quality(common_ids, factor=histone)
            #sys.exit(1)
            ids, ref, sym, X = self.rp_hdf_reader(factor=histone, iid=common_ids, gene=None)
        else:
            with h5py.File(self.c.get_dnase_bin) as dnase:
                common_ids = list(map(lambda x: x.decode('utf-8').replace('dnase_', ''), dnase['IDs']))
            # filter by quality some ENCODE3 phase data has 0 peaks!!! and even be selected by l1
            common_ids = self.filter_low_quality(common_ids, factor='DNase')
            rpi = self.rp_hdf_reader(factor='DNase', only_id=True)
            common_ids = list(np.intersect1d(common_ids, rpi))
            ids, ref, sym, X = self.rp_hdf_reader(factor='DNase', iid=common_ids, gene=None)

        # lisa fastq mode with in-house HDF5
        # combine with public dataset
        # [public, in-house]
        if other_h5 != None:
            with h5py.File(other_h5) as self_h5:
                self_rp = self_h5["RP"][...]
                self_ids = np.array([ i.decode('utf-8') for i in self_h5["IDs"][...] ])
            print(self_ids)
            ids = np.concatenate([ids, self_ids])
            X = np.hstack([X, self_rp])

        Y = np.zeros(len(ref))
        if len(np.intersect1d(ref, gene_set))>1:
            print('input refseq ...')
            Y[np.in1d(ref, gene_set)] = 1
        elif len(np.intersect1d(sym, gene_set))>1:
            print('input symbol ...')
            Y[np.in1d(sym, gene_set)] = 1
        else:
            raise BaseException

        map_dict = {}
        for i, (r, s) in enumerate(zip(ref, sym)):
            # get the first refseq for a symbol
            map_dict.setdefault(s, (i, r))

        # variable after removing replicates
        si = []
        ss = []
        rr = []
        Ym = []
        for s, (i, r) in map_dict.items():
            si.append(i)
            Ym.append(Y[i])
            rr.append(r)
            ss.append(s)
        Y = np.array(Ym, dtype=np.int32)
        X_unique = X[si]
        return np.array(ids), np.array(rr), np.array(ss), X, X_unique, Y, ref

    def rp_hdf_reader(self, factor, iid=None, gene=None, only_id=False):
        """
        c: Config object
        """
        h5 = self.c.get_rp(factor)
        with h5py.File(h5) as h5in:
            ids = h5in['IDs'][...]
            if only_id: return ids
            # be careful!!! to uppercase uniformly
            # python 3 bytes
            ref = np.array(list(map(lambda x: x.decode('utf-8').split(':')[-2].upper(), h5in['RefSeq'][...])))
            sym = np.array(list(map(lambda x: x.decode('utf-8').split(':')[-1].upper(), h5in['RefSeq'][...])))
            if gene == None:
                if iid == None:
                    return ids, ref, sym, h5in['RP'][...]
                elif isinstance(iid, int):
                    sample_index, = np.where(ids==str(iid))
                    return h5in['RP'][:, sample_index]
                elif isinstance(iid, list):
                    map_id = {}
                    for i, c in enumerate(ids):
                        map_id[c.decode('utf-8').split('_')[0]] = i

                    idx = np.array([ map_id[str(i)] for i in iid ])
                    iid = np.array(iid)
                    reindex = np.argsort(idx)
                    index = idx[reindex]
                    iid = iid[reindex]
                    mat = h5in['RP'][:, index]
                    return iid, ref, sym, mat
                raise TypeError
            elif isinstance(gene, list):
                if iid == None:
                    mat = np.zeros((len(gene), h5in['RP'].shape[1]))
                    for i, c in enumerate(gene):
                        gene_index, = np.where(ref == c)
                        mat[i,:] = h5in['RP'][gene_index,:][...]
                    return mat
                elif isinstance(iid, int):
                    mat = np.zeros(len(gene))
                    sample_index, = np.where(ids==str(iid))
                    for i, c in enumerate(gene):
                        gene_index, = np.where(ref == c)
                        mat[i] = h5in['RP'][gene_index, sample_index]
                    return mat
                elif isinstance(iid, list):
                    mat = np.zeros((len(gene), len(iid)))
                    for g, gc in enumerate(gene):
                        gene_index, = np.where(ref == gc)
                        for s, sc in enumerate(iid):
                            sample_index, = np.where(ids==str(sc))
                            mat[g, s] = h5in['RP'][gene_index, sample_index]
                    return mat
                raise TypeError
            elif isinstance(gene, str):
                return self.rp_hdf_reader(iid, [gene])[0]
            raise TypeError

    def rp_count_common_id(self, factor):
        ci = set(self.count_hdf_reader(factor, only_id=True).tolist())
        rpi = set(self.rp_hdf_reader(factor, only_id=True).tolist())
        return list(ci & rpi)

    def count_hdf_reader(self, factor, iid=None, index=None, only_id=False):
        """
        c: Config object
        """
        h5 = self.c.genome_count(factor)
        with h5py.File(h5) as h5in:
            ids = h5in['IDs'][...]
            if only_id: return ids
            count = h5in['OrderCount']
            if index==None:
                if iid == None:
                    print('warning: load too much data')
                    return count[...] # time consuming
                elif isinstance(iid, int):
                    sample_index, = np.where(ids==str(iid))
                    return count[:, sample_index][...]
                elif isinstance(iid, list):
                    map_ids = {}
                    for ii, i in enumerate(ids):
                        map_ids[i.decode('utf-8').split('_')[0]] = ii
                    idx = np.array([ map_ids[i] for i in iid ])
                    # hdf index must be increasing order
                    idx_sort = np.sort(idx)
                    mat = count[:, idx_sort]

                    # map back to original iid order
                    map_ids = {}
                    for ii, i in enumerate(idx_sort):
                        map_ids[i] = ii
                    idx = np.array([ map_ids[i] for i in idx ])
                    return mat[:, idx]
            elif isinstance(index, list):
                assert len(index) <= count.shape[0]+1
                if iid == None:
                    if len(index) == 1:
                        return count[index[0], :][...]
                    elif len(index) == 2:
                        assert index[0] < index[1]
                        return count[index[0]:index[1],:][...] # time consuming
                    elif len(index) > 2:
                        mat = np.zeros((len(index), count.shape[1]))
                        for i in index:
                            mat[i] = count[i,:][...]
                        return mat
                elif isinstance(iid, int):
                    sample_index, = np.where(ids==str(iid))
                    if len(index) == 1:
                        return count[index[0], sample_index]
                    elif len(index) == 2:
                        assert index[0] < index[1]
                        return count[index[0]:index[1], sample_index][...]
                    elif len(index) > 2:
                        mat = np.zeros(len(index))
                        for i, c in enumerate(index):
                            mat[i] = count[c, sample_index]
                        return mat
                elif isinstance(iid, list):
                    mat = np.zeros((len(index), len(iid)))
                    for i, sample in enumerate(iid):
                        sample_index, = np.where(ids==str(sample))
                        for j, region in enumerate(index):
                            mat[j,i] = count[region, sample_index]
                    return mat
            elif isinstance(index, int):
                return self.count_hdf_reader(iid, [index])

            raise TypeError

def get_chrom_len(c):
    chroms = {}
    with open(c.chrom) as inf:
        for line in inf:
            line = line.strip().split()
            chroms[line[0]] = line[1]
    return chroms

def get_clean_path(f):
    return os.path.abspath(os.path.expanduser(f))

def log_transform(x):
    pcount = 1
    xt = x + pcount
    #med = np.median(xt, 0)
    med = np.mean(xt, 0)
    return np.log2(xt) - np.log2(med)

def load_motif_meta(c):
    motifs = pd.read_table(c.get_motif_meta)
    return motifs.ix[:, ["id", "source", "dbd", "symbol", "refseq"]]

def load_chip_meta(c, chip):
    import codecs
    # UnicodeDecodeError: 'utf-8' codec can't decode byte 0x96 in position 432: invalid start byte
    # Python3 http://stackoverflow.com/questions/22216076/unicodedecodeerror-utf8-codec-cant-decode-byte-0xa5-in-position-0-invalid-s
    # for strange unicode string
    # str = unicode(str, errors='replace')
    # str = unicode(str, errors='ignore')
    # for file
    with codecs.open(c.tf_chipseq_meta, 'r', encoding='utf-8', errors='ignore') as inf:
        inf.readline()
        for line in inf:
            line = line.strip().split('\t')
            result = ''
            sep= ''
            if line[5] != 'NA':
               result += line[5]
               sep='-'
            else:
               sep=''
            if line[6] != 'NA':
               result += sep+line[6]
               sep='-'
            #else:
            #   sep=''
            #if not pd.isnull(line[7]):
            #   result += sep+line[7]
            chip[line[0]] = ','.join([result, line[8]]) # annotate cell and TF name

def site_binarize(c, bed):
    """ take an BED file, overlap with 1kb window,
    return a list of index 0-based """
    if not os.path.exists(bed): raise
    os.system("bedtools intersect -wa -u -a %s -b %s | cut -f 4 - > %s.index" % (c.genome_100bp_window, bed, bed))
    indexes = []
    with open("%s.index" % bed) as inf:
        for line in inf:
            r = int(line.strip())-1
            indexes.append(r)
    return indexes

def chrom_bin_boundary(c):
    chrom_bin = {}
    with open(c.genome_window) as inf:
        for line in inf:
            line = line.strip().split()
            chrom_bin[line[0]] = int(line[-1])-1
    return chrom_bin

def ks_test(x, y):
    ##one-side significant less x < y
    test = ks_2samp(x, y)
    d = test[0]
    p = test[1]/2
    return p

def ks_test_df(delta):
    """ a text table with two column names
    e.g.
           1 1 1 1 1 0 0 0 0 0 <-differential gene: 1, header 1
           A B C D E F G H I K <-refseq name, header 2
    motif1
    motif2
    motif3
    """
    header1 = pd.read_table(delta, nrows=1, header=None)
    delta_val = pd.read_table(delta, skiprows=1, header=0)
    pos = header1.iloc[-1].values==1
    result = delta_val.apply(lambda x: ks_test(x[~pos], x[pos]), axis=1)
    result.to_csv("%s.ks" % delta)

class Model(object):
    def __init__(self, X, X_unique, Y, refs, ids, symbols, refseqs, top=10, seed=999, mode="logit", penalty='l1'):
        """ building
        1. l1 regularized logistic regression or ( gene list )
        2. lasso linear regression ( gene fold change )
        with sklearn with gene rp and differential gene,
        set, and apply it to genome window count data
        used to predict TFBS """
        self.X = X_unique
        self.X_all = X
        self.Y = Y
        self.Y_count = None
        self.ids = ids
        self.symbols = symbols
        self.refseqs = refseqs # first refseqs
        self.all_refseq = refs # all refseqs
        self.bins = None
        self.top = top
        self.seed = seed

        self.mode = mode
        self.penalty = penalty

        if mode == "linear":
            self.model = Lasso(alpha=0.1, random_state=self.seed)
        elif mode == 'logit':
            self.model = LogisticRegression(penalty=penalty, tol=0.01, dual=False, C=1.0, random_state=self.seed)

        # selected features
        self.sid = None
        self.performance = None
        self.coefs = None
        self.params = None

    def renew_model(self, C=1, seed=999):
        self.seed = seed
        if self.mode == "linear":
            self.model = Lasso(alpha=C, random_state=self.seed)
        elif self.mode == "logit":
            self.model = LogisticRegression(penalty=self.penalty, tol=0.01, dual=False, C=C, random_state=self.seed)

    def _wrap_l2_fs(self):
        # really slow!!!
        #C = np.arange(1e-8, 1.0, 1e-8)
        #epsilon = 1e-8
        #parameters = {'C': C}
        #scorer = make_scorer(average_precision_score)

        #gs = GridSearchCV(self.model, parameters, cv=3, n_jobs=10, scoring=scorer)
        #gs.fit(self.X, self.Y)
        #self.params = gs.best_params_.get('C', None)
        # temporarily use a fixed C for testing only
        self.params = 0.1
        self.renew_model(self.params)
        self.model.fit(self.X, self.Y)

        top = np.argsort(np.abs(self.model.coef_[0]))[::-1][:self.top]
        return self.ids[top], self.X[:, top], self.X_all[:, top]

    def _wrap_fs(self):
        C = np.arange(1e-8, 1.0, 1e-8)
        high = len(C)-1
        low = 0
        epsilon = 1e-6
        while low <= high:
            mid = (low + high) // 2
            self.renew_model(C[mid])
            self.model.fit(self.X, self.Y)
            s = np.abs(self.model.coef_[0]) > epsilon
            snum = np.sum(s)
            if snum == self.top:
                print('feature ready', len(self.ids[s]))
                return self.ids[s], self.X[:, s], self.X_all[:, s]
            elif snum > self.top+2:
                # too many samples, prefer stronger regularization, lower C
                high = mid - 1
            elif snum < self.top-2:
                low = mid
            else:
                break
        print('feature ready', len(self.ids[s]))
        return self.ids[s], self.X[:, s], self.X_all[:, s]

    def _linear_wrap_fs(self):
        N = self.X.shape[1]

        def adjusted_rsq(y, y_pred):
            SS_Residual = sum((y-y_pred)**2)
            SS_Total = sum((y-np.mean(y))**2)
            r_squared = 1 - float(SS_Residual)/SS_Total
            adjusted_r_squared = 1 - (1-r_squared)*(len(y)-1)/(len(y)-N-1)
            return adjusted_r_squared

        C = np.arange(1e-5, 0.01, 1e-5)
        high = len(C)-1
        low = 0
        epsilon = 1e-8
        while low <= high:
            mid = (low + high) // 2
            self.renew_model(C[mid])
            self.model.fit(self.X, self.Y)
            s = np.abs(self.model.coef_) >= epsilon
            snum = np.sum(s)
            if snum == self.top:
                return self.ids[s], self.X[:, s]
            elif snum < self.top-20: # default top 30, search range 45 ~ 10
                high = mid - 1
            elif snum > self.top+5:
                low = mid + 1
            else:
                break
            print(snum, low, high, mid, C[mid])
            print(adjusted_rsq(self.Y, self.model.predict(self.X)))
            print (self.model.coef_)
        return self.ids[s], self.X[:, s]

    def up_sample_minor_class(self):
        rand = np.random.RandomState(self.seed)
        pos, = np.where(self.Y == 1)
        neg, = np.where(self.Y == 0)

        up_sample_idx = rand.choice(pos, len(neg) - len(pos), replace=True)
        up_sample_y = np.hstack((self.Y, self.Y[up_sample_idx]))
        up_sample_x = np.vstack((self.X, self.X[up_sample_idx]))
        # add gaussian noise
        up_sample_x += rand.randn(*up_sample_x.shape)

        self.X = up_sample_x
        self.Y = up_sample_y
        self.symbols = np.hstack((self.symbols, self.symbols[up_sample_idx]))
        self.refseqs = np.hstack((self.refseqs, self.refseqs[up_sample_idx]))

    def down_sample_major_class(self):
        rand = np.random.RandomState(self.seed)
        pos, = np.where(self.Y == 1)
        neg, = np.where(self.Y == 0)

        down_sample_idx = rand.choice(neg, 500, replace=False) # how many background genes
        down_sample_y = np.hstack((self.Y[pos], self.Y[down_sample_idx]))
        down_sample_x = np.vstack((self.X[pos], self.X[down_sample_idx]))
        self.X = down_sample_x
        self.Y = down_sample_y
        self.symbols = np.hstack((self.symbols[pos], self.symbols[down_sample_idx]))
        self.refseqs = np.hstack((self.refseqs[pos], self.refseqs[down_sample_idx]))

    def get_dnase_most_info_id(self):
        self.fit()
        return self.sid[np.argsort(np.abs(self.coefs))][::-1]

    def linear_fit(self):
        N = self.X.shape[1]

        def adjusted_rsq(y, y_pred):
            SS_Residual = sum((y-y_pred)**2)
            SS_Total = sum((y-np.mean(y))**2)
            r_squared = 1 - float(SS_Residual)/SS_Total
            adjusted_r_squared = 1 - (1-r_squared)*(len(y)-1)/(len(y)-N-1)
            return adjusted_r_squared

        self.sid, self.X = self._linear_wrap_fs()
        # self.model = LinearRegression(n_jobs=8)
        # self.model.fit(self.X, self.Y)

        parameters = {'alpha': np.linspace(1e-4,0.1,20)}
        scorer = make_scorer(r2_score)
        # scorer = make_scorer(adjusted_rsq)

        self.renew_model()
     
        gs = GridSearchCV(self.model, parameters, cv=3, n_jobs=10, scoring=scorer)
        gs.fit(self.X, self.Y)
        self.params = gs.best_params_.get('alpha', None)
        self.renew_model(self.params)
        self.model.fit(self.X, self.Y)

        self.coefs = self.model.coef_
        marge_rp = self.model.predict(self.X)
        self.performance = adjusted_rsq(self.Y, marge_rp)

        # fpr, tpr, thresholds = roc_curve(self.Y, marge_rp, pos_label=1)
        # self.tpr = tpr
        # self.fpr = fpr
        # self.marge_rp_dict = {}
        # for v, k in zip(marge_rp, self.symbols):
        #     self.marge_rp_dict[k.upper()] = v


    def fit(self):
        """ fit a logit regression model with gene rp and differential gene set
        """
        if self.penalty == 'l1':
            self.sid, self.X, self.X_all = self._wrap_fs()
        elif self.penalty == 'l2':
            self.sid, self.X, self.X_all = self._wrap_l2_fs()
        else:
            raise

        #self.up_sample_minor_class() # even worse
        parameters = {'C': np.linspace(1e2,1e4,100)} # larger C means smaller lambda, tend to preserve the specified top selected samples.
        #scorer = make_scorer(roc_auc_score)
        scorer = make_scorer(average_precision_score)
        self.renew_model()
        gs = GridSearchCV(self.model, parameters, cv=3, n_jobs=5, scoring=scorer)
        gs.fit(self.X, self.Y)
        self.params = gs.best_params_.get('C', None)
        self.renew_model(self.params)
        self.model.fit(self.X, self.Y)
        self.coefs, = self.model.coef_
        marge_rp_unique = self.model.predict_log_proba(self.X)[:,1]
        marge_rp = self.model.predict_log_proba(self.X_all)[:,1]

        self.performance = (roc_auc_score(self.Y, marge_rp_unique),
                            average_precision_score(self.Y, marge_rp_unique))

        fpr, tpr, thresholds = roc_curve(self.Y, marge_rp_unique, pos_label=1)
        self.tpr = tpr
        self.fpr = fpr
        self.marge_rp_dict = {}
        for v, k in zip(marge_rp, self.all_refseq):
            self.marge_rp_dict[k.upper()] = v

    def __getstate__(self):
        self.down_sample_major_class()
        return (self.X, self.Y, self.symbols, self.refseqs, self.sid, self.coefs, self.params, self.performance, self.marge_rp_dict, self.tpr, self.fpr)

    def __setstate__(self, state):
        self.X, self.Y, self.symbols, self.refseqs, self.sid, self.coefs, self.params, self.performance, self.marge_rp_dict, self.tpr, self.fpr = state

    def test(self, factor, h5, ct_Y, other_h5):
        """
        h5: a HDF object for loading read count
        ct_Y: a binary numpy array
           gold standard, chip-seq peak or motif site
        """
        if other_h5 != None:
            original_order_ids = self.sid
            # separate self.sid into 2 groups
            with h5py.File(other_h5) as store:
                self_ids = np.array([ i.decode("utf-8") for i in store["IDs"][...] ])

            # 1. selected in-house ids
            in_index, = np.where(np.in1d(original_order_ids, self_ids))
            in_house = original_order_ids[in_index]
            if len(in_house) > 0:
                with h5py.File(other_h5) as store:
                    map_id = {}
                    for i, c in enumerate(self_ids):
                        map_id[c] = i
                    idx = np.array([ map_id[str(i)] for i in in_house ])
                    # hdf index must be increasing order
                    idx_sort = np.sort(idx)
                    mat = store['OrderCount'][:, idx_sort]
                    # map back to original iid order
                    map_ids = {}
                    for ii, i in enumerate(idx_sort):
                        map_ids[i] = ii
                    idx = np.array([ map_ids[i] for i in idx ])
                    self_data = mat[:, idx]
                    
                    # 2. selected public ids
                    public_index, = np.where(np.in1d(original_order_ids, self_ids, invert=True))
                    public_data = original_order_ids[public_index]
                    public_ct = h5.count_hdf_reader(factor, iid=public_data.tolist())
                    ct = np.hstack([public_ct, self_data])
                    ids = np.concatenate([public_data, in_house])
            else:
                public_index, = np.where(np.in1d(original_order_ids, self_ids, invert=True))
                ids = original_order_ids[public_index]
                ct = h5.count_hdf_reader(factor, iid=ids.tolist())

            # mapping back to original order for coefficient match
            map_ids = {}
            for ii, i in enumerate(ids):
                map_ids[i] = ii
            index = np.array([ map_ids[i] for i in original_order_ids ])
            ct = ct[:, index]
        else:
            ct = h5.count_hdf_reader(factor, iid=self.sid.tolist())

        self.Y_count = ct

        if ct_Y != None:
            pred = np.dot(ct, self.coefs)
            self.ct_performance = (roc_auc_score(ct_Y, pred),
                                   average_precision_score(ct_Y, pred))
            return self.ct_performance
        return 'no chip-seq'

    def get_refseq_bin_info(self, ann_obj):
        """ get refseq TSS bin information """
        ann = {}
        # a: (GeneTSS, GeneBin)
        for a in ann_obj:   
            # TypeError: iter() returned non-iterator of type 
            # Python3 rename next to __next__
            # refseq => (Chrom, TSS, bin center, bin index)
            refseq = a[0].label.split(':')[0] # support cluster, uniform interface
            ann[refseq] = (a[0].chrom, a[0].start, a[1].get_center(), a[1].label)

        chrs= []; tss = []; idx = []; bin_centers =[]
        for r in self.refseqs:
            chrs.append(ann[r][0])
            idx.append(int(ann[r][3])-1)
            tss.append(int(ann[r][1]))
            bin_centers.append(int(ann[r][2]))
        return chrs, np.array(idx, np.int32), \
               np.array(tss, np.int32), \
               np.array(bin_centers, np.int32)

    def cluster_motif_score(self, motif_cutoffs_h5, cutoffs, prefix, ann_obj, chrom_bin, true_y, config, mapping):
        """
        region based delta RP
        cluster histone mark or DNase 1kb read count (self.Ycount )
        then compute 99% motif hits occurence in cluster and across cluster
        """
        near_gene_bin = []; boundary_mask = []
        diff = []
        geneset = ann_obj.get_gene_set()
        for a in ann_obj:
            ref, sym = a[0].label.split(':')
            if (ref in geneset) or (sym in geneset):
                diff.append(1)
            else:
                diff.append(0)
            near_gene_bin.append(int(a[1].label)-1)
            boundary_mask.append(chrom_bin[a[0].chrom])

        diff = np.array(diff) # gene
        near_gene_bin = np.array(near_gene_bin) # gene,
        boundary_mask = np.array(boundary_mask) # gene,

        idx_2d = []
        for i in range(-100, 100):
            idx_2d.append(near_gene_bin + i)
        idx_2d = np.vstack(idx_2d) # bin x gene
        diff = np.repeat(diff.reshape(1, len(diff)), 200, axis=0) # bin x gene

        idx_2d_right = (idx_2d - boundary_mask)<=0 # bin x gene
        idx_2d_left = idx_2d >= 0
        boundary_mask = idx_2d_left & idx_2d_right

        idx_2d = idx_2d * boundary_mask
        diff = diff * boundary_mask # mask out of boundary

        idx = idx_2d.T.reshape(-1) # (gene x bin, )
        diff = diff.T.reshape(-1)  # (gene x bin, )

        idx, uindex = np.unique(idx, return_index=True) # remove duplicates of bins
        diff = diff[uindex][1:] # the first element is the masked bin
        idx = idx[1:]

        norm = log_transform(self.Y_count[idx])
        if isinstance(true_y, np.ndarray):
            TFBS_near_gene = np.dot(norm, self.coefs)
            print('near gene', roc_auc_score(true_y[idx], TFBS_near_gene), average_precision_score(true_y[idx], TFBS_near_gene))

        # with 10 samples, do not need to adjust
        log_count = norm * np.abs(self.coefs) # so many examples, adjust by coefficient

        # filter regions to get the HPR (high plasticity region, not sure this help)
        # TODO: add TSS and CGI filter
        # np.sd
        #row_sd = np.std(norm, axis=1) 
        #filter_index, = np.where(row_sd >= np.percentile(row_sd,99))
        #log_count = log_count[filter_index]
   
        n_clusters = 10
        kmeans = MiniBatchKMeans(n_clusters=n_clusters, random_state=0).fit(log_count)

        # ChIP-seq peak density in clusters
        if isinstance(true_y, np.ndarray):
            tf_pos = np.asarray(true_y[idx], dtype=np.bool)
            #tf_pos = np.asarray(true_y[idx][filter_index], dtype=np.bool)

        #diff = np.asarray(diff, dtype=np.bool)
        labels_order = np.argsort(kmeans.labels_)

        # order by cluster labels
        
        if isinstance(true_y, np.ndarray):
            tf_pos = tf_pos[labels_order] # tf location
        #diff = diff[labels_order]     # differential genes
        log_count = log_count[labels_order]
        labels = kmeans.labels_[labels_order]
        cluster_stat = np.unique(labels, return_counts=True)

        prop_dict = {}

        cluster_dict = {}
        for cluster, ct in zip(*cluster_stat):
            cluster_dict[cluster] = ct

        centers = []
        for n in range(n_clusters):
            centers.append(np.mean(log_count[labels==n,:], axis=0))

        cluster_mean = pd.DataFrame(np.vstack(centers), columns=self.sid, index=['cluster%s' %s for s in range(n_clusters)])
        # cluster x features
        cluster_mean.to_csv(prefix + '_cluster_aver_signal.csv')

      
        if isinstance(true_y, np.ndarray):
            tf_stat = {}
            for i, j in zip(*np.unique(labels[tf_pos], return_counts=True)):
                tf_stat[i] = j
            # prop_dict['chipseq'] = [1.0*np.sum(tf_pos)/len(tf_pos)]
            prop_dict['chipseq'] = []
            for key in cluster_dict:
                if key in tf_stat:
                    prop_dict['chipseq'].append(tf_stat[key])
                else:
                    prop_dict['chipseq'].append(0)
    
        with h5py.File(config.tf_chipseq) as chip_occur:
            chip_ids = chip_occur['IDs']
            for mf in chip_ids:
                #print(mf)
                chip_region = chip_occur[mf][...] -1 # -1 to get 0-based index of peak window
                chip = np.zeros(mapping[-1]+1, dtype=np.int32) # 1kb window
                chip[mapping[chip_region]] = 1 # 1kb 0-1 vector

                is_s = np.asarray(chip, dtype=np.bool)
                is_s = is_s[idx] # near the gene chip-seq peak window
                #is_s = is_s[idx][filter_index] # get top 1% bin
                is_s = is_s[labels_order] # order by cluster label
                labels_chip = labels[is_s] # get the labels for the motif hit
                chip_stat = {}
                for i, j in zip(*np.unique(labels_chip, return_counts=True)):
                    chip_stat[i] = j
                prop_dict[mf] = []
                for key in cluster_dict:
                    prop_dict[mf].append(chip_stat.get(key,0))

        for cutoff, motif_cutoff in zip(cutoffs, motif_cutoffs_h5):
           with h5py.File(motif_cutoff) as motif_occur:
               # motif_region = motif_occur['IsMotifRegion%d' % cutoff] # actually bool or operation not faster than the index intersection method
               motif_ids = motif_occur['IDs']
               for mf in range(len(motif_ids)):
                   #print("%d\t%s"%(cutoff, motif_ids[mf]))
                   motif_region = motif_occur[motif_ids[mf]][...] # motif is 0-based regions 
                   motif = np.zeros(mapping[-1]+1, dtype=np.int32) # 1kb window
                   motif[mapping[motif_region]] = 1 # 1kb 0-1 vector
                   is_s = np.asarray(motif, dtype=np.bool) # use 100bp index to convert and map to 1kb
                   is_s = is_s[idx] # near the gene motif
                   #is_s = is_s[idx][filter_index] # near the gene motif
                   is_s = is_s[labels_order] # order by cluster label
                   labels_motif = labels[is_s] # get the labels for the motif hit
                   motif_stat = {}
                   for i, j in zip(*np.unique(labels_motif, return_counts=True)):
                       #print(i,j)
                       motif_stat[i] = j
                   #prop_dict[motif_ids[mf]] = [1.0*np.sum(is_s)/len(is_s)]
                   prop_dict[motif_ids[mf]] = []
                   for key in cluster_dict:
                       prop_dict[motif_ids[mf]].append(motif_stat.get(key,0))

        prop_dict['cluster_number'] = []
        for key in cluster_dict:
            prop_dict['cluster_number'].append(cluster_dict[key])

        # cluster x (chipseq + motif)
        #prop =  pd.DataFrame(prop_dict, index=['background']+['cluster%s' %s for s in range(10)])
        prop =  pd.DataFrame(prop_dict, index=['cluster%s' %s for s in cluster_dict])
        prop.to_csv(prefix + "_cluster_stat_kl_divergence.csv")

    def get_deletion_rp(self, motif_or_chip, chrom_bin, ann_obj, promoter=False, DNase=False):
        """ intermediate result for deleting motif or chip peak
        motif_or_chip: a binary array, 1 for TF region in 1kb window

        S: H3K27ac      e.g 20 sample x 2000 gene x 200 bin
        E: motif occurence matrix, e.g. 2000 gene x 200 bin
        W: adjuted weight matrix,  e.g. 2000 gene x 200 bin
        """
        chrs, idx, tss, bc = self.get_refseq_bin_info(ann_obj)
        boundary_mask = np.array([ chrom_bin[c] for c in chrs ])

        bc_2d = []; idx_2d = []
        for i in range(-100, 100):
            # bin center
            bc_2d.append(bc + i*ann_obj.weight.bl)
            # TSS bin index array
            idx_2d.append(idx + i)
        bc_2d = np.vstack(bc_2d).T
        idx_2d = np.vstack(idx_2d).T

        idx_2d_left = idx_2d >= 0
        idx_2d_right = (idx_2d.T - boundary_mask).T <= 0
        boundary_mask = idx_2d_left & idx_2d_right

        distances = (bc_2d.T - tss).T
        W = ann_obj.weight.balance_weight(distances)

        idx_2d = idx_2d * boundary_mask
        # mask element got the 0th element
        E = motif_or_chip[idx_2d]
        # E = np.logical_xor(E, np.ones_like(E))
        E = np.logical_not(E)
        # gene x bin x sample
        S = self.Y_count[idx_2d]
        # sample x gene x bin
        S = S.transpose(2, 0, 1)
        # another fix! multiply window size
        S = ann_obj.weight.bl * S

        # in case mask 0th bin is not 0, multiple mask again
        deletion_rp = S * E * W * boundary_mask
        # sample x gene
        deletion_rp = np.sum(deletion_rp, axis=2)
        # gene   x sample
        return deletion_rp.T

    def get_delta_rp(self, del_rp):
        """ Motif Delta RP: f(log_RP) - f(log_deletion_rp)
        """

        log_del_rp = log_transform(del_rp)
        delta = np.dot(self.X, self.coefs) - np.dot(log_del_rp, self.coefs)
        return delta

    def get_delta_rp_openmp(self, count_Y, motif_h5, chrom_bin, ann_obj, cutoff=95, DNase=False, window_map=None, dc_chip=True):
        """ reduce dimension from tensor4 to matrix for openmp """
        chrs, idx, tss, bc = self.get_refseq_bin_info(ann_obj)
        boundary_mask = np.array([ chrom_bin[c] for c in chrs ])
        bc_2d = []; idx_2d = []
        for i in range(-100, 100):
            # bin center
            bc_2d.append(bc + i*ann_obj.weight.bl)
            # TSS bin index array
            idx_2d.append(idx + i)

        bc_2d = np.vstack(bc_2d).T
        idx_2d = np.vstack(idx_2d).T

        idx_2d_left = idx_2d >= 0
        idx_2d_right = (idx_2d.T - boundary_mask).T <= 0
        boundary_mask = idx_2d_left & idx_2d_right

        distances = (bc_2d.T - tss).T
        W = ann_obj.weight.balance_weight(distances)
        # gene x bin
        idx_2d = idx_2d * boundary_mask
        # gene x bin x sample
        S = self.Y_count[idx_2d]
        # sample x gene x bin
        S = S.transpose(2, 0, 1)
        # another fix! multiply window size
        S = ann_obj.weight.bl * S

        if isinstance(DNase, np.ndarray):
            has_dnase = True
            if count_Y == None: # no chip-seq avaible as gold standard
                count_Y = None
            else:
                # chip-seq to be consistent with motif process
                # overlap chip-seq with DNase first
                # this is for the corresponding ChIP-seq of differential expression data in the same paper, 0-based already, look at scripts/lisa !! DNase 0-based
                count_Y_index = np.intersect1d(DNase, count_Y) # overlap with dnase at 100bp window 
                count_Y = np.zeros(window_map[-1]+1, dtype=np.int32)
                count_Y[window_map[count_Y_index]] = 1
        elif isinstance(DNase, bool):
            has_dnase = False
        else:
            raise

        marge_rp = np.dot(self.X, self.coefs)
        precompute = S * W * boundary_mask # nice!, precompute signal * weight * boundary masking

        # reshape 3-D to get 2-D array to openmp parallel operation
        precompute = \
            precompute.reshape((precompute.shape[0], -1)) # sample x (gene, bin)

        mode = theano.Mode(linker='cvm', optimizer='fast_run')
        theano.config.exception_verbosity = 'high'
        theano.config.openmp = True
        ##theano.config.blas.ldflags = ' -lopenblas '
        ##theano.config.openmp_elemwise_minsize = 5000 # minmal element number to use openmp
        T = theano.tensor
        f = theano.function
        x = T.imatrix('E') # each motif tensor
        precomp = theano.shared(precompute.astype(theano.config.floatX), name='precompute')
        r = theano.shared(marge_rp.astype(theano.config.floatX), name='MARGE RP')
        c = theano.shared(self.coefs.astype(theano.config.floatX), name='coefficients')

        # sample x (gene,bin)
        y = T.extra_ops.repeat(x, precompute.shape[0], axis=0)
        tensor_del = y * precomp # sample x (gene,bin) # for openmp broadcast
        tensor_del = T.reshape(tensor_del, (c.shape[0],r.shape[0],200)) # sample x gene x bin
        tensor_del = T.transpose(T.sum(tensor_del, axis=2), (1,0)) + T.constant(1) # one motif
        tensor_del_med = T.mean(tensor_del, axis=0)  # one motif
        log_tensor_del = T.log2(tensor_del) - T.log2(tensor_del_med)
        tensor_delta = r - T.dot(log_tensor_del, c)
        theano_delta_rp = f([x], tensor_delta, mode=mode)

        if count_Y != None:
            #chip = np.logical_not(count_Y[idx_2d]) # gene x bin
            chip = np.logical_not(count_Y[idx_2d].reshape(1, -1)) # 1x(gene,bin)
            chip = theano_delta_rp(chip)

        with h5py.File(motif_h5) as motif_occur:
            if dc_chip: # cistrome dc tf chip-seq
                IDs = motif_occur['IDs'][...]
                TFs = IDs
            else:  # cistrome motif
                IDs = motif_occur['IDs'][...]
                TFs = motif_occur['TFs'][...]
            motif_num = len(IDs)
            for mf in range(0, motif_num):
                if has_dnase:
                    if dc_chip: # cistrome dc tf chip-seq 1-based index from scripts/tf_100bp.py
                        motifp = np.zeros(window_map[-1] + 1, dtype=np.int32)
                        # intersect dnase 0-based with chip-seq peak 0-based at 100bp windows, and map back to 1kb bins
                        dnase_motif = np.intersect1d(DNase, motif_occur[IDs[mf]][...]-1) # 1-based index -1 to 0-based index
                        dnase_motif = window_map[dnase_motif]
                        motifp[dnase_motif] = 1
                    else: # cistrome motif 0-based index from scripts/motif_100bp.py through np.where
                        # intersect dnase 0-based with motif 0-based at 100bp windows, and map back to 1kb bins
                        motifp = np.zeros(window_map[-1] + 1, dtype=np.int32)
                        dnase_motif = np.intersect1d(DNase, motif_occur[IDs[mf]][...])
                        dnase_motif = window_map[dnase_motif]
                        motifp[dnase_motif] = 1
                else:
                    if dc_chip: # cistrome dc tf chip-seq 1-based index from scripts/tf_100bp.py
                        motifp = np.zeros(window_map[-1] + 1, dtype=np.int32)
                        motif_map = window_map[motif_occur[IDs[mf]][...]-1] # chip-seq 1-based window -1
                        motifp[motif_map] = 1
                    else:
                        motifp = np.zeros(window_map[-1] + 1, dtype=np.int32)
                        motif_map = window_map[motif_occur[IDs[mf]][...]] # motif 0-based window
                        motifp[motif_map] = 1

                E = motifp[idx_2d] # gene x bin
                E = E.reshape(1, -1) # openmp, 1x(gene,bin)

                batch_ids = IDs[mf]
                batch_tfs = TFs[mf]

                E = np.logical_not(E)

                delta = theano_delta_rp(E)
                if mf == 0:
                    if count_Y != None:
                        yield (['peak', batch_ids], ['chipseq', batch_tfs], [chip, delta])
                    else:
                        yield ([batch_ids], [batch_tfs], [delta])
                else:
                    yield ([batch_ids], [batch_tfs], [delta])

    def get_mini_batch_delta_rp(self, count_Y, motif_h5, chrom_bin, ann_obj, cutoff=95, DNase=False, window_map=None, theano_mode=True):
        """
        idx_2d: numpy array of 2d index for selecting H3K27ac count from self.Y_count
        self.Y_count => S
        motif_or_chip => E

        E: motif occurence matrix, e.g. 1000  motif x 1         x 2000 gene x 200 bin
        S: H3K27ac                 e.g               20(sample) x 2000 gene x 200 bin
        W: adjuted weight matrix,  e.g.                           2000 gene x 200 bin
        """
        chrs, idx, tss, bc = self.get_refseq_bin_info(ann_obj)
        boundary_mask = np.array([ chrom_bin[c] for c in chrs ])
        bc_2d = []; idx_2d = []
        for i in range(-100, 100):
            # bin center
            bc_2d.append(bc + i*ann_obj.weight.bl)
            # TSS bin index array
            idx_2d.append(idx + i)

        bc_2d = np.vstack(bc_2d).T
        idx_2d = np.vstack(idx_2d).T

        idx_2d_left = idx_2d >= 0
        idx_2d_right = (idx_2d.T - boundary_mask).T <= 0
        boundary_mask = idx_2d_left & idx_2d_right

        distances = (bc_2d.T - tss).T
        W = ann_obj.weight.balance_weight(distances)

        # gene x bin
        idx_2d = idx_2d * boundary_mask

        #if promoter:
        #    # mask TSS, +/- 1bin
        #    idx_2d[:, 99:101] = 0

        # gene x bin x sample
        S = self.Y_count[idx_2d]
        # sample x gene x bin
        S = S.transpose(2, 0, 1)
        # another fix! multiply window size
        S = ann_obj.weight.bl * S

        if isinstance(DNase, np.ndarray):
            has_dnase = True
            if count_Y == None: # no chip-seq avaible as gold standard
                count_Y = None
            else:
                # chip-seq to be consistent with motif process
                # overlap chip-seq with DNase first
                count_Y_index = np.intersect1d(DNase, count_Y)
                count_Y = np.zeros(3209513)
                count_Y[window_map[count_Y_index]] = 1
        elif isinstance(DNase, bool):
            has_dnase = False
        else:
            raise

        marge_rp = np.dot(self.X, self.coefs)
        precompute = S * W * boundary_mask # nice!, precompute signal * weight * boundary masking
        # manual specify theano mode, use c implementation with faster mode
        if theano_mode: # if gpu is availble, turn on gpu in .theanorc
            T = theano.tensor
            f = theano.function

            tensor4_broad2way = T.TensorType('int32', (False, True, False, False))
            x = tensor4_broad2way('E') # a batch of motifs tensor

            precomp = theano.shared(precompute.astype(theano.config.floatX), name='precompute')
            r = theano.shared(marge_rp.astype(theano.config.floatX), name='MARGE RP')
            c = theano.shared(self.coefs.astype(theano.config.floatX), name='coefficients')

            tensor_del = x * precomp
            tensor_del = T.transpose(T.sum(tensor_del, axis=3), (0,2,1)) + T.constant(1) # a batch of motif
            tensor_del_med = T.mean(tensor_del, axis=1) # no theano median ops, a batch of motif
            log_tensor_del = T.log2(tensor_del) - T.log2(tensor_del_med.dimshuffle(0, 'x', 1) )
            # gene x sample %*% sample, - motif x gene x sample %*% sample,
            tensor_delta = r - T.dot(log_tensor_del, c) # a batch of motif
            theano_delta_rp = f([x], tensor_delta)

        with h5py.File(motif_h5) as motif_occur:
            # h5py does not support two dimensional index
            # fetch all elements into 2d numpy array
            # use 2-d index to get 3-d motif occurence array

            # whole batch motif, numpy array method
            # tried for 1000 genes, 21 samples, 1061 motifs
            # consume too much memory, over 70GB, not plausible
            batch_size = 25
            IDs = motif_occur['IDs'][...]
            TFs = motif_occur['TFs'][...]
            motif_num = len(IDs)
            if not has_dnase: # extract 0-1 motif hit matrix for 1kb windows
                motifp = motif_occur['IsMotifRegion%d' % cutoff]
            for mf in range(0, motif_num, batch_size):
                if has_dnase: # intersect dnase with motif at 100bp windows, and map back to 1kb bins
                    if motif_num - mf < batch_size: # left motifs is smaller than batch_size
                        motifp = np.zeros((3209513, motif_num - mf)) # 3209513 bins number for 1kb windows
                        batch_upper = motif_num
                    else:
                        motifp = np.zeros((3209513, batch_size))
                        batch_upper = mf+batch_size
                    for ii, i in enumerate(range(mf, batch_upper)):
                        dnase_motif = np.intersect1d(DNase, motif_occur[IDs[i]][...])
                        dnase_motif = window_map[dnase_motif]
                        ##print dnase_motif.shape
                        motifp[dnase_motif, ii] = 1
                    E = motifp[idx_2d]
                else:
                    # try mini-batch of motifs instead
                    # gene x bin x mini-batch motif
                    E = motifp[:,mf:(mf+batch_size)][idx_2d]
                #print E.shape

                # mini-batch motif x gene x bin
                E = E.transpose(2, 0, 1)
                #print E.shape

                batch_ids = IDs[mf:(mf+batch_size)]
                batch_tfs = TFs[mf:(mf+batch_size)]

                if mf == 0 and count_Y != None: # insert ChIP-seq peak at 0th position
                #   1 x gene x bin
                    batch_ids = np.hstack(['peak', batch_ids])
                    batch_tfs = np.hstack(['chipseq', batch_tfs])
                    chip = count_Y[idx_2d][np.newaxis, :, :]
                    #(mini-batch motif+1) x gene x bin
                    E = np.concatenate([chip, E], axis=0)

                # motif x 1 x gene x bin
                E = E[:, np.newaxis, :, :]
                # mask element got the 0th element
                E = np.logical_not(E)
                # print E.shape

                # in case mask 0th bin is not 0, multiple mask again
                # actually 0th bin always be 0, because chromosome end
                # motif x sample x gene x bin
                if not theano_mode:
                    #deletion_rp = E * S * W * boundary_mask
                    deletion_rp = E * precompute

                    # motif x sample x gene ->tranpose-> motif x gene x sample
                    deletion_rp = np.sum(deletion_rp, axis=3).transpose(0, 2, 1)

                    pcount = 1
                    xt = deletion_rp + pcount
                    # motif x sample
                    #med = np.median(xt, axis=1)
                    med = np.mean(xt, axis=1) # use mean for theano compability

                    # motif x gene x sample - motif x 1 x sample
                    log_del_rp = np.log2(xt) - np.log2(med[:, np.newaxis, :])
                    # gene x sample %*% sample, - motif x gene x sample %*% sample,
                    delta = marge_rp - log_del_rp.dot(self.coefs)
                    # motif x gene
                else:
                    delta = theano_delta_rp(E)
                yield (batch_ids, batch_tfs, delta)


def add_node(node, parent):
    # boundary case, sometimes chip-seq peak beta cluster get Infinity distance
    if np.isinf(node.dist): 
        dist = 0.999
    else:
        dist = node.dist
    newNode = dict( node_id = node.id, children = [], length = dist )
    parent["children"].append( newNode )
    if node.left: add_node(node.left, newNode)
    if node.right: add_node(node.right, newNode)
