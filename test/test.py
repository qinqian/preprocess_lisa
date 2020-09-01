import unittest
from marge2_conf import Config
from marge2 import Annotation, Region
import numpy as np

class TestAnnotation(unittest.TestCase):
    def setUp(self):
        self.conf = Config('/data/home/qqin/MARGE/scripts/marge2_rawmotif.conf')
        self.bl = 1
        self.ann = Annotation(self.conf, '/data/home/qqin/12_data/siRNA/AR_up_genes.txt', '/data/home/qqin/MARGE/PhaseA_motif_deltarp/43060_treat.bw', self.bl)
        self.iid = 43060

    def test_getrp(self):
        self.assertEqual(self.conf.get('basics', 'meta'), '/data/home/qqin/12_data/margeK27ac.csv')
        h5 = HDF(self.conf)
        n = 0
        for a in self.ann:
            if a:
                n+=1
                if n>5: break
                rp = self.ann.get_rp(a[0])
                self.assertTrue(abs(rp[0] - h5.rp_hdf_reader(iid=self.iid, gene=a[0].label)) < 50)
                print rp[0], h5.rp_hdf_reader(iid=self.iid, gene=a[0].label)

    def test_deletion(self):
        h5 = HDF(self.conf)

        print self.ann.get_rp(Region('chr21', 32728040, 32728041, 'NM_203446'), Region('chr21', 32720000, 32732000, 12000))
        index = [3074291, 3074303]
        ct = h5.count_hdf_reader(iid=self.iid, index=index)
        val = self.ann.get_delta_rp(Region('chr21', 32728040, 32728041, 'NM_203446'), Region('chr21', 32720000, 32732000, 12))
        print 'original\tdelete one true\tdelete one estimation'
        est = 1000 * np.dot(val[0], ct) # from hdf5 count, val[1] from bigwig
        print "%s\t%s\t%s" %(h5.rp_hdf_reader(iid=self.iid, gene='NM_203446'), val[1], est)

        print 'index\tregion\tdelete one true\tdelete one estimation'
        step = 1
        s = 0; z = 0
        for i in range(3074291, 3074303, step):
            index = [i, i+step]
            i = i - 3074291
            val = self.ann.get_delta_rp(Region('chr21', 32728040, 32728041, 'NM_203446'), Region('chr21', 32720000+i*1000, 32720000+(i+step)*1000, step))
            ct = h5.count_hdf_reader(iid=self.iid, index=index)
            est = 1000 * val[0] * np.sum(ct) # from hdf5 count, val[1] from bigwig
            s += est; z+=val[1][0]
            print "%s\t%.50s\t%s\t%s" % (str(index), str((32720000+i*1000, 32720000+(i+step)*1000)), val[1][0], est[0])
        print s,z
    def test_deletion2(self):
        h5 = HDF(self.conf)

        print self.ann.get_rp(Region('chr3', 156674415, 156674416, 'NM_015508'), Region('chr3', 156670000, 156680000, 10000))
        index = [647821, 647831]
        ct = h5.count_hdf_reader(iid=self.iid, index=index)
        val = self.ann.get_delta_rp(Region('chr3', 156674415, 156674416, 'NM_015508'), Region('chr3', 156670000, 156680000, 10))
        print 'original\tdelete one true\tdelete one estimation'
        est = 1000 * np.dot(val[0], ct) # from hdf5 count, val[1] from bigwig
        print "%s\t%s\t%s" %(h5.rp_hdf_reader(iid=self.iid, gene='NM_015508'), val[1], est)

        print 'index\tregion\tdelete one true\tdelete one estimation'
        step = 1
        s = 0; z = 0
        for i in range(647821, 647832, step):
            index = [i, i+step]
            i = i - 647821
            val = self.ann.get_delta_rp(Region('chr3', 156674415, 156674416, 'NM_015508'), Region('chr3', 156670000+i*1000, 156670000+(i+step)*1000, step))
            ct = h5.count_hdf_reader(iid=self.iid, index=index)
            est = 1000 * val[0] * np.sum(ct) # from hdf5 count, val[1] from bigwig
            s += est; z+=val[1][0]
            print "%s\t%.50s\t%s\t%s" % (str(index), str((156670000+i*1000, 156670000+(i+step)*1000)), val[1][0], est[0])
        print s,z
    def test_deletion3(self):
        h5 = HDF(self.conf)

        gene = Region('chr12', 121274585, 121274586, 'NM_001270486')
        b = Region('chr12', 121290000, 121300000, 10000)
        b_est = Region('chr12', 121290000, 121300000, 10)
        start = 121290000
        index = [2221103, 2221113]

        print self.ann.get_rp(gene, b)

        ct = h5.count_hdf_reader(iid=self.iid, index=index)
        val = self.ann.get_delta_rp(gene, b_est)

        print 'original\tdelete one true\tdelete one estimation'
        est = 1000 * np.dot(val[0], ct) # from hdf5 count, val[1] from bigwig
        print "%s\t%s\t%s" %(h5.rp_hdf_reader(iid=self.iid, gene='NM_001270486'), val[1], est)

        print 'index\tregion\tdelete one true\tdelete one estimation'
        step = 1
        s = 0; z = 0
        for i in range(2221103, 2221113, step):
            index = [i, i+step]
            i = i - 2221103
            val = self.ann.get_delta_rp(gene, Region('chr12', start+i*1000, start+(i+step)*1000, step))
            ct = h5.count_hdf_reader(iid=self.iid, index=index)

            est = 1000 * val[0] * np.sum(ct) # from hdf5 count, val[1] from bigwig
            s += est; z+=val[1][0]
            print "%s\t%.50s\t%s\t%s" % (str(index), str((start+i*1000, start+(i+step)*1000)), val[1][0], est[0])
        print s,z


    def test_deletion2(self):
        h5 = HDF(self.conf)

        gene = Region('chr12', 121274585, 121274586, 'NM_001270486')
        b = Region('chr12', 121270000, 121280000, 10000)
        b_est = Region('chr12', 121270000, 121280000, 10)
        start = 121270000

        print self.ann.get_rp(gene, b)

        index = [2221083, 2221093]
        ct = h5.count_hdf_reader(iid=self.iid, index=index)
        val = self.ann.get_delta_rp(gene, b_est)

        print 'original\tdelete one true\tdelete one estimation'
        est = 1000 * np.dot(val[0], ct) # from hdf5 count, val[1] from bigwig
        print "%s\t%s\t%s" %(h5.rp_hdf_reader(iid=self.iid, gene='NM_001270486'), val[1], est)

        print 'index\tregion\tdelete one true\tdelete one estimation'
        step = 1
        s = 0; z = 0
        for i in range(2221083, 2221093, step):
            index = [i, i+step]
            i = i - 2221083
            val = self.ann.get_delta_rp(gene, Region('chr12', start+i*1000, start+(i+step)*1000, step))
            ct = h5.count_hdf_reader(iid=self.iid, index=index)

            est = 1000 * val[0] * np.sum(ct) # from hdf5 count, val[1] from bigwig
            s += est; z+=val[1][0]
            print "%s\t%.50s\t%s\t%s" % (str(index), str((start+i*1000, start+(i+step)*1000)), val[1][0], est[0])
        print s,z

    def test_hdf_reader(self):
        h5 = HDF(self.conf)
        self.assertEqual(h5.rp_hdf_reader(iid=None, gene=None)[3].shape, (52876, 4030))
        print h5.rp_hdf_reader(iid=[self.iid, 39148, 35378], gene=['NM_015508', 'NM_001077654']).shape
        print h5.rp_hdf_reader(iid=self.iid, gene=['NM_015508', 'NM_001077654']).shape
        print h5.rp_hdf_reader(iid=self.iid, gene='NM_001077654')
        print h5.rp_hdf_reader(iid=self.iid, gene=None).shape
        print h5.rp_hdf_reader(iid=None, gene='NM_001077654').shape
        print h5.count_hdf_reader(iid=None, index=[0,1,2,3]).shape
        print h5.count_hdf_reader(iid=None, index=[0,5]).shape
        #print h5.count_hdf_reader(iid=None, index=None).shape
        print h5.count_hdf_reader(iid=39148, index=0)
        print h5.count_hdf_reader(iid=[self.iid, 39148]).shape
        print h5.count_hdf_reader(iid=39148, index=[0,5]).shape
        print h5.count_hdf_reader(iid=[self.iid, 39148], index=[0,1,2,3]).shape
