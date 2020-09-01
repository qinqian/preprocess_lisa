from multiprocessing import Pool,Process,cpu_count
import os, glob
def funcs(i):
    if not os.path.exists("%s.100bp" % i):
        #os.system('seqpos2 -f hg38_window100bp_both10bp.fa -p {0} -o {0}.100bp'.format(i)) 
        os.system('seqpos2 -f mm10_window100bp_both10bp.fa -p {0} -o {0}.100bp'.format(i)) 
    return i

p = Pool(6)
results = []
multi = p.map_async(funcs, glob.glob('cistrome/*pwm'),  callback=results.append)
p.close()
p.join()
multi.wait()
if multi.successful():
    print 'done'
print results
