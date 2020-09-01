from multiprocessing import Pool,Process,cpu_count
import os, glob
import pickle

def funcs(pair):
    i, j = pair
    pi = os.path.basename(i).replace('.pwm','')
    pj = os.path.basename(j).replace('.pwm','')
    f = os.popen("java -cp ape-2.0.1.jar ru.autosome.macroape.EvalSimilarity %s %s" % (i,j))
    c = f.read()
    c = c.split('\n')
    sim = c[21].split()[-1]
    f.close()
    return (pi, pj, sim)

pairs = []
pwms = glob.glob('cistrome/*pwm')
l = len(pwms)
for i in range(l):
    for j in range(i+1,l):
        pairs.append((pwms[i], pwms[j]))

print(len(pairs))

p = Pool(6)
results = []
multi = p.map_async(funcs, pairs, callback=results.append)
p.close()
p.join()
multi.wait()
with open('ape_pairs.pkl', 'wb') as s:
    pickle.dump(results, s)

