import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import PCA
import sys
from random import random

BASES = ['A', 'C', 'G', 'T']
MER_LEN = 6

def i_to_mer(i):
    ret = [0]*MER_LEN
    for n in range(0, MER_LEN):
        ret[6-n-1] = BASES[(i >> n*2) & 3]
    return "".join(ret)

def mer_to_i(mer):
    i = BASES.index(mer[0])
    for n in range(1, len(mer)):
        i = (i << 2) | BASES.index(mer[n])
    return i

model_in = open(sys.argv[1])
seed_len = int(sys.argv[2])

kmers = list()
lv_means = list()
lv_stdvs = list()
sd_means = list()
sd_stdvs = list()
ig_lambdas = list()
weights = list()

model_in.readline()
for line in model_in:
    tabs = line.strip().split()
    kmer = tabs[0]
    lv_mean, lv_stdv, sd_mean, sd_stdv, ig_lambda, weight = map(float, tabs[1:])
    
    kmers.append(kmer)
    lv_means.append(lv_mean)
    lv_stdvs.append(lv_stdv)
    sd_means.append(sd_mean)
    sd_stdvs.append(sd_stdv)
    ig_lambdas.append(ig_lambda)
    weights.append(weight)

seed_coords = list()
colors = list()
xs = list()
ys = list()


sys.stdin.readline()
seq = sys.stdin.readline().strip()

mers = [mer_to_i(seq[e:e+MER_LEN]) for e in range(0, len(seq)-MER_LEN+1)]

seed_count = len(mers)-seed_len+1

print(seed_count)

for i in range(0, seed_count):
    seed = tuple(mers[i:i+seed_len])

    print ("%d\t%s" % (i, " ".join(([i_to_mer(s) for s in seed]) ) ))

    seed_coords.append( tuple([lv_means[e] for e in seed]) + 
                        tuple([lv_stdvs[e] for e in seed]) +
                        tuple([sd_means[e] for e in seed]))
    colors.append(i / (seed_count-1))

    #event = (lv_means[merid], lv_stdvs[merid], sd_means[merid])
    #xs.append(lv_means[merid])
    #ys.append(lv_stdvs[merid])
    #events.append(event)
    #colors.append(i / float(len(seq)-6))
    #print("%d\t%s\t%f\t%f\t%f" % ((merid, kmer) + event))


pca = PCA(np.array(seed_coords))

print(seed_count)

pca1 = [p[0] for p in pca.Y]
pca2 = [p[1] for p in pca.Y]

cmap = plt.get_cmap('inferno')
#plt.scatter(pca1, pca2, s=2)
#plt.scatter(pca1, pca2, cmap=cmap, c=colors, s=2)
plt.plot(pca1, pca2)

for i in range(0, seed_count):
    plt.annotate(
        str(i),
        xy=pca.Y[i][:2], xytext = (0, 0),
        textcoords='offset points', ha='right', va='bottom')


plt.show()
