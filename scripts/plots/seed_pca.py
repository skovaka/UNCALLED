import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import PCA
import sys

BASES = ['A', 'C', 'G', 'T']

def i_to_mer(i, length):
    ret = [0]*length
    for n in range(0, length):
        ret[length-n-1] = BASES[(i >> n*2) & 3]
    return "".join(ret)

def mer_to_i(mer):
    i = BASES.index(mer[0])
    for n in range(1, len(mer)):
        i = (i << 2) | BASES.index(mer[n])
    return i


model_in = open(sys.argv[1])
length = int(sys.argv[2])
frac = float(sys.argv[3])


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

seed_means = list()

print("Generating seeds...")
for i in range(0, 4**length, int(1/frac)):
    seed = i_to_mer(i, length)
    seed_means.append([lv_means[mer_to_i(seed[j:j+6])] for j in range(0, length-6+1)])


print("Doing PCA...")
pca = PCA(np.array(seed_means))

pca1 = [p[0] for p in pca.Y]
pca2 = [p[1] for p in pca.Y]

plt.scatter(pca1, pca2, s=1, c=["C1"]*len(pca1))
plt.show()
