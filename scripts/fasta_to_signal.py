import sys
from kmer_model import *

model = KmerModel(sys.argv[1])

if len(sys.argv) > 3:
    norm = (float(sys.argv[2]), float(sys.argv[3]))
else:
    norm = (1, 0)

rev = "rev" in sys.argv

seq = ""

for line in sys.stdin:
    if line[0] == ">":
        continue

    seq += line.strip()

if rev:
    seq = model.reverse_comp(seq)

i = 0
for k in model.seq_to_kmers(seq):
    print("%d\t%s\t%.2f\t%.2f\t%.2f\t%.2f" % ((i, model.i_to_kmer(k),) + model.get_norm_vals(k, norm)))
    i += 1
        
