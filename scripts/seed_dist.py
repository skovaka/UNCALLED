from read_sim import *
from kmer_model import *
import numpy as np
from itertools import chain

BASES = ['A', 'C', 'G', 'T']

def events_to_vector(events):
    #return np.array([mean for mean, stdv in events])
    return np.array(list(chain(*events)))

def i_to_kmer(i, k):
    ret = [0] * k
    for n in range(0, k):
        ret[k - n - 1] = BASES[(i >> n * 2) & 3]
    return "".join(ret)

model = KmerModel("../kmer_models/r9.2_180mv_250bps_6mer/template_median68pA.model")

seq1 = i_to_kmer(0, 13)
seq2 = i_to_kmer(3, 13)

print(seq1, seq2)

events1_1 = simulate_read(seq1, model)
events1_2 = simulate_read(seq1, model)

events2_1 = simulate_read(seq2, model)
events2_2 = simulate_read(seq2, model)

vec1_1 = events_to_vector(events1_1)
vec1_2 = events_to_vector(events1_2)

vec2_1 = events_to_vector(events2_1)
vec2_2 = events_to_vector(events2_2)

print(np.linalg.norm(vec1_1-vec1_2))
print(np.linalg.norm(vec2_1-vec2_2))

print(np.linalg.norm(vec1_1-vec2_1))
print(np.linalg.norm(vec1_1-vec2_2))

print(np.linalg.norm(vec1_2-vec2_1))
print(np.linalg.norm(vec1_2-vec2_2))
