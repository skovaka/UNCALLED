from read_sim import *
from kmer_model import *
import sys
import numpy as np

BASES = ['A', 'C', 'G', 'T']

def edit_dist(s1, s2):
    return sum([c1 != c2 for c1, c2 in zip(s1, s2)])

def get_mutants(model, seq, dist):
    pos = list(range(0, dist))

    out = set()

    while pos[-1] < len(seq):
        vals = [0]*dist
        
        done = False
        while not done:
            
            neighbor = ""
            prev = 0
            i = 0
            for p in pos:
                neighbor += seq[prev:p] + model.BASES[vals[i]]
                prev = p+1
                i += 1
            neighbor += seq[prev:]
            
            if not neighbor in out:
                out.add(neighbor)
                yield neighbor
            
            i = 0
            vals[i] += 1
            while vals[i] == len(model.BASES):
                vals[i] = 0
                i += 1
                
                if i < dist:
                    vals[i] += 1
                else:
                    done = True
                    break
        
        i = 0
        pos[i] += 1
        while i < dist-1 and pos[i] == pos[i+1]:
            pos[i] = i
            i += 1
            pos[i] += 1

def mutate_seq(seq, d):
    seq = list(seq)
    locs = np.random.choice(range(0, len(seq)), d, False)

    for i in locs:
        opts = list(BASES)
        opts.remove(seq[i])
        seq[i] = np.random.choice(opts)
        

    return "".join(seq)


def match_prob(model, seq1, seq2):

    kmers = model.seq_to_kmers(seq1)

    events = simulate_read(seq2, model)

    seed_prob = 0
    min_prob = 0
    for kmer, event in zip(kmers, events):
        ng, ig = model.event_match_probs(event, kmer, (1, 0))
        prob = np.log(ng * ig)
        seed_prob += prob
        if prob < min_prob:
            min_prob = prob

    return (seed_prob, min_prob)

    
if __name__ == "__main__":
    model = KmerModel(sys.argv[1])
    fasta = open(sys.argv[2])
    dist = int(sys.argv[3])
    N = int(sys.argv[4])

    event_thresh = -5.29
    seed_thresh  = -3.35

    seq = filter(lambda s: s[0] != ">", fasta.readlines())
    seq = "".join(map(lambda s: s.strip(), seq))

    event_count = len(seq) - 6 + 1

    event_pass = seed_pass = both_pass = 0

    for i in range(0, N):
        seq2 = mutate_seq(seq, dist)
        seed_prob, min_prob = match_prob(model, seq, seq2)

        good_seed = seed_prob >= seed_thresh * event_count
        good_event = min_prob >= event_thresh
        seed_pass += good_seed
        event_pass += good_event
        both_pass += good_seed and good_event

    print ("%.2f%% pass event thresh (%.2f)" % (100.0 * event_pass / N, event_thresh))
    print ("%.2f%% pass seed thresh (%.2f)" % (100.0 * seed_pass / N, seed_thresh))
    print ("%.2f%% pass both" % (100.0 * both_pass / N))



