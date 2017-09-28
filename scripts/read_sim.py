import sys
from kmer_model import KmerModel
from numpy.random import normal, wald


def simulate_read(seq, model):

    kmers = model.seq_to_kmers(seq)

    model_params = [model.get_norm_vals(k) for k in kmers]

    events = [(normal(lvm, lvs), wald(sdm, igl)) \
                for lvm, lvs, sdm, igl in model_params]

    return events

if __name__ == "__main__":

    model = KmerModel(sys.argv[1])

    seq = sys.stdin.readline().strip()

    events = simulate_read(seq, model)

    for mean, stdv in events:
        print("%.3f\t%.3f" % (mean, stdv))
