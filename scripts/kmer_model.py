import sys
import os
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

MODEL_FILES = [SCRIPT_DIR+'/../kmer_models/r9.2_180mv_250bps_6mer/template_median68pA.model']

class KmerModel:

    BASES = ['A', 'C', 'G', 'T']

    BASE_IDS = {'A': 0, 'a': 0,
                'C': 1, 'c': 1,
                'G': 2, 'g': 2,
                'T': 3, 't': 3,
                'N': -1, 'n': -1}

    BASE_COMPL = {'A': 'T', 'a': 't',
                  'C': 'G', 'c': 'g',
                  'G': 'C', 'g': 'c',
                  'T': 'A', 't': 'a'}

    def __init__(self, model_fname=MODEL_FILES[0]):

        model_file = open(model_fname)
        header = model_file.readline().strip().split()
        
        self.has_lambda = len(header) == 7
        self.ig_lambda = -1

        self.lv_means = list()
        self.lv_stdvs = list()
        self.sd_means = list()
        self.sd_stdvs = list()
        
        self.k = -1
        self.kmer_count = -1

        self.rev_ids = list()

        for line in model_file:
            tabs = line.strip().split()

            if self.k < 0:
                self.k = len(tabs[0])
                self.kmer_count = len(self.BASES)**self.k

            k_id = self.kmer_to_i(tabs[0])

            self.rev_ids.append(self.kmer_to_i(self.reverse_comp(tabs[0])))
            
            self.lv_means.append(float(tabs[1])) 
            self.lv_stdvs.append(float(tabs[2]))
            self.sd_means.append(float(tabs[3]))
            self.sd_stdvs.append(float(tabs[4]))

            if (self.has_lambda):
                if self.ig_lambda > 0:
                    if float(tabs[5]) != self.ig_lambda:
                        self.ig_lambda = -1
                        self.has_lambda = False
                else:
                    self.ig_lambda = float(tabs[5])

        self.model_mean = np.mean(self.lv_means)
        self.model_stdv = np.std(self.lv_means)

    def get_norm_params(self, events):
        evt_mean = np.mean([mean for mean, stdv, length in events])
        evt_stdv = np.std([mean for mean, stdv, length in events])

        scale = evt_stdv / self.model_stdv
        shift = evt_mean - (scale * self.model_mean)

        return (scale, shift)

    def get_norm_vals(self, i, norm_params):
        scale, shift = norm_params
        
        norm_mean = self.lv_means[i] * scale + shift;

        ig_lambda = self.ig_lambda
        if not self.has_lambda:
            ig_lambda = pow(self.sd_means[k_id], 3) / pow(self.sd_stdvs[k_id], 2);
        
        return (norm_mean, self.lv_stdvs[i], self.sd_means[i], ig_lambda)

    def event_match_probs(self, event, k_id, norm_params):
        evt_mean, evt_stdv = event

        lv_mean, lv_stdv, sd_mean, ig_lambda = self.get_norm_vals(k_id, norm_params)

        ng = np.exp(-pow(evt_mean - lv_mean, 2) / (2 * pow(lv_stdv, 2))) / np.sqrt(2 * np.pi * pow(lv_stdv, 2))
        ig = np.sqrt(ig_lambda / (2 * np.pi * pow(evt_stdv, 3))) * np.exp(-ig_lambda * pow(evt_stdv - sd_mean, 2) / (2 * evt_stdv * pow(sd_mean, 2)));

        return ng, ig;
    

    def i_to_kmer(self, i):
        ret = [0] * self.k
        for n in range(0, self.k):
            ret[self.k - n - 1] = self.BASES[(i >> n * 2) & 3]
        return "".join(ret)

    def kmer_to_i(self, kmer):
        i = self.BASE_IDS[kmer[0]]

        if i < 0:
            return None

        for n in range(1, self.k):

            if self.BASE_IDS[kmer[n]] < 0:
                return None

            i = (i << 2) | self.BASE_IDS[kmer[n]]

        return i

    def reverse_comp(self, seq):
        rev = ['']*len(seq)
        for i in range(0, len(seq)):
            rev[len(seq)-i-1] = self.BASE_COMPL.get(seq[i], 'N')
        return ''.join(rev)



if __name__ == "__main__":
    model = KmerModel()



