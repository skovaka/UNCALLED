import sys
import os


SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

MODEL_FILES = [SCRIPT_DIR+'/../kmer_models/r9.2_180mv_250bps_6mer/template_median68pA.model']

class KmerModel:

    BASES = ['A', 'C', 'G', 'T']

    BASE_IDS = {'A': 0, 'a': 0,
                'C': 1, 'c': 1,
                'G': 2, 'g': 2,
                'T': 3, 't': 3}

    BASE_COMPL = {'A': 'T', 'a': 't',
                  'C': 'G', 'c': 'g',
                  'G': 'C', 'g': 'c',
                  'T': 'A', 't': 'a'}

    def __init__(self, model_fname):
        sys.stderr.write("Reading %s\n" % model_fname)

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

    

    def i_to_kmer(self, i):
        ret = [0] * self.k
        for n in range(0, self.k):
            ret[self.k - n - 1] = self.BASES[(i >> n * 2) & 3]
        return "".join(ret)

    def kmer_to_i(self, kmer):
        i = self.BASE_IDS[kmer[0]]
        for n in range(1, self.k):
            i = (i << 2) | self.BASE_IDS[kmer[n]]
        return i

    def reverse_comp(self, seq):
        rev = ['']*len(seq)
        for i in range(0, len(seq)):
            rev[len(seq)-i-1] = self.BASE_COMPL.get(seq[i], 'N')
        return ''.join(rev)



if __name__ == "__main__":
    model = KmerModel(MODEL_FILES[0])



