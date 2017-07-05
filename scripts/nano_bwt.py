import sys
from math import ceil

MER_LEN = 6
BASES = ['A', 'C', 'G', 'T']

class PoreModel:

    em_means = [0]*(4**MER_LEN)
    em_stdevs = [0]*(4**MER_LEN)
    es_means = [0]*(4**MER_LEN)
    es_stdevs = [0]*(4**MER_LEN)
    
    @staticmethod
    def load_from_txt(filename):
        ##REPLACE WITH FAST5 PARSING
        file_in = open(filename)
        i = 0
        for line in file_in:
            fields = map(float, line.strip().split())
            PoreModel.em_means[i] = fields[0]
            PoreModel.em_stdevs[i] = fields[1]
            PoreModel.es_means[i] = fields[2]
            PoreModel.es_stdevs[i] = fields[3]
            i += 1

    @staticmethod
    def signal_cmp(i1, i2):

        if i1 >= len(PoreModel.es_means):
            if i2 >= len(PoreModel.es_means):
                return 0
            return 1
        elif i2 >= len(PoreModel.es_means):
            return -1

        if PoreModel.em_means[i1] < PoreModel.em_means[i2]:
            return -1
        elif PoreModel.em_means[i1] > PoreModel.em_means[i2]:
            return 1

        if PoreModel.es_means[i1] < PoreModel.es_means[i2]:
            return -1
        elif PoreModel.es_means[i1] > PoreModel.es_means[i2]:
            return 1

        if PoreModel.em_stdevs[i1] < PoreModel.em_stdevs[i2]:
            return -1
        elif PoreModel.em_stdevs[i1] > PoreModel.em_stdevs[i2]:
            return 1

        if PoreModel.es_stdevs[i1] < PoreModel.es_stdevs[i2]:
            return -1
        elif PoreModel.es_stdevs[i1] > PoreModel.es_stdevs[i2]:
            return 1

        return 0


class NanoBWT:
    
    def __init__(self, seq):
        self.signals = [4**MER_LEN]*(len(seq)-MER_LEN+2)
        
        #Note last value remains 4096 - equivalent to '$'
        for i in range(0, len(self.signals)-1):
            self.signals[i] = mer_to_i(seq[i:i+MER_LEN])

        rotations = sorted(range(0, len(self.signals)), cmp=self.rot_cmp)
        self.bwt = [self.signals[i-1] for i in rotations]
        
        for i in self.bwt:
            print i

    def rot_cmp(self, r1, r2):
        for i in range(0, len(self.signals)):
            c = PoreModel.signal_cmp(self.signals[r1+i], self.signals[r2+i])
            if c != 0:
                return c


def i_to_mer(i):
    ret = [0]*MER_LEN
    for n in range(0, MER_LEN):
        ret[MER_LEN-n-1] = BASES[(i >> n*2) & 3]
    return "".join(ret)

def mer_to_i(mer):
    i = BASES.index(mer[0])
    for n in range(1, MER_LEN):
        i = (i << 2) | BASES.index(mer[n])
    return i

def read_fasta(filename):
    fasta_in = open(filename)

    headers = list()
    seqs = list()

    header = None
    seq = list()

    for line in fasta_in:

        if line[0] == ">":

            if header != None:
                headers.append(header)
                seqs.append("".join(seq))

            seq = list()
            header = line.strip()[1:]

        else:
            seq.append(line.strip())

    headers.append(header)
    seqs.append("".join(seq))

    return dict(zip(headers, seqs))



if __name__ == "__main__":
    PoreModel.load_from_txt(sys.argv[1])
    fasta = read_fasta(sys.argv[2])

    for seq in fasta:
        seq =fasta[seq]
        signals = [4**MER_LEN]*(len(seq)-MER_LEN+2)
        
        #Note last value remains 4096 - equivalent to '$'
        for i in range(0, len(signals)-1):
            signals[i] = mer_to_i(seq[i:i+MER_LEN])

        nbt = NanoBWT(seq)
