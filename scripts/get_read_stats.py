import h5py
import sys
import numpy as np
import argparse
import os
from kmer_model import *

MATCH = 0
STAY = 1
SKIP = 2

LEN_FACTOR = 4000

if __name__ == "__main__":
    i = 0
    model = KmerModel(sys.argv[1])
    fast5_dir = sys.argv[2]
    seed_lens = map(int, sys.argv[3:])

    for filename in os.listdir(fast5_dir):
        path = os.path.join(os.path.abspath(fast5_dir), filename)

        hdf = h5py.File(path,'r')

        readnum = filename.split('_')[-2]
        readname='R'+readnum[1:4]+'_'+readnum[4:]

        #events = hdf['Analyses']['EventDetection_000']['Reads'][readname]['Events']
        reads = hdf['Analyses']['EventDetection_000']['Reads']
        events = hdf['Analyses']['Basecall_1D_000']['BaseCalled_template']['Events']

        bc_fields = events.dtype.fields
        
        event_info = list()

        match_lens = list()
        stay_lens = list()
        skip_lens = list()

        p_kmer = None
        p_length = None

        for e in events:
            kmer = e[4].decode('UTF-8')
            length = e[3]*LEN_FACTOR

            if not p_kmer or p_kmer[1:] == kmer[:-1]:
                event_info.append(MATCH)

                if p_length != None:
                    #match_lens.append( max((length, p_length)) )    
                    match_lens.append((p_length + length) / 2)    

            elif p_kmer == kmer:
                event_info.append(STAY)

                if p_length != None:
                    #stay_lens.append( max((length, p_length)) )    
                    stay_lens.append( (p_length + length) / 2 )    

            elif p_kmer[2:] == kmer[:-2]:

                s_kmer = p_kmer[1:] + kmer[-2]

                event_info.append(SKIP)

                if p_length != None:
                    #skip_lens.append( max((length, p_length)) )    
                    skip_lens.append((p_length + length) / 2)    

            p_kmer = kmer
            p_length = length

        skip_counts = [0]*len(seed_lens)

        for i in range(0, len(event_info)):

            for j in range(0, len(seed_lens)):
                if (i+seed_lens[j] <= len(event_info) and 
                    SKIP in event_info[i:i+seed_lens[j]]):

                    skip_counts[j] += 1


        avg_skip = sum(skip_lens) / len(skip_lens)
        std_skip = np.std(skip_lens)
        avg_stay = sum(stay_lens) / len(stay_lens)
        std_stay = np.std(stay_lens)
        avg_match = sum(match_lens) / len(match_lens)
        std_match = np.std(match_lens)

        sk_ratio = avg_skip - avg_match

        out = "%d\t%.2f (%.2f)\t%.2f (%.2f)\t%.2f (%.2f)\t%.2f" % (len(event_info), avg_match, std_match, avg_stay, std_stay, avg_skip, std_skip, avg_skip-avg_match) 
        for sc in skip_counts:
            out += "\t%.1f%%" % (100.0 * sc / len(event_info))

        print out

        sys.stdout.flush()

    






