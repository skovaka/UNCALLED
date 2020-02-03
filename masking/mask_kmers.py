#!/usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import argparse

def mask_seq(seq, kmer):
    ranges = list()
    
    k = len(kmer)
    n = 0

    useq = seq.upper()

    i = useq.find(args.kmer)
    while i >= 0:
        j = i + k
        if len(ranges) == 0 or i > ranges[-1][1]:
            ranges.append( (i,j) )
        else:
            #sys.stderr.write(str(ranges[-1]) + " ")
            ranges[-1] = (ranges[-1][0],j)
            #sys.stderr.write(str(ranges[-1]) + "\n")

        n += 1
        i = useq.find(args.kmer, i+1)

    masked = list()
    i = 0
    for j,k in ranges:
        masked.append(seq[i:j] + "N"*(k-j))
        i = k

    if len(ranges) > 0:
        masked.append(seq[ranges[-1][1]:])
        return "".join(masked), n


    return seq, 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("infile", type=str, default=None)
    parser.add_argument("-k", "--kmer", required=False, type=str, default=None)
    parser.add_argument("-o", "--outfile", required=False, type=str, default=None)
    args = parser.parse_args(sys.argv[1:])

    run = 0
    prev = None

    seq = list()

    masked = 0
    headers = list()
    masked_seqs = list()

    k=len(args.kmer)

    Ns = "N"*k


    for line in open(args.infile):
        if line[0] == ">":
            if len(seq) > 0:
                m,n = mask_seq("".join(seq), args.kmer)
                masked_seqs.append(m)
                masked += n
                seq = list()
            headers.append(line.strip())
        else:
            seq.append(line.strip())

    m,n = mask_seq("".join(seq), args.kmer)
    masked_seqs.append(m)
    masked += n

    sys.stderr.write("masked %d occurences of %s\n" % (masked, args.kmer))

    for h,s in zip(headers, masked_seqs):
        print(h)
        print(s)
