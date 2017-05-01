#!/usr/bin/python2.7

import sys
import random

if __name__ == "__main__":
    length = int(sys.argv[1])
    print("> random_seq_of_len_%s" % length)
    print("".join([random.choice(["A", "C", "G", "T"]) for _ in range(length)]))
