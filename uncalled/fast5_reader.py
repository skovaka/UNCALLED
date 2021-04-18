import sys
import time
import glob
import numpy as np
import pandas as pd
import os
from _uncalled import _Fast5Dict

def fast5_glob(root, recursive):
    suffix = "*.fast5"
    if recursive:
        suffix = os.path.join("**", suffix)
    return glob.iglob(os.path.join(root, suffix), recursive=recursive)

class Fast5Dict(_Fast5Dict):
    def __init__(self, seqsum, root_dir=".", recursive=False):

        fast5_paths = {
            os.path.basename(path) : path
            for path in fast5_glob(root_dir, recursive)
        }

        seqsum = pd.read_csv(seqsum, sep="\t", usecols=["filename", "read_id"])
        #seqsum.replace({"filename" : fast5_paths}, inplace=True)

        groups = seqsum.groupby("filename").groups
        groups = {
            fast5_paths[fast5] : list(seqsum.loc[rows, "read_id"])
            for fast5,rows in groups.items()
        }

        print(groups)

        _Fast5Dict.__init__(self, groups)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("seqsum", type=str)
    parser.add_argument("root", type=str)
    parser.add_argument("-r", "--recursive", action="store_true")
    args = parser.parse_args(sys.argv[1:])

    fast5s = Fast5Dict(args.seqsum, args.root, args.recursive)
