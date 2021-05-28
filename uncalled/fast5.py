import sys
import time
import glob
import numpy as np
import pandas as pd
import os
from _uncalled import _Fast5Dict, _Fast5Iter
import uncalled as unc
import re

def is_fast5(fname):
    return fname.endswith(".fast5")

def is_read_id(read_id):
    return re.match("[a-z0-9]+(-[a-z0-9])+", read_id) is not None

def parse_fast5_paths(fast5s, recursive):
    fast5_paths = list()

    def add_fast5(fname):
        if fname.startswith("#") or not is_fast5(fname):
            return False

        path = os.path.abspath(fname)
        if not os.path.isfile(path):
            sys.stderr.write("Warning: fast5 file \"%s\" does not exist\n" % fname)
            return False

        fast5_paths.append(path)
        return True

    for path in fast5s:
        path = path.strip()

        if not os.path.exists(path):
            sys.stderr.write("Error: \"%s\" does not exist\n" % path)
            sys.exit(1)

        isdir = os.path.isdir(path)

        #Recursive directory search 
        if isdir and recursive:
            for root, dirs, files in os.walk(path):
                for fname in files:
                    add_fast5(os.path.join(root, fname))

        #Non-recursive directory search 
        elif isdir and not recursive:
            for fname in os.listdir(path):
                add_fast5(os.path.join(path, fname))

        #Read fast5 name directly
        elif path.endswith(".fast5"):
            add_fast5(path)

        #Read fast5 filenames from text file
        else:
            with open(path) as infile:
                for line in infile:
                    add_fast5(line.strip())

    return fast5_paths

def parse_read_ids(reads):
    if reads is None:
        return []

    if os.path.exists(reads):
        with open(reads) as reads_in:
            return [line.split()[0] for line in reads_in]

    return reads.split(",")

class Fast5Reader:

    def __init__(self, fast5s=None, index=None, reads=None, recursive=None, conf=None):
        self.conf = unc.Config() if conf is None else conf
        self.prms = self.conf.fast5_reader

        if fast5s is None: 
            fast5s = self.prms.fast5_files

        if recursive is not None: self.prms.recursive = recursive

        self.prms.fast5_files = parse_fast5_paths(fast5s, recursive)

        if reads is not None: 
            self.prms.read_filter = parse_read_ids(reads)

        if index     is not None: self.prms.fast5_index = index

        self.indexed = len(self.prms.fast5_index) != 0

        if self.indexed:
            self._load_index()
            
    def _load_index(self):
        index = None
        names = None

        fast5_paths = {os.path.basename(path) : path for path in self.prms.fast5_files}

        with open(self.prms.fast5_index) as infile:
            head = infile.readline().split()
            if "filename" in head and "read_id" in head:
                names = head
                header = 0

            elif len(head) == 2:
                if is_fast5(head[0]) and is_read_id(head[1]):
                    names = ["filename", "read_id"]
                elif is_read_id(head[0]) and is_fast5(head[1]):
                    names = ["read_id", "filename"]
                header = None

            infile.seek(0)

            if names is not None:
                index = pd.read_csv(infile, sep="\t", names=names, header=header)

        if index is None:
            raise RuntimeError("Unable to read index \"%s\"" % self.prms.fast5_index)

        if len(self.prms.read_filter) > 0:
            index = index[index["read_id"].isin(self.prms.read_filter)]

        if len(index) == 0:
            sys.stderr.write("Warning: fast5 index is empty after applying read filter. Treating Fast5Reader as un-indexed (only iteration supported)\n")
            self.indexed = False
            return

        self.prms.read_filter = list(index["read_id"])

        groups = index.groupby("filename").groups
        groups = {
            fast5_paths[os.path.basename(fast5)] : list(index.loc[rows, "read_id"])
            for fast5,rows in groups.items()
        }


        self._dict = _Fast5Dict(groups, self.prms)
    
    def get_read_file(self, read_id):
        if not self.indexed:
            raise RuntimeError("Fast5 index is required to query fast5 filenames")
        return self._dict.get_read_file(read_id)

    def __getitem__(self, read_id):
        if not self.indexed:
            raise RuntimeError("Fast5 index is required for dict-like fast5 access (iteration still supported)")

        read = self._dict[read_id]
        if read.empty():
            raise KeyError(read_id)

        return read
    
    def __iter__(self):
        self._iter = _Fast5Iter(self.prms)
        return self

    def __next__(self):
        if self._iter.empty():
            raise StopIteration
        return self._iter.next_read()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("seqsum", type=str)
    parser.add_argument("root", type=str)
    parser.add_argument("-r", "--recursive", action="store_true")
    args = parser.parse_args(sys.argv[1:])

    fast5s = Fast5Dict(args.seqsum, args.root, args.recursive)

Opt = unc.config.Opt
FAST5_OPTS = (
    Opt(unc.config.FAST5_PARAM, "fast5_reader", nargs="+", type=str),
    Opt(("-r", "--recursive"), "fast5_reader", action="store_true"),
    Opt(("-l", "--read-filter"), "fast5_reader", type=parse_read_ids),
    Opt(("-x", "--fast5-index"), "fast5_reader"),
    Opt(("-n", "--max-reads"), "fast5_reader")
)
