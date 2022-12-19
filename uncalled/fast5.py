import sys
from time import time
import glob
import numpy as np
import pandas as pd
import os
from _uncalled import _Fast5Dict, _Fast5Iter, _Fast5Reader
import uncalled as unc
import re
from collections import defaultdict

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

    if isinstance(fast5s, str):
        fast5s = [fast5s]

    for path in fast5s:
        path = path.strip()

        if not os.path.exists(path):
            raise ValueError("Error: \"%s\" does not exist\n" % path)

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

    if isinstance(reads, str):
        if os.path.exists(reads):
            with open(reads) as reads_in:
                return [line.split()[0] for line in reads_in]
        else:
            return reads.split(",")

    return list(reads)

class Fast5Reader:

    #def __init__(self, fast5s=None, index=None, reads=None, recursive=None, conf=None):
    def __init__(self, read_index, conf=None):#fast5s=None, index=None, reads=None, recursive=None, conf=None):

        self.conf = unc.Config(conf)# if conf is None else conf

        #TODO probably not ll this
        conf_prms = self.conf.fast5_reader
        self.prms = _Fast5Iter.Params()
        self.prms.max_reads = conf_prms.max_reads
        self.prms.max_buffer = conf_prms.max_buffer
        self.prms.load_bc = conf_prms.load_bc
        self.prms.bc_group = conf_prms.bc_group

        self.indexed = read_index.indexed

        #TODO rewrite Fast5Dict to take read->filename, filename->path
        #maybe read_array, filename_array, path_array, all sorted by filename?
        if self.indexed:
            #self._dict = _Fast5Dict(read_index.get_fast5_dict(), self.prms)
            idx = read_index.read_files
            self._dict = _Fast5Dict(read_index.file_paths, list(idx["read_id"]), list(idx["filename"]), self.prms)
        else:
            self._dict = None

        #if reads is None: 
        #    reads = conf_prms.read_filter
        #self.prms.read_filter = parse_read_ids(reads)

        #if recursive is None: 
        #    recursive = conf_prms.recursive
        #self.prms.recursive = recursive

        #if fast5s is None: 
        #    fast5s = conf_prms.fast5_files
        #self.prms.fast5_files = parse_fast5_paths(fast5s, recursive)

        #if index is None:
        #    index = conf_prms.fast5_index
        ##self.prms.fast5_index = index

        #if isinstance(index, str): 
        #    if len(index) > 0:
        #        self._load_index_file(index)
        #elif isinstance(index, pd.DataFrame):
        #    self._load_index_df(index)
        #else:
        #    raise ValueError("Unknown fast5 index format")

    def _load_index_df(self, df):
        fast5_paths = {os.path.basename(path) : path for path in self.prms.fast5_files}

        idx = df.set_index("filename").sort_index()
        fast5_reads = dict()
        for fast5 in idx.index.unique():
            path = fast5_paths[os.path.basename(fast5)]
            reads = idx.loc[fast5, "read_id"]
            if isinstance(reads, str):
                reads = [reads]
            else:
                reads = list(reads)
            fast5_reads[path] = reads

        self._dict = _Fast5Dict(fast5_reads, self.prms)

        self.indexed = True
            
    def _load_index_file(self, filename):
        index = None
        names = None



        with open(filename) as infile:
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
                index = pd.read_csv(infile, sep="\t", names=names, header=header, usecols=["filename", "read_id"])

        if index is None:
            raise RuntimeError("Unable to read index \"%s\"" % self.prms.fast5_index)

        if len(self.prms.read_filter) > 0:
            index = index[index["read_id"].isin(self.prms.read_filter)]

        if len(index) == 0:
            sys.stderr.write("Warning: fast5 index is empty after applying read filter. Treating Fast5Reader as un-indexed (only iteration supported)\n")
            self.indexed = False
            return

        self.prms.read_filter = list(index["read_id"])

        self._load_index_df(index)

        self.indexed = True

    def get_read_filter(self):
        if len(self.prms.read_filter) > 0:
            return set(self.prms.read_filter)
        return None
    
    def get_read_file(self, read_id):
        if not self.indexed:
            raise RuntimeError("Fast5 index is required to query fast5 filenames")
        return self._dict.get_read_file(read_id)

    def __contains__(self, read_id):
        if not self.indexed:
            raise RuntimeError("Fast5 index is required for dict-like fast5 access (iteration still supported)")
        #TODO faster solution?
        return not self._dict[read_id].empty()

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
        try:
            while not self._iter.empty():
                return self._iter.next_read()
        except:
            sys.stderr.write("Failed to load fast5 file\n")
        raise StopIteration

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("seqsum", type=str)
    parser.add_argument("root", type=str)
    parser.add_argument("-r", "--recursive", action="store_true")
    args = parser.parse_args(sys.argv[1:])

    fast5s = Fast5Dict(args.seqsum, args.root, args.recursive)

