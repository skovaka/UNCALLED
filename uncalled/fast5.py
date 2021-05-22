import sys
import time
import glob
import numpy as np
import pandas as pd
import os
from _uncalled import _Fast5Dict, _Fast5Iter
import uncalled as unc
import re

def fast5_path(fname):
    if fname.startswith("#") or not fname.endswith("fast5"):
        return None

    path = os.path.abspath(fname)
    if not os.path.isfile(path):
        sys.stderr.write("Warning: skipping \"%s\" (not a fast5 file)\n" % fname)
        return None

    return path

def iter_fast5_fnames(fast5s, recursive):
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
                    fast5_name = fast5_path(os.path.join(root, fname))
                    if fast5_name is not None: yield fast5_name

        #Non-recursive directory search 
        elif isdir and not recursive:
            for fname in os.listdir(path):
                fast5_name = fast5_path(os.path.join(path, fname))
                if fast5_name is not None: yield fast5_name

        #Read fast5 name directly
        elif path.endswith(".fast5"):
            fast5_name = fast5_path(path)
            if fast5_name is not None: yield fast5_name

        #Read fast5 filenames from text file
        else:
            with open(path) as infile:
                for line in infile:
                    fast5_name = fast5_path(line.strip())
                    if fast5_name is not None: yield fast5_name

def iter_reads(reads):
    if type(reads) in (list, tuple):
        return reads

    if reads is None:
        return []

    if os.path.exists(reads):
        #with open(reads) as reads_in:
        return (l.strip() for l in open(reads))
        sys.stderr.write("Error: failed to open read filter\n")
        sys.exit(1)

    return reads.split(",")

def get_fast5_reader(fast5s, recursive, read_filter):
    reader = unc.Fast5Reader()

    for fast5 in iter_fast5_fnames(fast5s, recursive):
        reader.add_fast5(fast5)

    for read in iter_reads(read_filter):
        reader.add_read(read)

    return reader


def fast5_glob(root, recursive):
    suffix = "*.fast5"
    if recursive:
        suffix = os.path.join("**", suffix)
    return glob.iglob(os.path.join(root, suffix), recursive=recursive)

#TODO make fast5 parameters
class Fast5Reader:

    def __init__(self, fast5s=None, index=None, reads=None, recursive=None, conf=None):
        self.conf = unc.Conf() if conf is None else conf
        self.prms = self.conf.fast5_reader

        if fast5s is None: 
            fast5s = self.prms.fast5_files
        self.parse_fast5_paths(fast5s)

        if reads is not None: 
            self.prms.read_filter = self.parse_read_str(reads)

        if index     is not None: self.prms.fast5_index = index
        if recursive is not None: self.prms.recursive = recursive

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
                if self.is_fast5(head[0]) and self.is_read_id(head[1]):
                    names = ["filename", "read_id"]
                elif self.is_read_id(head[0]) and self.is_fast5(head[1]):
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

    @staticmethod
    def parse_read_str(reads):
        """Takes either a comma-separated string of read IDs or a path to a file containing read IDs and returns a list of read IDs"""
        if reads is None:
            return []

        if os.path.exists(reads):
            with open(reads) as reads_in:
                return [line.split()[0] for line in reads_in]

        return reads.split(",")
    
    def parse_fast5_paths(self, fast5s=None):
        if fast5s is None:
            fast5s = self.prms.fast5_files
        elif type(fast5s) == str:
            fast5s = fast5s.split(",")

        fast5_paths = list()

        def add_fast5(fname):
            if fname.startswith("#") or not self.is_fast5(fname):
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
            if isdir:
                if self.prms.recursive:
                    for root, dirs, files in os.walk(path):
                        for fname in files:
                            fast5_name = os.path.join(root, fname)
                            add_fast5(fast5_name)
                else:
                    for fname in os.listdir(path):
                        fast5_name = os.path.join(path, fname)
                        add_fast5(fast5_name)

            #Read fast5 name directly
            elif self.is_fast5(path):
                add_fast5(path)

            #Read fast5 filenames from text file
            else:
                with open(path) as infile:
                    for line in infile:
                        add_fast5(line.strip())

        self.prms.fast5_files = fast5_paths

    @staticmethod
    def is_fast5(fname):
        return fname.endswith(".fast5")

    @staticmethod
    def is_read_id(read_id):
        return re.match("[a-z0-9]+(-[a-z0-9])+", read_id) is not None

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


class Fast5Dict(_Fast5Dict):
    def __init__(self, root_dir=".", seqsum=None, recursive=False, conf=unc.Conf()):
        fast5_paths = {os.path.basename(path) : path for path in iter_fast5_fnames(root_dir, recursive)}

        #TODO clean this up with paramters
        #if type(root_dir) == list:
        #    for fname in root_dir:
        #        with open(fname) as infile:
        #            fast5_paths = {os.path.basename(path.strip()) : path.strip() for path in infile}
        #else:
        #    fast5_paths = {
        #        os.path.basename(path) : path
        #        for path in fast5_glob(root_dir, recursive)
        #    }

        seqsum = pd.read_csv(seqsum, sep="\t", usecols=["filename", "read_id"])

        groups = seqsum.groupby("filename").groups
        groups = {
            fast5_paths[fast5] : list(seqsum.loc[rows, "read_id"])
            for fast5,rows in groups.items()
        }

        _Fast5Dict.__init__(self, groups, conf.fast5_reader)

class Fast5Iter(_Fast5Iter):
    def __init__(self, params):
        _Fast5Iter.__init__(self, params)

    def __iter__(self):
        return self
    
    def __next__(self):
        if self.empty():
            raise StopIteration
        return self.next_read()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("seqsum", type=str)
    parser.add_argument("root", type=str)
    parser.add_argument("-r", "--recursive", action="store_true")
    args = parser.parse_args(sys.argv[1:])

    fast5s = Fast5Dict(args.seqsum, args.root, args.recursive)

Opt = unc.ArgParser.Opt
FAST5_OPTS = (
    Opt("fast5_files", "fast5_reader", nargs="+", type=str),
    Opt(("-r", "--recursive"), "fast5_reader", action="store_true"),
    Opt(("-l", "--read-filter"), "fast5_reader", type=Fast5Reader.parse_read_str),
    Opt(("-x", "--fast5-index"), "fast5_reader"),
    Opt(("-n", "--max-reads"), "fast5_reader")
)
