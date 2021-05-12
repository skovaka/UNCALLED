import sys
import time
import glob
import numpy as np
import pandas as pd
import os
from _uncalled import _Fast5Dict, _Fast5Iter
import uncalled as unc

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

        if fast5s is not None: self.parse_fast5_paths(fast5s)
        if index is not None: self.prms.fast5_index = index
        if reads is not None: self.prms.read_filter = self.parse_read_str(reads)
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
                if is_fast5(head[0]) and is_read(head[1]):
                    names = ["filename", "read_id"]
                elif is_read(head[0]) and is_fast5(head[1]):
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
            elif is_fast5(path):
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
