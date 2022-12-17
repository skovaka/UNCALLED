import sys
from time import time
import numpy as np
import pandas as pd
import os

def ReadIndex:
    def __init__(self, index_filename=None, file_paths=None, read_filter=None, file_suffix=".fast5"):
        self.file_suffix = file_suffix
        self.index_filename = index_filename
        self.file_paths = None

        self._load_filter(read_filter)

    def _load_filter(self, reads):
        if reads is None:
            self.read_filter = None

        elif isinstance(reads, str):
            if os.path.exists(reads):
                with open(reads) as reads_in:
                    self.read_filter = {line.split()[0] for line in reads_in}
            else:
                self.read_filter = set(reads.split(","))
        else:
            self.read_filter = set(reads)

    def _is_read_file(self, fname):
        return fname.endswith(self.file_suffix) and not fname.startswith("#")

    def load_index(self, filename=None):
        if filename is None:
            if self.index_filename is None:
                raise ValueError("No read index file specified")
            filename = self.index_filename

        read_files = None
        with open(filename) as infile:
            head = infile.readline().split()
            if "filename" in head and "read_id" in head:
                names = head
                header = 0

            elif len(head) == 2:
                if self._is_read_file(head[0]) and is_read_id(head[1]):
                    names = ["filename", "read_id"]
                elif is_read_id(head[0]) and self._is_read_file(head[1]):
                    names = ["read_id", "filename"]
                header = None

            else:
                raise RuntimeError("Could not find 'filename' and 'read_id' columns")

            infile.seek(0)
            read_files = pd.read_csv(infile, sep="\t", names=names, header=header, usecols=["read_id", "filename"])

        if read_files is None:
            raise RuntimeError("Unable to read index \"%s\"" % filename)
        
        self.read_files = read_files.set_index("read_id")["filename"]

        if self.read_filter is None:
            reads = self.read_files.index.intersection(self.read_filter)
            if len(reads) < len(self.read_filter):
                n = len(self.read_filter) - len(reads)
                sys.stderr.write("Warning: {n} read_ids misisng from read index\n")
            self.read_files = self.read_files.loc[reads]

    def load_paths(self, paths, recursive):
        self.file_paths = dict()

        def add_file(path):
            if not self._is_read_file(path):
                return False

            if not os.path.isfile(path):
                sys.stderr.write("Warning: fast5 file \"%s\" does not exist\n" % fname)
                return False

            #TODO is abspath nesissary?
            path,fname = os.path.split(os.path.abspath(path))

            file_paths[fname] = path

            return True

        if isinstance(paths, str):
            paths = [paths]

        for path in paths:
            path = path.strip()

            if not os.path.exists(path):
                raise ValueError("Error: \"%s\" does not exist\n" % path)

            isdir = os.path.isdir(path)

            #Recursive directory search 
            if isdir and recursive:
                for root, dirs, files in os.walk(path):
                    for fname in files:
                        add_file(os.path.join(root, fname))

            #Non-recursive directory search 
            elif isdir and not recursive:
                for fname in os.listdir(path):
                    add_file(os.path.join(path, fname))

            #Read fast5 name directly
            elif path.endswith(".fast5"):
                add_file(path)

            #Read fast5 filenames from text file
            else:
                with open(path) as infile:
                    for line in infile:
                        add_file(line.strip())

        return file_paths


