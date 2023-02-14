import sys
from time import time
import numpy as np
import pandas as pd
import os
import re
from _uncalled import ReadBufferBC

from ont_fast5_api.fast5_interface import get_fast5_file
import pyslow5

def is_read_id(read_id):
    return re.match("[a-z0-9]+(-[a-z0-9])+", read_id) is not None

#class ReaderBase:
#    def __init__(self, filename):

SUFFIXES = {"fast5", "slow5", "blow5"}

class ReadIndex:
    def __init__(self, file_paths=None, read_filter=None, index_filename=None, recursive=False):
        #self.suffix_re = re.compile(suffix_regex)
        self.index_filename = index_filename
        self.file_info = dict()
        self.read_files = None

        if file_paths is not None:
            self.load_paths(file_paths, recursive)

        self._load_filter(read_filter)
        self.load_index_file(index_filename)

        self.infile = None
        self.infile_name = None
        
    @property
    def indexed(self):
        return len(self.file_info) > 0
    
    def subset(self, read_ids):
        ret = ReadIndex(read_filter=read_ids)
        ret.load_index_df(self.read_files.reset_index())
        ret.file_paths = self.file_info
        return ret

    def load_index_file(self, fname=None):
        if fname is None or len(fname) == 0:
            if self.index_filename is None or len(self.index_filename) == 0:
                return
            fname = self.index_filename

        read_files = None
        with open(fname) as infile:
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

        self.load_index_df(read_files)

    def load_index_df(self, df):
        if len(df.columns) > 2:
            df = df[["read_id", "filename"]]

        if self.read_filter is not None:
            df = df[df["read_id"].isin(self.read_filter)]

        filenames = df.set_index("read_id")["filename"].str.rsplit("/", n=2, expand=True).iloc[:,-1]

        if self.read_files is None:
            self.read_files = filenames.sort_index()
        else:
            filenames = filenames[~filenames.index.isin(self.read_files.index)]
            self.read_files = pd.concat([self.read_files, filenames]).sort_index()

    def load_paths(self, paths, recursive):
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
                        self._add_read_file(os.path.join(root, fname))

            #Non-recursive directory search 
            elif isdir and not recursive:
                for fname in os.listdir(path):
                    self._add_read_file(os.path.join(path, fname))

            #Read fast5 name directly
            elif self._is_read_file(path):
                self._add_read_file(path)

            #Read fast5 filenames from text file
            else:
                with open(path) as infile:
                    for line in infile:
                        self._add_read_file(line.strip())

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
        return not fname.startswith("#") and fname.split(".")[-1] in SUFFIXES

    def get_read_file(self, read_id):
        filename = self.read_files.loc[read_id]
        return filename

    def _add_read_file(self, path):
        if not self._is_read_file(path):
            return False

        fname = os.path.basename(path)
        Reader = FILE_READERS[fname.split(".")[-1]]
        self.file_info[fname] = (Reader, path)

        if Reader == Slow5Reader:
            self._open(fname)
            self.load_index_df(pd.DataFrame({
                "read_id"  : self.infile.read_ids,
                "filename" : fname
            }))

        return True

    def _open(self, filename):
        if filename == self.infile_name:
            return self.infile
        if self.infile is not None:
            self.infile.close()
        self.infile_name = filename

        Reader,path = self.file_info[filename]
        self.infile = Reader(path)
        return self.infile

    def __getitem__(self, read_id):
        self._open(self.get_read_file(read_id))
        return self.infile[read_id]
        #return self._get(read_id, self.get_read_file(read_id))

    def close(self):
        if self.infile is not None:
            self.infile.close()

class Fast5Reader:
    #def __init__(self, file_paths=None, read_filter=None, index_filename=None, recursive=False):
    def __init__(self, filename):
        self.infile = get_fast5_file(filename, mode="r")

    def __getitem__(self, read_id):
        f5 = self.infile.get_read(read_id)
        return self._dict_to_read(f5)
    
    def _dict_to_read(self, f5):
        channel = f5.get_channel_info()
        attrs = f5.handle[f5.raw_dataset_group_name].attrs#["start_time"])
        read = ReadBufferBC(f5.read_id, channel["channel_number"], attrs["read_number"], attrs["start_time"], f5.get_raw_data(scale=True))
        #read.filename = filename

        bc_group = f5.get_latest_analysis("Basecall_1D")
        seg_group = f5.get_latest_analysis("Segmentation")
        if bc_group is not None and seg_group is not None:
            moves = np.array(f5.get_analysis_dataset(bc_group, "BaseCalled_template")["Move"])
            move_attr = f5.get_analysis_dataset(bc_group, "Summary")["basecall_1d_template"].attrs
            stride = move_attr["block_stride"]
            seg_attr = f5.get_analysis_dataset(seg_group, "Summary")["segmentation"].attrs
            template_start = seg_attr["first_sample_template"]

            read.set_moves(moves, template_start, stride)

        return read

    def __iter__(self):
        for r in self.infile:
            yield self._dict_to_read(r)

    def close(self):
        self.infile.close()

class Slow5Reader:
    def __init__(self, filename):
        self.infile = pyslow5.Open(filename, mode="r")
        self.read_ids,_ = self.infile.get_read_ids()

    def _dict_to_read(self, d):
        return ReadBufferBC(d["read_id"], 0, 0, 0, d["signal"])

    def __getitem__(self, read_id):
        r = self.infile.get_read(read_id, pA=True)
        return self._dict_to_read(r)

    def __iter__(self):
        for r in self.infile.seq_read():
            yield self._dict_to_read(r)

    def close(self):
        self.infile.close()

FILE_READERS = {
    "fast5" : Fast5Reader,
    "slow5" : Slow5Reader,
    "blow5" : Slow5Reader
}
