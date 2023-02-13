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

class ReadIndex:
    def __init__(self, file_paths=None, read_filter=None, index_filename=None, recursive=False, suffix_regex="\.fast5"):
        self.suffix_re = re.compile(suffix_regex)
        self.index_filename = index_filename
        self.file_paths = dict()
        self.read_files = None

        if file_paths is not None:
            self.load_paths(file_paths, recursive)

        self._load_filter(read_filter)
        self.load_index_file(index_filename)

        self.infile = None
        self.infile_name = None
        
    @property
    def indexed(self):
        return len(self.file_paths) > 0
    
    def subset(self, read_ids):
        ret = ReadIndex(read_filter=read_ids)
        ret.suffix_re = self.suffix_re
        ret.load_index_df(self.read_files.reset_index())
        ret.file_paths = self.file_paths
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
            elif self.suffix_re.search(fname):
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
        return self.suffix_re.search(fname) and not fname.startswith("#")

    def get_read_file(self, read_id):
        filename = self.read_files.loc[read_id]
        return self.file_paths[filename]

    def _add_read_file(self, path):
        if not self._is_read_file(path):
            return False

        fname = os.path.basename(path)
        self.file_paths[fname] = path

        return True

    def _open(self, filename):
        if filename == self.infile_name:
            return False
        elif self.infile is not None:
            self.infile.close()
        self.infile_name = filename
        return True

    def __iter__(self):
        prev = prev = None
        for read_id, filename in self.read_files.items():
            if prev != filename:
                path = self.file_paths[filename]
            yield self._get(read_id, path)

    def __getitem__(self, read_id):
        return self._get(read_id, self.get_read_file(read_id))

    def close(self):
        if self.infile is not None:
            self.infile.close()

class Fast5Reader(ReadIndex):
    def __init__(self, file_paths=None, read_filter=None, index_filename=None, recursive=False):
        ReadIndex.__init__(self, file_paths, read_filter, index_filename, recursive, "\.fast5")

    def _open(self, filename):
        if ReadIndex._open(self, filename):
            self.infile = get_fast5_file(self.infile_name, mode="r")

    def _get(self, read_id, filename):
        self._open(filename)
        f5 = self.infile.get_read(read_id)
        channel = f5.get_channel_info()
        attrs = f5.handle[f5.raw_dataset_group_name].attrs#["start_time"])
        read = ReadBufferBC(f5.read_id, channel["channel_number"], attrs["read_number"], attrs["start_time"], f5.get_raw_data(scale=True))
        read.filename = filename

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

class Slow5Reader(ReadIndex):
    def __init__(self, file_paths=None, read_filter=None, recursive=False):
        ReadIndex.__init__(self, file_paths, read_filter, None, recursive, "\.[bs]low5$")

        for fname,path in self.file_paths.items():
            self._open(path)
            read_ids,n = self.infile.get_read_ids()
            self.load_index_df(pd.DataFrame({"read_id" : read_ids, "filename" : fname}))

    def _add_read_file(self, path):
        if ReadIndex._add_read_file(self, path):
            self._open(path)
            read_ids,n = self.infile.get_read_ids()
            self.load_index_df(pd.DataFrame({"read_id" : read_ids, "filename" : path}))
            return True
        return False

    def _open(self, filename):
        if ReadIndex._open(self, filename):
            self.infile = pyslow5.Open(self.infile_name, mode="r")

    def _get(self, read_id, filename):
        self._open(filename)
        
        s5 = self.infile.get_read(read_id, pA=True)

        read = ReadBufferBC(s5["read_id"], 0, 0, 0, s5["signal"])
        
        return read

