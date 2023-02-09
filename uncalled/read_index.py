import sys
from time import time
import numpy as np
import pandas as pd
import os
import re
from ont_fast5_api.fast5_interface import get_fast5_file
from _uncalled import ReadBufferBC

def is_read_id(read_id):
    return re.match("[a-z0-9]+(-[a-z0-9])+", read_id) is not None

class ReadIndex:
    def __init__(self, index_filename=None, file_paths=None, read_filter=None, file_suffix=".fast5"):
        self.file_suffix = file_suffix
        self.index_filename = index_filename
        self.file_paths = list()
        self.file_paths_d = dict()
        self.read_files = None

        self._load_filter(read_filter)
        self.load_index_file(index_filename)
        
    @property
    def indexed(self):
        return len(self.file_paths) > 0
    
    def subset(self, read_ids):
        ret = ReadIndex(read_filter=read_ids, file_suffix=self.file_suffix)
        ret.load_index_df(self.read_files.reset_index())
        ret.file_paths = self.file_paths
        ret.file_paths_d = self.file_paths_d
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

        df = df.set_index("read_id").sort_index()

        if self.read_files is None:
            self.read_files = df
        else:
            self.read_files = pd.concat([self.read_files, df])


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
            elif path.endswith(self.file_suffix):
                self._add_read_file(path)

            #Read fast5 filenames from text file
            else:
                with open(path) as infile:
                    for line in infile:
                        self._add_read_file(line.strip())

    #def get_fast5_dict(self):
    #    fast5_reads = dict()
    #    idx = self.read_files.set_index("filename").sort_index()
    #    for fast5 in idx.index.unique():
    #        path = self.file_paths[os.path.basename(fast5)]
    #        reads = idx.loc[fast5, "read_id"]
    #        if isinstance(reads, str):
    #            reads = [reads]
    #        else:
    #            reads = list(reads)
    #        fast5_reads[path] = reads
    #    return fast5_reads

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

    def get_read_file(self, read_id):
        filename = self.read_files.loc[read_id, "filename"]
        return self.file_paths_d[filename]

    def _add_read_file(self, path):
        if not self._is_read_file(path):
            return False

        #TODO make sure fast5 reader doesn't crash on failure
        #if not os.path.isfile(path):
        #    sys.stderr.write("Warning: fast5 file \"%s\" does not exist\n" % fname)
        #    return False
        #TODO is abspath nesissary?

        fname = os.path.basename(path)
        self.file_paths_d[fname] = path

        self.file_paths.append(path)

        return True

class Fast5Reader(ReadIndex):
    #def __init__(self, index_filename=None, file_paths=None, read_filter=None):
    def __init__(self, ri, conf):
        self.read_files = ri.read_files
        self.file_suffix = ".fast5" 
        self.index_filename = ri.index_filename
        self.file_paths = ri.file_paths
        self.file_paths_d = ri.file_paths_d
        self.conf = conf

        #ReadIndex.__init__(self, index_filename, file_paths, read_filter, ".fast5")
        self.infile = None
        self.infile_name = None

    def _open(self, filename):
        if filename == self.infile_name:
            return
        if self.infile is not None:
            self.infile.close()
        self.infile_name = filename
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

    def __iter__(self):
        prev = prev = None
        for read_id, filename in self.read_files["filename"].items():
            if prev != filename:
                path = self.get_read_file(filename)
            yield self._get(read_id, path)

    def __getitem__(self, read_id):
        return self._get(read_id, self.get_read_file(read_id))


    def close(self):
        if self.infile is not None:
            self.infile.close()

