import sys
from time import time
import numpy as np
import pandas as pd
import os
import re
from _uncalled import ReadBuffer
from . import config
from uuid import UUID
from types import SimpleNamespace

from ont_fast5_api.fast5_interface import get_fast5_file
import pyslow5
import pod5

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

def is_read_id(read_id):
    return re.match("[a-z0-9]+(-[a-z0-9])+", read_id) is not None

#("paths", None, list, "Paths to fast5, slow5, or pod5 files, or to directories containing those files (optionally recursive)"),
#("read_filter", None, list, "List of read IDs to load, or file containing one read ID per line"),
#("read_index", None, str, "File containing a mapping of read IDs to filenames"),
#("recursive", None, bool, "Recursively search 'paths' for fast5, slow5, or pod5 files"),
#("read_count", None, int, "Maximum number of reads to load"),
#("load_signal", False, bool, "Must be set to true to load signal from FAST5/SLOW5/POD5"),

class ReadIndex:
    #def __init__(self, file_paths=None, read_filter=None, index_filename=None, recursive=False):
    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("read_index", *args, **kwargs)

        self.file_info = dict()
        self.read_files = None
        self.index_files = set()
        self.base_paths = set()

        self._default_model = None
        self.infile = None
        self.infile_name = None
        self.prev_read = None
        self._load_filter(self.prms.read_filter)

        self.load_index_file(self.prms.read_index)

        if self.prms.paths is not None and self.prms.load_signal:
            self.load_paths(self.prms.paths)

    @property
    def default_model(self):
        if self._default_model is None:
            try:
                self._open(next(iter(self.file_info.keys())))   
            except StopIteration:
                return None, None

            flowcell,kit = self.infile.get_run_info()
            self._default_model = (flowcell,kit)

        return self._default_model
        
    @property
    def indexed(self):
        return len(self.file_info) > 0
    
    def subset(self, read_ids):
        ret = ReadIndex(read_count=self.prms.read_count, load_signal=self.prms.load_signal, read_filter=read_ids)                        

        if self.read_files is not None:
            ret.load_index_df(self.read_files.reset_index())
        ret.file_info = self.file_info
        ret._default_model = self._default_model
        return ret

    def load_index_file(self, fname=None):
        if fname is None or len(fname) == 0:
            if self.prms.read_index is None or len(self.prms.read_index) == 0:
                return
            fname = self.prms.read_index

        if fname in self.index_files:
            return
        
        sys.stderr.write(f"Loading read index '{fname}'\n")

        self.index_files.add(fname)

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

        #filenames = df.set_index("read_id")["filename"].str.rsplit("/", n=2, expand=True).iloc[:,-1]
        filenames = df.set_index("read_id")["filename"].apply(os.path.basename)

        if self.read_files is None:
            self.read_files = filenames.sort_index()
        else:
            filenames = filenames[~filenames.index.isin(self.read_files.index)]
            self.read_files = pd.concat([self.read_files, filenames]).sort_index()

    def build_index(self, root, out_fname):
        files = list(self.iter_dir(root))
        if len(files) == 1: return

        sys.stderr.write(f"Writing read index to '{out_fname}'\n")
        index_dfs = list()

        for fname in files:
            infile = self._get_reader(fname)(fname)
            sys.stderr.write(f"Indexing '{fname}'\n")
            index_dfs.append(pd.DataFrame({
                "read_id"  : infile.read_ids,
                "filename" : fname
            }))
            infile.close()
        df = pd.concat(index_dfs)
        self.load_index_df(df)
        df.to_csv(out_fname, sep="\t", index=False)
            
    def iter_dir(self, path):
        if self.prms.recursive:
            for root, dirs, files in os.walk(path):
                for fname in files:
                    if self._is_read_file(fname):
                        yield os.path.join(root, fname)
        else:
            for fname in os.listdir(path):
                if self._is_read_file(fname):
                    yield os.path.join(path, fname)

    def load_paths(self, paths):
        if paths is None:                             
            paths = self.prms.paths                   
                                              
        if paths is None or not self.prms.load_signal:
            return                                    



        if isinstance(paths, str):
            paths = [paths]

        for path in paths:
            path = path.strip()

            if path in self.base_paths:
                continue
            self.base_paths.add(path)

            if not os.path.exists(path):
                sys.stderr.write("Error: \"%s\" does not exist\n" % path)
                continue

            isdir = os.path.isdir(path)

            if isdir and self.prms.read_index is None:
                fname = os.path.join(path, self.prms.default_read_index)
                if os.path.exists(fname):
                    self.load_index_file(fname)
                else:
                    self.build_index(path, fname)

            if isdir:
                for fname in self.iter_dir(path):
                    self._add_read_file(fname)
                #for root, dirs, files in os.walk(path):
                #    for fname in files:
                #        self._add_read_file(os.path.join(root, fname))

            #Non-recursive directory search 
            #elif isdir and not recursive:
            #    for fname in os.listdir(path):
            #        self._add_read_file(os.path.join(path, fname))

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
        if self.read_files is not None:
            filename = self.read_files.loc[read_id]
        elif len(self.file_info) == 1:
            filename, = self.file_info.keys()
        else:
            raise IndexError("Must provide read index TSV file for random access over multiple files")

        return filename
    
    def _get_reader(self, fname):
        return FILE_READERS[fname.split(".")[-1]]

    def _add_read_file(self, path):
        if not self._is_read_file(path):
            return False

        fname = os.path.basename(path)
        Reader = self._get_reader(fname)
        self.file_info[fname] = (Reader, path)

        #if Reader != Fast5Reader:
        #    self._open(fname)
        #    self.load_index_df(pd.DataFrame({
        #        "read_id"  : self.infile.read_ids,
        #        "filename" : fname
        #    }))

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

    def __contains__(self, read_id):
        if self.read_files is not None:
            return read_id in self.read_files.index
        elif len(self.file_info) == 1:
            filename, = self.file_info.keys()
            self._open(filename)
            return read_id in self.infile

    def __getitem__(self, read_id):
        if self.prev_read is not None and self.prev_read.id == read_id: 
            return self.prev_read
        read = None

        if self.has_signal:
            try:
                filename = self.get_read_file(read_id)
                self._open(filename)
                read = self.infile[read_id]
            except Exception as e:
                sys.stderr.write(f"Warning: failed to open read {read_id} ({e})\n")
                read = None

        if read is None:
            read = ReadBuffer(read_id, 0, 0, 0, [])

        self.prev_read = read
        return read

    @property
    def has_signal(self):
        return self.prms.load_signal and len(self.file_info) > 0
        #return self.read_files is not None
    
    def __iter__(self):
        for filename in self.file_info.keys():
            self._open(filename)
            for r in self.infile:
                yield r

    def get(self, read_id, default=None):
        if read_id in self:
            return self[read_id]
        return default

    def close(self):
        if self.infile is not None:
            self.infile.close()

class ReaderBase:
    def __init__(self):
        self._read_ids = None

    @property
    def read_ids(self):
        if self._read_ids is None:
            self._read_ids = self.get_read_ids()
        return self._read_ids

    def __contains__(self, read_id):
        return read_id in self.read_ids

class Fast5Reader(ReaderBase):
    #def __init__(self, file_paths=None, read_filter=None, index_filename=None, recursive=False):
    def __init__(self, filename):
        ReaderBase.__init__(self)
        self.infile = get_fast5_file(filename, mode="r")

    def __getitem__(self, read_id):
        f5 = self.infile.get_read(read_id)
        return self._dict_to_read(f5)

    def get_read_ids(self):
        return self.infile.get_read_ids()
    
    def _dict_to_read(self, f5):
        channel = f5.get_channel_info()
        attrs = f5.handle[f5.raw_dataset_group_name].attrs#["start_time"])
        read = ReadBuffer(f5.read_id, channel["channel_number"], attrs["read_number"], attrs["start_time"], f5.get_raw_data(scale=True))
        #read.filename = filename

        bc_group = f5.get_latest_analysis("Basecall_1D")
        seg_group = f5.get_latest_analysis("Segmentation")
        if bc_group is not None and seg_group is not None:
            try:
                moves = np.array(f5.get_analysis_dataset(bc_group, "BaseCalled_template")["Move"])
                move_attr = f5.get_analysis_dataset(bc_group, "Summary")["basecall_1d_template"].attrs
                stride = move_attr["block_stride"]
                seg_attr = f5.get_analysis_dataset(seg_group, "Summary")["segmentation"].attrs
                template_start = seg_attr["first_sample_template"]
                read.set_moves(moves, template_start, stride)
            except:
                pass

        return read

    def __iter__(self):
        for r in self.infile.get_reads():
            yield self._dict_to_read(r)

    def get_run_info(self):
        rid = self.infile.get_read_ids()[0]
        r = self.infile.get_read(rid)
        tags = r.get_context_tags()
        kit = tags.get("sequencing_kit", None)
        flowcell = tags.get("flowcell_type", None)
        if flowcell is None:
            tags = r.get_tracking_id()
            flowcell = tags.get("flow_cell_product_code", None)
        upper = lambda s: s.upper() if s is not None else None
        return upper(flowcell), upper(kit)

    def close(self):
        self.infile.close()

class Pod5Reader(ReaderBase):
    def __init__(self, filename):
        ReaderBase.__init__(self)
        self.infile = pod5.Reader(filename)

    def __contains__(self, read_id):
        return read_id in self.read_ids

    def get_run_info(self):
        info = self.infile.run_info_table.read_pandas()
        info = info.iloc[0]
        return info["flow_cell_product_code"].upper(), info["sequencing_kit"].upper()

    def get_read_ids(self):
        read_ids = self.infile.read_table.read_pandas()["read_id"]
        return read_ids.apply(lambda r: str(UUID(bytes=r)))

    def __getitem__(self, read_id):
        r = next(self.infile.reads(selection=[read_id]))
        c = r.calibration
        signal = (r.signal + c.offset) * c.scale
        return ReadBuffer(read_id, 0, 0, 0, signal)

    def __iter__(self):
        for r in self.infile.seq_read():
            yield self._dict_to_read(r)

    def close(self):
        self.infile.close()

class Slow5Reader(ReaderBase):
    #def __init__(self, file_paths=None, read_filter=None, index_filename=None, recursive=False):
    def __init__(self, filename):
        ReaderBase.__init__(self)
        self.infile = pyslow5.Open(filename, mode="r")

    def __contains__(self, read_id):
        return read_id in self.read_ids

    def get_run_info(self):
        flowcell = self.infile.get_header_value("flow_cell_product_code")
        kit = self.infile.get_header_value("sequencing_kit")
        return flowcell.upper(), kit.upper()

    def get_read_ids(self):
        return self.infile.get_read_ids()[0]

    def _dict_to_read(self, d):
        return ReadBuffer(d["read_id"], 0, 0, 0, d["signal"])

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
    "blow5" : Slow5Reader,
    "pod5"  : Pod5Reader
}

SUFFIXES = set(FILE_READERS.keys())
