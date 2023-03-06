import sys
import _uncalled
import os
from glob import glob
import numpy as np
import pandas as pd
from time import time

from ...pore_model import PoreModel
from ..aln_track import AlnTrack
from ...aln import LAYER_META, parse_layers
from . import TrackIO

#LAYERS = [("dtw",l) for l in ("kmer", "current", "stdv", "length_ms")]
LAYERS = [("seq","kmer"), ("dtw","current"), ("dtw","stdv"), ("dtw","length_ms")]

class ModelTrainer(TrackIO):
    FORMAT = "model"

    def __init__(self, filename, write, tracks, track_count):
        TrackIO.__init__(self, filename, write, tracks, track_count)

        self.tprms = self.conf.train

        self.row_dtype = [
            (name, np.dtype(dt.lower() if isinstance(dt,str) else dt))
            for (_,name),dt in LAYER_META.loc[LAYERS, "dtype"].items()
        ]
        self.itemsize = sum((d.itemsize for _,d in self.row_dtype))

        self.output = self.input = None

        if self.write_mode:
            self.init_write_mode()
        else:
            self.init_read_mode()

    def init_read_mode(self, load_index=True):
        self.input = open(self._filename("data"), "rb")

        if load_index:
            self.kmer_index = pd.read_csv(self._filename("index"), index_col="kmer", sep="\t")

    def _filename(self, name, itr=None):
        if itr is None:
            itr = self.iter
        return os.path.join(self.filename, f"it{self.iter}.{name}")

    def init_write_mode(self):
        TrackIO.init_write_mode(self)

        self.model = None

        if self.tprms.append and not self.prms.buffered:
            prev_models = glob(f"{self.filename}/it*.model.tsv")
            if len(prev_models) == 0:
                raise ValueError("--append can only be used with existing model training directory")

            max_itr = -1
            fname = None
            for m in prev_models:
                itr = int(os.path.basename(m).split(".")[0][2:])
                if itr > max_itr:
                    fname = m
                    max_itr = itr

            self.iter = max_itr + int(not self.tprms.skip_dtw)
            self.conf.pore_model.name = fname
            self.set_model(PoreModel(self.conf.pore_model))

        elif self.tprms.skip_dtw:
            self.iter = self.tprms.iterations
            self.kmer_counts = None

        else:
            os.makedirs(self.filename, exist_ok=True)
            self.conf.to_toml(os.path.join(self.filename, "conf.toml"))
            self.iter = 1
            self.kmer_counts = None


        #self.index_file = open(self._filename("index"), "w")

        self.buff_len = 0

        self.out_buffer = list()

    def set_model(self, model):
        self.model = model
        self.kmer_counts = pd.Series(0, index=self.model.KMERS)
        self.kmer_index = None #pd.DataFrame({"start" : 0}, index=self.model.KMERS)
        self.full_kmers = set()


    #def write_layers(self, track, groups):
        #if self.model is None:
        #    self.set_model(track.model)

        #track.calc_layers(LAYERS)
    def write_alignment(self, aln):
        
        #mask = track.layers[("cmp", "dist")] <= 1
        #dtw = track.layers.loc[mask, LAYERS]["dtw"].set_index("kmer", drop=True) #\
         #.sort_index() 

        df = aln.to_pandas(LAYERS, index=["seq.kmer"])

        mask = df["mvcmp","dist"] <= 1
        dtw = track.layers.loc[mask, LAYERS]["dtw"].set_index("kmer", drop=True) #\

        if self.prms.buffered:
            self.out_buffer.append(dtw)
        else:
            self.write_buffer([dtw])

    def write_buffer(self, out=[], force=False):
        #t = time()
        #for df in out:
        #    df = df.drop(self.full_kmers, errors="ignore")
        #    kc = df.index.value_counts()
        #    self.kmer_counts[kc.index] += kc.to_numpy()

        #    full = self.kmer_counts >= self.tprms.kmer_samples
        #    if np.any(full):
        #        self.full_kmers.update(self.kmer_counts.index[full])
        #        #self.kmer_counts = self.kmer_counts[~full]

        #    self.out_buffer.append(df.reset_index().to_records(index=False,column_dtypes=dict(self.row_dtype)))
        #    self.buff_len += len(df)

        if len(out) > 0:
            df = pd.concat(out) \
                   .sort_index() \
                   .drop(self.full_kmers, errors="ignore") 


            kc = df.index.value_counts()
            self.kmer_counts[kc.index] += kc.to_numpy()

            full = self.kmer_counts >= self.tprms.kmer_samples
            if np.any(full):
                self.full_kmers.update(self.kmer_counts.index[full])
                #self.kmer_counts = self.kmer_counts[~full]

            self.out_buffer.append(df.reset_index().to_records(index=False,column_dtypes=dict(self.row_dtype)))
            self.buff_len += len(df)

        if self.buff_len == 0 or (self.buff_len * self.itemsize < self.tprms.buffer_size*10**6 and not force):
            return

        if self.output is None:
            self.output = open(self._filename("data"), "wb")

        out = np.concatenate(self.out_buffer)
        out.sort(kind="mergesort")

        kmers, counts = np.unique(out["kmer"], return_counts=True)
        df = pd.DataFrame({"start" : 0, "length" : counts}, index=kmers)
        df["start"].iloc[1:] = counts.cumsum()[:-1]
        if self.kmer_index is None:
            self.kmer_index = df
            df.to_csv(self._filename("index"), sep="\t", index_label="kmer", mode="w")
            #self.index_file.flush()
        else:
            df["start"] += self.kmer_index.iloc[-1].sum()
            df.to_csv(self._filename("index"), sep="\t", index_label="kmer", header=False, mode="a")
            #self.index_file.flush()
            self.kmer_index = pd.concat([self.kmer_index, df])

        self.output.write(out.tobytes())

        self.out_buffer = list()
        self.buff_len = 0

    def is_full(self):
        return self.model is not None and len(self.full_kmers) == self.model.KMER_COUNT

    def next_model(self, load_index=False):
        self.close()
        self.init_read_mode(load_index=load_index)

        self.kmer_index.sort_index(inplace=True)

        t = time()
        model_rows = list()
        for kmer in self.kmer_index.index.unique():
            chunks = self.kmer_index.loc[[kmer]]
            rows = np.zeros(chunks["length"].sum(), dtype=self.row_dtype)

            i = 0
            for _,(start,length) in chunks.iterrows():
                self.input.seek(start*self.itemsize)
                rows[i:i+length] = np.fromfile(self.input, self.row_dtype, length)
                i += length

            if self.tprms.use_median:
                avg = np.median
            else:
                avg = np.mean

            k = self.model.kmer_to_str(kmer)

            model_rows.append((
                kmer,
                avg(rows["current"]),
                avg(rows["stdv"]),
                avg(rows["length_ms"]),
                np.std(rows["current"]),
                np.std(rows["stdv"]),
                np.std(rows["length_ms"]),
                len(rows)
            ))

        df = pd.DataFrame(model_rows, columns=["kmer", "mean", "stdv_mean", "dwell_mean", "stdv", "stdv_stdv", "dwell_stdv", "count"]).set_index("kmer").reindex(self.model.KMERS)

        subs_locs = np.array([0,self.model.K-1])

        #TODO plot pA change for each possible substitution for every kmer

        for kmer in df.index[df["mean"].isna()]:
            subs = list()
            for i in subs_locs:
                old = self.model.kmer_base(kmer, i)
                subs.append(self.model.set_kmer_base(kmer, i, [b for b in range(4) if b != old]))
            subs = np.unique(np.concatenate(subs))
            means = df.loc[subs, "mean"].dropna().sort_values()
            if len(means) > 0:
                mid = means.index[len(means)//2]
                df.loc[kmer] = df.loc[mid]
            else:
                df.loc[kmer] = self.model.to_df().loc[kmer]
            df.loc[kmer, "count"] = 0

        df["count"] = df["count"].astype(int)
            
        model_out = PoreModel(df=df, reverse=self.model.PRMS.reverse, complement=self.model.PRMS.complement, extra_cols=True)
        outfile = self._filename("model.tsv")
        model_out.to_tsv(outfile)

        self.set_model(PoreModel(outfile, reverse=self.model.PRMS.reverse, complement=self.model.PRMS.complement, extra_cols=True))

        self.iter += 1
        print("write", self.model.PRMS.reverse, self.model.PRMS.complement)
        return self.model #outfile#self._filename("model.tsv")
        #return outfile
        

    def close(self):
        if not self.prms.buffered:
            self.write_buffer(force=True)

        if self.output is not None:
            self.output.close()
            self.output = None

            #self.kmer_index.to_csv(, sep="\t", index_label="kmer")

        if self.input is not None:
            self.input.close()
            self.input = None
            #self.index_out.close()
            #self.index_out = None

