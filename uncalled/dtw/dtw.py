#!/usr/bin/env python3

import sys, os
import numpy as np
import argparse
from collections import defaultdict
import re
import time
from typing import NamedTuple
from matplotlib.colors import Normalize
import pandas as pd
import scipy.stats
import copy

from ..pafstats import parse_paf, PafEntry
from ..config import Config
from _uncalled import PORE_MODELS, BwaIndex, nt

class ReadAln:

    REFMIR_COL  = "refmir"
    REF_COL    = "ref"
    START_COL  = "start"
    LENGTH_COL = "length"
    MEAN_COL   = "mean"
    KMER_COL   = "kmer"

    def __init__(self, index, aln, df=None, is_rna=False, ref_bounds=None):
        if not isinstance(aln, PafEntry):
            raise RuntimeError("ReadAlns can only be initialized from PafEntrys currently")

        self.index = index
        self.read_id = aln.qr_name
        self.is_rna = is_rna

        self.clipped = False
    
        self.set_ref_bounds(aln, ref_bounds)

        if self.empty: 
            return

        self._init_mirror_coords()

        if df is not None:
            self.df = df[(df.index >= self.ref_start) & (df.index <= self.ref_end)]

            has_ref = self.df.index.name == self.REF_COL
            has_refmir = self.REFMIR_COL in self.df.columns

            #TODO check for required columns
            if not has_ref and not has_refmir:
                raise RuntimeError("ReadAln DataFrame must include a column named \"%s\" or \"%s\"" % (self.REF_COL, self.REFMIR_COL))
            
            if has_ref and not has_refmir:
                self.df[self.REFMIR_COL] = self.index.ref_to_refmir(self.ref_id, self.df.index, self.is_fwd, self.is_rna)

            elif not has_ref and has_refmir:
                self.df[REF_COL] = self.index.mirref_to_ref(self.df[REFMIR_COL])

            self.df.sort_values("refmir", inplace=True)
        

    def set_ref_bounds(self, aln, ref_bounds):
        if ref_bounds is None:
            self.ref_bounds = (aln.rf_name, aln.rf_st, aln.rf_en, aln.is_fwd)
        else:
            if aln.rf_st < ref_bounds[1]:
                ref_st = ref_bounds[1]
                clipped = True
            else:
                ref_st = aln.rf_st

            if aln.rf_en > ref_bounds[2]:
                ref_en = ref_bounds[2]
                clipped = True
            else:
                ref_en = aln.rf_en

            if ref_st > ref_en:
                self.empty = True
                return

            self.ref_bounds = (aln.rf_name, ref_st, ref_en, aln.is_fwd)

        self.ref_id = self.index.get_ref_id(self.ref_name)

        self.empty = False

    #TODO ReadAln stores RefCoords, handles all this conversion?

    def refmir_to_ref(self, refmir):
        return self.index.refmir_to_ref(refmir) 

    def ref_to_refmir(self, ref):
        return self.index.ref_to_refmir(self.ref_name, ref, ref, self.is_fwd, self.is_rna)[0]

    def ref_to_samp(self, ref):
        return self.refmir_to_samp(self.ref_to_refmir(ref))
        
    def refmir_to_samp(self, refmir):
        i = np.clip(self.df['refmir'].searchsorted(refmir), 0, len(self.df)-1)
        return self.df['sample'][i]

    def _init_mirror_coords(self):
        self.refmir_start, self.refmir_end = self.index.ref_to_refmir(*self.ref_bounds, self.is_rna)

    @property
    def ref_name(self):
        """The sequence name of the alignment reference coordinate"""
        return self.ref_bounds[0]

    @property
    def ref_start(self):
        """The start of the alignment reference coordinate"""
        return self.ref_bounds[1]

    @property
    def ref_end(self):
        """The end of the alignment reference coordinate"""
        return self.ref_bounds[2]

    @property
    def is_fwd(self):
        return self.ref_bounds[3]
    
    def sort_ref(self):
        self.df.sort_values(self.REF_COL, inplace=True)

    def sort_refmir(self):
        self.df.sort_values(self.REFMIR_COL, inplace=True)

    def get_samp_bounds(self):
        samp_min = self.df['start'].min()
        max_i = self.df['start'].argmax()
        samp_max = self.df['start'].iloc[max_i] + self.df['length'].iloc[max_i]
        return samp_min, samp_max
    
    #def set_bands(self, bands):

    def set_subevent_aln(self, aln, ref_mirrored=False, kmer_str=False, ref_col=REFMIR_COL, start_col=START_COL, length_col=LENGTH_COL, mean_col=MEAN_COL, kmer_col=KMER_COL):

        aln["cuml_mean"] = aln[length_col] * aln[mean_col]

        grp = aln.groupby(ref_col)

        if kmer_str:
            kmers = [nt.kmer_rev(nt.str_to_kmer(k,0)) for k in grp[kmer_col].first()]
        else:
            kmers = grp[kmer_col].first()

        if ref_mirrored:
            refmirs = grp[ref_col].first()
            refs = self.refmir_to_ref(refmirs)
        else:
            refs = grp[ref_col].first()

        lengths = grp[length_col].sum()

        self.df = pd.DataFrame({
            self.REF_COL    : refs,
            self.KMER_COL   : kmers,
            self.START_COL  : grp[start_col].min(),
            self.LENGTH_COL : lengths,
            self.MEAN_COL   : grp["cuml_mean"].sum() / lengths
        })
        
        if ref_mirrored:
            self.df[self.REFMIR_COL] = refmirs

        self.df = self.df.set_index(self.REF_COL).sort_values(self.REF_COL)
        

    def get_index_kmers(self, index, kmer_shift=4):
        """Returns the k-mer sequence at the alignment reference coordinates"""
        start = self.refmir_start - kmer_shift

        if start < 0:
            lpad = -start
            start = 0
        else:
            lpad = 0

        print(start, self.refmir_end, self.is_rna, "AHAHSDF")
        kmers = index.get_kmers(start, self.refmir_end, self.is_rna)
        return np.array(kmers)

        #return np.insert(kmers, 0, [0]*lpad)
    
class Track:
    CONF_FNAME = "conf.toml"
    ALN_DIR = "alns"
    ALN_SUFFIX = ".pkl"

    INDEX_FNAME = "filename_mapping.txt"
    INDEX_HEADER = "read_id\tfilename\taln_file"

    #TODO make fast5_file alias for filename in Fast5Reader
    #also, in fast5 reader make static fast5 filename parser

    WRITE_MODE = "w"
    READ_MODE = "r"
    MODES = {WRITE_MODE, READ_MODE}

    KMER_LAYER = 0
    PA_LAYER = 1
    DWELL_LAYER = 2
    PA_DIFF_LAYER = 3

    KS_LAYERS = [PA_LAYER, DWELL_LAYER]

    LAYER_META = [
        ("K-mer",              False),
        ("Current (pA)",       True),
        ("Dwell Time (ms/bp)", True),
        ("pA Difference",      True)
    ]

    def __init__(self, path, mode="r", ref_bounds=None, full_overlap=None, conf=None, overwrite=False, index=None, mm2s=None):
        self.path = path.strip("/")
        self.mode = mode
        self.index = index

        self.conf = Config(conf)

        if mode == self.WRITE_MODE:
            os.makedirs(self.aln_dir, exist_ok=overwrite)

        self.fname_mapping_file = open(self.fname_mapping_filename, mode)

        if mode == self.READ_MODE:
            self.conf.load_toml(self.config_filename)

            if len(self.conf.fast5_reader.fast5_index) == 0:
                self.conf.fast5_reader.fast5_index = self.fname_mapping_filename

            self._load_index()

        elif mode == self.WRITE_MODE:
            self.conf.to_toml(self.config_filename)
            self.fname_mapping_file.write(self.INDEX_HEADER + "\n")

        print(ref_bounds)
        #TODO make TrackParams, deal with kwargs
        if ref_bounds is None:
            ref_bounds = self.conf.align.ref_bounds

        #TODO arguments overload conf params
        if self.conf.align.mm2_paf is not None:
            read_filter = set(self.conf.fast5_reader.read_filter)
            self.mm2s = {p.qr_name : p
                     for p in parse_paf(
                        self.conf.align.mm2_paf,
                        ref_bounds=self.conf.align.ref_bounds,
                        full_overlap=self.conf.browser.full_overlap,
            )}

        #TODO static bwa_index parameters, or instance
        if self.index is None and len(self.conf.mapper.bwa_prefix) > 0:
            self.index = BwaIndex(self.conf.mapper.bwa_prefix, True)

        #TODO PROBABLY NEED TO RETHINK FWD/REV COMPL, BUT EITHER WAY CLEAN THIS UP
        model_name = self.conf.mapper.pore_model
        if model_name.endswith("_compl"):
            model_name = model_name[:-5]+"templ"
        self.model = PORE_MODELS[model_name]

        print(ref_bounds)

        if ref_bounds is not None:
            self.load_region(ref_bounds, full_overlap)

    @property
    def read_ids(self):
        return list(self.fname_mapping.index)

    def _load_index(self):
        self.fname_mapping = pd.read_csv(self.fname_mapping_file, sep="\t", index_col="read_id")

    def save_aln(self, aln, fast5_fname):
        if self.mode != "w":
            raise RuntimeError("Must be write mode to add read to track")

        aln_fname = self.aln_fname(aln.read_id)
        aln.df.sort_index().to_pickle(aln_fname)
        self.fname_mapping_file.write("\t".join([aln.read_id, fast5_fname, aln_fname]) + "\n")

    def load_aln(self, read_id, ref_bounds=None):
        mm2 = self.mm2s[read_id]
        df = pd.read_pickle(self.aln_fname(read_id)).sort_index()
        return ReadAln(self.index, mm2, df, is_rna=not self.conf.read_buffer.seq_fwd)

    #TODO parse mm2 every time to enable changing bounds
    #eventually use some kind of tabix-like indexing
    def load_region(self, ref_bounds=None, partial_overlap=False):

        #self.mat = TrackMatrix(self, ref_bounds, self.mm2s, conf=self.conf)

        self.ref_bounds = ref_bounds
        self.width = self.ref_end-self.ref_start
        self.height = None

        read_meta = defaultdict(list)

        read_rows = defaultdict(list)
        mask_rows = list()

        self.ref_coords = pd.Series(
                np.arange(self.width, dtype=int),
                index=pd.RangeIndex(self.ref_start, self.ref_end)
        )

        def _add_row(df, mm2_paf):

            dtw_roi = df.loc[self.ref_start:self.ref_end-1]

            xs = self.ref_coords.reindex(
                self.ref_coords.index.intersection(dtw_roi.index)
            )

            roi_mask = np.ones(self.width)#, dtype=bool)
            roi_mask[xs] = False
            mask_rows.append(roi_mask)

            pa_diffs = self.model.match_diff(dtw_roi['mean'], dtw_roi['kmer'])
            dwell = 1000 * dtw_roi['length'] / self.conf.read_buffer.sample_rate

            def _add_layer_row(layer, vals):
                row = np.zeros(self.width)
                row[xs] = vals
                read_rows[layer].append(row)

            _add_layer_row(self.KMER_LAYER, dtw_roi['kmer'])
            _add_layer_row(self.PA_LAYER, dtw_roi['mean'])
            _add_layer_row(self.PA_DIFF_LAYER, pa_diffs)
            _add_layer_row(self.DWELL_LAYER, dwell)

            read_meta['ref_start'].append(df.index.min())
            read_meta['id'].append(mm2_paf.qr_name)
            read_meta['fwd'].append(mm2_paf.is_fwd)

            #self.mm2s[mm2_paf.qr_name] = mm2_paf

        for read_id,read in self.fname_mapping.iterrows():
            mm2 = self.mm2s.get(read_id, None)
            if mm2 is None: 
                continue
            df = pd.read_pickle(read["aln_file"])
            df.sort_index(inplace=True)
            _add_row(df, mm2)

        self.reads = pd.DataFrame(read_meta) \
                     .sort_values(['fwd', 'ref_start'], ascending=[False, True])

        self.has_fwd = np.any(self.reads['fwd'])
        self.has_rev = not np.all(self.reads['fwd'])

        read_order = self.reads.index.to_numpy()

        layer_names = sorted(read_rows.keys())

        mat = np.stack([
            np.stack(read_rows[l])[read_order] for l in layer_names
        ])
        mask = np.stack([np.stack(mask_rows)[read_order].astype(bool) for _ in layer_names])

        self.mat = np.ma.masked_array(mat, mask=mask)

        self.height = len(self.reads)

        #self.layer_extrema = np.array([
        #    [np.min(layer), np.max(layer)]
        #    for layer in self.mat
        #])
        #TODO define get_soft_minmax(num_stdvs):

        self.norms = [Normalize(np.min(layer), np.max(layer)) for layer in self.mat]
        for l in [self.PA_LAYER, self.DWELL_LAYER]:
            layer = self.mat[l]
            self.norms[l].vmax = min(
                layer.max(),
                np.ma.median(layer) + 2 * layer.std()
            )

        return mat


    def close(self):
        self.fname_mapping_file.close()

    @property
    def config_filename(self):
        return os.path.join(self.path, self.CONF_FNAME)

    @property
    def fname_mapping_filename(self):
        return os.path.join(self.path, self.INDEX_FNAME)
    
    @property
    def aln_dir(self):
        return os.path.join(self.path, self.ALN_DIR)
    
    @property
    def read_count(self):
        return len(self.fname_mapping)
    
    def aln_fname(self, read_id):
        return os.path.join(self.aln_dir, read_id+self.ALN_SUFFIX)

    def sort(self, layer, ref):
        order = np.argsort(self.mat[layer,:,ref])
        self.mat = self.mat[:,order,:]
        self.reads = self.reads.iloc[order]

    @property
    def ref_name(self):
        return self.ref_bounds[0]

    @property
    def ref_start(self):
        return self.ref_bounds[1]

    @property
    def ref_end(self):
        return self.ref_bounds[2]

    def calc_ks(self, track_b):
        ks_stats = np.zeros((len(self.KS_LAYERS), self.width))

        for i,l in enumerate(self.KS_LAYERS):
            for rf in range(self.width):
                a = self[l,:,rf]
                b = track_b[l,:,rf]
                ks_stats[i,rf] = scipy.stats.ks_2samp(a,b,mode="asymp")[0]

        return ks_stats

    def __getitem__(self, i):
        return self.mat.__getitem__(i)
    
    @property
    def kmer(self):
        return self.mat[self.KMER_LAYER]

    @property
    def pa(self):
        return self.mat[self.PA_LAYER]
    
    @property
    def dwell(self):
        return self.mat[self.DWELL_LAYER]
    
    @property
    def pa_diff(self):
        return self.mat[self.PA_DIFF_LAYER]



def ref_coords(coord_str):
    spl = coord_str.split(":")
    ch = spl[0]
    st,en = spl[1].split("-")

    coord = (ch, int(st), int(en))

    if len(spl) == 2:
        return coord
    else:
        return coord + (spl[2] == "+",)

class BcFast5Aln(ReadAln):
    BCE_K = 4
    CIG_OPS_STR = "MIDNSHP=X"
    CIG_RE = re.compile("(\d+)(["+CIG_OPS_STR+"])")
    CIG_OPS = set(CIG_OPS_STR)
    CIG_INCR_ALL = {'M','=', 'X'}
    CIG_INCR_RD = CIG_INCR_ALL | {'I','S'}
    CIG_INCR_RF = CIG_INCR_ALL | {'D','N'}

    SUB = 0
    INS = 1
    DEL = 2
    ERR_TYPES = [SUB, INS, DEL]
    ERR_MARKS = ['o', 'P', '_']
    ERR_SIZES = [100, 150, 150]
    ERR_WIDTHS = [0,0,5]

    def __init__(self, index, read, paf, ref_bounds=None):
        self.seq_fwd = read.conf.read_buffer.seq_fwd #TODO just store is_rna
        ReadAln.__init__(self, index, paf, is_rna=not self.seq_fwd, ref_bounds=ref_bounds)
        if self.empty: return

        self.refgap_bps = list()
        self.sub_bps = list()
        self.ins_bps = list()
        self.del_bps = list()
        self.err_bps = None

        self.flip_ref = paf.is_fwd != self.seq_fwd

        self.empty = (
                not read.f5.bc_loaded or 
                (not self.parse_cs(paf) and
                 not self.parse_cigar(paf))
        )
        if self.empty: 
            return

        #TODO make c++ build this 
        moves = np.array(read.f5.moves, bool)
        bce_qrs = np.cumsum(read.f5.moves)
        bce_samps = read.f5.template_start + np.arange(len(bce_qrs)) * read.f5.bce_stride

        samp_bps = pd.DataFrame({
            'sample' : bce_samps,#[moves],
            'bp'     : np.cumsum(read.f5.moves),#[moves],
        })

        self.df = samp_bps.join(self.bp_refmir_aln, on='bp').dropna()
        self.df.reset_index(inplace=True, drop=True)

        if self.err_bps is not None:
            self.errs = samp_bps.join(self.err_bps.set_index('bp'), on='bp').dropna()
            self.errs.reset_index(inplace=True, drop=True)
        else:
            self.errs = None

        self.ref_gaps = self.df[self.df['bp'].isin(self.refgap_bps)].index

        self.subs = self.df[self.df['bp'].isin(self.sub_bps)].index
        self.inss = self.df[self.df['bp'].isin(self.ins_bps)].index
        self.dels = self.df[self.df['bp'].isin(self.del_bps)].index

        self.empty = len(self.df) == 0
        if self.empty: 
            return

        self.y_min = -paf.rf_en if self.flip_ref else paf.rf_st
        self.y_max = self.y_min + self.df['refmir'].max()

    def parse_cs(self, paf):
        cs = paf.tags.get('cs', (None,)*2)[0]
        if cs is None: return False

        sys.stderr.write("Loading cs tag\n")

        #TODO rename to general cig/cs
        bp_refmir_aln = list()
        err_bps = list()

        if self.seq_fwd:
            qr_i = paf.qr_st
            #rf_i = paf.rf_st
        else:
            qr_i = paf.qr_len - paf.qr_en 
            #rf_i = -paf.rf_en+1

        mr_i = self.refmir_start

        cs_ops = re.findall("(=|:|\*|\+|-|~)([A-Za-z0-9]+)", cs)

        if paf.is_fwd != self.seq_fwd:
            cs_ops = reversed(cs_ops)

        for op in cs_ops:
            c = op[0]
            if c in {'=',':'}:
                l = len(op[1]) if c == '=' else int(op[1])
                bp_refmir_aln += zip(range(qr_i, qr_i+l), range(mr_i, mr_i+l))
                qr_i += l
                mr_i += l

            elif c == '*':
                self.sub_bps.append(qr_i)
                bp_refmir_aln.append((qr_i,mr_i))
                err_bps.append( (qr_i,mr_i,self.SUB) )
                qr_i += 1
                mr_i += 1

            elif c == '-':
                self.ins_bps.append(qr_i)
                err_bps.append( (qr_i,mr_i,self.DEL) )
                l = len(op[1])
                mr_i += l

            elif c == '+':
                self.del_bps.append(qr_i)
                err_bps.append( (qr_i,mr_i,self.INS) )

                l = len(op[1])
                qr_i += l

            elif c == '~':
                l = int(op[1][2:-2])
                self.refgap_bps.append(qr_i)
                mr_i += l

            else:
                print("UNIMPLEMENTED ", op)

        self.bp_refmir_aln = pd.DataFrame(bp_refmir_aln, columns=["bp","refmir"], dtype='Int64')
        self.bp_refmir_aln.set_index("bp", inplace=True)

        #TODO type shouldn't have to be 64 bit
        self.err_bps = pd.DataFrame(err_bps, columns=["bp","refmir","type"], dtype='Int64')

        return True        


    def parse_cigar(self, paf):
        cig = paf.tags.get('cg', (None,)*2)[0]
        if cig is None: return False

        bp_refmir_aln = list()#defaultdict(list)
        self.refgap_bps = list()

        #mr_i = self.refmir_start
        if self.seq_fwd:
            qr_i = paf.qr_st
        else:
            qr_i = paf.qr_len - paf.qr_en 

        if self.flip_ref:
            mr_i = self.ref_to_refmir(paf.rf_en)
        else:
            mr_i = self.ref_to_refmir(paf.rf_st)

        mr_bounds = range(self.refmir_start, self.refmir_end)

        cig_ops = self.CIG_RE.findall(cig)

        if paf.is_fwd != self.seq_fwd:
            cig_ops = list(reversed(cig_ops))

        for l,c in cig_ops:
            l = int(l)
            incr_qr = c in self.CIG_INCR_RD
            incr_rf = c in self.CIG_INCR_RF
            qr_j = qr_i + (l if incr_qr else 1)
            mr_j = mr_i + (l if incr_rf else 1)

            if c == "M":
                for qr, mr in zip(range(qr_i, qr_j), range(mr_i, mr_j)):
                    if mr in mr_bounds:
                        bp_refmir_aln.append((qr,mr))
                #bp_refmir_aln += zip(range(qr_i, qr_j), range(mr_i, mr_j))
            elif c == "N":
                if mr_i in mr_bounds:
                    bp_refmir_aln.append((qr_i,mr))

            if incr_qr:
                qr_i = qr_j 

            if incr_rf:
                mr_i = mr_j 

        self.bp_refmir_aln = pd.DataFrame(bp_refmir_aln, columns=["bp","refmir"], dtype='Int64')
        self.bp_refmir_aln.set_index("bp", inplace=True)

        return True

    def get_xy(self, i):
        df = self.df.loc[i]
        return (df['sample'], df['refmir']-0.5)

    def plot_scatter(self, ax, real_start=False, samp_min=None, samp_max=None):
        if samp_min is None: samp_min = 0
        if samp_max is None: samp_max = self.df['sample'].max()
        i = (self.df['sample'] >= samp_min) & (self.df['sample'] <= samp_max)

        return ax.scatter(self.df['sample'][i], self.df['refmir'][i], color='orange', zorder=2,s=20)

    def plot_step(self, ax, real_start=False, samp_min=None, samp_max=None):
        i = (self.df['sample'] >= samp_min) & (self.df['sample'] <= samp_max)

        ret = ax.step(self.df['sample'][i], self.df['refmir'][i], color='orange', zorder=1, where='post')

        if self.errs is not None:
            for t in self.ERR_TYPES:
                e = self.errs[self.errs['type'] == t]
                ax.scatter(
                    e['sample'], e['refmir'], 
                    color='red', zorder=3, 
                    s=self.ERR_SIZES[t],
                    linewidth=self.ERR_WIDTHS[t],
                    marker=self.ERR_MARKS[t]
                )

        return ret

