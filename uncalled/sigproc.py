"""Stores datastructures to compute and store normalized events"""
import sys, os, argparse, re
from collections import defaultdict
import pandas as pd
import numpy as np

from .config import ParamGroup, Config
from .pafstats import parse_paf
from _uncalled import EventDetector, EventProfiler, Normalizer, PORE_MODELS

#TODO refactor into SignalProcessor -> ProcRead
#eventually turn into C++

class ProcReadParams(ParamGroup):
    _name = "proc_read"
ProcReadParams._def_params(
    ("detect_events", True, bool, "Will detect events if True"),
    ("profile", False, bool, "Will mask events within low stdv windows if True"),
    ("normalizer", None, None, ""),
)
#Config._EXTRA_GROUPS["proc_read"] = ProcReadParams

class ProcRead:

    #TODO pass conf instad of event_detector/profiler
    def __init__(self, fast5_read, unc_meta=None, conf=None):
        self.f5 = fast5_read
        self.id = self.f5.id

        self.conf = conf if conf is not None else Config()

        #TODO probably need to rethink fwd/rev compl, but either way clean this up
        model_name = self.conf.mapper.pore_model
        if model_name.endswith("_compl"):
            model_name = model_name[:-5]+"templ"
        self.model = PORE_MODELS[model_name]

        self.prms = conf.proc_read
        self.has_events = self.prms.detect_events
        #self.profile_events = profile_events

        if unc_meta is not None:
            self.init_meta(unc_meta)
        else:
            self.init_f5()

        #self.init_bc_mm2(mm2_paf)

    def init_meta(self, meta):
        #self.signal = np.array(meta.proc_signal)
        self.df = pd.DataFrame(meta.events)
        self.mask_idxs = meta.masked_event_idxs
        pass

    #TODO breakup into more functions
    #implement a lot more
    def init_f5(self):
        #TODO pass np array
        signal = np.array(self.f5.signal)

        if self.prms.detect_events:
            evdt = EventDetector(self.conf.event_detector)
            raw_events = evdt.get_events(signal)
            self.store_events(raw_events)
            signal = self.df['mean']
        else:
            self.df = pd.DataFrame(zip(np.arange(len(signal)), signal), columns=['start', 'mean'])
            raw_events = None

        #TODO make profiler work for raw
        #TODO also make profiler return 2D array
        if self.prms.profile and raw_events is not None:

            
            self.df.insert(self.df.shape[1], "prof_mean", 0)
            self.df.insert(self.df.shape[1], "prof_stdv", 0)
            self.df.insert(self.df.shape[1], "mask", False)

            profiler = EventProfiler(self.conf.event_profiler)
            i = 0
            prof = None
            for e in raw_events:
                profiler.add_event(e)
                if profiler.is_full():
                    prof = profiler.get_prof()
                    self.df.loc[i, "prof_mean"] = prof.win_mean
                    self.df.loc[i, "prof_stdv"] = prof.win_stdv
                    self.df.loc[i, "mask"] = prof.mask
                    i += 1

            for j in range(i, len(self.df)):
                self.df.loc[j, "prof_mean"] = prof.win_mean
                self.df.loc[j, "prof_stdv"] = prof.win_stdv
                self.df.loc[j, "mask"] = prof.mask

            self.mask_idxs = np.arange(len(self.df))[self.df['mask'].to_numpy(bool)]

            signal = signal[self.df['mask']]
        else:
            self.df['mask'] = np.ones(len(self.df), bool)

        norm_params = list()

        #if conf.normalizer.len
        if self.prms.normalizer is None and self.conf.normalizer.len > 0:
            model = PORE_MODELS[self.conf.mapper.pore_model]
            self.conf.normalizer.tgt_mean = model.get_means_mean()
            self.conf.normalizer.tgt_stdv = model.get_means_stdv() #TODO set when pore model is set
            self.prms.normalizer = Normalizer(self.conf.normalizer)


        #Full-read normalization
        if self.prms.normalizer is None:
            self.df.insert(self.df.shape[1], "norm_sig", 0)
            self.scale = self.model.get_means_stdv() / np.std(signal)
            self.shift = self.model.get_means_mean() - self.scale * np.mean(signal)
            self.df.loc[self.df['mask'], "norm_sig"] = self.scale * signal + self.shift
            norm_params.append( (0, len(self.f5), self.scale, self.shift) )
        
        #Pre-defined shift/scale
        elif (isinstance(self.prms.normalizer, (tuple,list,np.ndarray)) and len(norm) == 2):
            self.df.insert(self.df.shape[1], "norm_sig", 0)
            self.scale,self.shift = norm
            self.signal = self.scale * signal + self.shift
            norm_params.append( (0, len(self.f5), self.scale, self.shift) )

        elif type(self.prms.normalizer) == Normalizer:
            self.df.insert(self.df.shape[1], "norm_sig", 0)
            norm = self.prms.normalizer
            chunk_len = int(self.conf.chunk_time * self.conf.sample_rate)
            i = 0
            norm_shift = list()
            norm_scale = list()
            norm_sig = list()
            while i < len(self.df):
                if self.has_events:
                    j = self.df['start'].searchsorted(self.df.loc[i, 'start'] + chunk_len)
                else:
                    j = min(i + chunk_len, len(self.df))
                idxs = np.arange(i,j)[self.df['mask'][i:j]]
                added = norm.push(signal[idxs])

                norm_params.append( (self.df['start'][i], self.df['start'][j-1]+self.df['length'][j-1], norm.get_scale(), norm.get_shift()) )

                while not norm.empty():
                    norm_sig.append(norm.pop())
                    i += 1

            self.df.loc[self.df['mask'], "norm_sig"] = norm_sig

        elif isinstance(self.prms.normalizer, pd.DataFrame):
            #TODO
            #linear regression from DTW
            #maybe also allow list of tuples, two lists, 2D array
            #chunked linreg? can test with read segments, maybe build in
            pass

        self.norm_params = pd.DataFrame(norm_params, columns=["start", "end", "scale", "shift"])
        self.norm_params.start.astype(int, False)
        self.norm_params.end.astype(int, False)

    def store_events(self, events):
        self.df = pd.DataFrame(
            [ [e.start, e.length, e.mean, e.stdv] for e in events ],
            columns = ['start','length','mean','stdv']
        )

    def sample_range(self, start, end):
        return self.df.loc[(self.df['start'] >= start) & (self.df['start'] <= end) & self.df['mask']]

    def get_norm_signal(self, samp_min, samp_max):
        ret = np.zeros(samp_max - samp_min)

        i = self.norm_params["end"].searchsorted(samp_min)

        while i < len(self.norm_params):
            n = self.norm_params.iloc[i]

            st = int(n["start"])
            if st >= samp_max: break

            en = int(n["end"])-1

            st = max(st, samp_min)
            en = min(en, samp_max)

            raw = self.f5.signal[st:en]
            ret[st-samp_min:en-samp_min] = (n["scale"] * raw) + n["shift"]
            i += 1

        return ret


    def plot_events(self, ax, samp_min=None, samp_max=None):
        if samp_min is not None and samp_max is not None:
            i = (self.df['start'] >= samp_min) & (self.df['start'] <= samp_max)
        else:
            i = np.arange(len(self.df))
        ax.step(self.df['start'][i], self.df['mean'][i], where='post', color='black', zorder=1)

        if 'prof_mean' in self.df:
            prof_mins = self.df['prof_mean'][i] - self.df['prof_stdv'][i]*3
            prof_maxs = self.df['prof_mean'][i] + self.df['prof_stdv'][i]*3
            ax.fill_between(
                self.df['start'][i], 
                prof_mins,
                prof_maxs,
                self.df['mask'][i],
                color='blue', 
                alpha=0.5,
                zorder=3
            )
            ax.fill_between(
                self.df['start'][i], 
                prof_mins, 
                prof_maxs, 
                ~self.df['mask'][i],
                color='red', 
                alpha=0.5,
                zorder=3
            )

            #TODO optionally highlight events (for mapping+other?)

    def event_to_samp(self, evt, mask=True):
        i = self.mask_idxs[evt] if mask else evt
        s = self.df.iloc[i]["start"]
        if hasattr(s,"to_numpy"):
            return s.to_numpy()
        return s

    def samp_to_event(self, evt):
        return self.df['start'].searchsorted(evt)

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description="Calculates enrichment from even/odd read until")
    parser.add_argument("fast5s", nargs="+", type=str)
    parser.add_argument("-r", "--recursive", action="store_true")
    parser.add_argument("-l", "--read-filter", required=False, type=str, default=None)
    parser.add_argument("-c", "--start-chunk", required=False, type=int, default=0)
    parser.add_argument("-C", "--max-chunks", required=False, type=int, default=10)
    parser.add_argument("-m", "--mm2-paf", required=False, type=str, default=None)
    args = parser.parse_args()

    plt.style.use(['seaborn'])

    conf = Config()
    conf.read_buffer.start_chunk = args.start_chunk
    conf.read_buffer.max_chunks = args.max_chunks
    conf.fast5_reader.load_bc = True
    conf.export_static()

    mm2s = defaultdict(list)
    if args.mm2_paf is not None:
        for p in parse_paf(args.mm2_paf):
            mm2s[p.qr_name].append(p)
    else:
        mm2s = dict()

    fast5s = Fast5Processor(
        args.fast5s, 
        args.read_filter, 
        args.recursive, 
        conf
    )

    for read in fast5s:
        if read.id in mm2s:
            fig,(ax1,ax2) = plt.subplots(2,1)
            for mm2 in mm2s[read.id]:
                aln = BcRawAln(read, mm2)
                aln.plot_step(ax2)
            ax1.get_shared_x_axes().join(ax1, ax2)
        else:
            fig,ax1 = plt.subplots()

        read.plot_events(ax1)
        print(read.id)
        plt.show()
        plt.close()

