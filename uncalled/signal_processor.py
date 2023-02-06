"""Stores datastructures to compute and store normalized events"""
import sys, os, argparse, re
from collections import defaultdict
import pandas as pd
import numpy as np

from .config import ParamGroup, Config
from .pafstats import parse_paf
from . import EventDetector, EventProfiler, Normalizer
from .pore_model import PoreModel
import _uncalled

class ProcessedRead(_uncalled._ProcessedRead):
    def __init__(self, read):
        if isinstance(read, pd.DataFrame):
            _uncalled._ProcessedRead.__init__(self)
            if not "stdv" in read.columns:
                read["stdv"] = 0
            read = read[["mean","stdv","start","length"]].to_records(index=False)
            self.set_events(read)
        else:
            _uncalled._ProcessedRead.__init__(self, read)

    def set_events(self, df):
        if isinstance(df, pd.DataFrame):
            if not "stdv" in df.columns:
                df["stdv"] = 0
            df = df[["mean","stdv","start","length"]].to_records(index=False)
        _uncalled._ProcessedRead.set_events(self, df)

    def event_bounds(self, samp_start, samp_end):
        start = np.searchsorted(self.events["start"], samp_start)
        end = np.searchsorted(self.events["start"], samp_end)
        if start > 0: start -= 1
        return start, end
        #if start > len(self.events): return start,start
        #if start > 0 && self.events["start"][start]

    def sample_range(self, start, end):
        #mask = (self.events["start"] >= start) & (self.events["start"] <= end)
        #return self.to_df()[mask]
        evt_st,evt_en = self.event_bounds(start, end)
        return self.to_df()[evt_st:evt_en]

    def get_norm_signal(self, samp_min=0, samp_max=None):
        n = self.norm[0]
        if samp_max is None:
            samp_max = len(self.signal)
        norm = self.signal.to_numpy()[samp_min:samp_max]
        return norm

        #ret = np.zeros(int(samp_max - samp_min))
        #i = self.norm[0]["end"].searchsorted(samp_min)

        #while i < len(self.norm):
        #    n = self.norm[i]

        #    st = int(n["start"])
        #    if st >= samp_max: break

        #    en = int(n["end"])-1

        #    st = max(st, samp_min)
        #    en = min(en, samp_max)

        #    raw = self.signal[st:en]
        #    ret[st-samp_min:en-samp_min] = (n["scale"] * raw) + n["shift"]
        #    i += 1
        #return ret

    def to_df(self):
        return pd.DataFrame(self.events)

class SignalProcessor:
    def __init__(self, model, conf):

        self.InstanceClass = getattr(_uncalled, f"SignalProcessorK{model.K}", None)
        if self.InstanceClass is None:
            raise ValueError(f"Invalid k-mer length {model.K}")

        self.model = model

        if isinstance(model, PoreModel):
            model = model.instance

        self.instance = self.InstanceClass(model, conf.event_detector, conf.normalizer)

    def __getattr__(self, name):
        return self.instance.__getattribute__(name)

    def process(self, read, normalize=False):
        return ProcessedRead(self.instance.process(read, normalize))
