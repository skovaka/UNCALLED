"""Stores datastructures to compute and store normalized events"""
import sys, os, argparse, re
from collections import defaultdict
import pandas as pd
import numpy as np

from .config import ParamGroup, Config
from .pafstats import parse_paf
from . import PoreModel, EventDetector, EventProfiler, Normalizer
from _uncalled import SignalProcessorK5, SignalProcessorK10, _ProcessedRead

class ProcessedRead(_ProcessedRead):
    def __init__(self, read, raw):
        _ProcessedRead.__init__(self, read)
        self.signal = raw.signal

    def sample_range(self, start, end):
        mask = (self.events["start"] >= start) & (self.events["start"] <= end)
        return self.to_df()[mask]

    def get_norm_signal(self, samp_min, samp_max):
        n = self.norm[0]
        print(self.norm)
        return self.signal[samp_min:samp_max] * n["scale"] + n["shift"]

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

KMER_CLASSES = {
    5  : SignalProcessorK5,
    10 : SignalProcessorK10,
}

class SignalProcessor:
    def __init__(self, model, conf):
        self.InstanceClass = KMER_CLASSES[model.K]

        if isinstance(model, PoreModel):
            model = model.instance

        print(model)

        self.instance = self.InstanceClass(model, conf.event_detector)

    def __getattr__(self, name):
        return self.instance.__getattribute__(name)

    def process(self, read):
        return ProcessedRead(self.instance.process(read), read)
