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
    def sample_range(self, start, end):
        mask = (self.events["start"] >= start) & (self.events["start"] <= end)
        return self.to_df()[mask]

    def to_df(self):
        return pd.DataFrame(self.events)

KMER_CLASSES = {
    5  : SignalProcessorK5,
    10 : SignalProcessorK10,
}

class SignalProcessor:
    def __init__(self, model, prms):
        self.InstanceClass = KMER_CLASSES[model.K]

        if isinstance(model, PoreModel):
            model = model.instance

        self.instance = self.InstanceClass(model, prms)

    def __getattr__(self, name):
        return self.instance.__getattribute__(name)

    def process(self, read):
        return ProcessedRead(self.instance.process(read))
