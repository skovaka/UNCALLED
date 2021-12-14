"""Stores datastructures to compute and store normalized events"""
import sys, os, argparse, re
from collections import defaultdict
import pandas as pd
import numpy as np

from .config import ParamGroup, Config
from .pafstats import parse_paf
from . import PoreModel, EventDetector, EventProfiler, Normalizer
from _uncalled import _SignalProcessor, _ProcessedRead

class ProcessedRead(_ProcessedRead):
    def sample_range(self, start, end):
        mask = (self.raw_events["start"] >= start) & (self.raw_events["start"] <= end)
        return self.to_df()[mask]

    def to_df(self):
        df = pd.DataFrame(self.raw_events[["start", "length", "stdv"]])
        df["mean"] = self.norm_events
        return df

class SignalProcessor(_SignalProcessor):
    def process(self, read):
        return ProcessedRead(_SignalProcessor.process(self, read))
