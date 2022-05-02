"""Stores datastructures to compute and store normalized events"""
import sys, os, argparse, re
from collections import defaultdict
import pandas as pd
import numpy as np

from .config import ParamGroup, Config
from .pafstats import parse_paf
from . import PoreModel, EventDetector, EventProfiler, Normalizer
from _uncalled import SignalProcessorK5, _ProcessedRead

class ProcessedRead(_ProcessedRead):
    def sample_range(self, start, end):
        mask = (self.events["start"] >= start) & (self.events["start"] <= end)
        return self.to_df()[mask]

    def to_df(self):
        return pd.DataFrame(self.events)

class SignalProcessor(SignalProcessorK5):
    def process(self, read):
        return ProcessedRead(SignalProcessorK5.process(self, read))
