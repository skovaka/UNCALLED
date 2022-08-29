import numpy as np
import pandas as pd
import sys
from ..aln_track import AlnTrack
from ...signal_processor import ProcessedRead
from . import TrackIO
import _uncalled

class Eventalign(TrackIO):
    FORMAT = "eventalign"

    def __init__(self, conf, mode):
        filename = conf.tracks.io.eventalign_in if mode == "r" else conf.tracks.io.eventalign_out
        TrackIO.__init__(self, filename, conf, mode)

        if self.filename == "-":
            self.out = sys.stdout
        else:
            self.out = open(self.filename, mode)

        self._header = True

        if self.write_mode:
            self.init_write_mode()
        else:
            self.init_read_mode()

    def init_write_mode(self):
        TrackIO.init_write_mode(self)

        flags = set(self.prms.eventalign_flags)
        self.write_read_name = "print-read-names" in flags
        self.write_signal_index = "signal-index" in flags
        self.write_samples = "samples" in flags

        header = ["contig", "position", "reference_kmer"]

        if self.write_read_name:
            header.append("read_name")
        else:
            header.append("read_index")

        header += ["strand", "event_index", "event_level_mean", "event_stdv", "event_length", "model_kmer", "model_mean", "model_stdv", "standardized_level"]

        if self.write_signal_index:
            header += ["start_idx", "end_idx"]

        if self.write_samples:
            header += ["samples"]

        self.out.write("\t".join(header) + "\n")

    def init_read_mode(self):
        if len(self.track_names) != 1:
            raise ValueError("Can only read eventalign TSV into a single track")
        name = self.track_names[0]
        self.init_track(None, name, name, self.conf.to_toml())
        #t = AlnTrack(self, None, name, name, self.conf)
        #self.tracks.append(t)

    #def write_dtw_events(self, track, events):
    #def write_layers(self, df, read, index=["pac","aln_id"]):
    def write_layers(self, track, groups):
        if "dtw" not in groups: 
            return

        df = track.layers["dtw"].dropna(how="all").sort_index()

        pacs = df.index.get_level_values(0)
        events = df.droplevel(1)#.set_index(track.coords.pac_to_ref(pacs))

        model = track.model
        kmers = events["kmer"]

        std_level = (events["current"] - model.model_mean) / model.model_stdv

        event_counts = events["events"]
        event_index = (event_counts.cumsum() - event_counts.iloc[0]).astype(int)
        if True: #TODO check for flipped ref
            event_index = event_index.max() - event_index
        event_index += 1

        evts = events.rename(columns={"current" : "mean", "current_stdv" : "stdv"})
        if self.read is not None:
            self.read.set_events(evts)
            read = self.read
        else:
            read = ProcessedRead(evts)

        if self.write_samples:
            signal = read.get_norm_signal()
            if len(signal) == 0:
                raise RuntimeError("Failed to output read signal")
        else:
            signal = []

        if self.write_read_name:
            read_id = track.alignments["read_id"].iloc[0]
        else:
            read_id = str(self.prev_aln_id)

        eventalign = _uncalled.write_eventalign_K5(
            self.conf, model.instance, read_id, track.coords.fwd, read,
            track.coords.ref_name, events.index-2, 
            self.write_signal_index,
            kmers, 
            event_index, #TODO properly rep skips?
            std_level, signal) #TODO compute internally?

        self.out.write(eventalign)

    def write_alignment(self, alns):
        pass

    def init_fast5(self, fast5):
        pass

    def init_read(self, read_id, fast5_id):
        pass

    def close(self):
        self.out.close()
