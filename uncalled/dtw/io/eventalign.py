import numpy as np
import pandas as pd
import sys
from ..aln_track import AlnTrack
from . import TrackIO
import _uncalled

class Eventalign(TrackIO):
    FORMAT = "eventalign"

    def __init__(self, conf, model, mode):
        filename = conf.tracks.io.eventalign_in if mode == "r" else conf.tracks.io.eventalign_out
        TrackIO.__init__(self, filename, conf, model, mode)

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
            header.append("read_idx")

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
        t = AlnTrack(self, None, name, name, self.conf, self.model)
        self.tracks.append(t)

    #def write_dtw_events(self, track, events):
    def write_layers(self, df, read, index=["pac","aln_id"]):
        for group in df.columns.levels[0]:
            if group == "dtw":
                break
            return

        track = self.tracks[0]
        df = df.sort_index()

        pacs = df.index.get_level_values(0)
        events = df[group].set_index(track.coords.pac_to_ref(pacs))

        model = track.model
        kmers = events["kmer"]

        std_level = (events["current"] - model.model_mean) / model.model_stdv

        evts = events.rename(columns={"current" : "mean", "current_stdv" : "stdv"})
        read.set_events(evts)
        #read = ProcessedRead(evts)

        if self.write_samples:
            signal = read.get_norm_signal()
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
            np.arange(len(events))[::-1]+1, #TODO properly rep skips?
            std_level, signal) #TODO compute internally?

        self.out.write(eventalign)

    def write_layers_old(self, df, index=["pac","aln_id"]):
        for group in df.columns.levels[0]:
            if group == "dtw":
                break
            return

        track = self.tracks[0]
        df = df.sort_index()

        pacs = df.index.get_level_values(0)
        events = df[group].set_index(track.coords.pac_to_ref(pacs))

        contig = track.coords.ref_name

        model = track.model
        kmers = events["kmer"]

        std_level = (events["current"] - model.model_mean) / model.model_stdv

        if self.conf.is_rna:
            model_kmers = model.kmer_rev(kmers)
        else:
            model_kmers = kmers

        if track.coords.fwd:
            ref_kmers = model_kmers
        else:
            ref_kmers = model.kmer_comp(model_kmers)

        #read_id = track.alignments.iloc[0]["read_id"]
        #sys.stderr.write(f"{self.prev_aln_id}\t{read_id}\n")

        #https://github.com/jts/nanopolish/issues/655

        stdvs = events["current_stdv"] if "current_stdv" in events else pd.NA

        #i32 read_idx, bool fwd,
        #ProcessedRead read, 
        #py::array_t<i64> pac_np,
        #py::array_t<typename ModelType::kmer_t> model_kmer_np,
        #py::array_t<i32> event_index_np,
        #py::array_t<float> std_level_np) {


        eventalign = pd.DataFrame(
            data = {
                "contig" : track.coords.ref_name,
                "position" : events.index-2,
                "reference_kmer" : model.kmer_to_str(ref_kmers),
                "read_index" : self.prev_aln_id,
                "strand" : "t",
                "event_index" : pd.RangeIndex(0,len(events))[::-1]+1,
                "event_level_mean" : events["current"],
                "event_stdv" : stdvs,
                "event_length" : events["length"] / track.conf.read_buffer.sample_rate,
                "model_kmer" : model.kmer_to_str(model_kmers),
                "model_mean" : model.means[kmers],
                "model_stdv" : model.stdvs[kmers],
                "standardized_level" : std_level,
                "start_idx" : events["start"],
                "end_idx" : events["start"] + events["length"],
            }, index = events.index)#.sort_values("position")

        #print(eventalign[:20])
        #sys.stdout.flush()


        eventalign.to_csv(
            self.out, sep="\t",
            header=self._header,
            float_format="%.5f",
            index=False)

        self._header = False

    def write_alignment(self, alns):
        pass

    def init_fast5(self, fast5):
        pass

    def init_read(self, read_id, fast5_id):
        pass

    def close(self):
        self.out.close()
