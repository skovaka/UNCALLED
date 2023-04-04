import numpy as np
import pandas as pd
import sys
from ..aln_track import AlnTrack
from ...aln import AlnDF
from ...signal_processor import ProcessedRead
from ...index import str_to_coord, RefCoord
from ...pore_model import PoreModel
from . import TrackIO
import _uncalled

class Eventalign(TrackIO):
    FORMAT = "eventalign"

    def __init__(self, filename, write, tracks, track_count):
        TrackIO.__init__(self, filename, write, tracks, track_count)

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

        self.header = ["contig", "position", "reference_kmer"]

        if self.write_read_name:
            self.header.append("read_name")
        else:
            self.header.append("read_index")

        self.header += ["strand", "event_index", "event_level_mean", "event_stdv", "event_length", "model_kmer", "model_mean", "model_stdv", "standardized_level"]

        if self.write_signal_index:
            self.header += ["start_idx", "end_idx"]

        if self.write_samples:
            self.header += ["samples"]

        TrackIO._init_output(self, self.prms.buffered)

        if not self.prms.buffered:
            self.output.write("\t".join(self.header) + "\n")

    def init_read_mode(self):
        name = self.filename

        self.output = None
        
        if self.conf.pore_model.name == "r94_dna":
            self.conf.pore_model.name = "r9.4_dna_450bps_6mer_npl"

        self.model = None

        self.init_track(name, name, self.conf)

    #def write_layers(self, track, groups):
    def write_alignment(self, aln):
        events = aln.to_pandas(["seq.kmer", "dtw"], ["seq.pos"]).sort_index().droplevel(0, axis=1)

        model = self.tracks.model
        kmers = events["kmer"]

        std_level = (events["current"] - model.model_mean) / model.model_stdv

        if "events" in events:
            event_counts = events["events"]
            event_index = (event_counts.cumsum() - event_counts.iloc[0]).astype(int)
            if not aln.fwd: #TODO check for flipped ref
                event_index = event_index.max() - event_index
            event_index += 1
        else:
            event_index = np.arange(len(events))

        evts = events.rename(columns={"current" : "mean", "current_sd" : "stdv"})
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
            read_id = aln.read_id #track.alignments["read_id"].iloc[0]
        else:
            read_id = str(aln.id)
        self.next_aln_id()
        
        #self.writer = getattr(_uncalled, f"write_eventalign_K{model.K}")
        if isinstance(model.instance, _uncalled.PoreModelU16):
            self.writer = _uncalled.write_eventalign_U16
        elif isinstance(model.instance, _uncalled.PoreModelU32):
            self.writer = _uncalled.write_eventalign_U32
        else:
            raise ValueError(f"Unknown PoreModel type: {model.instance}")

        eventalign = self.writer(
            self.conf, model.instance, read_id, aln.seq.fwd, read,
            aln.seq.name, events.index-2, 
            self.write_signal_index, kmers, 
            event_index, #TODO properly rep skips?
            std_level, signal) #TODO compute internally?

        self._set_output(eventalign)

    def iter_alns(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, full_overlap=None, ref_index=None):

        print("HERE")
        #read_filter = set(self.conf.read_index.read_filter)

        csv_iter = pd.read_csv(
            self.filename, sep="\t", chunksize=10000,
            usecols=["read_name","contig","position", "event_index",
                     "start_idx","event_level_mean","event_stdv","model_mean",
                     "event_length","strand","model_kmer"])
        print("HEREereerr")

        if self.model is None:
            self.model = self.tracks.model

        kmer_trim = ref_index.trim

        aln_id = 1

        print("a")
        def iter_layers(events, aln_id):
            groups = events.groupby(["contig", "read_name"])
            for (contig,read_id), df in groups:

                if self.read_filter is not None and read_id not in self.read_filter:

                    continue

                #df.drop(df.index[df["model_mean"] == 0], inplace=True)

                start = df["position"].min()
                end = df["position"].max()+1
                
                sam = None
                for read_sam in self.tracks.bam_in.iter_read(read_id):
                    if read_sam.reference_name == contig and start >= read_sam.reference_start and end < read_sam.reference_end:
                        sam = read_sam
                        break

                print("b")
                sys.stdout.flush();
                aln = self.tracks.bam_in.sam_to_aln(sam, load_moves=False)
                #aln = self.tracks.init_alignment(self.next_aln_id(), read_id, sam.reference_id, coords, sam)
                print("c")
                sys.stdout.flush();
                
                fwd = int( (df["event_index"].iloc[0] < df["event_index"].iloc[-1]) == (df["position"].iloc[0] < df["position"].iloc[-1]))

                pos = df["position"].to_numpy()
                df["mpos"] = self.tracks.index.pos_to_mpos(pos, fwd, self.conf.is_rna)-kmer_trim[0]
                df["mean_cml"]  = df["event_length"] * df["event_level_mean"]
                df["stdv_cml"]  = df["event_length"] * df["event_stdv"]

                grp = df.groupby("mpos")

                lengths = grp["event_length"].sum()
                df = pd.DataFrame({
                    "start" : grp["start_idx"].min(),
                    "length" : lengths,
                    "mean" : self.model.pa_to_norm(grp["mean_cml"].sum() / lengths),
                    "stdv" : self.model.pa_sd_to_norm(grp["stdv_cml"].sum() / lengths)
                }).set_index(grp["mpos"].min())
                print("d")

                #coords = _uncalled.IntervalIndexI64([(df.index.min()-kmer_trim[0], df.index.max()+1+kmer_trim[1])])
                coords = _uncalled.IntervalIndexI64([(df.index.min(), df.index.max()+1)])
                df = df.reindex(coords.expand())
                df["length"].fillna(-1, inplace=True)
                aln = self.tracks.init_alignment(self.next_aln_id(), read_id, sam.reference_id, coords, sam)
                dtw = AlnDF(aln.seq, df["start"], df["length"], df["mean"], df["stdv"])
                aln.set_dtw(dtw)
                yield aln

        leftover = pd.DataFrame()

        for events in csv_iter:
            events['event_length'] = np.round(events['event_length'] * self.model.sample_rate).astype(int)
            events['sum'] = events['event_level_mean'] * events['event_length']

            events = pd.concat([leftover, events])

            i = events["read_name"] == events["read_name"].iloc[-1]
            leftover = events[i]
            events = events[~i]

            for aln in iter_layers(events, aln_id):
                yield aln#s,layers

        if len(leftover) > 0:
            for alns,layers in iter_layers(leftover, aln_id):
                yield aln#s,layers


    def init_fast5(self, fast5):
        pass

    def init_read(self, read_id, fast5_id):
        pass

    def close(self):
        if not self.prms.buffered and self.output is not None:
            self.output.close()
