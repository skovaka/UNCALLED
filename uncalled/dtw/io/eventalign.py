import numpy as np
import pandas as pd
import sys
from ..aln_track import AlnTrack
from ...signal_processor import ProcessedRead
from ...index import str_to_coord, RefCoord
from ...pore_model import PoreModel
from . import TrackIO
import _uncalled

class Eventalign(TrackIO):
    FORMAT = "eventalign"

    def __init__(self, filename, write, tracks):
        TrackIO.__init__(self, filename, write, tracks)

        if self.filename == "-":
            self.out = sys.stdout
        else:
            self.out = open(self.filename, "w" if write else "r")

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
        #if len(self.track_names) != 1:
        #    raise ValueError("Can only read eventalign TSV into a single track")
        #name = self.track_names[0]
        name = self.filename
        
        if self.conf.pore_model.name == "r94_dna":
            self.conf.pore_model.name = "r9.4_dna_450bps_6mer_npl"

        self.init_track(1, name, name, self.conf)

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

        if "events" in events:
            event_counts = events["events"]
            event_index = (event_counts.cumsum() - event_counts.iloc[0]).astype(int)
            if track.all_rev: #TODO check for flipped ref
                event_index = event_index.max() - event_index
            event_index += 1
        else:
            event_index = np.arange(len(events))

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
        
        writer = self.writer = getattr(_uncalled, f"write_eventalign_K{model.K}")

        eventalign = writer(
            self.conf, model.instance, read_id, track.coords.fwd, read,
            track.coords.ref_name, events.index-2, 
            self.write_signal_index,
            kmers, 
            event_index, #TODO properly rep skips?
            std_level, signal) #TODO compute internally?

        self.out.write(eventalign)

    #def init_fast5(self, filename):

    #    row = self.cur.execute("SELECT id FROM fast5 WHERE filename = ?", (filename,)).fetchone()
    #    if row is not None:
    #        return row[0]
    #        
    #    self.cur.execute("INSERT INTO fast5 (filename) VALUES (?)", (filename,))
    #    fast5_id = self.cur.lastrowid
    #    #self.con.commit()

    #    return fast5_id

    #def init_read(self, read_id, fast5_id):
    #    self.cur.execute("INSERT OR IGNORE INTO read VALUES (?,?)", (read_id, fast5_id))

    def iter_alns(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, full_overlap=None, ref_index=None):

        read_filter = set(self.conf.fast5_reader.read_filter)

        sample_rate = self.conf.read_buffer.sample_rate

        csv_iter = pd.read_csv(
            self.filename, sep="\t", chunksize=10000,
            usecols=["read_name","contig","position", "event_index",
                     "start_idx","event_level_mean","model_mean",
                     "event_length","strand","model_kmer"])

        model = PoreModel(self.conf.pore_model)

        kmer_trim = ref_index.trim

        aln_id = 1

        def iter_layers(events, aln_id):
            groups = events.groupby(["contig", "read_name"])
            for (contig,read_id), df in groups:

                if len(read_filter) > 0 and read_id not in read_filter:
                    continue

                df.drop(df.index[df["model_mean"] == 0], inplace=True)

                start = df["position"].min()
                end = df["position"].max()+1
                
                fwd = int( (df["event_index"].iloc[0] < df["event_index"].iloc[-1]) == (df["position"].iloc[0] < df["position"].iloc[-1]))

                kmers = model.str_to_kmer(df["model_kmer"])
                if self.conf.is_rna:
                    kmers = model.kmer_rev(kmers)
                else:
                    fwd = int(df["event_index"].iloc[0] < df["event_index"].iloc[-1])

                ref_coord = RefCoord(contig, start, end, fwd)
                coords = ref_index.get_coord_space(ref_coord, self.conf.is_rna, kmer_trim=True)


                df["kmer"] = kmers

                pacs = coords.ref_to_pac(df["position"].to_numpy()+kmer_trim[0])

                layers = df.rename(columns={
                        "start_idx" : "start",
                        "event_length" : "length",
                        "event_level_mean" : "current",
                     }).set_index(pd.MultiIndex.from_product(
                        [[fwd], pacs, [aln_id]], names=("fwd","pac","aln_id")))

                samp_start = layers["start"].min()
                e = layers["start"].argmax()
                samp_end = layers["start"].iloc[e] + layers["length"].iloc[e]

                layers = pd.concat({"dtw" : layers}, names=("group","layer"), axis=1).sort_index()

                alns = pd.DataFrame({
                        "id" : [aln_id],
                        "track_id" : [1],
                        "read_id" : [read_id],
                        "ref_name" : [coords.ref_name],
                        "ref_start" : [coords.refs.start],
                        "ref_end" : [coords.refs.stop],
                        "fwd" :     [coords.fwd],
                        "samp_start" : [samp_start],
                        "samp_end" : [samp_end],
                        "tags" : [""]}).set_index("id")

                aln_id += 1

                yield alns,layers

        leftover = pd.DataFrame()

        for events in csv_iter:
            events['event_length'] = np.round(events['event_length'] * sample_rate).astype(int)
            events['sum'] = events['event_level_mean'] * events['event_length']

            events = pd.concat([leftover, events])

            i = events["read_name"] == events["read_name"].iloc[-1]
            leftover = events[i]
            events = events[~i]

            for alns,layers in iter_layers(events, aln_id):
                yield alns,layers

        if len(leftover) > 0:
            for alns,layers in iter_layers(leftover, aln_id):
                yield alns,layers


    def write_alignment(self, alns):
        pass

    def init_fast5(self, fast5):
        pass

    def init_read(self, read_id, fast5_id):
        pass

    def close(self):
        self.out.close()
