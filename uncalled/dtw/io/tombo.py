import numpy as np
import pandas as pd
import sys, os
from ..aln_track import AlnTrack
from ...signal_processor import ProcessedRead
from ...index import str_to_coord, RefCoord
from ...pore_model import PoreModel
from ont_fast5_api.fast5_interface import get_fast5_file
from . import TrackIO
import _uncalled
from ...aln import AlnDF

class Tombo(TrackIO):
    FORMAT = "tombo"

    def __init__(self, filename, write, tracks, track_count):
        TrackIO.__init__(self, filename, write, tracks, track_count)

        self._header = True

        self.fast5_in = None
        self.read_index = tracks.read_index
        self.read_id_in = None

        if self.write_mode:
            self.init_write_mode()
        else:
            self.init_read_mode()

    def init_write_mode(self):
        raise RuntimeError("Writing in tombo format not supported")

    def init_read_mode(self):
        name = self.filename
        
        if self.conf.pore_model.name == "r94_rna":
            self.conf.pore_model.name = "r94_rna_tombo"

        self.conf.read_index.paths = self.prms.tombo_in

        self.init_track(name, name, self.conf)

    #def init_fast5(self, fast5):
    #    self.
    #    pass

    def iter_alns(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, full_overlap=None, ref_index=None):

        
        read_filter = self.read_index.read_filter #f5reader.get_read_filter()
        paths = self.read_index.file_info.values() #f5reader.prms.paths

        if self.conf.read_index.read_count is not None and len(paths) > self.conf.read_index.read_count:
            paths = paths[:self.conf.read_index.read_count]

        aln_id = 1

        for _,fast5_fname in paths:
            fast5_basename = os.path.basename(fast5_fname)

            try:
                fast5 = get_fast5_file(fast5_fname, mode="r")
            except:
                sys.stderr.write(f"Unable to open \"{fast5_fname}\". Skipping.\n")
                continue

            read, = fast5.get_reads()

            if not (read_filter is None or read.read_id in read_filter): 
                continue

            if not 'BaseCalled_template' in read.handle['Analyses']['RawGenomeCorrected_000']:
                #TODO debug logs
                continue

            handle = read.handle['Analyses']['RawGenomeCorrected_000']['BaseCalled_template']
            attrs = handle.attrs
            if attrs["status"] != "success":
                #TODO debug logs
                continue

            is_rna = handle.attrs["rna"]
            if is_rna != self.conf.is_rna:
                raise RuntimeError("Reads appear to be RNA but --rna not specified")

            self.fast5_in = fast5_fname
            read_id = read.read_id

            aln_attrs = dict(handle["Alignment"].attrs)

            chrom = aln_attrs["mapped_chrom"]

            fwd = aln_attrs["mapped_strand"] == "+"
            if fwd:
                start = aln_attrs["mapped_start"]-1
                end = aln_attrs["mapped_end"]+3
            else:
                start = aln_attrs["mapped_start"]-3
                end = aln_attrs["mapped_end"]+1

            if start < 0:
                clip = -start
                start = 0
            else:
                clip = 0

            sam = None
            for read_sam in self.tracks.bam_in.iter_read(read_id):
                if read_sam.reference_name == chrom and start >= read_sam.reference_start and end < read_sam.reference_end:
                    sam = read_sam
                    break

            #ref_bounds = RefCoord(aln_attrs["mapped_chrom"],start, end,fwd)

            sig_fwd = fwd != is_rna

            #coords = ref_index.get_coord_space(ref_bounds, is_rna=is_rna, load_kmers=True, kmer_trim=True)
            #aln_id,_ = tracks.init_alignment(read.read_id, fast5_fname, coords)
            if sig_fwd:
                coords = _uncalled.IntervalIndexI64([(start+2, end-2)])
                step = 1
            else:
                coords = _uncalled.IntervalIndexI64([(-end+2, -start-2)])
                step = -1

            tombo_events = np.array(handle["Events"])[clip:]

            #if not ref_bounds.fwd:
            #    tombo_events = tombo_events[::-1]
                
            tombo_start = handle["Events"].attrs["read_start_rel_to_raw"]
            
            raw_len = len(read.get_raw_data())
            starts = tombo_events["start"]


            lengths = tombo_events["length"]
            currents = tombo_events["norm_mean"]

            if is_rna:
                starts = raw_len - tombo_start - starts - tombo_events["length"]

            #if is_rna == ref_bounds.fwd:
            #refs = coords.refs[2:-2]
            #pacs = coords.pacs[2:-2] #pos_to_pac(df["position"].to_numpy()+kmer_trim[0])
            #kmers = coords.ref_kmers.droplevel(0).loc[refs]

            #idx = pd.MultiIndex.from_product(
            #        [[self.in_id], [fwd], pacs, [aln_id]], names=("track.id","seq.fwd","seq.pac","aln.id"))

            df = pd.DataFrame({
                    "mpos" : np.array(coords.expand()),#[2:-2]
                    "start"  : starts[::step],
                    "length" : lengths[::step],
                    "current"   : currents[::step],
                    #"kmer" : kmers.to_numpy()
                 })#, index=idx)

            df["length"].fillna(-1, inplace=True)
            aln = self.tracks.init_alignment(self.next_aln_id(), read_id, read_sam.reference_id, coords, read_sam)

            #print(len(aln.seq), len(df["start"]), len(df["length"]), len(df["current"]))
            #dtw = AlnDF(aln.seq, df["start"], df["length"], df["current"], None) #df["stdv"])
            dtw = AlnDF(aln.seq, starts[::step], lengths[::step], currents[::step]) #df["stdv"])
            aln.set_dtw(dtw)

            #layers = pd.concat({"dtw" : layers}, names=("group","layer"), axis=1).sort_index()

            #alns = pd.DataFrame({
            #        "id" : [aln_id],
            #        "track.id" : [self.in_id],
            #        "read_id" : [read.read_id],
            #        "ref_name" : [coords.ref_name],
            #        "ref_start" : [coords.refs.start],
            #        "ref_end" : [coords.refs.stop],
            #        "seq.fwd" :     [coords.fwd],
            #        "samp_start" : [samp_start],
            #        "samp_end" : [samp_end],
            #        "tags" : [""]}).set_index(["track.id", "id"])

            aln_id += 1

            yield aln#s,layers
            
    def write_alignment(self, alns):
        pass


    def init_fast5(self, fast5):
        pass

    def init_read(self, read_id, fast5_id):
        pass

    def close(self):
        pass
