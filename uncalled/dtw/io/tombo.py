import numpy as np
import pandas as pd
import sys, os
from ..aln_track import AlnTrack
from ...signal_processor import ProcessedRead
from ...index import str_to_coord, RefCoord
from ...pore_model import PoreModel
from ...fast5 import Fast5Reader
from ont_fast5_api.fast5_interface import get_fast5_file
from . import TrackIO
import _uncalled

class Tombo(TrackIO):
    FORMAT = "tombo"

    def __init__(self, filename, write, tracks):
        TrackIO.__init__(self, filename, write, tracks)

        self._header = True

        self.fast5_in = None
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

        self.conf.fast5_reader.fast5_files = self.prms.tombo_in

        self.init_track(1, name, name, self.conf)

    #def init_fast5(self, fast5):
    #    self.
    #    pass

    def iter_alns(self, layers, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, full_overlap=None, ref_index=None):

        
        f5reader = Fast5Reader(conf=self.conf)
        read_filter = f5reader.get_read_filter()
        fast5_files = f5reader.prms.fast5_files

        if self.conf.fast5_reader.max_reads > 0 and len(fast5_files) > self.conf.fast5_reader.max_reads:
            fast5_files = fast5_files[:self.conf.fast5_reader.max_reads]

        aln_id = 1

        for fast5_fname in fast5_files:
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
            self.read_id_in = read.read_id

            aln_attrs = dict(handle["Alignment"].attrs)

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

            ref_bounds = RefCoord(aln_attrs["mapped_chrom"],start, end,fwd)

            sig_fwd = ref_bounds.fwd != is_rna

            coords = ref_index.get_coord_space(ref_bounds, is_rna=is_rna, load_kmers=True, kmer_trim=True)
            #aln_id,_ = tracks.init_alignment(read.read_id, fast5_fname, coords)

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
            refs = coords.refs[2:-2]
            pacs = coords.pacs[2:-2] #ref_to_pac(df["position"].to_numpy()+kmer_trim[0])

            kmers = coords.ref_kmers.droplevel(0).loc[refs]

            idx = pd.MultiIndex.from_product(
                    [[fwd], pacs, [aln_id]], names=("fwd","pac","aln_id"))

            layers = pd.DataFrame({
                    "start"  : starts,
                    "length" : lengths,
                    "current"   : currents,
                    "kmer" : kmers.to_numpy()
                 }, index=idx)

            samp_start = layers["start"].min()
            e = layers["start"].argmax()
            samp_end = layers["start"].iloc[e] + layers["length"].iloc[e]

            layers = pd.concat({"dtw" : layers}, names=("group","layer"), axis=1).sort_index()

            alns = pd.DataFrame({
                    "id" : [aln_id],
                    "track_id" : [1],
                    "read_id" : [read.read_id],
                    "ref_name" : [coords.ref_name],
                    "ref_start" : [coords.refs.start],
                    "ref_end" : [coords.refs.stop],
                    "fwd" :     [coords.fwd],
                    "samp_start" : [samp_start],
                    "samp_end" : [samp_end],
                    "tags" : [""]}).set_index("id")

            aln_id += 1

            yield alns,layers
            
    def write_alignment(self, alns):
        pass


    def init_fast5(self, fast5):
        pass

    def init_read(self, read_id, fast5_id):
        pass

    def close(self):
        pass
