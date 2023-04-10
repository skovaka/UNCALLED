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
        
        self.conf.pore_model.name = os.path.join("tombo", self.conf.pore_model.name)

        self.conf.read_index.paths = self.prms.tombo_in

        self.track_in = self.init_track(name, name, self.conf)

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

            if sig_fwd:
                mpos = pd.RangeIndex(start+2, end-2)
                step = 1
            else:
                mpos = pd.RangeIndex(-end+2, -start-2)
                step = -1

            tombo_events = np.array(handle["Events"])[clip:]

            tombo_start = handle["Events"].attrs["read_start_rel_to_raw"]
            
            raw_len = len(read.get_raw_data())
            starts = tombo_events["start"]


            lengths = tombo_events["length"]
            currents = self.track_in.model.pa_to_norm(tombo_events["norm_mean"])

            if is_rna:
                starts = raw_len - tombo_start - starts - tombo_events["length"]

            df = pd.DataFrame({
                    "mpos" : mpos,#[2:-2]
                    "start"  : starts[::step],
                    "length" : lengths[::step],
                    "current"   : currents[::step],
                    #"kmer" : kmers.to_numpy()
                 })#, index=idx)

            coords = RefCoord(read_sam.reference_name, start, end, fwd)

            #lengths.fillna(-1, inplace=True)
            aln = self.tracks.init_alignment(self.track_in.name, self.next_aln_id(), read_id, read_sam.reference_id, coords, read_sam)

            dtw = AlnDF(aln.seq, np.array(starts[::step]), np.array(lengths[::step]), np.array(currents[::step])) #df["stdv"])

            aln.set_dtw(dtw)

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
