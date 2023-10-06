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
import os
import json

class M6anet(TrackIO):
    FORMAT = "eventalign"

    def __init__(self, filename, write, tracks, track_count):
        TrackIO.__init__(self, filename, write, tracks, track_count)

        self._header = True

        if self.write_mode:
            self.init_write_mode()
        else:
            self.init_read_mode()

    @property
    def info_file(self):
        return os.path.join(self.filename, "data.info")

    @property
    def json_file(self):
        return os.path.join(self.filename, "data.json")

    def init_write_mode(self):
        TrackIO.init_write_mode(self)

        try:
            os.mkdir(self.filename)
        except FileExistsError:
            pass

        if not self.prms.buffered:
            self.json_out = open(self.json_file, "w")
            self.info_out = open(self.info_file, "w")
            self.info_out.write("transcript_id,transcript_position,start,end,n_reads\n")
        else:
            self.out_buffer = list()
        self.offset = 0 

        self.init_stats()

        self.model = self.track_out.model

    def init_stats(self):
        self.transcript_id = None
        self.prev_pos = None
        self.ref_stats = dict()
        self.ref_kmers = dict()

    def write_transcript(self):
        for pos in sorted(self.ref_stats.keys()):
            data = self.ref_stats[pos]
            vals, = data.values()
            n_reads = len(vals)
            if n_reads < 20: continue
            d = {self.transcript_id : {str(pos) : data}}
            json_line = json.dumps(d,separators=(',',':')) + "\n"
            end = self.offset + len(json_line)

            if self.prms.buffered:
                self.out_buffer.append((self.transcript_id, pos, n_reads, json_line))
            else:
                self.json_out.write(json_line)
                self.info_out.write(f"{self.transcript_id},{pos},{self.offset},{end},{n_reads}\n")

            self.offset = end
        self.init_stats()

    def write_buffer(self, lines):
        for tid,pos,n_reads,json_line in lines:
            self.json_out.write(json_line)
            start = self.offset
            end = self.offset + len(json_line)
            self.info_out.write(f"{tid},{pos},{self.offset},{end},{n_reads}\n")
            self.offset = end

    def write_alignment(self, aln):
        if self.transcript_id != aln.seq.name and self.transcript_id != None:
            self.write_transcript()

        drach_mask = self.model.pattern_mask("HCARD", aln.seq.kmer)
        drach_idx = np.flatnonzero(drach_mask)

        na_mask = aln.dtw.na_mask

        coords = aln.seq.pos

        for i in drach_idx:
            if i == 0 or i == len(aln.seq)-1 or not np.all(na_mask[i-1:i+2]): continue
            pos = coords[i]


            if pos in self.ref_stats:
                (kmer,pos_stats), = self.ref_stats[pos].items()
            else:
                lbase = self.model.BASES[self.model.kmer_base(aln.seq.kmer[i-1],0)]
                rbase = self.model.BASES[self.model.kmer_base(aln.seq.kmer[i+1],self.model.K-1)]
                kmer = (lbase + self.model.kmer_to_str(aln.seq.kmer[i]) + rbase)[::-1]
                pos_stats = list()
                self.ref_stats[pos] = {kmer : pos_stats}

            data = np.zeros(10, dtype=float)
            j = 0
            for k in reversed(range(i-1,i+2)):
                dwell = aln.dtw.samples.lengths[k] / self.model.sample_rate
                data[j] = dwell
                data[j+1] = self.model.norm_to_pa_sd(aln.dtw.current_sd[k])
                data[j+2] = self.model.norm_to_pa(aln.dtw.current[k])
                #data[j] = round(dwell,8)
                #data[j+1] = round(self.model.norm_to_pa_sd(aln.dtw.current_sd[k]),5)
                #data[j+2] = round(self.model.norm_to_pa(aln.dtw.current[k]),1)
                j += 3
            data[j] = aln.id

            pos_stats.append(data.tolist())
        
        self.transcript_id = aln.seq.name
        #print(drach_mask)


    def close(self):
        if len(self.ref_stats) > 0:
            self.write_transcript()

        if not self.prms.buffered:
            self.json_out.close()
            self.info_out.close()
