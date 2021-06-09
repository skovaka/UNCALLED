"""Converts tombo or nanopolish alignments to uncalled DTW track format"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd
from ont_fast5_api.fast5_interface import get_fast5_file

from .dtw import ReadAln, Track, ref_coords
from ..config import Config, ArgParser, ParamGroup, Opt
from ..fast5 import Fast5Reader, FAST5_OPTS
from ..index import BWA_OPTS
from _uncalled import nt, PORE_MODELS

import progressbar as progbar

CONVERT_OPTS = BWA_OPTS + FAST5_OPTS + (
    Opt(("-m", "--mm2-paf"), "align", required=True),
    Opt("--rna", fn="set_r94_rna"),
    Opt(("-R", "--ref-bounds"), "align", type=ref_coords),
    Opt(("-f", "--force-overwrite"), action="store_true"),
    Opt(("-o", "--out-path"), "align", required=True),
)


#Maps is_rna bool to model
MODELS = [
    PORE_MODELS["r94_dna_templ"], 
    PORE_MODELS["r94_rna_templ"]
]

#TODO make argparser accept functions as subcommands, rather than class wtih main?
class nanopolish:
    """Convert from nanopolish eventalign TSV to uncalled DTW track"""
    OPTS = CONVERT_OPTS + (Opt("eventalign_tsv", type=str, default=None),)

    @staticmethod
    def main(conf):

        f5reader = Fast5Reader(conf=conf)
        conf.fast5_reader.load_bc = True

        sample_rate = conf.read_buffer.sample_rate
        print(sample_rate)

        sys.stderr.write("Parsing TSV\n")
        alns = pd.read_csv(
                conf.eventalign_tsv, sep="\t", 
                usecols=["read_name", "contig","position","start_idx","event_level_mean","event_length","reference_kmer","model_kmer"]
                )#.rename(columns={"start_idx" : "start", "event_length" : "length"})

        alns['event_length'] = np.round(alns['event_length'] * sample_rate).astype(int)
        alns['sum'] = alns['event_level_mean'] * alns['event_length']

        sys.stderr.write("Grouping\n")

        read_groups = alns.groupby(["contig", "read_name"])

        sys.stderr.write("Writing alignments\n")

        track = Track(conf.align.out_path, "w", conf, conf.force_overwrite)

        for (contig, read_id), rows in read_groups.groups.items():
            mm2 =  track.mm2s[read_id]
            if contig != mm2.rf_name: continue

            aln = ReadAln(track.index, mm2, is_rna=conf.is_rna)
            df = alns.iloc[rows]

            aln.set_subevent_aln(df, kmer_str=True, ref_col="position", start_col="start_idx", length_col="event_length", mean_col="event_level_mean", kmer_col="reference_kmer")

            fast5_name = os.path.basename(f5reader.get_read_file(read_id))

            track.save_aln(aln, fast5_name)

class tombo:

    OPTS = CONVERT_OPTS 

    @staticmethod
    def main(conf):
        """Convert from Tombo resquiggled fast5s to uncalled DTW track"""
        f5reader = Fast5Reader(conf=conf)
        conf.fast5_reader.load_bc = True

        track = Track(conf.align.out_path, "w", conf, conf.force_overwrite)
        
        fast5_files = f5reader.prms.fast5_files
        
        pbar = progbar.ProgressBar(
                widgets=[progbar.Percentage(), progbar.Bar(), progbar.ETA()], 
                maxval=len(fast5_files)).start()
        pbar_count = 0

        for fast5_name in fast5_files:

            fast5_basename = os.path.basename(fast5_name)

            with get_fast5_file(fast5_name, mode="r") as fast5:
                pbar_count += 1

                read, = fast5.get_reads()

                if read.read_id not in track.mm2s: continue

                if not 'BaseCalled_template' in read.handle['Analyses']['RawGenomeCorrected_000']:
                    #TODO debug logs
                    continue

                handle = read.handle['Analyses']['RawGenomeCorrected_000']['BaseCalled_template']
                attrs = handle.attrs
                if attrs["status"] != "success":
                    #TODO debug logs
                    continue

                is_rna = handle.attrs["rna"]
                if is_rna != conf.is_rna:
                    raise RuntimeError("Reads appear to be RNA but --rna not specified")

                aln_attrs = dict(handle["Alignment"].attrs)
                ch = aln_attrs["mapped_chrom"]
                st = aln_attrs["mapped_start"]
                en = aln_attrs["mapped_end"]

                fwd = aln_attrs["mapped_strand"] == "+"
                sig_fwd = (fwd != is_rna)

                tombo_events = pd.DataFrame(np.array(handle["Events"]))
                tombo_start = handle["Events"].attrs["read_start_rel_to_raw"]
                
                raw_len = len(read.get_raw_data())
                samps = tombo_events["start"]
                refs = np.arange(st, en)

                K = 5
                shift = 1
                bases = tombo_events["base"].to_numpy().astype(str)
                kmers = [nt.str_to_kmer("".join(bases[i:i+K]), 0) for i in range(len(bases)-K+1)]
                kmers = shift*[kmers[0]] + kmers + (K-shift-1)*[kmers[-1]]
                #TODO load kmers from genome
                
                signal = tombo_events["norm_mean"]
                scale = MODELS[is_rna].get_means_stdv() / np.std(signal)
                shift = MODELS[is_rna].get_means_mean() - scale * np.mean(signal)
                signal = scale * signal + shift

                if not sig_fwd:
                    samps = raw_len - tombo_start - samps - tombo_events["length"]
                    kmers = nt.kmer_rev(kmers)

                mm2 =  track.mm2s[read.read_id]
                aln = ReadAln(track.index, mm2, is_rna=is_rna)

                aln.df = pd.DataFrame({
                    "ref"    : refs,
                    "start"  : samps,
                    "kmer"   : kmers,
                    "length" : tombo_events["length"],
                    "mean"   : signal
                }).set_index("ref")

                track.save_aln(aln, fast5_basename)

                pbar.update(pbar_count)

        track.close()

        pbar.finish()

SUBCMDS = [nanopolish, tombo]
