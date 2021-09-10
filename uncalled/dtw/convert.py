"""Import DTW alignments produced by other tools"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd
from ont_fast5_api.fast5_interface import get_fast5_file
import scipy.stats

from . import ReadAln, AlnTrack, RefCoord
from ..config import Config, ParamGroup
from ..argparse import ArgParser, Opt
from ..fast5 import Fast5Reader, FAST5_OPTS
from ..index import BWA_OPTS
from .. import nt, PoreModel

import progressbar as progbar

CONVERT_OPTS = (Opt("index_prefix", "track"),) + FAST5_OPTS + (
    Opt(("-m", "--mm2-paf"), "dtw", required=True),
    Opt("--rna", fn="set_r94_rna"),
    Opt(("-R", "--ref-bounds"), "track", type=RefCoord),
    Opt(("-f", "--overwrite"), "track", action="store_true"),
    Opt(("-o", "--out-path"), "track", "path", required=True),
)

NANOPOLISH_OPTS = CONVERT_OPTS + (Opt("eventalign_tsv", type=str, default=None),)
def nanopolish(conf):
    """Convert from nanopolish eventalign TSV to uncalled DTW track"""

    f5reader = Fast5Reader(conf=conf)
    conf.fast5_reader.load_bc = True

    sample_rate = conf.read_buffer.sample_rate

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

    track = AlnTrack(mode="w", conf=conf)

    for (contig, read_id), rows in read_groups.groups.items():
        #mm2 =  track.mm2s[read_id]
        #if contig != mm2.rf_name: continue

        aln = ReadAln(track.index, mm2, is_rna=conf.is_rna)
        df = alns.iloc[rows]

        aln.set_subevent_aln(df, kmer_str=True, ref_col="position", start_col="start_idx", length_col="event_length", mean_col="event_level_mean", kmer_col="reference_kmer")

        fast5_name = os.path.basename(f5reader.get_read_file(read_id))

        track.save_aln(aln, fast5_name)

TOMBO_OPTS = CONVERT_OPTS 
def tombo(conf):
    """Convert from Tombo resquiggled fast5s to uncalled DTW track"""
    f5reader = Fast5Reader(conf=conf)
    conf.fast5_reader.load_bc = True
    conf.pore_model.name = "r94_rna_tombo"

    track = AlnTrack(mode="w", conf=conf)
    
    fast5_files = f5reader.prms.fast5_files

    if conf.fast5_reader.max_reads > 0 and len(fast5_files) > conf.fast5_reader.max_reads:
        fast5_files = fast5_files[:conf.fast5_reader.max_reads]
    
    pbar = progbar.ProgressBar(
            widgets=[progbar.Percentage(), progbar.Bar(), progbar.ETA()], 
            maxval=len(fast5_files)).start()
    pbar_count = 0

    model = PoreModel(conf.pore_model)

    for fast5_name in fast5_files:

        fast5_basename = os.path.basename(fast5_name)

        with get_fast5_file(fast5_name, mode="r") as fast5:
            pbar_count += 1

            read, = fast5.get_reads()

            #if read.read_id not in track.mm2s: continue

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

            #aln = ReadAln(track.index, mm2, is_rna=is_rna)

            aln_attrs = dict(handle["Alignment"].attrs)
            ref_bounds = RefCoord(
                aln_attrs["mapped_chrom"],
                aln_attrs["mapped_start"]-1,
                aln_attrs["mapped_end"]+3,
                aln_attrs["mapped_strand"] == "+")
            sig_fwd = (ref_bounds.fwd != is_rna)

            track.init_read_aln(read.read_id, ref_bounds)

            tombo_events = pd.DataFrame(np.array(handle["Events"]))
            tombo_start = handle["Events"].attrs["read_start_rel_to_raw"]
            
            raw_len = len(read.get_raw_data())
            starts = tombo_events["start"]

            #shift = 1
            #bases = tombo_events["base"].to_numpy().astype(str)
            #kmers_old = [nt.str_to_kmer("".join(bases[i:i+K]), 0) for i in range(len(bases)-K+1)]
            #kmers_old = shift*[kmers_old[0]] + kmers_old + (K-shift-1)*[kmers_old[-1]]
            kmers = track.load_aln_kmers(store=False)

            currents = tombo_events["norm_mean"]
            lengths = tombo_events["length"]

            if not sig_fwd:
                starts = np.flip(raw_len - tombo_start - starts - tombo_events["length"] - 1)
                currents = np.flip(currents)
                lengths = np.flip(lengths)

                #kmers_old = np.flip(nt.kmer_rev(kmers_old))

            lr = scipy.stats.linregress(currents, model[kmers])
            currents = lr.slope * currents + lr.intercept

            track.read_aln.set_aln(pd.DataFrame({
                "mrefs"    : track.read_aln.mrefs,
                "start"  : starts,
                "length" : lengths,
                "current"   : currents
            }).set_index("mrefs"))

            track.save_read(fast5_basename)

            pbar.update(pbar_count)

    track.close()

    pbar.finish()

SUBCMDS = [
    (nanopolish, NANOPOLISH_OPTS), 
    (tombo, TOMBO_OPTS)
]
