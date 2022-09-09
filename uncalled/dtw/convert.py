"""Import DTW alignments produced by other tools into a database

subcommand options:
nanopolish  Convert alignments produced by nanopolish eventalign
tombo       Convert alignments produced by tomobo resquiggle"""

import sys, os
import numpy as np
import argparse
import pandas as pd
from ont_fast5_api.fast5_interface import get_fast5_file
import scipy.stats

from . import Tracks
from ..fast5 import Fast5Reader
from ..index import str_to_coord, RefCoord
from ..pore_model import PoreModel

import progressbar as progbar


def new(conf):
    """New convert interface"""
    conf.tracks.layers = ["dtw"]
    tracks = Tracks(conf=conf)

    for read_id, read in tracks.iter_reads():
        aln = read.alns[0]
        read.init_alignment(read_id, read.fast5s.get_read_file(read_id), aln.coords)
        
        dtw = aln.layers["dtw"].droplevel(1)#.set_index(aln.coords.ref_to_mref(aln.layer_refs, aln.all_fwd), drop=True)
        read.add_layers("dtw", dtw)
        read.write_alignment()

    tracks.close()

def nanopolish(conf):
    """Convert from nanopolish eventalign TSV to uncalled DTW track"""

    f5reader = Fast5Reader(conf=conf)
    conf.fast5_reader.load_bc = True

    read_filter = f5reader.get_read_filter()

    sample_rate = conf.read_buffer.sample_rate

    sys.stderr.write("Parsing TSV\n")
    csv_iter = pd.read_csv(
        conf.eventalign_tsv, sep="\t", chunksize=10000,
        usecols=["read_name","contig","position",
                 "start_idx","event_level_mean",
                 "event_length","strand","model_kmer"])

    tracks = Tracks(conf=conf)
    model = tracks[tracks.output_track].model

    def add_alns(events):
        groups = events.groupby(["contig", "read_name"])
        for (contig,read_id), df in groups:

            if not (read_filter is None or read_id in read_filter):
                continue

            start = df["position"].min()
            end = df["position"].max()+1
            fwd = df["strand"].iloc[0] == "t"

            #start = aln_attrs["mapped_start"]-2
            #end = aln_attrs["mapped_end"]+2
            if start < 0:
                clip = -start
                start = 0
            else:
                clip = 0

            ref_coord = RefCoord(contig, start, end, fwd)
            coords = tracks.index.get_coord_space(ref_coord, conf.is_rna, kmer_trim=True)

            df["mref"] = coords.ref_to_mref(df["position"].to_numpy()+2)

            kmers = model.str_to_kmer(df["model_kmer"])
            if conf.is_rna:
                kmers = model.kmer_rev(kmers)
            df["kmer"] = kmers

            df = df.rename(columns={
                    "start_idx" : "start",
                    "event_length" : "length",
                    "event_level_mean" : "current",
                 }).set_index("mref")

            fast5_name = f5reader.get_read_file(read_id)

            tracks.init_alignment(read_id, fast5_name, coords, {"dtw" : df})
            tracks.write_alignment()
            print(read_id)

    leftover = pd.DataFrame()

    for events in csv_iter:
        events['event_length'] = np.round(events['event_length'] * sample_rate).astype(int)
        events['sum'] = events['event_level_mean'] * events['event_length']

        events = pd.concat([leftover, events])

        i = events["read_name"] == events["read_name"].iloc[-1]
        leftover = events[i]
        events = events[~i]

        add_alns(events)

    if len(leftover) > 0:
        add_alns(leftover)

    tracks.close()

def tombo(conf):
    """Convert from Tombo resquiggled fast5s to uncalled DTW track"""
    f5reader = Fast5Reader(conf=conf)
    conf.fast5_reader.load_bc = True
    conf.pore_model.name = "r94_rna_tombo"
    read_filter = f5reader.get_read_filter()

    tracks = Tracks(conf=conf)
    
    fast5_files = f5reader.prms.fast5_files

    if conf.fast5_reader.max_reads > 0 and len(fast5_files) > conf.fast5_reader.max_reads:
        fast5_files = fast5_files[:conf.fast5_reader.max_reads]
    
    pbar = progbar.ProgressBar(
            widgets=[progbar.Percentage(), progbar.Bar(), progbar.ETA()], 
            maxval=len(fast5_files)).start()
    read_count = 0

    model = PoreModel(conf.pore_model)

    for fast5_name in fast5_files:

        fast5_basename = os.path.basename(fast5_name)

        try:
            fast5 = get_fast5_file(fast5_name, mode="r")
        except:
            sys.stderr.write(f"Unable to open \"{fast5_name}\". Skipping.")
            continue

        read_count += 1

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
        if is_rna != conf.is_rna:
            raise RuntimeError("Reads appear to be RNA but --rna not specified")

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

        ref_bounds = RefCoord(
            aln_attrs["mapped_chrom"],
            start, end,
            aln_attrs["mapped_strand"] == "+")

        sig_fwd = ref_bounds.fwd != is_rna

        coords = tracks.index.get_coord_space(ref_bounds, is_rna=is_rna, load_kmers=True, kmer_trim=True)
        aln_id,_ = tracks.init_alignment(read.read_id, fast5_name, coords)

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
        kmers = coords.ref_kmers.droplevel(0).loc[refs]
        #print(list(model.kmer_to_str(coords.kmers)[:10]))
        #print(list(tombo_events["base"][:15]))
        #print(len(coords.refs), len(coords.mrefs), len(coords.kmers), len(starts))
        #print(coords.refs)

        #TODO store scaling factors in pore model
        #currents = currents * 10.868760552593136 + 91.25486108714513

        df = pd.DataFrame({
                #"ref" : coords.refs,
                "start"  : starts,
                "length" : lengths,
                "current"   : currents,
                "kmer" : kmers
             }, index=refs)#.set_index("refs")

        tracks.add_layers("dtw", df, aln_id=aln_id)
        tracks.write_alignment()

        #track.save_read(fast5_basename)

        pbar.update(read_count)

    tracks.close()

    pbar.finish()

    if read_count == 0:
        sys.stderr.write("Warning: no reads were found (try including the --recursive option?)\n")

