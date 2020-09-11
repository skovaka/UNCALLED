#!/usr/bin/env python

# MIT License
#
# Copyright (c) 2018 Sam Kovaka <skovaka@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import sys                         
import os
import argparse
import numpy as np
import uncalled as unc

def add_index_opts(p, conf):
    p.add_argument("fasta_filename", type=str, help="FASTA file to index")
    p.add_argument("-o", "--bwa-prefix", default=None, type=str, help="Index output prefix. Will use input fasta filename by default")
    p.add_argument("-s", "--max-sample-dist", default=100, type=int, help="Maximum average sampling distance between reference self-alignments.")
    p.add_argument("--min-samples", default=50000, type=int, help="Minimum number of self-alignments to produce (approximate, due to deterministically random start locations)")
    p.add_argument("--max-samples", default=1000000, type=int, help="Maximum number of self-alignments to produce (approximate, due to deterministically random start locations)")
    p.add_argument("-k", "--kmer-len", default=5, type=int, help="Model k-mer length")
    p.add_argument("-1", "--matchpr1", default=0.6334, type=float, help="Minimum event match probability")
    p.add_argument("-2", "--matchpr2", default=0.9838, type=float, help="Maximum event match probability")
    p.add_argument("-f", "--pathlen-percentile", default=0.05, type=float, help="")
    p.add_argument("-m", "--max-replen", default=100, type=int, help="")
    p.add_argument("--probs", default=None, type=str, help="Find parameters with specified target probabilites (comma separated)")
    p.add_argument("--speeds", default=None, type=str, help="Find parameters with specified speed coefficents (comma separated)")

def add_bwa_opt(p, conf):
    p.add_argument("bwa_prefix", type=str, help="BWA prefix to mapping to. Must be processed by \"uncalled index\".")
    p.add_argument("-p", "--idx-preset", default=conf.idx_preset, type=str, help="Mapping mode")


def add_ru_opts(p, conf):
    #TODO: selectively enrich or deplete refs in index
    p.add_argument("-c", "--max-chunks", default=conf.max_chunks, required=True, type=int, help="Will give up on a read after this many chunks have been processed. Only has effect when --unblock is set")
    p.add_argument("--chunk-time", required=False, type=float, default=1, help="Length of chunks in seconds")

    modes = p.add_mutually_exclusive_group(required=True)
    modes.add_argument("-D", "--deplete", action='store_const', const=mapping.RealtimeMode.DEPLETE, dest='realtime_mode', help="Will eject reads that align to index")
    modes.add_argument("-E", "--enrich", action='store_const',  const=mapping.RealtimeMode.ENRICH, dest='realtime_mode', help="Will eject reads that don't align to index")

    active = p.add_mutually_exclusive_group()
    active.add_argument("--full", action='store_const', const=mapping.ActiveChs.FULL, dest='active_chs', help="Will monitor all pores if set (default)")
    active.add_argument("--even", action='store_const', const=mapping.ActiveChs.EVEN, dest='active_chs', help="Will only monitor even pores if set")
    active.add_argument("--odd", action='store_const', const=mapping.ActiveChs.ODD, dest='active_chs', help="Will only monitor odd pores if set")

def add_sim_opts(p, conf):
    p.add_argument("fast5s", nargs='+', type=str, help="Reads to mapping. Can be a directory which will be recursively searched for all files with the \".fast5\" extension, a text file containing one fast5 filename per line, or a comma-separated list of fast5 file names.")
    p.add_argument("-r", "--recursive", action="store_true")
    p.add_argument("--ctl-seqsum", required=True, type=str, help="")
    p.add_argument("--unc-seqsum", required=True, type=str, help="")
    p.add_argument("--unc-paf", required=True, type=str, help="")
    p.add_argument("--sim-speed", required=False, default=1.0, type=float, help="")

def add_realtime_opts(p, conf):
    p.add_argument('--host', default=conf.host, help='MinKNOW server host.')
    p.add_argument('--port', type=int, default=conf.port, help='MinKNOW server port.')
    p.add_argument('--duration', type=float, default=conf.duration, help='Duration to map real-time run in hours. Should be slightly longer than specified runtime to add wiggle room.')

def add_list_ports_opts(p, conf):
    p.add_argument('--log-dir', default='/var/log/MinKNOW', help='Directory to find MinKNOW log files')

def add_fast5_opts(p, conf):
    p.add_argument("fast5s", nargs='+', type=str, help="Reads to mapping. Can be a directory which will be recursively searched for all files with the \".fast5\" extension, a text file containing one fast5 filename per line, or a comma-separated list of fast5 file names.")
    p.add_argument("-r", "--recursive", action="store_true")

    p.add_argument("-l", "--read-list", default=None, type=str, help="Only map reads listed in this file")
    p.add_argument("-n", "--max-reads", type=int, default=None, help="Maximum number of reads to map")

#TODO get defautls from conf
def add_map_opts(p, conf):
    p.add_argument("-t", "--threads", default=conf.threads, type=int, help="Number of threads to use for mapping")
    p.add_argument("--num-channels", default=conf.num_channels, type=int, help="Number of channels used in sequencing. If provided will use unique mapper for each channel. Useful for streaming normalization simulation.")
    p.add_argument("-e", "--max-events", default=conf.max_events, type=int, help="Will give up on a read after this many events have been processed")

def get_parser(conf):
    parser = argparse.ArgumentParser(description="Rapidly maps raw nanopore signal to DNA references", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sp = parser.add_subparsers(dest="subcmd")

    index_parser = sp.add_parser("index", help="Calculates reference-specific parameters needed to map to a given a BWA-index.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #add_bwa_opt(index_parser)
    add_index_opts(index_parser, conf)

    map_parser = sp.add_parser("map", help="Map fast5 files to a BWA index that has been processed by \"uncalled index\"", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_bwa_opt(map_parser, conf)
    add_fast5_opts(map_parser, conf)
    add_map_opts(map_parser, conf)

    rt_parser = sp.add_parser("realtime", help="Perform real-time targeted sequencing", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_bwa_opt(rt_parser, conf)
    add_map_opts(rt_parser, conf)
    add_ru_opts(rt_parser, conf)
    add_realtime_opts(rt_parser, conf)

    sim_parser = sp.add_parser("sim", help="Simulate real-time target sequencing", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_bwa_opt(sim_parser, conf)
    add_sim_opts(sim_parser, conf)
    add_map_opts(sim_parser, conf)
    add_ru_opts(sim_parser, conf)

    ps_parser = sp.add_parser("pafstats", help="Computes speed and accuracy of UNCALLED mappings. Given an UNCALLED PAF file, will compute mean/median BP mapped per second, number of BP required to map each read, and total number of milliseconds to map each read. Can also optionally compute accuracy with respect to reference alignments, for example output by minimap2.")
    pafstats.add_opts(ps_parser)

    lp_parser = sp.add_parser("list-ports", help="List the port of all MinION devices detected in the current MinKNOW session", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_list_ports_opts(lp_parser, conf)

    return parser

def load_args(args, conf):
    for a, v in vars(args).items():
        if v == None:
            #sys.stderr.write("%s IS NONE\n" % a)
            continue

        if (not a.startswith("_")):
            if hasattr(conf, a):
                setattr(conf, a, v)
            #else:
            #    sys.stderr.write("%s\t%s\n" % (a, str(v)))
