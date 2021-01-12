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

def get_parser(conf):
    parser = argparse.ArgumentParser(
            description="Rapidly maps raw nanopore signal to DNA references", 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    sp = parser.add_subparsers(dest="subcmd")

    index_parser = sp.add_parser(
            "index", 
            help="Builds the UNCALLED index of a FASTA reference", 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_index_opts(index_parser, conf)

    map_parser = sp.add_parser(
            "map", 
            help="Map fast5 files to a DNA reference", 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_bwa_opt(map_parser, conf)
    add_fast5_opts(map_parser, conf)
    add_map_opts(map_parser, conf)

    rt_parser = sp.add_parser(
            "realtime", 
            help="Perform real-time targeted (ReadUntil) sequencing",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_bwa_opt(rt_parser, conf)
    add_map_opts(rt_parser, conf)
    add_ru_opts(rt_parser, conf)
    add_realtime_opts(rt_parser, conf)

    sim_parser = sp.add_parser(
            "sim", 
            help="Simulate real-time target sequencing.", 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    add_bwa_opt(sim_parser, conf)
    add_sim_opts(sim_parser, conf)
    add_map_opts(sim_parser, conf)
    add_ru_opts(sim_parser, conf)

    ps_parser = sp.add_parser(
            "pafstats", 
            help="Computes speed and accuracy of UNCALLED mappings.", #Given an UNCALLED PAF file, will compute mean/median BP mapped per second, number of BP required to map each read, and total number of milliseconds to map each read. Can also optionally compute accuracy with respect to reference alignments, for example output by minimap2.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    unc.pafstats.add_opts(ps_parser)
    #TODO move here

    #lp_parser = sp.add_parser("list-ports", help="List the port of all MinION devices detected in the current MinKNOW session", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #add_list_ports_opts(lp_parser, conf)

    return parser

def add_index_opts(p, conf):
    p.add_argument(
            "fasta_filename", 
            type=str, 
            help="FASTA file to index"
    )
    p.add_argument(
            "-o", "--bwa-prefix", 
            type=str, default=None, 
            help="Index output prefix. Will use input fasta filename by default"
    )
    p.add_argument(
            "-s", "--max-sample-dist", 
            type=int, default=100, 
            help="Maximum average sampling distance between reference alignments."
    )
    p.add_argument(
            "--min-samples", 
            type=int, default=50000, 
            help="Minimum number of alignments to produce (approximate, due to deterministically random start locations)"
    )
    p.add_argument(
            "--max-samples", 
            type=int, default=1000000, 
            help="Maximum number of alignments to produce (approximate, due to deterministically random start locations)"
    )
    p.add_argument(
            "-k", "--kmer-len", 
            type=int, default=5,
            help="Model k-mer length"
    )
    p.add_argument(
            "-1", "--matchpr1", 
            default=0.6334, type=float, 
            help="Minimum event match probability"
    )
    p.add_argument(
            "-2", "--matchpr2", 
            type=float, default=0.9838, 
            help="Maximum event match probability"
    )
    p.add_argument(
            "-f", "--pathlen-percentile", 
            type=float, default=0.05, 
            help=""
    )
    p.add_argument(
            "-m", "--max-replen", 
            type=int, default=100, 
            help=""
    )
    p.add_argument(
            "--probs", 
            type=str, default=None, 
            help="Find parameters with specified target probabilites (comma separated)"
    )
    p.add_argument(
            "--speeds", 
            type=str, default=None, 
            help="Find parameters with specified speed coefficents (comma separated)"
    )

def add_bwa_opt(p, conf):
    p.add_argument(
            "bwa_prefix", 
            type=str, 
            help="BWA prefix to mapping to. Must be processed by \"uncalled index\"."
    )
    p.add_argument(
            "-p", "--idx-preset", 
            type=str, default=conf.idx_preset, 
            help="Mapping mode"
    )

def add_ru_opts(p, conf):
    #TODO: selectively enrich or deplete refs in index
    p.add_argument(
            "-c", "--max-chunks", 
            type=int, default=conf.max_chunks, required=True, 
            help="Will give up on a read after this many chunks have been processed. Only has effect when --unblock is set"
    )
    p.add_argument(
            "--chunk-time", 
            type=float, default=1, required=False, 
            help="Length of chunks in seconds"
    )

    modes = p.add_mutually_exclusive_group(required=True)
    modes.add_argument(
            "-D", "--deplete", action='store_const', 
            const=unc.RealtimePool.DEPLETE, dest='realtime_mode', 
            help="Will eject reads that align to index"
    )
    modes.add_argument(
            "-E", "--enrich", action='store_const',
            const=unc.RealtimePool.ENRICH, dest='realtime_mode', 
            help="Will eject reads that don't align to index"
    )

    active = p.add_mutually_exclusive_group()
    active.add_argument(
            "--full", action='store_const', 
            const=unc.RealtimePool.FULL, dest='active_chs',
            help="Will monitor all pores if set (default)"
    )
    active.add_argument(
            "--even", action='store_const', 
            const=unc.RealtimePool.EVEN, dest='active_chs', 
            help="Will only monitor even pores if set"
    )
    active.add_argument(
            "--odd", action='store_const', 
            const=unc.RealtimePool.ODD, dest='active_chs', 
            help="Will only monitor odd pores if set")

def add_sim_opts(p, conf):
    p.add_argument(
            "fast5s", nargs='+', type=str, 
            help="Reads to unc. Can be a directory which will be recursively searched for all files with the \".fast5\" extension, a text file containing one fast5 filename per line, or a comma-separated list of fast5 file names."
    )
    p.add_argument(
            "-r", "--recursive", 
            action="store_true"
    )
    p.add_argument(
            "--ctl-seqsum", 
            type=str, required=True, 
            help=""
    )
    p.add_argument(
            "--unc-seqsum", 
            type=str, required=True, 
            help=""
    )
    p.add_argument(
            "--unc-paf", 
            type=str, required=True, 
            help=""
    )
    p.add_argument(
            "--sim-speed", 
            type=float, default=conf.sim_speed, 
            help="")

def add_realtime_opts(p, conf):
    p.add_argument(
            "--host", 
            type=str, default=conf.host, 
            help=unc.Conf.host.__doc__
    )
    p.add_argument(
            "--port", 
            type=int, default=conf.port, 
            help=unc.Conf.port.__doc__
    )
    p.add_argument(
            "--duration", 
            type=float, default=conf.duration, 
            help=unc.Conf.duration.__doc__
    )

#def add_list_ports_opts(p, conf):
#    p.add_argument(
#            "--log-dir", 
#            type=str, default="/var/log/MinKNOW", 
#            help="Directory to find MinKNOW log files"
#    )

def add_fast5_opts(p, conf):
    p.add_argument(
            "fast5s", nargs='+', type=str, 
            help="Reads to map. Can be a directory which will be recursively searched for all files with the \".fast5\" extension, a text file containing one fast5 filename per line, or a comma-separated list of fast5 file names."
    )
    p.add_argument(
            "-r", "--recursive", 
            action="store_true"
    )
    p.add_argument(
            "-l", "--read-list", 
            type=str, default=None, 
            help=unc.Conf.read_list.__doc__
    )
    p.add_argument(
            "-n", "--max-reads", 
            type=int, default=None, 
            help=unc.Conf.max_reads.__doc__
    )

#TODO get defautls from conf
def add_map_opts(p, conf):
    p.add_argument(
            "-t", "--threads", 
            type=int, default=conf.threads, 
            help="Number of threads to use for mapping"
    )
    p.add_argument(
            "--num-channels", 
            type=int, default=conf.num_channels, 
            help="Number of channels used in sequencing. If provided will use unique mapper for each channel. Useful for streaming normalization simulation."
    )
    p.add_argument(
            "-e", "--max-events", 
            type=int, default=conf.max_events, 
            help="Will give up on a read after this many events have been processed"
    )

def load_conf(argv):
    conf = unc.Conf()
    parser = get_parser(conf)
    args = parser.parse_args()

    for a, v in vars(args).items():
        if v == None:
            continue

        if (not a.startswith("_")):
            if hasattr(conf, a):
                setattr(conf, a, v)

    return parser, conf, args
