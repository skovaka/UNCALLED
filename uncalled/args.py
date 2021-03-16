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

class ArgParser:

    DEF_SUBCMDS = {
        "index" : (
            "Builds the UNCALLED index of a FASTA reference", 
            ["index"]
        ),
        "map" : (
            "Map fast5 files to a DNA reference",
            ["bwa", "fast5", "map"]
        ),
        "realtime" : (
            "Perform real-time targeted (ReadUntil) sequencing",
            ["bwa", "map", "ru", "realtime"],
        ),
        "sim" : (
            "Simulate real-time targeted sequencing.", 
            ["bwa", "sim", "map", "ru"],
        ),
        #"pafstats" : (
        #    "Computes speed and accuracy of UNCALLED mappings.",
        #    ["pafstats"]
        #)
    }

    OPTS_FORMAT = "add_%s_opts"

    def __init__(self, 
            argv=sys.argv[1:],
            subcmds=DEF_SUBCMDS, 
            conf=unc.Conf(),
            description="Rapidly maps raw nanopore signal to DNA references"):

        self.conf = conf
        print("CONF ", conf)

        self.parser = argparse.ArgumentParser(
                description="Rapidly maps raw nanopore signal to DNA references", 
                formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        subparsers = self.parser.add_subparsers(dest="subcmd")

        for subcmd, (desc, opts) in subcmds.items():
            sp = subparsers.add_parser(
                    subcmd, help=desc, 
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
            for opt_group in opts:
                opts_fn = getattr(self, self.OPTS_FORMAT % opt_group)
                opts_fn(sp)

        self.parse_args(argv)


    def parse_args(self, argv):
        args = self.parser.parse_args(argv)

        if hasattr(args, "conf"):
            if getattr(args, "conf") is not None:
                self.conf.load_toml(args.conf)
            delattr(args, "conf")

        if getattr(args,"rna",False):
            self.conf.set_r94_rna()
            delattr(args, "rna")

        for a, v in vars(args).items():

            if not a.startswith("_"):
                #Store argument in C++ Conf object if applicable present
                if v is not None and hasattr(self.conf, a):
                    setattr(self.conf, a, v)

                #Otherwise store within python
                else:
                    setattr(self, a, v) 


    def add_index_opts(self, p):
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

    def add_bwa_opts(self, p):
        p.add_argument(
                "bwa_prefix", 
                type=str, 
                help="BWA prefix to mapping to. Must be processed by \"uncalled index\"."
        )
        p.add_argument(
                "-p", "--idx-preset", 
                type=str, default=self.conf.idx_preset, 
                help="Mapping mode"
        )

    def add_ru_opts(self, p):
        #TODO: selectively enrich or deplete refs in index

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

    def add_sim_opts(self, p):
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
                type=float, default=self.conf.sim_speed, 
                help="")

    def add_realtime_opts(self, p):
        p.add_argument(
                "--host", 
                type=str, default=self.conf.host, 
                help=unc.Conf.host.__doc__
        )
        p.add_argument(
                "--port", 
                type=int, default=self.conf.port, 
                help=unc.Conf.port.__doc__
        )
        p.add_argument(
                "--duration", 
                type=float, default=self.conf.duration, 
                help=unc.Conf.duration.__doc__
        )

    #def add_list_ports_opts(p, conf):
    #    p.add_argument(
    #            "--log-dir", 
    #            type=str, default="/var/log/MinKNOW", 
    #            help="Directory to find MinKNOW log files"
    #    )

    def add_fast5_opts(self, p):
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
                type=str, default=None#, 
                #help=unc.Conf.read_list.__doc__
        )
        p.add_argument(
                "-n", "--max-reads", 
                type=int, default=None, 
                help=unc.Conf.max_reads.__doc__
        )

    #TODO get defautls from conf
    def add_map_opts(self, p):
        p.add_argument(
                "-t", "--threads", 
                type=int, default=self.conf.threads, 
                help="Number of threads to use for mapping"
        )
        p.add_argument(
                "--num-channels", 
                type=int, default=self.conf.num_channels, 
                help="Number of channels used in sequencing. If provided will use unique mapper for each channel. Useful for streaming normalization simulation."
        )
        p.add_argument(
                "-e", "--max-events", 
                type=int, default=self.conf.max_events, 
                help="Will give up on a read after this many events have been processed"
        )
        p.add_argument(
                "-c", "--max-chunks", 
                type=int, default=self.conf.max_chunks, required=True, 
                help="Will give up on a read after this many chunks have been processed. Only has effect when --unblock is set"
        )
        p.add_argument(
                "--chunk-time", 
                type=float, default=1, required=False, 
                help="Length of chunks in seconds"
        )
        p.add_argument(
                "--conf", 
                type=str, default=None, required=False, 
                help="Config file"
        )
        p.add_argument(
                "--rna", 
                action="store_true",
                help="Will use RNA parameters if set"
        )

        #TODO move to different parser set
        p.add_argument(
                "--mm2", 
                type=str, default=None, required=False, 
                help="Minimap2 PAF file for comparison"
        )


