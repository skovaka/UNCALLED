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

    BWA_OPTS = (
        ("bwa_prefix", {
            "type" : str, 
            "help" : "BWA prefix to mapping to. Must be processed by \"uncalled index\"."
        }),
        ("-p", "--idx-preset", {
            "type" : str, 
            "default" : unc.Conf().idx_preset, 
            "help" : "Mapping mode"
        }),
    )

    FAST5_OPTS = (
        ("fast5s", {
            "nargs" : '+', 
            "type" : str, 
            "help" : "Reads to map. Can be a directory which will be recursively searched for all files with the \".fast5\" extension, a text file containing one fast5 filename per line, or a comma-separated list of fast5 file names."
        }),
        ("-r", "--recursive", {
            "action" : "store_true"
        }),
        ("-l", "--read-list", {
            "type" : str, 
            "default" : None, 
            #"help" : unc.Conf.read_list.__doc__
        }),
        ("-n", "--max-reads", {
            "type" : int,
            "default" : None, 
            "help" : unc.Conf.max_reads.__doc__
        }),
    )

    MAP_OPTS = (
        ("-t", "--threads", {
            "type" : int, 
            "default" : unc.Conf().threads, 
            "help" : "Number of threads to use for mapping"
        }),
        ("--num-channels", {
            "type" : int, 
            "default" : unc.Conf().num_channels, 
            "help" : "Number of channels used in sequencing. If provided will use unique mapper for each channel. Useful for streaming normalization simulation."
        }),
        ("-e", "--max-events", {
            "type" : int, 
            "default" : unc.Conf().max_events, 
            "help" : "Will give up on a read after this many events have been processed"
        }),
        ("-c", "--max-chunks", {
            "type" : int, 
            "default" : unc.Conf().max_chunks, 
            "required" : True, 
            "help" : "Will give up on a read after this many chunks have been processed. Only has effect when --unblock is set"
        }),
        ("--chunk-time", {
            "type" : float, 
            "default" : 1, 
            "required" : False, 
            "help" : "Length of chunks in seconds"
        }),
        ("--conf", {
            "type" : str, 
            "default" : None, 
            "required" : False, 
            "help" : "Config file"
        }),
        ("--rna", {
            "action" : "store_true",
            "help" : "Will use RNA parameters if set"
        }),

        #TODO move to different parser set
        ("-o", "--out-prefix", {
            "type" : str, 
            "default" : None, 
            "required" : False, 
            "help" : "Output prefix"
        }),
        ("--mm2", {
            "type" : str, 
            "default" : None, 
            "required" : False, 
            "help" : "Minimap2 PAF file for comparison"
        }),
    ) #end MAP_OPTS

    SIM_OPTS = (
        ("fast5s", {
            "nargs" : '+', 
            "type" : str, 
            "help" : "Reads to unc. Can be a directory which will be recursively searched for all files with the \".fast5\" extension, a text file containing one fast5 filename per line, or a comma-separated list of fast5 file names."
        }),
        ("-r", "--recursive", {
            "action" : "store_true"
        }),
        ("--ctl-seqsum", {
            "type" : str, 
            "required" : True, 
            "help" : ""
        }),
        ("--unc-seqsum", {
            "type" : str, 
            "required" : True, 
            "help" : ""
        }),
        ("--unc-paf", {
            "type" : str, 
            "required" : True, 
            "help" : ""
        }),
        ("--sim-speed", {
            "type" : float, 
            "default" : unc.Conf().sim_speed, 
            "help" : ""
        }),
    )

    REALTIME_OPTS = (
        ("--host", {
            "type" : str, 
            "default" : unc.Conf().host, 
            "help" : unc.Conf.host.__doc__
        }),
        ("--port", {
            "type" : int, 
            "default" : unc.Conf().port, 
            "help" : unc.Conf.port.__doc__
        }),
        ("--duration", {
            "type" : float, 
            "default" : unc.Conf().duration, 
            "help" : unc.Conf.duration.__doc__
        }),
    )

    INDEX_OPTS = (
        ("fasta_filename", {
            "type" : str, 
            "help" : "FASTA file to index"
        }),
        ("-o", "--bwa-prefix", {
            "type" : str, 
            "default" : None, 
            "help" : "Index output prefix. Will use input fasta filename by default"
        }),
        ("-s", "--max-sample-dist", {
            "type" : int, 
            "default" : 100, 
            "help" : "Maximum average sampling distance between reference alignments."
        }),
        ("--min-samples", {
            "type" : int, 
            "default" : 50000, 
            "help" : "Minimum number of alignments to produce (approximate, due to deterministically random start locations}),"
        }),
        ("--max-samples", {
            "type" : int, 
            "default" : 1000000, 
            "help" : "Maximum number of alignments to produce (approximate, due to deterministically random start locations}),"
        }),
        ("-k", "--kmer-len", {
            "type" : int, 
            "default" : 5,
            "help" : "Model k-mer length"
        }),
        ("-1", "--matchpr1", {
            "type" : float, 
            "default" : 0.6334, 
            "help" : "Minimum event match probability"
        }),
        ("-2", "--matchpr2", {
            "type" : float, 
            "default" : 0.9838, 
            "help" : "Maximum event match probability"
        }),
        ("-f", "--pathlen-percentile", {
            "type" : float, 
            "default" : 0.05, 
            "help" : ""
        }),
        ("-m", "--max-replen", {
            "type" : int, 
            "default" : 100, 
            "help" : ""
        }),
        ("--probs", {
            "type" : str, 
            "default" : None, 
            "help" : "Find parameters with specified target probabilites (comma separated}),"
        }),
        ("--speeds", {
            "type" : str, 
            "default" : None, 
            "help" : "Find parameters with specified speed coefficents (comma separated}),"
        }),
    ) #end INDEX_OPTS


    PAFSTATS_OPTS = (
        ("infile",  {
            "type" : str, 
            "help" : "PAF file output by UNCALLED"
        }),
        ("-n", "--max-reads", {
            "required" : False, 
            "type" : int, 
            "default" : None, 
            "help" : "Will only look at first n reads if specified"
        }),
        ("-r", "--ref-paf", {
            "required" : False, 
            "type" : str, 
            "default" : None, 
            "help" : "Reference PAF file. Will output percent true/false positives/negatives with respect to reference. Reads not mapped in reference PAF will be classified as NA."
        }),
        ("-a", "--annotate", {
            "action" : 'store_true', 
            "help" : "Should be used with --ref-paf. Will output an annotated version of the input with T/P F/P specified in an 'rf' tag"
        }),
    ) #end pafstats opts

    DEF_SUBCMDS = {
        "index" : (
            "Builds the UNCALLED index of a FASTA reference", 
            INDEX_OPTS
            #["index"]
        ),
        "map" : (
            "Map fast5 files to a DNA reference",
            BWA_OPTS + FAST5_OPTS + MAP_OPTS
            #["bwa", "fast5", "map"]
        ),
        "realtime" : (
            "Perform real-time targeted (ReadUntil) sequencing",
            BWA_OPTS + MAP_OPTS + REALTIME_OPTS #TODO ADD RU OPTS!!
            #["bwa", "map", "ru", "realtime"],
        ),
        "sim" : (
            "Simulate real-time targeted sequencing.", 
            BWA_OPTS + SIM_OPTS + MAP_OPTS #TODO ADD RU OPTS!!
            #["bwa", "sim", "map", "ru"],
        ),
        "pafstats" : (
            "Computes speed and accuracy of UNCALLED mappings.",
            PAFSTATS_OPTS
            #["pafstats"]
        )
    }

    OPTS_FORMAT = "add_%s_opts"

    def __init__(self, 
            argv=sys.argv[1:],
            subcmds=DEF_SUBCMDS, 
            conf=unc.Conf(),
            description="Rapidly maps raw nanopore signal to DNA references"):

        self.conf = conf
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
            for opt in opts:
                args = opt[:-1]
                kwargs = opt[-1]
                arg = sp.add_argument(*args, **kwargs)

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
                if v is not None or not hasattr(self.conf, a): # and hasattr(self.conf, a):
                    setattr(self.conf, a, v)

                #Otherwise store within python
                #else:
                #    setattr(self, a, v) 

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

