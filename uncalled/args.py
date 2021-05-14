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
import toml
import inspect



class ArgParser:

    class Groups:
        GLOBAL    = ""
        MAPPER    = "mapper"
        READ      = "read_buffer"
        NORM      = "normalizer"
        EVENT     = "event_detector"
        PROFILER  = "event_profiler"
        SEED      = "seed_tracker"
        FAST5     = "fast5_reader"
        REALTIME  = "realtime"
        SIMULATOR = "simulator"
        INDEX     = "index"

    class Opt:
        def __init__(self, args, group_name=None, param=None, **kwargs):
            self.args = args if type(args) == tuple else (args,)
            self.group_name = group_name

            self.param = param
            if self.param is None:
                for arg in self.args:
                    if arg.startswith("--") or not arg.startswith("-"):
                        self.param = arg.strip("-").replace("-","_")
                        break
                if self.param is None:
                    sys.stderr.write("Error: must specify parameter name for flag \"%s\"\n" % self.args)
                    sys.exit(0)

            self.extra_kw = kwargs

        def get_dest(self, conf):
            if self.group_name is None:
                return (conf, self.param)
            return (conf.get_group(self.group_name), self.param)

        def get_kwargs(self, conf):
            if self.group_name is None:
                return self.extra_kw

            group = conf.get_group(self.group_name)

            if self.param is not None:
                param_name = self.param
            else:
                param_name = flag_to_var(self.args[-1])

            default = getattr(group, param_name)
            doc = getattr(type(group), param_name).__doc__

            if hasattr(group, "_types"):
                _type = getattr(group, "_types")[param_name]
            else:
                _type = type(default)

            kwargs = {"help" : doc}

            if _type != bool:
                kwargs["type"] = _type
                kwargs["default"] = default

            kwargs.update(self.extra_kw)

            return kwargs

        def set_val(self, conf, val):
            if self.group is None or len(self.group) == 0:
                group = conf
            else:
                group = getattr(conf, self.group)

            setattr(group, self.param, val)


    BWA_OPTS = (
        Opt("bwa_prefix", Groups.MAPPER),
        Opt(("-p", "--idx-preset"), Groups.MAPPER),
    )

    FAST5_OPTS = (
        Opt("fast5_files", Groups.FAST5, nargs="+", type=str),
        Opt(("-r", "--recursive"), Groups.FAST5, action="store_true"),
        Opt(("-l", "--read-filter"), Groups.FAST5, type=unc.Fast5Reader.parse_read_str),
        Opt(("-x", "--fast5-index"), Groups.FAST5),
        Opt(("-n", "--max-reads"), Groups.FAST5)
    )

    MAP_OPTS = (
        Opt(("-t", "--threads"), Groups.GLOBAL),
        Opt("--num-channels", Groups.READ),
        Opt(("-e", "--max-events"), Groups.MAPPER),
        Opt(("-c", "--max-chunks"), Groups.READ),
        Opt("--chunk-time", Groups.READ),
        Opt("--conf", 
            type = str, 
            default = None, 
            required = False, 
            help = "Config file"
        ),
        Opt("--rna", 
            action = "store_true",
            help = "Will use RNA parameters if set"
        ),

        #TODO move to different parser set
        Opt(("-o", "--out-prefix"), 
            type = str, 
            default = None, 
            required = False, 
            help = "Output prefix"
        ),
        Opt("--mm2", 
            type = str, 
            default = None, 
            required = False, 
            help = "Minimap2 PAF file for comparison"
        ),
    ) #end MAP_OPTS

    SIM_OPTS = (
        Opt("fast5s", 
            nargs = '+', 
            type = str, 
            help = "Reads to unc. Can be a directory which will be recursively searched for all files with the \".fast5\" extension, a text file containing one fast5 filename per line, or a comma-separated list of fast5 file names."
        ),
        Opt(("-r", "--recursive"), 
            action = "store_true"
        ),
        Opt("--ctl-seqsum", Groups.SIMULATOR),
        Opt("--unc-seqsum", Groups.SIMULATOR),
        Opt("--unc-paf", Groups.SIMULATOR),
        Opt("--sim-speed", Groups.SIMULATOR),
    )

    REALTIME_OPTS = (
        Opt("--host", Groups.REALTIME),
        Opt("--port", Groups.REALTIME),
        Opt("--duration", Groups.REALTIME),
    )

    INDEX_OPTS = (
        Opt("fasta_filename", Groups.INDEX),
        Opt(("-o", "--bwa-prefix"), Groups.INDEX),
        Opt(("-s", "--max-sample-dist"), Groups.INDEX),
        Opt("--min-samples", Groups.INDEX),
        Opt("--max-samples", Groups.INDEX),
        Opt(("-k", "--kmer-len"), Groups.INDEX),
        Opt(("-1", "--matchpr1"), Groups.INDEX),
        Opt(("-2", "--matchpr2"), Groups.INDEX),
        Opt(("-f", "--pathlen-percentile"), Groups.INDEX),
        Opt(("-m", "--max-replen"), Groups.INDEX),
        Opt("--probs", Groups.INDEX),
        Opt("--speeds", Groups.INDEX),
    ) #end INDEX_OPTS


    PAFSTATS_OPTS = (
        Opt("infile",  
            type = str, 
            help = "PAF file output by UNCALLED"
        ),
        Opt(("-n", "--max-reads"), 
            required = False, 
            type = int, 
            default = None, 
            help = "Will only look at first n reads if specified"
        ),
        Opt(("-r", "--ref-paf"), 
            required = False, 
            type = str, 
            default = None, 
            help = "Reference PAF file. Will output percent true/false positives/negatives with respect to reference. Reads not mapped in reference PAF will be classified as NA."
        ),
        Opt(("-a", "--annotate"), 
            action = 'store_true', 
            help = "Should be used with --ref-paf. Will output an annotated version of the input with T/P F/P specified in an 'rf' tag"
        ),
    ) #end pafstats opts

    DEFAULT_SUBCMDS = {
        "index" : (
            "Builds the UNCALLED index of a FASTA reference", 
            INDEX_OPTS
        ),
        "map" : (
            "Map fast5 files to a DNA reference",
            BWA_OPTS + FAST5_OPTS + MAP_OPTS
        ),
        "realtime" : (
            "Perform real-time targeted (ReadUntil) sequencing",
            BWA_OPTS + MAP_OPTS + REALTIME_OPTS #TODO ADD RU OPTS!!
        ),
        "sim" : (
            "Simulate real-time targeted sequencing.", 
            BWA_OPTS + SIM_OPTS + MAP_OPTS #TODO ADD RU OPTS!!
        ),
        "pafstats" : (
            "Computes speed and accuracy of UNCALLED mappings.",
            PAFSTATS_OPTS
        )
    }

    OPTS_FORMAT = "add_%s_opts"

    def __init__(self, 
            argv=sys.argv[1:],
            subcmds=DEFAULT_SUBCMDS, 
            conf=None,
            description="Rapidly maps raw nanopore signal to DNA references"):

        if conf is None:
            self.conf = unc.Conf()
        else:
            self.conf = conf

        self.parser = argparse.ArgumentParser(
                description="Rapidly maps raw nanopore signal to DNA references", 
                formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        subparsers = self.parser.add_subparsers(dest="subcmd")
        
        self.dests = dict()

        for subcmd, (desc, opts) in subcmds.items():
            sp = subparsers.add_parser(
                    subcmd, help=desc, 
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )
            for opt in opts:
                arg = sp.add_argument(*opt.args, **opt.get_kwargs(self.conf))
                self.dests[arg.dest] = opt.get_dest(self.conf)

        self.parse_args(argv)

    def parse_args(self, argv):
        args = self.parser.parse_args(argv)

        if hasattr(args, "conf"):
            if getattr(args, "conf") is not None:
                self.conf.load_toml(args.conf)

        if getattr(args,"rna",False):
            self.conf.set_r94_rna()

        for name, value in vars(args).items():

            if not name.startswith("_"):
                if name in self.dests:
                    group, param = self.dests[name]
                else: #TODO only needed for subcommand, which I should rethink
                    group = self.conf
                    param = name

                if value is not None or not hasattr(group, param):
                    setattr(group, param, value)
                #if value is not None or not hasattr(self.conf, name):
                #    setattr(self.conf, name, value)


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

