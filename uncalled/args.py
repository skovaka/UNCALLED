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

def flag_to_var(flag):
    return flag.strip("-").replace("-","_")

def is_public(name):
    print(name)
    return not name.startswith("_")

TOML_TYPES = {int, float, str, list}

def param_valid(name, val):
    return (not name.startswith("_") and
            type(val) in TOML_TYPES and
            (not hasattr(val, "__len__") or len(val) > 0))

def conf_to_toml(conf, filename=None):
    out = dict()

    for param in conf._GLOBAL_PARAMS:
        val = getattr(conf,param)
        if param_valid(param, val):
            out[param] = val

    for param, val in vars(conf).items():
        if param_valid(param, val):
            out[param] = val

    for group_name in conf._PARAM_GROUPS:
        group = getattr(conf, group_name)
        out[group_name] = dict()
        for param in dir(group):
            val = getattr(group,param)
            if param_valid(param, val):
                out[group_name][param] = val

    if filename is None:
        return toml.dumps(out)
    else:
        return toml.dump(out, filename)

def toml_to_conf(filename, conf):
    toml_dict = toml.load(filename)

    for name, val in toml_dict.items():
        if isinstance(val, dict):
            group = getattr(conf, name, None)
            if group is None:
                group = dict()
                setattr(conf, name, group)
            for param, param_val in val.items():
                group[param] = param_val
        else:
            setattr(conf, name, val)
    
    return conf

class ArgParser:

    class Groups:
        GLOBAL   = unc.Conf
        MAPPER   = unc.Conf.mapper
        READ     = unc.Conf.read_buffer
        NORM     = unc.Conf.normalizer
        EVENT    = unc.Conf.event_detector
        PROFILER = unc.Conf.event_profiler
        SEED     = unc.Conf.seed_tracker
        FAST5    = unc.Conf.fast5_reader
        REALTIME = unc.Conf.realtime

    class Opt:
        def __init__(self, args, group=None, kw={}, param=None):
            self.args = args if type(args) == tuple else (args,)
            self.group = group
            self.kw = kw
            self.param = param

    BWA_OPTS = (
        Opt("bwa_prefix", Groups.MAPPER),
        Opt(("-p", "--idx-preset"), Groups.MAPPER),
    )

    FAST5_OPTS = (
        Opt("fast5s", kw={
            "nargs" : '+', 
            "type" : str, 
            "help" : "Reads to map. Can be a directory which will be searched for all files with the \".fast5\" extension, a text file containing one fast5 filename per line, or a comma-separated list of fast5 file names."
        }),
        Opt(("-r", "--recursive"), kw={
            "action" : "store_true",
            "help" : "Will perform recursive directory search for fast5 files if specified",
        }),
        Opt(("-l", "--read-list"), kw={
            "type" : str, 
            "default" : None, 
            "help" : "List of read IDs, either comma separated or path to text file containing one read ID per line. Will only map these reads if specified"
        }),
        Opt(("-n", "--max-reads"), Groups.FAST5)
    )

    MAP_OPTS = (
        Opt(("-t", "--threads"), Groups.GLOBAL),
        Opt("--num-channels", Groups.READ),
        Opt(("-e", "--max-events"), Groups.MAPPER),
        Opt(("-c", "--max-chunks"), Groups.READ),
        Opt("--chunk-time", Groups.READ),
        Opt("--conf", kw={
            "type" : str, 
            "default" : None, 
            "required" : False, 
            "help" : "Config file"
        }),
        Opt("--rna", kw={
            "action" : "store_true",
            "help" : "Will use RNA parameters if set"
        }),

        #TODO move to different parser set
        Opt(("-o", "--out-prefix"), kw={
            "type" : str, 
            "default" : None, 
            "required" : False, 
            "help" : "Output prefix"
        }),
        Opt("--mm2", kw={
            "type" : str, 
            "default" : None, 
            "required" : False, 
            "help" : "Minimap2 PAF file for comparison"
        }),
    ) #end MAP_OPTS

    SIM_OPTS = (
        Opt("fast5s", kw={
            "nargs" : '+', 
            "type" : str, 
            "help" : "Reads to unc. Can be a directory which will be recursively searched for all files with the \".fast5\" extension, a text file containing one fast5 filename per line, or a comma-separated list of fast5 file names."
        }),
        Opt(("-r", "--recursive"), kw={
            "action" : "store_true"
        }),
        Opt("--ctl-seqsum", kw={
            "type" : str, 
            "required" : True, 
            "help" : ""
        }),
        Opt("--unc-seqsum", kw={
            "type" : str, 
            "required" : True, 
            "help" : ""
        }),
        Opt("--unc-paf", kw={
            "type" : str, 
            "required" : True, 
            "help" : ""
        }),
        Opt("--sim-speed", kw={
            "type" : float, 
            "default" : unc.Conf().sim_speed, 
            "help" : ""
        }),
    )

    REALTIME_OPTS = (
        Opt("--host", kw={
            "type" : str, 
            "default" : unc.Conf().host, 
            "help" : unc.Conf.host.__doc__
        }),
        Opt("--port", kw={
            "type" : int, 
            "default" : unc.Conf().port, 
            "help" : unc.Conf.port.__doc__
        }),
        Opt("--duration", kw={
            "type" : float, 
            "default" : unc.Conf().duration, 
            "help" : unc.Conf.duration.__doc__
        }),
    )

    INDEX_OPTS = (
        Opt("fasta_filename", kw={
            "type" : str, 
            "help" : "FASTA file to index"
        }),
        Opt(("-o", "--bwa-prefix"), kw={
            "type" : str, 
            "default" : None, 
            "help" : "Index output prefix. Will use input fasta filename by default"
        }),
        Opt(("-s", "--max-sample-dist"), kw={
            "type" : int, 
            "default" : 100, 
            "help" : "Maximum average sampling distance between reference alignments."
        }),
        Opt("--min-samples", kw={
            "type" : int, 
            "default" : 50000, 
            "help" : "Minimum number of alignments to produce (approximate, due to deterministically random start locations}),"
        }),
        Opt("--max-samples", kw={
            "type" : int, 
            "default" : 1000000, 
            "help" : "Maximum number of alignments to produce (approximate, due to deterministically random start locations}),"
        }),
        Opt(("-k", "--kmer-len"), kw={
            "type" : int, 
            "default" : 5,
            "help" : "Model k-mer length"
        }),
        Opt(("-1", "--matchpr1"), kw={
            "type" : float, 
            "default" : 0.6334, 
            "help" : "Minimum event match probability"
        }),
        Opt(("-2", "--matchpr2"), kw={
            "type" : float, 
            "default" : 0.9838, 
            "help" : "Maximum event match probability"
        }),
        Opt(("-f", "--pathlen-percentile"), kw={
            "type" : float, 
            "default" : 0.05, 
            "help" : ""
        }),
        Opt(("-m", "--max-replen"), kw={
            "type" : int, 
            "default" : 100, 
            "help" : ""
        }),
        Opt("--probs", kw={
            "type" : str, 
            "default" : None, 
            "help" : "Find parameters with specified target probabilites (comma separated}),"
        }),
        Opt("--speeds", kw={
            "type" : str, 
            "default" : None, 
            "help" : "Find parameters with specified speed coefficents (comma separated}),"
        }),
    ) #end INDEX_OPTS


    PAFSTATS_OPTS = (
        Opt("infile",  kw={
            "type" : str, 
            "help" : "PAF file output by UNCALLED"
        }),
        Opt(("-n", "--max-reads"), kw={
            "required" : False, 
            "type" : int, 
            "default" : None, 
            "help" : "Will only look at first n reads if specified"
        }),
        Opt(("-r", "--ref-paf"), kw={
            "required" : False, 
            "type" : str, 
            "default" : None, 
            "help" : "Reference PAF file. Will output percent true/false positives/negatives with respect to reference. Reads not mapped in reference PAF will be classified as NA."
        }),
        Opt(("-a", "--annotate"), kw={
            "action" : 'store_true', 
            "help" : "Should be used with --ref-paf. Will output an annotated version of the input with T/P F/P specified in an 'rf' tag"
        }),
    ) #end pafstats opts

    DEF_SUBCMDS = {
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
                #flags = opt[0]
                #if type(flags) == str:
                #    flags = (flags,)

                if opt.group is None:
                    kwargs = opt.kw
                else:
                    if opt.group != self.Groups.GLOBAL:
                        group = opt.group.__get__(self.conf, unc.Conf)
                    else:
                        group = self.conf

                    if opt.param is not None:
                        param_name = opt.param
                    else:
                        param_name = flag_to_var(opt.args[-1])

                    default = getattr(group, param_name)
                    doc = getattr(type(group), param_name).__doc__
                    kwargs = {
                        "default" : default,
                        "type"    : type(default),
                        "help"    : doc
                    }

                arg = sp.add_argument(*opt.args, **kwargs)

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
                if v is not None or not hasattr(self.conf, a):
                    setattr(self.conf, a, v)


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

