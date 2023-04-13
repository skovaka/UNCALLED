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
import argparse

from .config import Config

from . import __title__, __version__, __summary__

FAST5_PARAM = "paths"
CONFIG_PARAM = "config"
SPECIAL_PARAMS = {FAST5_PARAM, CONFIG_PARAM}

def comma_split(s):
    return s.split(",")
 
def ref_coords(coord_str):             
    spl = coord_str.split(":")         
    ch = spl[0]                        
    st,en = spl[1].split("-")          
                                       
    coord = (ch, int(st), int(en))     
                                       
    if len(spl) == 2:                  
        return coord                   
    else:                              
        return coord + (spl[2] == "+",)

class Formatter(argparse.ArgumentDefaultsHelpFormatter,argparse.RawDescriptionHelpFormatter):
    pass

class ArgParser:
    def __init__(self, 
            subcmds, 
            desc,
            config=None):

        if config is None:
            self.config = Config()
        else:
            self.config = config

        self.dests = dict()
        self.fns = dict()

        self.parser = argparse.ArgumentParser(
                description=desc, 
               #formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                formatter_class=argparse.RawDescriptionHelpFormatter,
                prog=__title__,
                usage="%(prog)s [subcommand] [-h] [-v]"
        )

        self.parser.add_argument("--cprof", type=str, default=None, help=argparse.SUPPRESS)
        self.parser.add_argument("-v", "--version", action="version", version=__version__)

        self._add_subcmds(self.parser, subcmds, None)

    def _add_subcmds(self, parser, subcmds, desc=None):
        subparsers = parser.add_subparsers(title="Subcommands", description=desc, help=argparse.SUPPRESS
        )
        for cmd,(module,doc,opts) in subcmds.items():

            #if isinstance(subcmd, tuple):
            #    subcmd, opts = subcmd
            #    main_func = subcmd
            #else:
            #    opts = getattr(subcmd, "OPTS", None)
            #    main_func = getattr(subcmd, "main", None)

            #subcmd_name = subcmd.__name__.split(".")[-1].strip("_")

            #if main_func is not None:
            #    desc = main_func.__doc__
            #else:
            #    desc = subcmd.__doc__

            sp = subparsers.add_parser(
                cmd, prog=" ".join([__title__, cmd]), 
                description=doc, formatter_class=Formatter
            )

            #if main_func is not None:
            if module is not None:
                sp.set_defaults(_cmd=(module, cmd))
                for opt in opts:
                    if type(opt) is Opt:
                        self._add_opt(cmd, sp, opt)
                    elif type(opt) is MutexOpts:
                        self._add_mutex_opts(cmd, sp, opt)

            elif isinstance(opts, dict):
                self._add_subcmds(sp, opts, None)

            else:
                raise RuntimeError(f"Invalid subcommand definition: \"{cmd}\"")

    def _add_opt(self, subcmd, parser, opt):
        arg = parser.add_argument(*opt.args, **opt.get_kwargs(self.config))

        if opt.has_dest:
            name = arg.dest
            dest = opt.get_dest(self.config)

            if name in self.dests and self.dests[name] != dest:
                raise RuntimeError("Conflicting ArgParser dests with name \"%s\"" % arg.dest)
            self.dests[(subcmd,name)] = dest

        #TODO don't need this?
        if opt.has_fn:
            self.fns[(subcmd,arg.dest)] = opt.fn

    def _add_mutex_opts(self, subcmd, parser, mutex):
        group = parser.add_mutually_exclusive_group()
        for opt in mutex.opts:
            self._add_opt(subcmd, group, opt)

    def parse_args(self, argv=sys.argv[1:]):
        self.config.args = " ".join(argv)

        args = self.parser.parse_args(argv)

        module, cmd = getattr(args, "_cmd", (None,None))

        fns =  getattr(args, "_fns", None)
        if fns is not None:
            for fn in fns:
                getattr(self.config, fn)()

        #TODO use functions with parameters for special opts
        config_toml = getattr(args, CONFIG_PARAM, None)
        if config_toml is not None:
            self.config.load_toml(config_toml)

        for name, value in vars(args).items():

            if not name.startswith("_") or name in SPECIAL_PARAMS:
                if (cmd,name) in self.dests:
                    group, param = self.dests[(cmd,name)]
                else: 
                    group = self.config
                    param = name

                if value is not None or not hasattr(group, param):
                    setattr(group, param, value)
                #if value is not None or not hasattr(self.config, name):
                #    setattr(self.config, name, value)

        #fast5s = getattr(args, FAST5_PARAM, None)
        #if fast5s is not None:
        #    self.config.read_index.paths = unc.fast5.parse_fast5_paths(fast5s, self.config.read_index.recursive)

        return module, cmd, self.config
    
    def print_help(self):
        self.parser.print_help()

class Opt:
    def __init__(self, args, group_name=None, param=None, fn=None, help_suffix=None, **kwargs):
        self.args = args if type(args) == tuple else (args,)

        self.group_name = group_name
        self.param = param
        self.fn = fn
        self.help_suffix = help_suffix
        
        self.has_fn = fn is not None
        self.has_dest = group_name is not None

        if self.has_fn and self.has_dest:
            sys.stderr.write("ArgParser Opt fn and dest cannot both be set\n")
            sys.exit(1)

        if self.param is None:
            for arg in self.args:
                if arg.startswith("--") or not arg.startswith("-"):
                    self.param = arg.strip("-").replace("-","_")
                    break
            if self.param is None:
                sys.stderr.write("Error: must specify parameter name for flag \"%s\"\n" % self.args)
                sys.exit(0)

        self.extra_kw = kwargs

    def get_dest(self, config):
        if self.group_name is None:
            return (config, self.param)
        return (config.get_group(self.group_name), self.param)

    def get_kwargs(self, config):
        if self.has_fn:
            kw = {'dest' : '_fns', 'action' : 'append_const', 'const' : self.fn}
            kw.update(self.extra_kw)
            return kw

        elif self.group_name is None:
            return self.extra_kw

        group = config.get_group(self.group_name)

        if self.param is not None:
            param_name = self.param
        else:
            param_name = flag_to_var(self.args[-1])

        if not hasattr(group, param_name):
            raise ValueError(f"Parameter {self.group_name}.{param_name} does not exist")

        default = getattr(group, param_name)

        if hasattr(group, "_types"):
            _type = getattr(group, "_types")[param_name]
        else:
            _type = type(default)

        doc = getattr(type(group), param_name).__doc__
        if self.help_suffix != None:
            doc = doc + self.help_suffix
        kwargs = {"help" : doc}

        if _type != bool:
            kwargs["type"] = _type
            kwargs["default"] = default

        kwargs.update(self.extra_kw)

        return kwargs

    def set_val(self, config, val):
        if self.group is None or len(self.group) == 0:
            group = config
        else:
            group = getattr(config, self.group)

        setattr(group, self.param, val)

class MutexOpts:
    def __init__(self, dest, opts, **kwargs):
        self.dest = dest
        self.opts = opts
        self.kwargs = kwargs

class Subcmd:
    def __init__(self, name, desc, opts, fn=None):
        self.name = name
        self.desc = desc
        self.opts = opts
        self.fn = fn

        self.has_opts = fn is not None

        if not self.has_opts and not all([type(o) is Subcmd for o in opts]):
            raise TypeError("Must specify function for subcommand \"%s\" unless all options are of type Subcmd" % name)
        elif self.has_opts and not all([type(o) in [Opt, MutexOpts] for o in opts]):
            raise TypeError("Cannot specify function for subcommand \"%s\" unless all options are of type Opt" % name)

