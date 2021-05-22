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
from _uncalled import _Conf
from collections import namedtuple

TOML_TYPES = {int, float, str, list, bool}

def param_valid(name, val):
    return (not name.startswith("_") and
            type(val) in TOML_TYPES and
            (not hasattr(val, "__len__") or len(val) > 0))

#TODO make this a factory function or whatever
Param = namedtuple("Param", ["name", "default", "type", "doc"])
class ParamGroup:
    _types = dict()
    def __init__(self):
        self._values = dict()

    def set(self, **kwargs):
        for arg, val in kwargs.items():
            if hasattr(self, arg):
                setattr(self, arg, val)
            else:
                raise RuntimeError("Unknown kwarg \"%s\" for %s" % (arg, type(self)))

    def from_kw(self, **kwargs):
        self.set(**kwargs)
        return self

    @classmethod
    def _def_params(_class, *params):
        for p in params:
            if type(p) != Param:
                p = Param._make(p)
            _class._def_param_property(p)

    @classmethod
    #def _prm(_class, name, default, type, docstr):
    def _def_param_property(_class, p):
        def getter(self):
            return self._values.get(p.name, p.default)

        def setter(self, value):
            self._values[p.name] = value

        setattr(_class, p.name, property(getter, setter, doc=p.doc))
        _class._types[p.name] = p.type


class Conf(_Conf):
    _EXTRA_GROUPS = dict()

    def __init__(self, toml=None):
        _Conf.__init__(self)

        for name,cls in self._EXTRA_GROUPS.items():
            setattr(self, name, cls())
        
        if toml is not None:
            self.load_toml(toml)

    def to_toml(self, filename=None):
        out = dict()

        for param in self._GLOBAL_PARAMS:
            val = getattr(self,param)
            if param_valid(param, val):
                out[param] = val

        for param, val in vars(self).items():
            if param_valid(param, val):
                out[param] = val

        groups = self._PARAM_GROUPS + list(self._EXTRA_GROUPS.keys())

        for group_name in groups:
            group = getattr(self, group_name)
            out[group_name] = dict()
            for param in dir(group):
                val = getattr(group,param)
                if param_valid(param, val):
                    out[group_name][param] = val

        if filename is None:
            return toml.dumps(out)
        else:
            with open(filename, "w") as fout:
                toml.dump(out, fout)

    def load_toml(self, filename):
        toml_dict = toml.load(filename)

        for name, val in toml_dict.items():
            if isinstance(val, dict):
                group = getattr(self, name, None)
                if group is None:
                    setattr(self, name, val)
                else:
                    for param, value in val.items():
                        setattr(group, param, value)
            else:
                setattr(self, name, val)

    def get_group(self, group_name):
        """Returns the group based on the group name"""

        if group_name is None or len(group_name) == 0:
            return self

        return getattr(self, group_name)

class ArgParser:
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


    def __init__(self, 
            argv=sys.argv[1:],
            subcmds=None, 
            conf=None,
            description="Rapidly maps raw nanopore signal to DNA references"):

        if conf is None:
            self.conf = unc.Conf()
        else:
            self.conf = conf

        self.parser = argparse.ArgumentParser(
                description=description, 
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
        self.conf.args = " ".join(argv)

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

