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
import copy
from _uncalled import _Conf
from collections import namedtuple

from . import __title__, __version__, __summary__

FAST5_PARAM = "fast5_files"
CONFIG_PARAM = "config_toml"
SPECIAL_PARAMS = {FAST5_PARAM, CONFIG_PARAM}

TOML_TYPES = {int, float, str, list, bool}

#TODO make this a factory function or whatever
Param = namedtuple("Param", ["name", "default", "type", "doc"])
class ParamGroup:
    _name = None

    def __init__(self):
        self._values = dict()

    def set(self, name, val):
        """Sets parameter with automatic type conversion"""
        if name not in self._types:
            raise KeyError(name)

        type_ = self._types[name]
        if not (type_ is None or isinstance(val, type_)):
            val = type_(val)

        self._values[name] = val
    
    @property
    def count(self):
        return len(self._types)

    @classmethod
    def _def_params(_class, *params, ignore_toml={}):
        _class._order = list()
        _class._types = dict()

        for p in params:
            _class._def_param_property(p)

        _class._ignore_toml = ignore_toml

        Config._EXTRA_GROUPS[_class._name] = _class 

    @classmethod
    #def _prm(_class, name, default, type, docstr):
    def _def_param_property(_class, p):
        if type(p) != Param:
            p = Param._make(p)

        def getter(self):
            return self._values.get(p.name, p.default)

        def setter(self, value):
            type_ = self._types[p.name]
            if not (type_ is None or isinstance(value, type_)):
                value = type_(value)
            self._values[p.name] = value

        setattr(_class, p.name, property(getter, setter, doc=p.doc))
        _class._types[p.name] = p.type
        _class._order.append(p.name)

class Config(_Conf):
    _EXTRA_GROUPS = dict()

    def __init__(self, conf=None, toml=None):
        _Conf.__init__(self)

        for name,cls in self._EXTRA_GROUPS.items():
            setattr(self, name, cls())

        if conf is not None:
            self.load_config(conf)
        
        if toml is not None:
            self.load_toml(toml)


    def _param_writable(self, name, val, group=None):
        return (not self.is_default(name, group) and
                not name.startswith("_") and
                type(val) in TOML_TYPES and
                (not hasattr(val, "__len__") or len(val) > 0))

    def load_config(self, other, ignore_defaults=True):
        for param in self._GLOBAL_PARAMS:
            if not (ignore_defaults and other.is_default(param)):
                setattr(self, param, copy.copy(getattr(other, param)))

        for param, val in vars(other).items():
            if not (isinstance(val, ParamGroup) or (ignore_defaults and other.is_default(param))):
                setattr(self, param, copy.copy(val))

        groups = other._PARAM_GROUPS + list(self._EXTRA_GROUPS.keys())

        for group_name in groups:
            ogroup = getattr(other, group_name)
            sgroup = getattr(self, group_name)
            for param in dir(ogroup):
                if not (param.startswith("_") or (ignore_defaults and other.is_default(param, group_name))):
                    setattr(sgroup, param, getattr(ogroup, param))

    def to_toml(self, filename=None):

        def fmt(val):
            if isinstance(val, str) and os.path.exists(val):
                return os.path.abspath(val)
            return val
        
        out = dict()

        for param in self._GLOBAL_PARAMS:
            val = getattr(self,param)
            if self._param_writable(param, val):
                out[param] = fmt(val)

        for param, val in vars(self).items():
            if self._param_writable(param, val):
                out[param] = fmt(val)

        groups = self._PARAM_GROUPS + list(self._EXTRA_GROUPS.keys())

        for group_name in groups:
            group = getattr(self, group_name)
            vals = dict()
            for param in dir(group):
                val = getattr(group,param)
                if self._param_writable(param, val, group_name):
                    vals[param] = fmt(val)
            if len(vals) > 0:
                out[group_name] = fmt(vals)

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
                    if self.is_default(name):
                        setattr(self, name, val)
                else:
                    for param, value in val.items():
                        if not hasattr(group, param):
                            raise ValueError("Unrecognized parameter in TOML: %s.%s" % (group, param))
                        if self.is_default(param, name):
                            setattr(group, param, value)
            else:
                if self.is_default(name):
                    setattr(self, name, val)

    def get_group(self, group_name):
        """Returns the group based on the group name"""

        if group_name is None or len(group_name) == 0:
            return self

        return getattr(self, group_name)
    
    @property
    def is_rna(self):
        return not self.read_buffer.seq_fwd

    def is_default(self, param, group=None):
        if group is None:
            sg = self
            dg = _DEFAULTS

            if not hasattr(dg, param) and hasattr(sg, param):
                return False

        else:
            sg = getattr(self, group, None)
            dg = getattr(_DEFAULTS, group, None)

        if sg is None and dg is None:
            return True
        elif sg is None or dg is None:
            return False


        return getattr(sg, param, None) == getattr(dg, param, None)

CONF_KW = "conf"

def _init_group(name, *args, **kwargs):
    conf = Config(kwargs.get(CONF_KW, rc))

    if not hasattr(conf, name):
        raise ValueError("Invalid parameter group: " + str(name))

    params = getattr(conf, name)

    if len(args) > params.count:
        raise ValueError("Too many arguments for " + name)

    arg_params = set()
    
    for i,val in enumerate(args):
        arg = params._order[i]
        arg_params.add(arg)
        setattr(params, arg, val)

    for arg, val in kwargs.items():
        if arg in arg_params:
            raise ValueError("Conflicting *arg and **kwarg values for %s.%s" % (name, param))

        if hasattr(params, arg):
            setattr(params, arg, val)

        elif arg != CONF_KW:
            raise ValueError("Unknown kwarg \"%s\" for %s parameters" % (arg, name))

    return conf, params

class ArgParser:

    def __init__(self, 
            subcmds=None, 
            desc="Rapidly maps raw nanopore signal to DNA references",
            config=None):

        if config is None:
            self.config = unc.Config()
        else:
            self.config = config

        self.dests = dict()
        self.fns = dict()

        self.parser = argparse.ArgumentParser(
                description=desc, 
                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                prog=__title__
        )

        self.parser.add_argument("-v", "--version", action="version", version=__version__)

        self._add_subcmds(self.parser, subcmds)

    def _add_subcmds(self, parser, subcmds):
        subparsers = parser.add_subparsers(title="subcommands")
        for subcmd in subcmds:

            if isinstance(subcmd, tuple):
                subcmd, opts = subcmd
                main_func = subcmd
            else:
                opts = getattr(subcmd, "OPTS", None)
                main_func = getattr(subcmd, "main", None)

            subcmd_name = subcmd.__name__.split(".")[-1]

            if main_func is not None:
                desc = main_func.__doc__
            else:
                desc = subcmd.__doc__

            sp = subparsers.add_parser(
                subcmd_name, help=desc, 
                formatter_class=argparse.ArgumentDefaultsHelpFormatter
            )

            if main_func is not None:
                sp.set_defaults(_cmd=(subcmd_name, main_func))
                for opt in opts:
                    if type(opt) is Opt:
                        self._add_opt(subcmd_name, sp, opt)
                    elif type(opt) is MutexOpts:
                        self._add_mutex_opts(subcmd_name, sp, opt)

            elif hasattr(subcmd, "SUBCMDS"):
                self._add_subcmds(sp, subcmd.SUBCMDS)

            else:
                raise RuntimeError("Subcommand module \"%s\" does not contain \"main\" function or \"SUBCMDS\" list" % subcmd.__name__)

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

        cmd_name, cmd = getattr(args, "_cmd", None)

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
                if (cmd_name,name) in self.dests:
                    group, param = self.dests[(cmd_name,name)]
                else: 
                    group = self.config
                    param = name

                if value is not None or not hasattr(group, param):
                    setattr(group, param, value)
                #if value is not None or not hasattr(self.config, name):
                #    setattr(self.config, name, value)

        fast5s = getattr(args, FAST5_PARAM, None)
        if fast5s is not None:
            self.config.fast5_reader.fast5_files = unc.fast5.parse_fast5_paths(fast5s, self.config.fast5_reader.recursive)

        return cmd, self.config
    
    def print_help(self):
        self.parser.print_help()

class Opt:
    def __init__(self, args, group_name=None, param=None, fn=None, **kwargs):
        self.args = args if type(args) == tuple else (args,)

        self.group_name = group_name
        self.param = param
        self.fn = fn
        
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

