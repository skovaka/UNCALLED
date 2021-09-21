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
            self.load_toml(text=toml)


    def _param_writable(self, name, val, group=None):
        if group is not None: 
            if name in getattr(getattr(self, group), "_ignore_toml", {}):
                return False
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

    def load_toml(self, filename=None, text=None):
        if filename is not None:
            toml_dict = toml.load(filename)
        elif text is not None:
            toml_dict = toml.loads(text)

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

