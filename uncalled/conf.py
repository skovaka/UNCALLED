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
