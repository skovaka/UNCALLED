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
import toml
import copy
from collections import namedtuple

from _uncalled import _Conf

TOML_TYPES = {int, float, str, list, bool}

#TODO make this a factory function or whatever
Param = namedtuple("Param", ["name", "default", "type", "doc"])
class ParamGroup:
    _name = None

    def __init__(self, init=None):
        self._values = dict()

    def _set(self, name, val):
        """Sets parameter with automatic type conversion"""
        if name not in self._types:
            raise KeyError(name)

        self._values[name] = self._convert(name, val)

    def _convert(self, name, val):
        type_ = self._types[name]
        if not type_ is None:
            if issubclass(type_, ParamGroup) and val is None:
                val = type_()

            elif val is not None and not isinstance(val, type_):
                val = type_(val)

        return val
    
    @property
    def count(self):
        return len(self._types)

    @classmethod
    def _def_params(_class, *params, ignore_toml={}, config_add=True):
        _class._order = list()
        _class._types = dict()

        for p in params:
            _class._def_param_property(p)

        _class._ignore_toml = ignore_toml

        if config_add:
            Config._EXTRA_GROUPS[_class._name] = _class 

    @classmethod
    #def _prm(_class, name, default, type, docstr):
    def _def_param_property(_class, p):
        if type(p) != Param:
            p = Param._make(p)

        _class._types[p.name] = p.type
        _class._order.append(p.name)

        def getter(self):
            if not p.name in self._values:
                self._values[p.name] = self._convert(p.name, p.default)
            return self._values[p.name]

        def setter(self, value):
            self._set(p.name, value)
            #type_ = self._types[p.name]

            #if not (type_ is None or isinstance(value, type_)):
            #    value = type_(value)

            #self._values[p.name] = value

        setattr(_class, p.name, property(getter, setter, doc=p.doc))

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


    def _param_writable(self, name, val, group=None, ignore=True):
        if group is not None: 
            if ignore and (name in getattr(self.get_group(group), "_ignore_toml", {}) 
                or name == "read_filter"): #TODO not great - should wrap Fast5Params in ParamGroup
                return False
        #print(group, name, self.is_default(name, group), not hasattr(val, "__len__") or len(val) > 0)
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

        def _load_group(group):
            ogroup = other.get_group(group)
            sgroup = self.get_group(group)
            for param in dir(ogroup):
                if not (param.startswith("_") or (ignore_defaults and other.is_default(param, group))):
                    val = getattr(ogroup, param)
                    if isinstance(val, ParamGroup):
                        _load_group(".".join([group,param]))
                    else:
                        setattr(sgroup, param, val)

        for group_name in groups:
            ogroup = getattr(other, group_name)
            sgroup = getattr(self, group_name)
            
            _load_group(group_name)

            #for param in dir(ogroup):
            #    if not (param.startswith("_") or (ignore_defaults and other.is_default(param, group_name))):
            #        setattr(sgroup, param, getattr(ogroup, param))


    def to_toml(self, filename=None, force_all=False):

        def fmt(val):
            if isinstance(val, str) and os.path.exists(val):
                return os.path.abspath(val)
            return val
        
        out = dict()

        for param in self._GLOBAL_PARAMS:
            val = getattr(self,param)
            if self._param_writable(param, val, ignore=not force_all):
                out[param] = fmt(val)

        for param, val in vars(self).items():
            if self._param_writable(param, val, ignore=not force_all):
                out[param] = fmt(val)

        groups = self._PARAM_GROUPS + list(self._EXTRA_GROUPS.keys())

        def write_group(group_name, out):
            group = self.get_group(group_name)

            vals = dict()
            for param in dir(group):
                val = getattr(group,param)
                if self._param_writable(param, val, group_name, ignore=not force_all):
                    vals[param] = fmt(val)
                elif isinstance(val, ParamGroup):
                    #out[param] = dict()
                    write_group(f"{group_name}.{param}", vals)

            if len(vals) > 0:
                out[group_name.split(".")[-1]] = fmt(vals)

        for group_name in groups:
            write_group(group_name, out)

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


        def load_group(name, vals):
            group = self.get_group(name)
            if group is None:
                if self.is_default(name):
                    setattr(self, name, vals)
            else:
                basename = name.split(".")[0]
                for param, value in vals.items():
                    if not hasattr(group, param):
                        sys.stderr.write(f"Unrecognized config parameter: \"{name}.{param}\". Skipping.\n")
                    if isinstance(value, dict):
                        load_group(f"{basename}.{param}", value)
                    elif self.is_default(param, name):
                        setattr(group, param, value)

        for name, val in toml_dict.items():
            if isinstance(val, dict):
                load_group(name, val)
            else:
                if self.is_default(name):
                    setattr(self, name, val)

    #support pickling via toml
    def __setstate__(self, toml):
        self.__init__(toml=toml)
        #self.load_toml(text=toml)

    def __getstate__(self):
        return self.to_toml(force_all=True)

    def get_group(self, group_name):
        """Returns the group based on the group name"""

        if group_name is None or len(group_name) == 0:
            return self
        
        if not isinstance(group_name, str):
            return group_name

        group = self
        for name in group_name.split("."):
            if not hasattr(group, name):
                raise ValueError("Unknown Config group: " + group_name)

            group = getattr(group, name)

        return group
    
    @property
    def is_rna(self):
        return not self.read_buffer.seq_fwd

    def is_default(self, param, group=None):
        sg = self.get_group(group)
        dg = _DEFAULTS.get_group(group)

        #if group is None:
        #    sg = self
        #    dg = _DEFAULTS

        #    if not hasattr(dg, param) and hasattr(sg, param):
        #        return False

        #else:
        #    sg = getattr(self, group, None)
        #    dg = getattr(_DEFAULTS, group, None)

        #if sg is None and dg is None:
        #    return True
        #elif sg is None or dg is None:
        #    return False

        return getattr(sg, param, None) == getattr(dg, param, None)

CONF_KW = "conf"

def _init_group(name, *args, copy_conf=True, **kwargs):
    if not copy_conf and CONF_KW in kwargs:
        conf = kwargs[CONF_KW]
    else:
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

