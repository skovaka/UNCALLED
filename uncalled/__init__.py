from .__about__ import (
    __title__,
    __version__, 
    __summary__, 
    __author__, 
    __email__,  
    __uri__,   
    __license__,  
    __copyright__,
)

from _uncalled import *

from _uncalled import _RefCoord

def str_to_coord(coord_str):
    spl = coord_str.split(":")         
    #name = spl[0]
    i = 1
    while True:
        try:
            st,en = map(int, spl[i].split("-"))
            break
        except ValueError:
            i += 1
    name = ":".join(spl[:i])
    
    if i < len(spl):
        return RefCoord(name,st,en)
    else:                              
        fwd = spl[i] == "+"
        return RefCoord(name,st,en,fwd)

class RefCoord(_RefCoord):
    def __init__(self, *args):
        if len(args) == 1:
            args = args[0]
            if isinstance(args, str):
                args = self._str_to_tuple(args)
            if not isinstance(args, tuple):
                raise ValueError(f"Invalid RefCoord: {args}")
        _RefCoord.__init__(self, *args)

    def _str_to_tuple(self, coord_str):
        spl = coord_str.split(":")         
        name = spl[0]                        
        st,en = map(int, spl[1].split("-"))
        
        if len(spl) == 2:                  
            return (name,st,en)
        elif len(spl) == 3 and spl[2] in {"+","-"}:
            fwd = spl[2] == "+"
            return (name,st,en,fwd)
        return None

from . import params

from .config import Config

config._DEFAULTS = Config()
config.rc = Config()

from .pore_model import PoreModel
from .index import load_index
