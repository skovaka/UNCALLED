from time import time

t = time()

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
print("a", time() - t)
t = time()

from _uncalled import *
print("b", time() - t)
t = time()

from . import params

from .config import Config

print("c", time() - t)
t = time()

from .argparse import ArgParser

print("d", time() - t)
t = time()

from .fast5 import Fast5Reader

print("e", time() - t)
t = time()

from .index import RefIndex
print("f", time() - t)
t = time()

from .pore_model import PoreModel
print("g", time() - t)
t = time()

#from . import dtw

print("h", time() - t)
t = time()

#from . import vis
print("i", time() - t)
t = time()
#from . import stats 
print("j", time() - t)
t = time()

config._DEFAULTS = Config()
config.rc = Config()
print("k", time() - t)
t = time()


#SUBCMDS = [
#    dtw, dotplot, sigplot, browser, convert, 
#    (compare, COMPARE_OPTS), 
#    (method_compare, METHOD_COMPARE_OPTS),
#    (refstats, REFSTATS_OPTS),
#]
