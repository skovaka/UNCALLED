
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

from .rt import ReadUntilClient

from .config import Config
from .argparse import ArgParser
from .fast5 import Fast5Reader
from .index import RefIndex

from .pore_model import PoreModel
from .sigproc import ProcRead

from . import dtw
from . import vis
#from . import 

config._DEFAULTS = Config()
config.rc = Config()


#SUBCMDS = [
#    dtw, dotplot, sigplot, browser, convert, 
#    (compare, COMPARE_OPTS), 
#    (method_compare, METHOD_COMPARE_OPTS),
#    (refstats, REFSTATS_OPTS),
#]
