#from . import dtw, convert
from collections import namedtuple

from .bcaln import Bcaln
from .tracks import Tracks
from .aln_track import LAYERS
from .dtw import Fast5Processor #TODO move this to main module (eventually sigproc, then to C++)
