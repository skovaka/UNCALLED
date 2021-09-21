#from . import dtw, convert
from collections import namedtuple
LayerMeta = namedtuple("LayerMeta", ["type", "label"])

LAYER_META = {
    "ref"     : LayerMeta(int, "Reference Coordinate"),
    "start"   : LayerMeta(int, "Sample Start"),
    "length"  : LayerMeta(int, "Sample Length"),
    "current" : LayerMeta(float, "Mean Current (pA)"),
    "kmer"    : LayerMeta(int, "Reference K-mer"),
    "mref"  : LayerMeta(int, "Mirrored Packed Ref. Coord."),
    "aln_id"      : LayerMeta(int, "Alignment ID"),
}

class RefCoord:
    def __init__(self, name=None, start=None, end=None, fwd=None):
        self.fwd = fwd
        if start is None and end is None:
            if isinstance(name ,str):
                self._init_str(name)
            elif isinstance(name, tuple):
                self._init_tuple(name)

        elif start is None:
            raise ValueError("RefCoords must include a start coordinate")
        
        else:
            self.name = name
            self.start = start
            self.end = end

    def union(self, other):
        if self.name != other.name or max(self.start, other.start) > min(self.end, other.end):
            return None
        return RefCoord(self.name, min(self.start, other.start), max(self.end, other.end), self.fwd)

    def intersect(self, other):
        start = max(self.start, other.start)
        end = min(self.end, other.end)
        if self.name != other.name or start > end:
            return None
        return RefCoord(self.name, start, end, self.fwd)

    def _init_str(self, coord_str):
        spl = coord_str.split(":")
        self.name = spl[0]
        coords = tuple(map(int, spl[1].split("-")))

        if len(coords) == 2:
            self.start, self.end = coords
        elif len(coords) == 1:
            self.start, = coords
            self.end = None
        else:
            raise ValueError("RefCoords must contain one or two coordinate")

        if len(spl) == 3:
            self.fwd = spl[2] == "+"
        else:
            self.fwd = None

    def _init_tuple(self, coords):
        self.name = coords[0]
        self.start = coords[1]
        if len(coords) > 2:
            self.end = coords[2]
        if len(coords) > 3:
            self.end = coords[3]

    def __repr__(self):
        s = "%s:%d" % (self.name, self.start)
        if self.end is not None:
            s += "-%d" % self.end
        if self.fwd is not None:
            s += " (%s)" % ("+" if self.fwd else "-")
        return s

from .bcaln import Bcaln
from .track_io import TrackIO
from .dtw import Fast5Processor #TODO move this to main module (eventually sigproc, then to C++)
