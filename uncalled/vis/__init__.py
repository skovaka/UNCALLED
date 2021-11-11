#from .mpl import browser
from . import trackplot, dotplot, sigplot, browser
from .trackplot import Trackplot
from .dotplot import Dotplot
from .. import config

class VisParams(config.ParamGroup):
    _name = "vis"
VisParams._def_params(
    ("track_colors", ["#AA0DFE", "#1CA71C", "#4676FF", "#d90000"], list, "Track Colors"),
    ("base_colors", ["#80ff80", "#6b93ff", "#ffe543", "#ff8080"], list, "Colors for each base (A,C,G,T/U)"), 
)
