from . import aggregation, shapes, stainings
from ..patches import Patches2D  # TODO: remove import in sopa>=2.0.0
from .aggregation import Aggregator, overlay_segmentation
from .stainings import StainingSegmentation
from .tissue import tissue_segmentation
from .methods import cellpose, baysor
from . import methods
