from . import _aggregation, _stainings, shapes
from ._aggregation import Aggregator, overlay_segmentation  # deprecated
from ._stainings import StainingSegmentation
from ._tissue import tissue
from .methods import cellpose, baysor, custom_staining_based, comseg
from . import methods
