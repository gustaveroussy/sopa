from . import _aggregation, _stainings, shapes
from .resolve import solve_conflicts, combine
from ._aggregation import Aggregator, overlay_segmentation  # deprecated in this module
from ._stainings import StainingSegmentation
from ._tissue import tissue, shapes_bounding_box
from .methods import cellpose, baysor, custom_staining_based, comseg, stardist, proseg
from . import methods


class Patches2D:
    def __init__(self, *args, **kwargs):
        raise NameError(
            "sopa.segmentation.Patches2D is deprecated and will be removed in sopa==2.1.0. Use `sopa.make_image_patches` or `sopa.make_transcript_patches` instead. For instance:\n"
            "    - sopa.make_image_patches(sdata, patch_width=1200, patch_overlap=50)\n"
            "    - sopa.make_transcript_patches(sdata, patch_width=1200, patch_overlap=50)\n"
            "See the migration guide to sopa 2.0.0 for more details: https://github.com/gustaveroussy/sopa/discussions/138"
        )
