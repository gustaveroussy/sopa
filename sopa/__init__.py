import importlib.metadata
import logging
import sys
from ._logging import configure_logger

__version__ = importlib.metadata.version("sopa")

log = logging.getLogger("sopa")
configure_logger(log)

if not any(f"--{option}" in sys.argv for option in ["version", "help"]):  # no import for cli helpers
    from spatialdata import read_zarr  # will set `dataframe.query-planning` to False

    from ._settings import settings
    from . import utils
    from . import shapes
    from . import io
    from . import spatial
    from . import segmentation
    from .aggregation import aggregate, overlay_segmentation
    from .patches import make_transcript_patches, make_image_patches
    from .utils import get_spatial_image, get_spatial_element, to_intrinsic, get_boundaries
