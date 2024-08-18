import importlib.metadata
import logging
import sys
from ._logging import configure_logger

__version__ = importlib.metadata.version("sopa")

log = logging.getLogger("sopa")
configure_logger(log)

if "--help" not in sys.argv:
    from . import utils
    from . import io
    from . import segmentation

    from .segmentation import tissue_segmentation
    from ._sdata import get_spatial_image, get_spatial_element, to_intrinsic
