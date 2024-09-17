from __future__ import annotations

import logging
from pathlib import Path

from spatialdata import SpatialData

from .utils import _general_tif_directory_reader

log = logging.getLogger(__name__)


def macsima(path: Path, **kwargs: int) -> SpatialData:
    """Read MACSIMA data as a `SpatialData` object

    Notes:
        For all dulicated name, their index will be added in brackets after, for instance you may find `DAPI (1)`.

    Args:
        path: Path to the directory containing the MACSIMA `.tif` images
        kwargs: Kwargs for the `_general_tif_directory_reader`

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    return _general_tif_directory_reader(path, **kwargs)
