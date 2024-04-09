# Readers for multiplex-imaging technologies
# In the future, we will completely rely on spatialdata-io (when all these functions exist)

from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd
from spatialdata import SpatialData

from .utils import _deduplicate_names, _general_tif_directory_reader

log = logging.getLogger(__name__)


def macsima(path: Path, **kwargs: int) -> SpatialData:
    """Read MACSIMA data as a `SpatialData` object

    Notes:
        For all dulicated name, their index will be added in brackets after, for instance you will often find `DAPI (000)` to indicate the DAPI channel of index `000`

    Args:
        path: Path to the directory containing the MACSIMA `.tif` images
        kwargs: Kwargs for `_general_tif_directory_reader`

    Returns:
        A `SpatialData` object with a 2D-image of shape `(C, Y, X)`
    """
    return _general_tif_directory_reader(
        path, files_to_channels=_get_channel_names_macsima, **kwargs
    )


def _parse_name_macsima(file):
    index = file.name[2:5] if file.name[0] == "C" else file.name[:3]
    match = re.search(r"_A-(.*?)_C-", file.name)
    if match:
        antibody = match.group(1)
        channel = re.search(r"_C-(.*?)\.tif", file.name).group(1)
        uid = f"{channel}-{index}"
    else:
        antibody = re.search(r"_A-(.*?)\.tif", file.name).group(1)
        uid = index
    return [antibody, uid]


def _get_channel_names_macsima(files):
    df_antibodies = pd.DataFrame([_parse_name_macsima(file) for file in files])
    return _deduplicate_names(df_antibodies)
