import logging
import re
from pathlib import Path

from spatialdata import SpatialData

from .utils import _deduplicate_names, _general_tif_directory_reader

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
    files = list(Path(path).glob("*.tif"))

    if any("A-" in file.name for file in files):  # non-ome.tif format
        return _general_tif_directory_reader(path, files_to_channels=_get_channel_names_macsima, **kwargs)

    return _general_tif_directory_reader(path, **kwargs)


def _parse_name_macsima(file):
    match = re.search(r"_A-(.*?)_C-", file.name)
    if match:
        antibody = match.group(1)
    else:
        antibody = re.search(r"_A-(.*?)\.tif", file.name).group(1)
    return antibody


def _get_channel_names_macsima(files):
    return _deduplicate_names([_parse_name_macsima(file) for file in files])
