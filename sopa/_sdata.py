import logging

import geopandas as gpd
import pandas as pd
from multiscale_spatial_image import MultiscaleSpatialImage
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import SpatialElement
from spatialdata.transformations import Identity, get_transformation, set_transformation

from ._constants import SopaKeys

log = logging.getLogger(__name__)


def _try_get_boundaries(
    sdata: SpatialData, shapes_key: str, return_key: bool
) -> gpd.GeoDataFrame | None:
    if shapes_key in sdata.shapes:
        return (shapes_key, sdata[shapes_key]) if return_key else sdata[shapes_key]


def get_boundaries(
    sdata: SpatialData, return_key: bool = False, warn: bool = False
) -> gpd.GeoDataFrame | tuple[str, gpd.GeoDataFrame] | None:
    for shapes_key in [SopaKeys.BAYSOR_BOUNDARIES, SopaKeys.CELLPOSE_BOUNDARIES]:
        res = _try_get_boundaries(sdata, shapes_key, return_key)
        if res is not None:
            return res

    error_message = "sdata object has no cellpose boundaries and no baysor boundaries. Consider running segmentation first."

    if not warn:
        raise ValueError(error_message)

    log.warn(error_message)
    return (None, None) if return_key else None


def get_intrinsic_cs(
    sdata: SpatialData, element: SpatialElement | str, name: str | None = None
) -> str:
    if name is None:
        name = f"_{element if isinstance(element, str) else id(element)}_intrinsic"

    if isinstance(element, str):
        element = sdata[element]

    for cs, transform in get_transformation(element, get_all=True).items():
        if isinstance(transform, Identity):
            return cs

    set_transformation(element, Identity(), name)
    return name


def to_intrinsic(
    sdata: SpatialData, element: SpatialElement | str, element_cs: SpatialElement | str
):
    if isinstance(element, str):
        element = sdata[element]
    cs = get_intrinsic_cs(sdata, element_cs)
    return sdata.transform_element_to_coordinate_system(element, cs)


def get_key(sdata: SpatialData, attr: str, key: str | None = None):
    if key is not None:
        return key

    elements = getattr(sdata, attr)

    if not len(elements):
        return None

    assert (
        len(elements) == 1
    ), f"Trying to get an element key of `sdata.{attr}`, but it contains multiple values and no dict key was provided"

    return next(iter(elements.keys()))


def get_element(sdata: SpatialData, attr: str, key: str | None = None):
    key = get_key(sdata, attr, key)
    return sdata[key] if key is not None else None


def get_item(sdata: SpatialData, attr: str, key: str | None = None):
    key = get_key(sdata, attr, key)
    return key, sdata[key] if key is not None else None


def get_intensities(sdata: SpatialData) -> pd.DataFrame | None:
    if not sdata.table.uns[SopaKeys.UNS_KEY][SopaKeys.UNS_HAS_INTENSITIES]:
        return None

    if sdata.table.uns[SopaKeys.UNS_KEY][SopaKeys.UNS_HAS_TRANSCRIPTS]:
        return sdata.table.obsm[SopaKeys.INTENSITIES_OBSM]

    return sdata.table.to_df()


def get_spatial_image(
    sdata: SpatialData, key: str | None = None, return_key: bool = False
) -> SpatialImage | tuple[str, SpatialImage]:
    key = get_key(sdata, "images", key)

    assert key is not None, "One image in `sdata.images` is required"

    image = sdata.images[key]
    if isinstance(image, MultiscaleSpatialImage):
        image = SpatialImage(next(iter(image["scale0"].values())))

    if return_key:
        return key, image
    return image
