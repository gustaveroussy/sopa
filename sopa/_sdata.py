import logging

import geopandas as gpd
import xarray as xr
from multiscale_spatial_image import MultiscaleSpatialImage
from spatialdata import SpatialData
from spatialdata.transformations import Identity, get_transformation, set_transformation

from ._constants import SopaKeys

log = logging.getLogger(__name__)


def _try_get_boundaries(
    sdata: SpatialData, shapes_key: str, return_key: bool
) -> gpd.GeoDataFrame | None:
    if hasattr(sdata.shapes, shapes_key):
        return (shapes_key, sdata[shapes_key]) if return_key else sdata[shapes_key]


def get_boundaries(sdata: SpatialData, return_key: bool = False) -> gpd.GeoDataFrame | None:
    for shapes_key in [SopaKeys.BAYSOR_BOUNDARIES, SopaKeys.CELLPOSE_BOUNDARIES]:
        res = _try_get_boundaries(sdata, shapes_key, return_key)
        if res is not None:
            return res

    log.warn("sdata object has no cellpose boundaries and no baysor boundaries")


def get_intrinsic_cs(sdata: SpatialData, element_name: str) -> str:
    for cs, transform in get_transformation(sdata[element_name], get_all=True).items():
        if isinstance(transform, Identity):
            return cs

    cs = f"_{element_name}_intrinsic"
    set_transformation(sdata[element_name], Identity(), cs)
    return cs


def to_intrinsic(sdata: SpatialData, element, element_name_cs: str):
    if isinstance(element, str):
        element = sdata[element]
    cs = get_intrinsic_cs(sdata, element_name_cs)
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


def get_spatial_image(sdata: SpatialData, key: str | None = None) -> tuple[str, xr.DataArray]:
    key = get_key(sdata, "images", key)

    assert key is not None, "One image in `sdata.images` is required"

    image = sdata.images[key]  # TODO: switch axes for c,y,x

    if isinstance(image, MultiscaleSpatialImage):
        return key, image["scale0"][key]

    return key, image
