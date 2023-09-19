import xarray as xr
from multiscale_spatial_image import MultiscaleSpatialImage
from spatialdata import SpatialData


def _get_key(sdata: SpatialData, attr: str, key: str | None = None):
    spatial_elements = getattr(sdata, attr)

    if key is not None:
        return key

    if not len(spatial_elements):
        return None

    assert (
        len(spatial_elements) == 1
    ), f"Trying to get an element of sdata.{attr}, but it contains multiple values and no dict key was provided"

    return next(iter(spatial_elements.keys()))


def _get_value(sdata: SpatialData, attr: str, key: str | None = None):
    spatial_elements = getattr(sdata, attr)

    if key is not None:
        return spatial_elements[key]

    if not len(spatial_elements):
        return None

    assert (
        len(spatial_elements) == 1
    ), f"Trying to get an element of sdata.{attr}, but it contains multiple values and no dict key was provided"

    return next(iter(spatial_elements.values()))


def _get_spatial_image(sdata: SpatialData, key: str | None = None) -> tuple[str, xr.DataArray]:
    if key is None:
        key = _get_key(sdata, "images")

    image = sdata.images[key]

    if isinstance(image, MultiscaleSpatialImage):
        return key, image["scale0"]

    return key, image
