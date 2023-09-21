import xarray as xr
from multiscale_spatial_image import MultiscaleSpatialImage
from spatialdata import SpatialData


def _get_key(sdata: SpatialData, attr: str, key: str | None = None):
    if key is not None:
        return key

    elements = getattr(sdata, attr)

    if not len(elements):
        return None

    assert (
        len(elements) == 1
    ), f"Trying to get an element key of sdata.{attr}, but it contains multiple values and no dict key was provided"

    return next(iter(elements.keys()))


def _get_element(sdata: SpatialData, attr: str, key: str | None = None):
    key = _get_key(sdata, attr, key)
    return getattr(sdata, attr)[key]


def _get_spatial_image(sdata: SpatialData, key: str | None = None) -> tuple[str, xr.DataArray]:
    key = _get_key(sdata, "images", key)

    image = sdata.images[key]  # TODO: switch axes for c,y,x

    if isinstance(image, MultiscaleSpatialImage):
        return key, image["scale0"]

    return key, image
