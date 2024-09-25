from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Iterator

import geopandas as gpd
import pandas as pd
from anndata import AnnData
from datatree import DataTree
from spatialdata import SpatialData
from spatialdata.models import SpatialElement
from spatialdata.transformations import Identity, get_transformation, set_transformation
from xarray import DataArray

from ._constants import SopaAttrs, SopaFiles, SopaKeys

log = logging.getLogger(__name__)


def get_boundaries(
    sdata: SpatialData, return_key: bool = False, warn: bool = False
) -> gpd.GeoDataFrame | tuple[str, gpd.GeoDataFrame] | None:
    """Gets the baysor boundaries or cellpose boundaries of a SpatialData object after running Sopa

    Args:
        sdata: A SpatialData object
        return_key: Whether to return the key of the shapes or not.
        warn: If `True`, prints a warning if no boundary is found. Else, raises an error.

    Returns:
        A `GeoDataFrame` containing the boundaries, or a tuple `(shapes_key, geo_df)`
    """
    VALID_BOUNDARIES = [
        SopaKeys.BAYSOR_BOUNDARIES,
        SopaKeys.COMSEG_BOUNDARIES,
        SopaKeys.CELLPOSE_BOUNDARIES,
    ]
    for shapes_key in VALID_BOUNDARIES:
        res = _try_get_boundaries(sdata, shapes_key, return_key)
        if res is not None:
            return res

    error_message = "sdata object has no valid segmentation boundary. Consider running Sopa segmentation first."

    if not warn:
        raise ValueError(error_message)

    log.warn(error_message)
    return (None, None) if return_key else None


def _try_get_boundaries(sdata: SpatialData, shapes_key: str, return_key: bool) -> gpd.GeoDataFrame | None:
    """Try to get a cell boundaries for a given `shapes_key`"""
    if shapes_key in sdata.shapes:
        return (shapes_key, sdata[shapes_key]) if return_key else sdata[shapes_key]


def get_intrinsic_cs(sdata: SpatialData, element: SpatialElement | str, name: str | None = None) -> str:
    """Gets the name of the intrinsic coordinate system of an element

    Args:
        sdata: A SpatialData object
        element: `SpatialElement`, or its key
        name: Name to provide to the intrinsic coordinate system if not existing. By default, uses the element id.

    Returns:
        Name of the intrinsic coordinate system
    """
    if name is None:
        name = f"_{element if isinstance(element, str) else id(element)}_intrinsic"

    if isinstance(element, str):
        element = sdata[element]

    for cs, transform in get_transformation(element, get_all=True).items():
        if isinstance(transform, Identity):
            return cs

    set_transformation(element, Identity(), name)
    return name


def to_intrinsic(sdata: SpatialData, element: SpatialElement | str, element_cs: SpatialElement | str) -> SpatialElement:
    """Transforms a `SpatialElement` into the intrinsic coordinate system of another `SpatialElement`

    Args:
        sdata: A SpatialData object
        element: `SpatialElement` to transform, or its key
        element_cs: `SpatialElement` of the target coordinate system, or its key

    Returns:
        The `SpatialElement` after transformation in the target coordinate system
    """
    if isinstance(element, str):
        element = sdata[element]
    cs = get_intrinsic_cs(sdata, element_cs)
    return sdata.transform_element_to_coordinate_system(element, cs)


def get_intensities(sdata: SpatialData) -> pd.DataFrame | None:
    """Gets the intensity dataframe of shape `n_obs x n_channels`"""
    assert SopaKeys.TABLE in sdata.tables, f"No '{SopaKeys.TABLE}' found in sdata.tables"

    adata = sdata.tables[SopaKeys.TABLE]

    if not adata.uns[SopaKeys.UNS_KEY][SopaKeys.UNS_HAS_INTENSITIES]:
        return None

    if adata.uns[SopaKeys.UNS_KEY][SopaKeys.UNS_HAS_TRANSCRIPTS]:
        return adata.obsm[SopaKeys.INTENSITIES_OBSM]

    return adata.to_df()


def iter_scales(image: DataTree) -> Iterator[DataArray]:
    """Iterates through all the scales of a `DataTree`

    Args:
        image: a `DataTree`

    Yields:
        Each scale (as a `DataArray`)
    """
    assert isinstance(image, DataTree), f"Multiscale iteration is reserved for type DataTree. Found {type(image)}"

    for scale in image:
        yield next(iter(image[scale].values()))


def get_spatial_element(
    element_dict: dict[str, SpatialElement],
    key: str | None = None,
    valid_attr: str | None = None,
    return_key: bool = False,
    as_spatial_image: bool = False,
) -> SpatialElement | tuple[str, SpatialElement]:
    """Gets an element from a SpatialData object.

    Args:
        sdata: SpatialData object.
        key: Optional element key. If `None`, returns the only element (if only one), or tries to find an element with `valid_attr`.
        return_key: Whether to also return the key of the element.
        valid_attr: Attribute that the element must have to be considered valid.
        as_spatial_image: Whether to return the element as a `SpatialImage` (if it is a `DataTree`)

    Returns:
        If `return_key` is False, only the element is returned, else a tuple `(element_key, element)`
    """
    assert len(element_dict), "No spatial element was found in the dict."

    if key is not None:
        return _return_element(element_dict, key, return_key, as_spatial_image)

    if len(element_dict) == 1:
        key = next(iter(element_dict.keys()))

        assert valid_attr is None or _get_spatialdata_attrs(element_dict[key]).get(
            valid_attr, True
        ), f"Element {key} is not valid for the attribute {valid_attr}."

        return _return_element(element_dict, key, return_key, as_spatial_image)

    assert valid_attr is not None, "Multiple elements found. Provide an element key."

    keys = [key for key, element in element_dict.items() if _get_spatialdata_attrs(element).get(valid_attr)]

    assert len(keys) > 0, f"No element with the attribute {valid_attr}. Provide an element key."
    assert len(keys) == 1, f"Multiple valid elements found: {keys}. Provide an element key."

    return _return_element(element_dict, keys[0], return_key, as_spatial_image)


def _get_spatialdata_attrs(element: SpatialElement) -> dict[str, Any]:
    if isinstance(element, DataTree):
        element = next(iter(element["scale0"].values()))
    return element.attrs.get("spatialdata_attrs", {})


def _update_spatialdata_attrs(element: SpatialElement, attrs: dict):
    if isinstance(element, DataTree):
        for image_scale in iter_scales(element):
            _update_spatialdata_attrs(image_scale, attrs)
        return

    old_attrs = element.uns if isinstance(element, AnnData) else element.attrs

    if "spatialdata_attrs" not in old_attrs:
        old_attrs["spatialdata_attrs"] = {}

    old_attrs["spatialdata_attrs"].update(attrs)


def get_spatial_image(
    sdata: SpatialData,
    key: str | None = None,
    return_key: bool = False,
    valid_attr: str = SopaAttrs.CELL_SEGMENTATION,
) -> DataArray | tuple[str, DataArray]:
    """Gets a DataArray from a SpatialData object (if the image has multiple scale, the `scale0` is returned)

    Args:
        sdata: SpatialData object.
        key: Optional image key. If `None`, returns the only image (if only one), or tries to find an image with `valid_attr`.
        return_key: Whether to also return the key of the image.
        valid_attr: Attribute that the image must have to be considered valid.

    Returns:
        If `return_key` is False, only the image is returned, else a tuple `(image_key, image)`
    """
    return get_spatial_element(
        sdata.images,
        key=key,
        valid_attr=valid_attr,
        return_key=return_key,
        as_spatial_image=True,
    )


def _return_element(
    element_dict: dict[str, SpatialElement], key: str, return_key: bool, as_spatial_image: bool
) -> SpatialElement | tuple[str, SpatialElement]:
    element = element_dict[key]

    if as_spatial_image and isinstance(element, DataTree):
        element = next(iter(element["scale0"].values()))

    return (key, element) if return_key else element


def get_cache_dir(sdata: SpatialData) -> Path:
    assert sdata.is_backed(), "SpatialData not saved on-disk. Save the object, or provide a cache directory."

    return sdata.path / SopaFiles.SOPA_CACHE_DIR
