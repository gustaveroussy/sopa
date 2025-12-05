import pandas as pd
import pytest
from spatialdata.models import PointsModel
from xarray import DataArray, DataTree

import sopa
from sopa.constants import ATTRS_KEY, SopaAttrs


def test_make_toy_dataset():
    assert sopa.io.toy_dataset(length=512) is not None
    assert sopa.io.toy_dataset(as_output=True, length=512) is not None


def test_make_blobs():
    assert sopa.io.blobs(length=512) is not None


def test_cache():
    sdata = sopa.io.toy_dataset(length=100)
    cache_dir = sopa.utils.get_cache_dir(sdata)

    assert cache_dir.exists()
    assert cache_dir.is_dir()

    sopa.utils.delete_cache(sdata)
    assert not cache_dir.exists()

    sdata.write("_test_cache.zarr")

    new_cache_dir = sopa.utils.get_cache_dir(sdata)

    assert new_cache_dir.exists()
    assert new_cache_dir.is_dir()
    assert new_cache_dir != cache_dir

    sopa.utils.delete_cache(sdata)

    assert not new_cache_dir.exists()

    import shutil

    shutil.rmtree("_test_cache.zarr")


def test_sdata_attrs_points():
    sdata = sopa.io.toy_dataset(length=100, add_second_points_key=True)

    points_key, points = sopa.utils.get_spatial_element(
        sdata.points, key=sdata.attrs.get(SopaAttrs.TRANSCRIPTS), return_key=True
    )

    assert points_key == "transcripts"

    feature_key = sopa.utils.get_feature_key(points)

    assert feature_key == "genes"

    del sdata[points_key].attrs[ATTRS_KEY]["feature_key"]

    with pytest.raises(ValueError):
        sopa.utils.get_feature_key(sdata[points_key], raise_error=True)

    del sdata.points[points_key]

    with pytest.raises(AssertionError):
        sopa.utils.get_spatial_element(sdata.points, key=sdata.attrs.get(SopaAttrs.TRANSCRIPTS))

    sopa.utils.get_spatial_element(sdata.points)

    with pytest.raises(ValueError):
        sopa.utils.get_feature_key(sdata["misc"], raise_error=True)


def test_sdata_attrs_shapes():
    sdata = sopa.io.toy_dataset(length=100)

    with pytest.raises(ValueError):
        sopa.utils.get_boundaries(sdata)

    sdata["cellpose_boundaries"] = sdata.shapes["cells"]

    sopa.utils.get_boundaries(sdata)

    sdata["baysor_boundaries"] = sdata.shapes["cells"]

    shapes_key, _ = sopa.utils.get_boundaries(sdata, return_key=True)

    assert shapes_key == "baysor_boundaries"

    del sdata.shapes["baysor_boundaries"]
    del sdata.shapes["cellpose_boundaries"]

    sdata.attrs[SopaAttrs.BOUNDARIES] = "cells"

    shapes_key, _ = sopa.utils.get_boundaries(sdata, return_key=True)

    assert shapes_key == "cells"


def test_sdata_attrs_images():
    sdata = sopa.io.toy_dataset(length=100)

    image_key, image = sopa.utils.get_spatial_image(sdata, return_key=True)

    assert image_key == "image"
    assert isinstance(image, DataArray)

    image_key, image = sopa.utils.get_spatial_element(
        sdata.images, return_key=True, key=sdata.attrs[SopaAttrs.CELL_SEGMENTATION]
    )

    assert image_key == "image"
    assert isinstance(image, DataArray)

    with pytest.raises(AssertionError):
        sopa.utils.get_spatial_element(sdata.images)

    image_key, image = sopa.utils.get_spatial_image(sdata, return_key=True, valid_attr=SopaAttrs.TISSUE_SEGMENTATION)

    assert image_key == "he_image"
    assert isinstance(image, DataArray)

    image_key, image = sopa.utils.get_spatial_element(
        sdata.images, return_key=True, key=sdata.attrs[SopaAttrs.TISSUE_SEGMENTATION], as_spatial_image=False
    )

    assert image_key == "he_image"
    assert isinstance(image, DataTree)


def test_add_spatial_element():
    sdata = sopa.io.toy_dataset(length=100)

    points = PointsModel.parse(pd.DataFrame({"x": [1, 2], "y": [1, 2]}))

    sopa.utils.add_spatial_element(sdata, "points_test", points)

    assert "points_test" in sdata.points

    with pytest.raises(AssertionError):
        sopa.utils.add_spatial_element(sdata, "points_test", points, overwrite=False)

    sopa.utils.add_spatial_element(sdata, "points_test", points)

    del sdata.points["points_test"]

    sdata.write("_test_add_spatial_element.zarr")

    sopa.settings.auto_save_on_disk = False

    sopa.utils.add_spatial_element(sdata, "points_test", points)

    assert not (sdata.path / "points" / "points_test").is_dir()

    sopa.settings.auto_save_on_disk = True

    with pytest.raises(AssertionError):
        sopa.utils.add_spatial_element(sdata, "points_test", points, overwrite=False)

    sopa.utils.add_spatial_element(sdata, "points_test", points)

    assert (sdata.path / "points" / "points_test").is_dir()

    import shutil

    shutil.rmtree("_test_add_spatial_element.zarr")
