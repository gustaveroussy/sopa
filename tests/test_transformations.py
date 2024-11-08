import numpy as np
import pandas as pd
import pytest
import spatialdata
from spatialdata import SpatialData
from spatialdata.models import PointsModel
from spatialdata.transformations import (
    Identity,
    Scale,
    get_transformation_between_coordinate_systems,
)

import sopa


@pytest.fixture
def sdata() -> SpatialData:
    sdata = sopa.io.toy_dataset(length=200)
    return sdata


def test_no_transform(sdata: SpatialData):
    image_key, image = sopa.utils.get_spatial_image(sdata, return_key=True)

    transformed_image = sopa.utils.to_intrinsic(sdata, image, image)
    assert transformed_image is image

    transformed_image = sopa.utils.to_intrinsic(sdata, image, image_key)
    assert transformed_image is image

    transformed_image = sopa.utils.to_intrinsic(sdata, image_key, image)
    assert transformed_image is image

    transformed_image = sopa.utils.to_intrinsic(sdata, image_key, image_key)
    assert transformed_image is image

    transformed_image = spatialdata.transform(
        image, transformation=Scale([2, 2], axes=["y", "x"]), maintain_positioning=True
    )

    assert transformed_image is not image


def test_elements_outside_of_sdata(sdata):
    # shouldn't raise an error, as they are in the same coordinate system
    sopa.utils.to_intrinsic(SpatialData(), sdata["image"], sdata["he_image"])

    # shouldn't raise an error, as they have a shared coordinate system transform
    sopa.utils.to_intrinsic(SpatialData(), sdata["transcripts"], sdata["he_image"])

    del sdata["transcripts"].attrs["transform"]["global"]

    # now it should
    with pytest.raises(ValueError):
        sopa.utils.to_intrinsic(SpatialData(), sdata["transcripts"], sdata["he_image"])


def test_transform_multiple_shortest_path():
    points1 = PointsModel.parse(
        pd.DataFrame({"x": [1, 2], "y": [1, 2]}),
        transformations={"global": Identity(), "cs": Scale([2, 2], ["y", "x"])},
    )
    points2 = PointsModel.parse(
        pd.DataFrame({"x": [4, 3], "y": [3, 1]}),
        transformations={"global": Scale([0.5, 0.5], ["y", "x"]), "cs": Identity()},
    )

    sdata = SpatialData(points={"points1": points1, "points2": points2})

    with pytest.raises(RuntimeError):
        get_transformation_between_coordinate_systems(sdata, points1, points2)

    df = sopa.utils.to_intrinsic(sdata, "points1", "points2").compute()

    expected_coordinates = np.array([[2, 2], [4, 4]])

    assert (df.values == expected_coordinates).all()
