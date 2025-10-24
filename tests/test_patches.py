import shutil

import dask
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import spatialdata
from shapely import Point
from spatialdata import SpatialData
from spatialdata.models import Image2DModel, PointsModel, ShapesModel

import sopa
from sopa._constants import SopaFiles, SopaKeys
from sopa.patches._patches import Patches1D, Patches2D
from sopa.patches._transcripts import _unassigned_to_zero
from sopa.segmentation._transcripts import _check_transcript_patches

dask.config.set({"dataframe.query-planning": False})
import dask.dataframe as dd  # noqa: E402


@pytest.fixture
def sdata() -> SpatialData:
    sdata = sopa.io.toy_dataset(length=512, cell_density=1e-3)
    return sdata


def test_patchify_image(sdata: SpatialData):
    sopa.make_image_patches(sdata, 300, 100)
    assert len(sdata[SopaKeys.PATCHES]) == 9

    sopa.make_image_patches(sdata, 512, 0)
    assert len(sdata[SopaKeys.PATCHES]) == 1


def test_patchify_inside_tissue_roi(sdata: SpatialData):
    sopa.make_image_patches(sdata, 80, 0)
    assert len(sdata[SopaKeys.PATCHES]) == 49

    sopa.segmentation.tissue(sdata)

    sopa.make_image_patches(sdata, 80, 0)
    assert len(sdata[SopaKeys.PATCHES]) == 42  # inside the tissue ROI

    del sdata.shapes[SopaKeys.ROI]


def test_pathify_transcripts():
    # reproduce issue #214

    xmax = np.float32(11475.797)
    xmin = np.float32(2.671875)

    Patches1D(xmin=xmin, xmax=xmax, patch_width=300.0, patch_overlap=30.0, tight=True, int_coords=False)


def test_patchify_baysor(sdata: SpatialData):
    with pytest.raises(AssertionError):
        sopa.utils.get_transcripts_patches_dirs(sdata)

    sopa.make_transcript_patches(sdata, 30, 10)
    assert len(sdata[SopaKeys.TRANSCRIPTS_PATCHES]) == 9
    assert len(sopa.utils.get_transcripts_patches_dirs(sdata)) == 9

    sopa.make_transcript_patches(sdata, 52, 0)
    assert len(sdata[SopaKeys.TRANSCRIPTS_PATCHES]) == 1
    assert len(sopa.utils.get_transcripts_patches_dirs(sdata)) == 1

    sopa.utils.delete_cache(sdata)


def test_patchify_baysor_inside_tissue_roi(sdata: SpatialData):
    sopa.make_transcript_patches(sdata, 5, 0, min_points_per_patch=0)

    assert len(sdata[SopaKeys.TRANSCRIPTS_PATCHES]) == 121

    sopa.segmentation.tissue(sdata)

    sopa.make_transcript_patches(sdata, 5, 0, min_points_per_patch=0)

    assert len(sdata[SopaKeys.TRANSCRIPTS_PATCHES]) == 107  # inside the tissue ROI

    sopa.utils.delete_cache(sdata)


def test_patches_inference_clustering():
    sdata = sopa.io.toy_dataset(length=200)

    with pytest.raises(AssertionError):  # two images, and no image_key provided
        sopa.patches.compute_embeddings(sdata, "dummy", 50)

    sopa.patches.compute_embeddings(sdata, "dummy", 50, image_key="he_image")

    assert sdata["dummy_embeddings"].shape == (4, 3)

    sopa.patches.compute_embeddings(sdata, "dummy", 25, level=-1, image_key="he_image")

    assert sdata["dummy_embeddings"].shape == (1, 3)

    sopa.patches.compute_embeddings(sdata, "dummy", 13, patch_overlap=3, level=-1, image_key="he_image")

    assert sdata["dummy_embeddings"].shape == (9, 3)

    sdata["he_image"].attrs["backend"] = "tiffslide"
    sdata["he_image"].attrs["metadata"] = {
        "level_downsamples": [1, 2, 4],
        "properties": {"tiffslide.objective-power": 50},
    }

    sopa.patches.compute_embeddings(sdata, "dummy", 10, magnification=10, image_key="he_image")
    assert sdata["dummy_embeddings"].shape == (9, 3)

    sopa.patches.compute_embeddings(sdata, "dummy", 11, magnification=10, image_key="he_image")
    assert sdata["dummy_embeddings"].shape == (4, 3)

    sopa.patches.compute_embeddings(sdata, "dummy", 50, magnification=100, image_key="he_image")
    assert sdata["dummy_embeddings"].shape == (16, 3)

    sopa.patches.cluster_embeddings(sdata, "dummy_embeddings")

    assert "cluster" in sdata["dummy_embeddings"].obs

    # Testing the xarray slicing for the read region
    sopa.settings.native_read_region = False
    sopa.patches.compute_embeddings(sdata, "dummy", 10, magnification=10, image_key="he_image")
    assert sdata["dummy_embeddings"].shape == (9, 3)


def test_gene_exlude_pattern():
    _default_gene_exclude_pattern = sopa.settings.gene_exclude_pattern

    num_to_exclude = 5

    gene_names = [
        pd.NA,  # should be excluded
        np.nan,  # should be excluded
        "bLaNk",  # should be excluded
        "SystemControl",  # should be excluded
        "negcontrol_gene",  # should be excluded
        "gene1",
        "gene2",
        "gene3",
        "CDNAN",  # should not be excluded although it contains 'nan'
    ]

    df = pd.DataFrame({
        "x": np.zeros(len(gene_names)),
        "y": np.zeros(len(gene_names)),
        "z": np.zeros(len(gene_names)),
        "genes": gene_names,
    })
    points = {
        "transcripts": PointsModel.parse(df, feature_key="genes"),
    }

    sdata = SpatialData(points=points)

    sopa.make_transcript_patches(sdata, min_points_per_patch=0)

    df = pd.read_csv(
        sopa.utils.get_cache_dir(sdata) / SopaFiles.TRANSCRIPT_CACHE_DIR / "0" / SopaFiles.TRANSCRIPTS_FILE
    )

    assert len(df) == len(sdata["transcripts"]) - num_to_exclude

    sopa.settings.gene_exclude_pattern = None

    sopa.make_transcript_patches(sdata, min_points_per_patch=0)

    df = pd.read_csv(
        sopa.utils.get_cache_dir(sdata) / SopaFiles.TRANSCRIPT_CACHE_DIR / "0" / SopaFiles.TRANSCRIPTS_FILE
    )

    assert len(df) == len(sdata["transcripts"])

    sopa.settings.gene_exclude_pattern = _default_gene_exclude_pattern


def test_patches_with_and_without_centroids():
    gdf = ShapesModel.parse(gpd.GeoDataFrame(geometry=[Point(100, 100).buffer(80)]))
    im = Image2DModel.parse(np.zeros((1, 200, 200), dtype=np.uint8))

    sdata = SpatialData(images={"im": im}, shapes={"gdf": gdf})

    patches = Patches2D(sdata, element="im", patch_width=67, patch_overlap=0, roi_key="gdf")
    assert len(patches) == 9
    assert patches.geo_df.index[0] == 0

    patches = Patches2D(sdata, element="im", patch_width=67, patch_overlap=0, roi_key="gdf", use_roi_centroids=True)
    assert len(patches) == 1
    assert patches.geo_df.index[0] == 0


def test_unassigned_to_zero():
    df = pd.DataFrame({"col1": ["a", "b", "zero", "b"], "col2": [-1, 0, 1, 2], "col3": [0, 1, 2, 10]})
    df = dd.from_pandas(df, npartitions=1)

    assert (_unassigned_to_zero(df["col1"], "zero").compute() == [1, 2, 0, 2]).all()
    assert (_unassigned_to_zero(df["col2"], -1).compute() == [0, 9223372036854775806, 1, 2]).all()
    assert (_unassigned_to_zero(df["col3"], 0).compute() == [0, 1, 2, 10]).all()


def test_move_sdata_transcript_cache():
    sdata = sopa.io.toy_dataset()

    sopa.segmentation.tissue(sdata)
    sdata.write("test.zarr")

    sopa.make_transcript_patches(sdata, patch_width=70, patch_overlap=0)

    shutil.move("test.zarr", "test_moved.zarr")

    sdata = spatialdata.read_zarr("test_moved.zarr")
    _check_transcript_patches(sdata)  # the cache should still be detected

    shutil.rmtree("test_moved.zarr")
