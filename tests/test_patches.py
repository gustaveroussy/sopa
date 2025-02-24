import numpy as np
import pandas as pd
import pytest
from spatialdata import SpatialData

import sopa
from sopa._constants import SopaFiles, SopaKeys
from sopa.patches._patches import Patches1D


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

    assert len(sdata[SopaKeys.TRANSCRIPTS_PATCHES]) == 115  # inside the tissue ROI

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

    sdata["he_image"].attrs["backend"] = "test"
    sdata["he_image"].attrs["metadata"] = {
        "level_downsamples": [1, 2, 4],
        "properties": {"test.objective-power": 50},
    }

    sopa.patches.compute_embeddings(sdata, "dummy", 10, magnification=10, image_key="he_image")
    assert sdata["dummy_embeddings"].shape == (9, 3)

    sopa.patches.compute_embeddings(sdata, "dummy", 11, magnification=10, image_key="he_image")
    assert sdata["dummy_embeddings"].shape == (4, 3)

    sopa.patches.compute_embeddings(sdata, "dummy", 50, magnification=100, image_key="he_image")
    assert sdata["dummy_embeddings"].shape == (4, 3)

    sopa.patches.cluster_embeddings(sdata, "dummy_embeddings")

    assert "cluster" in sdata["dummy_embeddings"].obs


def test_gene_exlude_pattern():
    _default_gene_exclude_pattern = sopa.settings.gene_exclude_pattern

    sdata = sopa.io.toy_dataset(add_nan_gene_name=True, length=1000)

    sopa.make_transcript_patches(sdata)

    df = pd.read_csv(
        sopa.utils.get_cache_dir(sdata) / SopaFiles.TRANSCRIPT_CACHE_DIR / "0" / SopaFiles.TRANSCRIPTS_FILE
    )

    assert len(df) == len(sdata["transcripts"]) - 1

    sopa.settings.gene_exclude_pattern = None

    sopa.make_transcript_patches(sdata)

    df = pd.read_csv(
        sopa.utils.get_cache_dir(sdata) / SopaFiles.TRANSCRIPT_CACHE_DIR / "0" / SopaFiles.TRANSCRIPTS_FILE
    )

    assert len(df) == len(sdata["transcripts"])

    sopa.settings.gene_exclude_pattern = _default_gene_exclude_pattern
