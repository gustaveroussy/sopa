import pytest
from spatialdata import SpatialData

import sopa
from sopa._constants import SopaKeys


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


def test_patchify_baysor(sdata: SpatialData):
    with pytest.raises(AssertionError):
        sopa.utils.get_transcripts_patches_dirs(sdata)

    sopa.make_transcript_patches(sdata, 30, 10)
    assert len(sdata[SopaKeys.TRANSCRIPT_PATCHES]) == 9
    assert len(sopa.utils.get_transcripts_patches_dirs(sdata)) == 9

    sopa.make_transcript_patches(sdata, 52, 0)
    assert len(sdata[SopaKeys.TRANSCRIPT_PATCHES]) == 1
    assert len(sopa.utils.get_transcripts_patches_dirs(sdata)) == 1

    sopa.utils.delete_cache(sdata)


def test_patchify_baysor_inside_tissue_roi(sdata: SpatialData):
    sopa.make_transcript_patches(sdata, 5, 0, min_points_per_patch=0)

    assert len(sdata[SopaKeys.TRANSCRIPT_PATCHES]) == 121

    sopa.segmentation.tissue(sdata)

    sopa.make_transcript_patches(sdata, 5, 0, min_points_per_patch=0)

    assert len(sdata[SopaKeys.TRANSCRIPT_PATCHES]) == 115  # inside the tissue ROI

    sopa.utils.delete_cache(sdata)
