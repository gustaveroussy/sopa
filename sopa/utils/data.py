import dask.dataframe as dd
import numpy as np
import pandas as pd
from spatialdata import SpatialData
from spatialdata.datasets import BlobsDataset


def blobs(
    *_,
    length: int = 1_024,
    n_points: int = 10_000,
    c_coords=["DAPI", "CD3", "EPCAM", "CD20"],
    **kwargs,
) -> SpatialData:
    _blobs = BlobsDataset(
        length=length, n_points=n_points, c_coords=c_coords, n_channels=len(c_coords), **kwargs
    )

    image = _blobs._image_blobs(
        _blobs.transformations,
        _blobs.length,
        _blobs.n_channels,
        _blobs.c_coords,
    )
    image.data = (image.data * 255).astype(np.uint8)

    points = _blobs._points_blobs(_blobs.transformations, _blobs.length, _blobs.n_points)
    genes = pd.Series(np.random.choice(list("abcdef"), size=len(points))).astype("category")
    points["genes"] = dd.from_pandas(genes, npartitions=points.npartitions)

    return SpatialData(images={"blob_image": image}, points={"blob_transcripts": points})
