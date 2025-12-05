import dask
import dask.array as da
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import spatialdata
import xarray as xr
from anndata import AnnData
from scipy.sparse import csr_matrix
from shapely import Point, Polygon, box
from spatialdata import SpatialData
from spatialdata.models import PointsModel, ShapesModel, TableModel

import sopa
from sopa.aggregation.channels import _aggregate_channels_aligned
from sopa.aggregation.transcripts import _count_transcripts_aligned
from sopa.constants import SopaKeys

dask.config.set({"dataframe.query-planning": False})
import dask.dataframe as dd  # noqa: E402


def test_aggregate_channels_aligned():
    image = np.random.randint(1, 10, size=(3, 8, 16))
    arr = da.from_array(image, chunks=(1, 8, 8))
    xarr = xr.DataArray(arr, dims=["c", "y", "x"])

    cell_size = 4
    cell_start = [(0, 0), (6, 2), (9, 3)]

    # One cell is on the first block, one is overlapping on both blocks, and one is on the last block
    cells = [box(x, y, x + cell_size - 1, y + cell_size - 1) for x, y in cell_start]

    mean_intensities = _aggregate_channels_aligned(xarr, cells, "average")
    min_intensities = _aggregate_channels_aligned(xarr, cells, "min")
    max_intensities = _aggregate_channels_aligned(xarr, cells, "max")

    true_mean_intensities = np.stack([
        image[:, y : y + cell_size, x : x + cell_size].mean(axis=(1, 2)) for x, y in cell_start
    ])
    true_min_intensities = np.stack([
        image[:, y : y + cell_size, x : x + cell_size].min(axis=(1, 2)) for x, y in cell_start
    ])
    true_max_intensities = np.stack([
        image[:, y : y + cell_size, x : x + cell_size].max(axis=(1, 2)) for x, y in cell_start
    ])

    assert (mean_intensities == true_mean_intensities).all()
    assert (min_intensities == true_min_intensities).all()
    assert (max_intensities == true_max_intensities).all()


def test_count_transcripts():
    df_pandas = pd.DataFrame({
        "x": [1, 2, 3, 7, 11, 1, 3, 2, 2, 1],
        "y": [1, 1, 2, 8, 0, 2, 4, 3, 3, 1],
        "gene": ["a", "a", "b", "c", "a", "c", "gene_control", "b", "b", "blank"],
    })
    df_pandas["gene"] = df_pandas["gene"].astype(object)
    df_pandas.loc[0, "gene"] = np.nan

    points = dd.from_pandas(df_pandas, npartitions=2)
    polygons = [
        Polygon(((1, 2), (3, 2), (3, 4), (1, 4))),
        Polygon(((0, 0), (2, 0), (2, 2), (0, 2))),
        Polygon(((0, 0), (3, 0), (3, 3), (0, 3))),
    ]

    gdf = gpd.GeoDataFrame(geometry=polygons)

    adata = _count_transcripts_aligned(gdf, points, "gene")
    expected = np.array([[0, 3, 1], [1, 0, 1], [1, 3, 1]])

    assert (adata.var_names == ["a", "b", "c"]).all()
    assert (adata.X.toarray() == expected).all()

    adata = _count_transcripts_aligned(gdf, points, "gene", only_excluded=True)
    expected = np.array([[0, 1], [1, 0], [1, 0]])

    assert (adata.var_names == ["blank", "gene_control"]).all()
    assert (adata.X.toarray() == expected).all()


def test_count_transcripts_with_low_quality():
    df_pandas = pd.DataFrame({
        "x": [1, 2, 3, 1],
        "y": [1, 1, 2, 2],
        "gene": ["a", "a", "b", "b"],
        SopaKeys.LOW_QUALITY_TRANSCRIPT_KEY: [True, False, False, False],
    })
    df_pandas["gene"] = df_pandas["gene"].astype(object)

    points = dd.from_pandas(df_pandas, npartitions=2)
    polygons = [box(-1, -1, 10, 10)]

    gdf = gpd.GeoDataFrame(geometry=polygons)

    adata = _count_transcripts_aligned(gdf, points, "gene")

    assert (adata.var_names == ["a", "b"]).all()
    assert (adata.X.toarray() == np.array([[1, 2]])).all()


def test_aggregate_bins():
    bins = [
        box(0, 0, 1, 1),
        box(1, 0, 2, 1),
        box(0, 1, 1, 2),
        box(5, 5, 6, 6),
        box(6, 5, 7, 6),
        box(5, 6, 6, 7),
        box(6, 6, 7, 7),
    ]

    gdf_bins = gpd.GeoDataFrame(geometry=bins)
    gdf_bins.index = [f"bin{i}" for i in range(7)]
    gdf_bins = ShapesModel.parse(gdf_bins)

    cells = [
        Polygon([(-1, -1), (0.5, 0.5), (1, -1)]),  # touches 0
        box(-1, -1, 4, 4),  # touches 0, 1, 2
        Polygon([(4, 4), (4, 9), (9, 4)]),  # touches 3, 4, 5, 6
    ]

    gdf_cells = gpd.GeoDataFrame(geometry=cells)
    gdf_cells = ShapesModel.parse(gdf_cells)

    counts = np.array([
        [3, 2, 1],
        [0, 1, 0],
        [10, 10, 15],
        [0, 0, 12],
        [2, 2, 2],
        [0, 4, 0],
        [32, 1, 2],
    ])

    expected_counts = np.array([
        [3, 2, 1],
        [13, 13, 16],
        [34, 7, 16],
    ])

    adata = AnnData(counts)
    adata.var_names = ["gene1", "gene2", "gene3"]

    adata.obs["bin_id"] = gdf_bins.index
    adata.obs["bin_key"] = "bins_2um"

    adata = TableModel.parse(adata, region="bins_2um", instance_key="bin_id", region_key="bin_key")

    sdata = SpatialData(shapes={"bins_2um": gdf_bins, "cells": gdf_cells}, tables={"table": adata})

    adata_aggr = sopa.aggregation.aggregate_bins(sdata, "cells", "table")

    assert list(adata_aggr.var_names) == ["gene1", "gene2", "gene3"]

    assert (expected_counts == adata_aggr.X).all()

    sdata.tables["table"].X = csr_matrix(sdata.tables["table"].X)

    adata_aggr = sopa.aggregation.aggregate_bins(sdata, "cells", "table")

    assert isinstance(adata_aggr.obsm["bins_assignments"], csr_matrix)

    assert (adata_aggr.obsm["bins_assignments"].nonzero()[1] == [0, 0, 1, 2, 3, 4, 5, 6]).all()

    assert list(adata_aggr.var_names) == ["gene1", "gene2", "gene3"]

    assert (adata_aggr.X.toarray() == expected_counts).all()


@pytest.mark.parametrize("as_sparse", [True, False])
def test_aggregate_bins_no_overlap(as_sparse: bool):
    gdf_bins = gpd.GeoDataFrame(geometry=[box(i, j, i + 0.9, j + 0.9) for i in range(4) for j in range(4)])
    gdf_bins.index = [f"bin{i}" for i in range(len(gdf_bins))]
    gdf_bins = ShapesModel.parse(gdf_bins)

    counts = np.array([
        [
            (i < 2) and (j < 2),
            1,
            (i >= 2 and j >= 2),
            (i < 2 and j >= 2),
        ]
        for i in range(4)
        for j in range(4)
    ]).astype(int)

    if as_sparse:
        counts = csr_matrix(counts)

    adata = AnnData(counts)
    adata.var_names = ["gene1", "gene2", "gene3", "gene4"]

    adata.obs["bin_id"] = gdf_bins.index
    adata.obs["bin_key"] = "bins_2um"

    gdf_cells = gpd.GeoDataFrame(geometry=[Point(0.7, 1.5), Point(1.3, 3), Point(2.5, 3.2)])
    gdf_cells.geometry = gdf_cells.buffer(0.6)
    gdf_cells = ShapesModel.parse(gdf_cells)

    adata = TableModel.parse(adata, region="bins_2um", instance_key="bin_id", region_key="bin_key")

    sdata = SpatialData(shapes={"bins_2um": gdf_bins, "cells": gdf_cells}, tables={"table": adata})

    adata_aggr1 = sopa.aggregation.aggregate_bins(sdata, "cells", "table", no_overlap=False, expand_radius_ratio=1.5)
    assert isinstance(adata_aggr1.obsm["bins_assignments"], csr_matrix)

    adata_aggr2 = sopa.aggregation.aggregate_bins(sdata, "cells", "table", no_overlap=True, expand_radius_ratio=1.5)
    assert isinstance(adata_aggr2.obsm["bins_assignments"], csr_matrix)

    X1 = adata_aggr1.X.toarray() if as_sparse else adata_aggr1.X
    X2 = adata_aggr2.X.toarray() if as_sparse else adata_aggr2.X

    assert (X2 <= X1).all()

    assert (np.array([[4, 6, 0, 2], [0, 6, 2, 4], [0, 6, 4, 2]]) == X1).all()

    assert (np.array([[4, 4, 0, 0], [0, 4, 0, 4], [0, 4, 4, 0]]) == X2).all()


def test_overlay():
    np.random.seed(0)
    x, y = np.meshgrid(np.arange(6), np.arange(6))
    x = x.flatten()
    y = y.flatten()
    genes = np.random.choice(["gene_a", "gene_b", "gene_c"], x.shape[0])
    points = PointsModel.parse(pd.DataFrame({"x": x, "y": y, "gene": genes}), feature_key="gene")

    gdf = gpd.GeoDataFrame(
        geometry=[
            box(0, 0, 1, 1),
            box(1, 1, 2, 2),
            box(2, 2, 4, 4),
            box(1, 2.5, 1.5, 5),
        ]
    )
    gdf.geometry = gdf.geometry.buffer(0.1)
    gdf = ShapesModel.parse(gdf)

    gdf2 = gpd.GeoDataFrame(
        geometry=[
            box(0, 0, 1, 5),
            box(3, 3, 5, 5),
        ]
    )
    gdf2.geometry = gdf2.geometry.buffer(0.1)
    gdf2 = ShapesModel.parse(gdf2)

    sdata = spatialdata.SpatialData(shapes={"cellpose_boundaries": gdf, "selection": gdf2}, points={"points": points})

    sopa.aggregate(sdata, aggregate_channels=False)

    assert (
        sdata["table"].X.toarray()
        == np.array([
            [2, 1, 1],
            [1, 0, 3],
            [5, 2, 2],
            [0, 1, 2],
        ])
    ).all()

    sopa.overlay_segmentation(sdata, "selection", area_ratio_threshold=0.5)

    assert (
        sdata["table"].X.toarray()
        == np.array([
            [1, 0, 3],
            [5, 2, 2],
            [0, 1, 2],
            [2, 6, 4],
            [3, 5, 1],
        ])
    ).all()
