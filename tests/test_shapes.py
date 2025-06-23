import dask.array as da
import geopandas as gpd
import numpy as np
import xarray as xr
from geopandas import GeoDataFrame
from shapely import MultiPolygon, Point, Polygon, box

from sopa.aggregation.channels import _aggregate_channels_aligned
from sopa.shapes import ensure_polygon, expand_radius, remove_overlap, to_valid_polygons

gdf_squares = gpd.GeoDataFrame(
    {"color": [0, 1, 2, 3, 4]},
    geometry=[
        box(0, 0, 0.8, 0.8),
        box(1, 1, 2, 2),
        box(0.5, 0.5, 1.5, 1.5),
        box(2.5, 2.5, 3.4, 3.4),
        box(1.4, 0.5, 2, 1.1),
    ],
)

gdf_circles = gpd.GeoDataFrame(
    {"color": [0, 1, 2]},
    geometry=[Point(0, 0).buffer(1), Point(0, 1.75).buffer(1), Point(0.5, 1).buffer(0.3)],
)


def test_remove_overlap():
    for gdf in [gdf_squares, gdf_circles]:
        res = remove_overlap(gdf)

        assert res.union_all().area / gdf.union_all().area > 0.999  # should conserve the original area

        res.geometry = res.geometry.buffer(-1e-3)  # to avoid polygons touching on single points

        assert len(gpd.sjoin(res, res)) == len(res)  # should not have overlaps

        res = expand_radius(gdf, 0.5, no_overlap=True)  # should not raise an error

        res.geometry = res.geometry.buffer(-1e-3)

        assert len(gpd.sjoin(res, res)) == len(res)


def test_remove_overlap_empty():
    gdf = GeoDataFrame(geometry=[box(-1, -1, 1, 1), Point(0, 0).buffer(0.5)])  # second shape is indluced in the first

    res = remove_overlap(gdf)

    assert res.geometry[1].is_empty

    image = np.random.randint(1, 10, size=(3, 8, 16))
    arr = da.from_array(image, chunks=(1, 8, 8))
    xarr = xr.DataArray(arr, dims=["c", "y", "x"])

    # we can still run aggregation on the empty shape
    X = _aggregate_channels_aligned(xarr, res.geometry, mode="average")
    assert (X[1] == 0).all()  # should be all zeros for the empty shape


def test_to_valid_polygons():
    gdf = GeoDataFrame(geometry=[Point(0, 0)])
    assert to_valid_polygons(gdf).empty

    mp = MultiPolygon([box(2, 2, 3, 3), box(-2, -2, 1, 1)])

    assert ensure_polygon(mp).area == 9
    assert ensure_polygon(mp, simple_polygon=False) == mp

    ext = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
    interior = [(0.1, 0.1), (0.1, 0.9), (0.9, 0.9), (0.9, 0.1), (0.1, 0.1)]
    polygon = Polygon(ext, [interior])

    assert ensure_polygon(polygon, False) == polygon
    assert ensure_polygon(polygon, True).area == 1
