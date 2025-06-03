import geopandas as gpd
from shapely.geometry import Point, box

from sopa.segmentation.shapes import expand_radius, remove_overlap

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
