import geopandas as gpd
from shapely.geometry import box

from sopa.segmentation.shapes import expand_radius


def test_expand_no_overlap():
    gdf = gpd.GeoDataFrame(
        geometry=[
            box(0, 0, 0.8, 0.8),
            box(1, 1, 2, 2),
            box(0.5, 0.5, 1.5, 1.5),
            box(2.5, 2.5, 3.4, 3.4),
            box(1.4, 0.5, 2, 1.1),
        ],
    )
    gdf["color"] = [0, 1, 2, 3, 4]  # for debugging only

    res = expand_radius(gdf, 0.5, no_overlap=True)
    res.geometry = res.geometry.buffer(-1e-5)  # to avoid polygons touching on single points

    assert len(gpd.sjoin(res, res)) == len(res)
