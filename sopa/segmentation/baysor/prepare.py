from pathlib import Path

import shapely
from shapely.geometry import Point, Polygon
from spatialdata import SpatialData

from ..._constants import ROI
from ...utils.tiling import Tiles2D
from ...utils.utils import _get_element


def patchify_transcripts(
    sdata: SpatialData,
    patches_dir: str,
    tile_width: float,
    tile_overlap: float,
    points_key: str | None = None,
):
    df = _get_element(sdata, "points", points_key)

    patches_dir = Path(patches_dir)

    tiles = Tiles2D.from_dataframe(df, tile_width, tile_overlap)

    polygon = None
    if ROI.KEY in sdata.shapes:
        geo_df = sdata.shapes[ROI.KEY]
        polygon: Polygon = sdata.transform_element_to_coordinate_system(geo_df, df).geometry[0]

    for i, bounds in enumerate(tiles):
        bbox = tiles.polygon(bounds)
        if polygon is not None and not polygon.intersects(bbox):
            continue

        patch_dir: Path = patches_dir / i
        patch_dir.mkdir(parents=True, exist_ok=True)

        tx0, ty0, tx1, ty1 = bounds
        where = (df.x >= tx0) & (df.x <= tx1) & (df.y >= ty0) & (df.y <= ty1)

        if polygon is None or polygon.contains(bbox):
            df[where].to_csv(patch_dir / "transcripts.csv", single_file=True)
        else:
            sub_df = df[where].compute()

            points = [Point(row) for row in df[["x", "y"]].values]
            tree = shapely.STRtree(points)
            indices = tree.query(polygon, predicate="intersects")

            sub_df.iloc[indices].to_csv(patch_dir / "transcripts.csv")
