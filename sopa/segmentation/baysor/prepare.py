import logging
from pathlib import Path

import shapely
from shapely.geometry import Point, Polygon
from spatialdata import SpatialData
from tqdm import tqdm

from ..._constants import ROI
from ...utils.tiling import Tiles2D
from ...utils.utils import _get_element

log = logging.getLogger(__name__)


def patchify_transcripts(
    sdata: SpatialData,
    patch_attrs_file: str,
    tile_width: float,
    tile_overlap: float,
    points_key: str | None = None,
):
    df = _get_element(sdata, "points", points_key)

    patches_dir = Path(patch_attrs_file).parent

    tiles = Tiles2D.from_dataframe(df, tile_width, tile_overlap)

    polygon = None
    if ROI.KEY in sdata.shapes:
        geo_df = sdata.shapes[ROI.KEY]
        polygon: Polygon = sdata.transform_element_to_coordinate_system(geo_df, df).geometry[0]

    patches_paths = []

    for i, bounds in tqdm(enumerate(tiles)):
        bbox = tiles.polygon(bounds)
        if polygon is not None and not polygon.intersects(bbox):
            continue

        patch_dir = (patches_dir / str(i)).absolute()
        patch_dir.mkdir(parents=True, exist_ok=True)
        patch_path = patch_dir / "transcripts.csv"

        patches_paths.append(patch_dir)

        tx0, ty0, tx1, ty1 = bounds
        where = (df.x >= tx0) & (df.x <= tx1) & (df.y >= ty0) & (df.y <= ty1)

        if polygon is None or polygon.contains(bbox):
            df[where].to_csv(patch_path, single_file=True)
        else:
            sub_df = df[where].compute()

            points = [Point(row) for row in sub_df[["x", "y"]].values]
            tree = shapely.STRtree(points)
            indices = tree.query(polygon, predicate="intersects")

            sub_df.iloc[indices].to_csv(patch_path)

    with open(patch_attrs_file, "w") as f:
        f.write("\n".join(map(str, patches_paths)))

    log.info(f"Patch informations saved in file {patch_attrs_file}")
