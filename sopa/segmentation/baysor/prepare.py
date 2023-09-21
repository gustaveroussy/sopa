from pathlib import Path

from spatialdata import SpatialData

from ...utils.tiling import Tiles2D


def patchify_transcripts(
    sdata: SpatialData, patches_dir: str, points_key: str, tile_width: float, tile_overlap: float
):
    df = sdata.points[points_key]

    patches_dir = Path(patches_dir)

    tiles = Tiles2D.from_dataframe(df, tile_width, tile_overlap)

    for i, (tx0, ty0, tx1, ty1) in enumerate(tiles):
        patch_dir: Path = patches_dir / i
        patch_dir.mkdir(parents=True, exist_ok=True)
        where = (df.x >= tx0) & (df.x <= tx1) & (df.y >= ty0) & (df.y <= ty1)
        df[where].to_csv(patch_dir / "transcripts.csv", single_file=True)
