from spatialdata import SpatialData

from ...utils.tiling import Tiles2D


def patchify_transcripts(sdata: SpatialData, points_key: str):
    df = sdata.points[points_key]

    xmin, ymin = df.x.min().compute(), df.y.min().compute()
    xmax, ymax = df.x.max().compute(), df.y.max().compute()

    tiles = Tiles2D(500, xmax, ymax, ymin=ymin, xmin=xmin, tight=True)

    for i, (tx0, ty0, tx1, ty1) in enumerate(tiles):
        where = (df.x >= tx0) & (df.x <= tx1) & (df.y >= ty0) & (df.y <= ty1)
        df[where].to_csv(f"tiles/tile_{i}.csv", single_file=True)
