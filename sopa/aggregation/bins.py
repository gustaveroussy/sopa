import logging

import geopandas as gpd
import numpy as np
from anndata import AnnData
from scipy.sparse import coo_matrix
from spatialdata import SpatialData

from ..segmentation.shapes import expand_radius
from ..utils import to_intrinsic

log = logging.getLogger(__name__)


def aggregate_bins(
    sdata: SpatialData,
    shapes_key: str,
    bins_key: str,
    expand_radius_ratio: float = 0,
) -> AnnData:
    """Aggregate bins (for instance, from Visium HD data) into cells.

    Args:
        sdata: The `SpatialData` object
        shapes_key: Key of the shapes containing the cell boundaries
        bins_key: Key of the table containing the bin-by-gene counts
        expand_radius_ratio: Cells polygons will be expanded by `expand_radius_ratio * mean_radius`. This help better aggregate bins from the cytoplasm.

    Returns:
        An `AnnData` object of shape with the cell-by-gene count matrix
    """
    bins_table: AnnData = sdata.tables[bins_key]

    bins_shapes_key = sdata.get_annotated_regions(bins_table)
    bins_shapes_key = bins_shapes_key[0] if isinstance(bins_shapes_key, list) else bins_shapes_key
    bins = sdata.shapes[bins_shapes_key].loc[sdata.get_instance_key_column(bins_table).values]
    bins = gpd.GeoDataFrame(geometry=bins.centroid.values)  # bins as points

    cells = to_intrinsic(sdata, shapes_key, bins_shapes_key).reset_index(drop=True)
    cells = expand_radius(cells, expand_radius_ratio)

    bin_within_cell = gpd.sjoin(bins, cells)

    indices_matrix = coo_matrix(
        (np.full(len(bin_within_cell), 1), (bin_within_cell["index_right"], bin_within_cell.index)),
        shape=(len(cells), len(bins)),
    )

    adata = AnnData(indices_matrix @ bins_table.X, obs=cells[[]], var=bins_table.var)
    adata.obsm["spatial"] = np.stack([cells.centroid.x, cells.centroid.y], axis=1)
    return adata
