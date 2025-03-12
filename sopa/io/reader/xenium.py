import logging
from pathlib import Path

from spatialdata import SpatialData

from ..._constants import ATTRS_KEY, SopaAttrs, SopaKeys
from ...utils import ensure_string_channel_names
from .utils import _default_image_kwargs

log = logging.getLogger(__name__)


def xenium(
    path: str | Path,
    image_models_kwargs: dict | None = None,
    imread_kwargs: dict | None = None,
    cells_boundaries: int = False,
    cells_table: int = False,
    nucleus_labels: int = False,
    cells_labels: int = False,
    cells_as_circles: int = False,
    nucleus_boundaries: int = False,
    qv_threshold: int | None = None,
    **kwargs: int,
) -> SpatialData:
    """Read Xenium data as a `SpatialData` object. For more information, refer to [spatialdata-io](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.xenium.html).

    This function reads the following files:
        - `transcripts.parquet`: transcripts locations and names
        - `experiment.xenium`: metadata file
        - `morphology_focus.ome.tif`: morphology image (or a directory, for recent versions of the Xenium)


    Args:
        path: Path to the Xenium directory containing all the experiment files
        image_models_kwargs: Keyword arguments passed to `spatialdata.models.Image2DModel`.
        imread_kwargs: Keyword arguments passed to `dask_image.imread.imread`.
        cells_boundaries: Whether to read cell boundaries
        cells_table: Whether to read cell table
        nucleus_labels: Whether to read nucleus labels
        cells_labels: Whether to read cell labels
        cells_as_circles: Whether to read cells as circles
        nucleus_boundaries: Whether to read nucleus boundaries
        qv_threshold: Whether to add a "low_quality_transcript" column to transcripts. Transcripts with a QV value below `qv_threshold` will not be used during segmentation.
        kwargs: Additional keyword arguments passed to `spatialdata_io.xenium

    Returns:
        A `SpatialData` object representing the Xenium experiment
    """
    from spatialdata_io.readers.xenium import xenium as xenium_spatialdata_io

    image_models_kwargs, imread_kwargs = _default_image_kwargs(image_models_kwargs, imread_kwargs)

    sdata: SpatialData = xenium_spatialdata_io(
        path,
        cells_table=cells_table,
        nucleus_labels=nucleus_labels,
        cells_labels=cells_labels,
        cells_as_circles=cells_as_circles,
        nucleus_boundaries=nucleus_boundaries,
        cells_boundaries=cells_boundaries,
        image_models_kwargs=image_models_kwargs,
        imread_kwargs=imread_kwargs,
        **kwargs,
    )

    if "table" in sdata.tables:
        sdata["table"].uns[ATTRS_KEY]["region"] = "cell_boundaries"
        sdata["table"].obs["region"] = "cell_boundaries"
        sdata["table"].obs["region"] = sdata["table"].obs["region"].astype("category")

    ensure_string_channel_names(sdata)

    ### Add Sopa attributes to detect the spatial elements
    if "morphology_focus" in sdata.images:
        sdata.attrs[SopaAttrs.CELL_SEGMENTATION] = "morphology_focus"

    if "he_image" in sdata.images:
        sdata.attrs[SopaAttrs.TISSUE_SEGMENTATION] = "he_image"

    if "transcripts" in sdata.points:
        sdata.attrs[SopaAttrs.TRANSCRIPTS] = "transcripts"
        if qv_threshold:
            sdata.points["transcripts"][SopaKeys.LOW_QUALITY_TRANSCRIPT_KEY] = (
                sdata.points["transcripts"].qv < qv_threshold
            )

        if qv_threshold is not None:
            assert "qv" in sdata.points["transcripts"].columns, "QV column not found in `sdata['transcripts']`"

            sdata.points["transcripts"][SopaKeys.LOW_QUALITY_TRANSCRIPT_KEY] = (
                sdata.points["transcripts"]["qv"] < qv_threshold
            )

    sdata.attrs[SopaAttrs.XENIUM_OUTPUT_PATH] = str(Path(path).resolve())

    return sdata
