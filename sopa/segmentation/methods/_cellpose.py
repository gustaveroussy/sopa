from __future__ import annotations

from typing import Callable

import numpy as np
from spatialdata import SpatialData

from ..._constants import SopaKeys
from ._staining import custom_staining_based


def cellpose_patch(
    diameter: float,
    channels: list[str],
    model_type: str = "cyto3",
    pretrained_model: str | bool = False,
    cellpose_model_kwargs: dict | None = None,
    **cellpose_eval_kwargs: int,
) -> Callable:
    """Creation of a callable that runs Cellpose segmentation on a patch

    Args:
        diameter: Cellpose diameter parameter
        channels: List of channel names
        model_type: Cellpose model type
        pretrained_model: Path to the pretrained model to be loaded, or `False`
        cellpose_model_kwargs: Kwargs to be provided to the `cellpose.models.CellposeModel` object
        **cellpose_eval_kwargs: Kwargs to be provided to `model.eval` (where `model` is a `cellpose.models.CellposeModel` object)

    Returns:
        A `callable` whose input is an image of shape `(C, Y, X)` and output is a cell mask of shape `(Y, X)`. Each mask value `>0` represent a unique cell ID
    """
    try:
        from cellpose import models
    except ImportError:
        raise ImportError(
            "To use cellpose, you need its corresponding sopa extra: `pip install 'sopa[cellpose]'` (normal mode) or `pip install -e '.[cellpose]'` (if using snakemake)"
        )

    cellpose_model_kwargs = cellpose_model_kwargs or {}

    if pretrained_model:
        model = models.CellposeModel(pretrained_model=pretrained_model, **cellpose_model_kwargs)
    else:
        model = models.Cellpose(model_type=model_type, **cellpose_model_kwargs)

    if isinstance(channels, str) or len(channels) == 1:
        channels = [0, 0]  # gray scale
    elif len(channels) == 2:
        channels = [1, 2]
    else:
        raise ValueError(f"Provide 1 or 2 channels. Found {len(channels)}")

    def _(patch: np.ndarray):
        mask, *_ = model.eval(patch, diameter=diameter, channels=channels, **cellpose_eval_kwargs)
        return mask

    return _


def cellpose(
    sdata: SpatialData,
    channels: list[str] | str,
    diameter: int,
    image_key: str | None = None,
    min_area: int | None = None,
    delete_cache: bool = True,
    recover: bool = False,
    flow_threshold: float = 2,
    cellprob_threshold: float = -6,
    cellpose_model_kwargs: dict | None = None,
    **cellpose_eval_kwargs: int,
):
    """Run Cellpose segmentation on a SpatialData object, and add a GeoDataFrame containing the cell boundaries.
    The segmentation is run on each patch, for memory efficiency.

    Args:
        sdata: A `SpatialData` object
        channels: Name of the channels to be used for segmentation (or list of channel names).
        diameter: The Cellpose parameter for the expected cell diameter (in pixel).
        image_key: Name of the image in `sdata` to be used for segmentation.
        min_area: Minimum area of a cell to be considered. By default, it is calculated based on the `diameter` parameter.
        delete_cache: Whether to delete the cache after segmentation.
        recover: If `True`, recover the cache from a failed segmentation, and continue.
        flow_threshold: Cellpose `flow_threshold` parameter.
        cellprob_threshold: Cellpose `cellprob_threshold` parameter.
        cellpose_model_kwargs: Dictionary of kwargs to be provided to the `cellpose.models.CellposeModel` object.
        **cellpose_eval_kwargs: Kwargs to be provided to `model.eval` (where `model` is a `cellpose.models.CellposeModel` object)
    """
    channels = channels if isinstance(channels, list) else [channels]

    method = cellpose_patch(
        diameter=diameter,
        channels=channels,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        cellpose_model_kwargs=cellpose_model_kwargs,
        **cellpose_eval_kwargs,
    )

    if min_area is None:
        min_area = (diameter / 2) ** 2  # by default, about 15% of the "normal cell" area

    custom_staining_based(
        sdata,
        method,
        channels,
        image_key=image_key,
        min_area=min_area,
        delete_cache=delete_cache,
        recover=recover,
        cache_dir_name=SopaKeys.CELLPOSE_BOUNDARIES,
        key_added=SopaKeys.CELLPOSE_BOUNDARIES,
    )
