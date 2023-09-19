import numpy as np
from cellpose import models


def cellpose_patch(
    diameter: float, channels: list[int], model_type: str = "cyto2", **cellpose_kwargs
):
    model = models.Cellpose(model_type=model_type)

    def _(patch: np.ndarray):
        mask, *_ = model.eval(patch, diameter=diameter, channels=channels, **cellpose_kwargs)
        return mask

    return _


def run_cellpose(channels):
    assert len(channels) in [1, 2], f"Provide one or two channel names. Found {len(channels)}"
    channels = [0, 0] if len(channels) == 1 else [1, 2]
    ...
