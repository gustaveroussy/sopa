import numpy as np
from cellpose import models


def cellpose_patch(
    diameter: float, channels: list[int], model_type: str = "cyto2", **cellpose_kwargs
):
    model = models.Cellpose(model_type=model_type)
    channels = list(range(len(channels)))

    def _(patch: np.ndarray):
        mask, *_ = model.eval(patch, diameter=diameter, channels=channels, **cellpose_kwargs)
        return mask

    return _
