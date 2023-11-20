import numpy as np
from cellpose import models


def cellpose_patch(
    diameter: float, channels: list[str], model_type: str = "cyto2", **cellpose_kwargs
):
    model = models.Cellpose(model_type=model_type)

    if isinstance(channels, str) or len(channels) == 1:
        channels = [0, 0]  # gray scale
    elif len(channels) == 2:
        channels = [1, 2]
    else:
        raise ValueError(f"Provide 1 or 2 channels. Found {len(channels)}")

    def _(patch: np.ndarray):
        mask, *_ = model.eval(patch, diameter=diameter, channels=channels, **cellpose_kwargs)
        return mask

    return _
