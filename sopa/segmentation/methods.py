from typing import Callable, List, Tuple, Union

import numpy as np


def cellpose_patch(
    diameter: float, channels: list[str], model_type: str = "cyto2", **cellpose_kwargs: int
) -> Callable:
    """Creation of a callable that runs Cellpose segmentation on a patch

    Args:
        diameter: Cellpose diameter parameter
        channels: List of channel names
        model_type: Cellpose model type
        **cellpose_kwargs: Kwargs to be provided to `model.eval` (where `model` is a `cellpose.models.Cellpose` object)

    Returns:
        A `callable` whose input is an image of shape `(C, Y, X)` and output is a cell mask of shape `(Y, X)`. Each mask value `>0` represent a unique cell ID
    """
    try:
        from cellpose import models
    except ImportError:
        raise ImportError(
            "To use cellpose, you need its corresponding sopa extra: `pip install 'sopa[cellpose]'`"
        )

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


def cellvit_patch(model_path: str, magnification: int, device: str = "cpu") -> Callable:

    import torch

    from .cellvit import load_cellvit

    model, inference_transforms, get_cell_predictions_with_tokens = load_cellvit(
        model_path, magnification, device
    )

    TYPE_NUCLEI_DICT = {
        1: "Neoplastic",
        2: "Inflammatory",
        3: "Connective",
        4: "Dead",
        5: "Epithelial",
    }

    def _(patch: np.ndarray):
        """
        predictions output format:
                * tissue_types: Raw tissue type prediction. Shape: (B, num_tissue_classes)
                * nuclei_binary_map: Raw binary cell segmentation predictions. Shape: (B, 2, H, W)
                * hv_map: Binary HV Map predictions. Shape: (B, 2, H, W)
                * nuclei_type_map: Raw binary nuclei type predictions. Shape: (B, num_nuclei_classes, H, W)
                * [Optional, if retrieve tokens]: tokens
                * [Optional, if regression loss]:
                * regression_map: Regression map for binary prediction. Shape: (B, 2, H, W)
        """
        patch = inference_transforms(patch).unsqueeze(0).to(device)
        with torch.no_grad():
            predictions = model(patch)
        instance_types, tokens = get_cell_predictions_with_tokens(
            model, predictions, magnification=magnification
        )

        return


def dummy_method(**method_kwargs):
    """A method builder builder (i.e. it returns a segmentation function).
    Kwargs can be provided and used in the below function"""

    def segmentation_function(image: np.ndarray) -> np.ndarray:
        """A dummy example of a custom segmentation method
        that creates one cell (with a padding of 10 pixels).

        Args:
            image: An image of shape `(C, Y, X)`

        Returns:
            A mask of shape `(Y, X)` containing one cell
        """
        mask = np.zeros(image.shape[1:], dtype=int)

        # one cell, corresponding to value 1
        mask[10:-10, 10:-10] = 1  # squared shaped

        return mask

    return segmentation_function
