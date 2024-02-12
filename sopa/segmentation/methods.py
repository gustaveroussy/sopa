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


TYPE_NUCLEI_DICT = {
    1: "Neoplastic",
    2: "Inflammatory",
    3: "Connective",
    4: "Dead",
    5: "Epithelial",
}


def cellvit_patch(model_path: str, magnification: int, device: str) -> Callable:

    import torch
    import torch.nn.functional as F
    import torchvision.transforms as T

    from sopa.segmentation.cellvit.models import (
        CellViT,
        CellViT256,
        CellViT256Shared,
        CellViTSAM,
        CellViTSAMShared,
        CellViTShared,
    )
    from sopa.segmentation.cellvit.utils import unflatten_dict

    model_checkpoint = torch.load(model_path, map_location="cpu")

    def get_model(
        run_conf: str,
    ) -> Union[
        CellViT, CellViTShared, CellViT256, CellViT256Shared, CellViTSAM, CellViTSAMShared,
    ]:
        """ Get the model from the run configuration"""
        model_type = run_conf["arch"]
        implemented_models = [
            "CellViT",
            "CellViTShared",
            "CellViT256",
            "CellViT256Shared",
            "CellViTSAM",
            "CellViTSAMShared",
        ]
        if model_type not in implemented_models:
            raise NotImplementedError(
                f"Unknown model type. Please select one of {implemented_models}"
            )
        if model_type in ["CellViT", "CellViTShared"]:
            if model_type == "CellViT":
                model_class = CellViT
            elif model_type == "CellViTShared":
                model_class = CellViTShared
            model = model_class(
                num_nuclei_classes=run_conf["data"]["num_nuclei_classes"],
                num_tissue_classes=run_conf["data"]["num_tissue_classes"],
                embed_dim=run_conf["model"]["embed_dim"],
                input_channels=run_conf["model"].get("input_channels", 3),
                depth=run_conf["model"]["depth"],
                num_heads=run_conf["model"]["num_heads"],
                extract_layers=run_conf["model"]["extract_layers"],
                regression_loss=run_conf["model"].get("regression_loss", False),
            )

        elif model_type in ["CellViT256", "CellViT256Shared"]:
            if model_type == "CellViT256":
                model_class = CellViT256
            elif model_type == "CellViTVIT256Shared":
                model_class = CellViT256Shared
            model = model_class(
                model256_path=None,
                num_nuclei_classes=run_conf["data"]["num_nuclei_classes"],
                num_tissue_classes=run_conf["data"]["num_tissue_classes"],
                regression_loss=run_conf["model"].get("regression_loss", False),
            )
        elif model_type in ["CellViTSAM", "CellViTSAMShared"]:
            if model_type == "CellViTSAM":
                model_class = CellViTSAM
            elif model_type == "CellViTSAMShared":
                model_class = CellViTSAMShared
            model = model_class(
                model_path=None,
                num_nuclei_classes=run_conf["data"]["num_nuclei_classes"],
                num_tissue_classes=run_conf["data"]["num_tissue_classes"],
                vit_structure=run_conf["model"]["backbone"],
                regression_loss=run_conf["model"].get("regression_loss", False),
            )
        return model

    def load_inference_transforms(run_config):
        """Load the inference transformations from the run_configuration"""

        transform_settings = run_conf["transformations"]
        if "normalize" in transform_settings:
            mean = transform_settings["normalize"].get("mean", (0.5, 0.5, 0.5))
            std = transform_settings["normalize"].get("std", (0.5, 0.5, 0.5))
        else:
            mean = (0.5, 0.5, 0.5)
            std = (0.5, 0.5, 0.5)
        inference_transforms = T.Compose([T.ToTensor(), T.Normalize(mean=mean, std=std)])
        return inference_transforms

    def get_cell_predictions_with_tokens(
        predictions: dict, magnification: int = 40
    ) -> Tuple[List[dict], torch.Tensor]:
        """Take the raw predictions, apply softmax and calculate type instances

        Args:
            predictions (dict): Network predictions with tokens. Keys:
            magnification (int, optional): WSI magnification. Defaults to 40.

        Returns:
            Tuple[List[dict], torch.Tensor]:
                * List[dict]: List with a dictionary for each batch element with cell seg results
                    Contains bbox, contour, 2D-position, type and type_prob for each cell
                * List[dict]: Network tokens on cpu device with shape (batch_size, num_tokens_h, num_tokens_w, embd_dim)
        """
        predictions["nuclei_binary_map"] = F.softmax(
            predictions["nuclei_binary_map"], dim=1
        )  # shape: (batch_size, 2, H, W)
        predictions["nuclei_type_map"] = F.softmax(
            predictions["nuclei_type_map"], dim=1
        )  # shape: (batch_size, num_nuclei_classes, H, W)
        # get the instance types
        (_, instance_types,) = model.calculate_instance_map(
            predictions, magnification=magnification
        )

        tokens = predictions["tokens"].to("cpu")

        return instance_types, tokens

    if magnification not in model_path:
        raise ValueError(f"Model path does not contain the magnification level.")

    # unpack checkpoint
    run_conf = unflatten_dict(model_checkpoint["config"], ".")
    model = get_model(run_conf)
    model.load_state_dict(model_checkpoint["model_state_dict"])
    model.eval()
    model.to(device)
    inference_transforms = load_inference_transforms(run_conf)

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
            predictions, magnification=magnification
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
