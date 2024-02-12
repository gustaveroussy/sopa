import torch
import torch.nn.functional as F
import torchvision.transforms as T

from .models import (
    CellViT,
    CellViT256,
    CellViT256Shared,
    CellViTSAM,
    CellViTSAMShared,
    CellViTShared,
)
from .utils import unflatten_dict


def get_model(
    run_conf: str,
) -> CellViT | CellViTShared | CellViT256 | CellViT256Shared | CellViTSAM | CellViTSAMShared:
    """Get the model from the run configuration"""
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
        raise NotImplementedError(f"Unknown model type. Please select one of {implemented_models}")
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


def load_inference_transforms(run_conf: dict):
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
    model, predictions: dict, magnification: int = 40
) -> tuple[list[dict], torch.Tensor]:
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
    (
        _,
        instance_types,
    ) = model.calculate_instance_map(predictions, magnification=magnification)

    tokens = predictions["tokens"].to("cpu")

    return instance_types, tokens


def load_cellvit(model_path: str, magnification: int, device: str):
    model_checkpoint = torch.load(model_path, map_location="cpu")

    if magnification not in model_path:
        raise ValueError(f"Model path does not contain the magnification level.")

    # unpack checkpoint
    run_conf = unflatten_dict(model_checkpoint["config"], ".")
    model = get_model(run_conf)
    model.load_state_dict(model_checkpoint["model_state_dict"])
    model.eval()
    model.to(device)
    inference_transforms = load_inference_transforms(run_conf)

    return model, inference_transforms
