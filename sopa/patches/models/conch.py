from torch import nn


class ConchFeatures(nn.Module):
    def __init__(self):
        super().__init__()
        self.conch, self.transform = _import_conch()
        self.conch.eval()

    def forward(self, x):
        return self.conch(self.transform(x))


def _import_conch():
    try:
        import einops_exts  # noqa: F401
        import timm  # noqa: F401
        from torchvision import transforms
        from transformers import AutoModel
    except ImportError:
        raise ImportError("To use CONCH, please install the CONCH dependencies: transformers, einops_exts, timm.")

    titan = AutoModel.from_pretrained("MahmoodLab/TITAN", trust_remote_code=True)
    conch, _ = titan.return_conch()
    transform = transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225))

    return conch, transform
