from torch import nn


class HOPTIMUSFeatures(nn.Module):
    def __init__(self):
        super().__init__()

        try:
            import timm
        except ImportError:
            raise ImportError(
                "Using the hoptimus0 model for inference requires installing the timm dependency, e.g. via `pip install timm`"
            )

        self.model = timm.create_model("hf_hub:bioptimus/H-optimus-0", pretrained=True)

    def forward(self, x):
        from torchvision import transforms

        transform = transforms.Compose(
            [
                transforms.Normalize(mean=(0.707223, 0.578729, 0.703617), std=(0.211883, 0.230117, 0.177517)),
            ]
        )

        return self.model(transform(x))
