import timm
from torch import nn
from torchvision import transforms


class HOPTIMUSFeatures(nn.Module):
    def __init__(self):
        super().__init__()
        self.model = timm.create_model("hf_hub:bioptimus/H-optimus-0", pretrained=True)

    def forward(self, x):
        transform = transforms.Compose(
            [
                transforms.Normalize(mean=(0.707223, 0.578729, 0.703617), std=(0.211883, 0.230117, 0.177517)),
            ]
        )

        return self.model(transform(x))
