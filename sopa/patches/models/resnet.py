from torch import nn


class Resnet50Features(nn.Module):
    def __init__(self):
        super().__init__()

        try:
            import timm
        except ImportError:
            raise ImportError(
                "Using the 'resnet50' model for inference requires installing the timm dependency, e.g. via `pip install timm`"
            )
        self.model = timm.create_model("resnet50", pretrained=True, num_classes=0)
        self.model.eval()

        from torchvision import transforms

        self.transform = transforms.Compose([
            transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225))
        ])

    def __call__(self, x):
        return self.model(self.transform(x))
