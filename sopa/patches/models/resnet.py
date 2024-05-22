from torch import nn


class Resnet50Features(nn.Module):
    def __init__(self):
        super().__init__()

        from torchvision import transforms
        from torchvision.models import resnet50

        resnet = resnet50(weights="IMAGENET1K_V2")

        self.features = nn.Sequential(
            resnet.conv1,
            resnet.bn1,
            resnet.relu,
            resnet.maxpool,
            resnet.layer1,
            resnet.layer2,
            resnet.layer3,
            nn.AdaptiveAvgPool2d(1),
            nn.Flatten(start_dim=1),
        )

        self.t = transforms.Compose(
            [transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225))]
        )

    def __call__(self, x):
        x = self.features(self.t(x))
        return x
