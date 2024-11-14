from torch import nn


class DummyFeatures(nn.Module):
    def __init__(self):
        super().__init__()

    def __call__(self, x):
        """(B, C, Y, X) -> (B, C)"""
        return x[:, :, 0, 0]
