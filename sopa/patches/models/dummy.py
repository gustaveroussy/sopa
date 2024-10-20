import torch
from torch import nn


class DummyFeatures(nn.Module):
    def __init__(self):
        super().__init__()

    def __call__(self, x):
        return torch.randn(len(x), 8)
