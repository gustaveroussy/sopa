from .conch import ConchFeatures
from .dinov2 import DINOv2Features
from .dummy import DummyFeatures
from .histo_ssl import HistoSSLFeatures
from .hoptimus0 import HOPTIMUSFeatures
from .resnet import Resnet50Features

__all__ = [
    "ConchFeatures",
    "DINOv2Features",
    "DummyFeatures",
    "HOPTIMUSFeatures",
    "HistoSSLFeatures",
    "Resnet50Features",
]

available_models = {
    "conch": ConchFeatures,
    "dinov2": DINOv2Features,
    "dummy": DummyFeatures,
    "hoptimus0": HOPTIMUSFeatures,
    "histo_ssl": HistoSSLFeatures,
    "resnet50": Resnet50Features,
}
