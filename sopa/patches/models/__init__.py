from .resnet import Resnet50Features
from .histo_ssl import HistoSSLFeatures
from .dinov2 import DINOv2Features
from .dummy import DummyFeatures

__all__ = ["Resnet50Features", "HistoSSLFeatures", "DINOv2Features", "DummyFeatures"]

available_models = {
    "resnet50": Resnet50Features,
    "histo_ssl": HistoSSLFeatures,
    "dinov2": DINOv2Features,
    "dummy": DummyFeatures,
}
