from __future__ import annotations
from typing import TYPE_CHECKING, Any

from .bins import aggregate_bins
from .channels import aggregate_channels
from .transcripts import count_transcripts
from .aggregation import aggregate, Aggregator
from .overlay import overlay_segmentation

if TYPE_CHECKING:
    from .nimbus import nimbus_aggregation as nimbus

__all__ = [
    "Aggregator",
    "aggregate",
    "aggregate_bins",
    "aggregate_channels",
    "count_transcripts",
    "nimbus",
    "overlay_segmentation",
]


def __getattr__(name: str) -> Any:
    if name == "nimbus":
        try:
            from .nimbus import nimbus_aggregation
        except ModuleNotFoundError as e:
            if e.name and e.name.startswith("nimbus"):
                raise ModuleNotFoundError(
                    "You need to install `nimbus_inference` to use `nimbus`."
                    "You can install it via `pip install Nimbus-Inference`."
                ) from e
            raise
        return nimbus_aggregation
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
