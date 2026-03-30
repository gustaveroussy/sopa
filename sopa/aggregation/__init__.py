from .bins import aggregate_bins
from .channels import aggregate_channels
from .transcripts import count_transcripts
from .aggregation import aggregate, Aggregator
from .overlay import overlay_segmentation


def nimbus(*args, **kwargs):
    """Lazily import Nimbus aggregation."""
    try:
        from .nimbus import nimbus_aggregation as _nimbus_aggregation
    except ModuleNotFoundError as e:
        if e.name and e.name.startswith("nimbus_inference"):
            raise ModuleNotFoundError(
                "Nimbus support requires the optional `nimbus-inference` dependency. "
                "Install a compatible Nimbus stack for your environment before calling `nimbus_aggregation`."
            ) from e
        raise

    return _nimbus_aggregation(*args, **kwargs)
