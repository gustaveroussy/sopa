from __future__ import annotations

import logging
import warnings

from ..aggregation import Aggregator as _Aggregator

log = logging.getLogger(__name__)


def overlay_segmentation(*args, **kwargs):
    from .. import overlay_segmentation as _overlay_segmentation

    warnings.warn(
        "overlay_segmentation is deprecated, use `sopa.overlay_segmentation` instead", DeprecationWarning, stacklevel=2
    )
    _overlay_segmentation(*args, **kwargs)


class Aggregator(_Aggregator):
    def __init__(self, *args, **kwargs):
        warnings.warn("Aggregator is deprecated, use `sopa.aggregate` instead", DeprecationWarning, stacklevel=2)
        super().__init__(*args, **kwargs)
