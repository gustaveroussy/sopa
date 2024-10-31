from __future__ import annotations

import logging

from ..aggregation import Aggregator as _Aggregator

log = logging.getLogger(__name__)


def overlay_segmentation(*args, **kwargs):
    from .. import overlay_segmentation as _overlay_segmentation

    log.warning("overlay_segmentation is deprecated, use `sopa.overlay_segmentation` instead")
    _overlay_segmentation(*args, **kwargs)


class Aggregator(_Aggregator):
    def __init__(self, *args, **kwargs):
        log.warning("Aggregator is deprecated, use `sopa.aggregate` instead")
        super().__init__(*args, **kwargs)
