import logging
from functools import partial
from pathlib import Path

from spatialdata import SpatialData

from ... import settings
from ..._constants import SopaAttrs, SopaFiles, SopaKeys
from ...utils import (
    delete_transcripts_patches_dirs,
    get_feature_key,
    get_transcripts_patches_dirs,
)
from .._transcripts import _check_transcript_patches, resolve

log = logging.getLogger(__name__)


def proseg(
    sdata: SpatialData,
    input_technology: str = "merscope",
    config: dict | str | None = None,
    key_added: str = SopaKeys.PROSEG_BOUNDARIES,
):
    _check_transcript_patches(sdata)

    prior_shapes_key = None
    if SopaKeys.PRIOR_SHAPES_KEY in sdata.shapes[SopaKeys.TRANSCRIPTS_PATCHES]:
        prior_shapes_key = sdata.shapes[SopaKeys.TRANSCRIPTS_PATCHES][
            SopaKeys.PRIOR_SHAPES_KEY
        ].iloc[0]

    if config is None or not len(config):
        config = _get_default_config(sdata, input_technology, prior_shapes_key)

    proseg_command = _get_proseg_command(config)

    proseg_run = Proseg(
        proseg_command,
    )

    patch_dir = get_transcripts_patches_dirs(sdata)
    if len(patch_dir) > 1:
        raise IndexError(
            "Proseg is fast enough to work on a single patch. Support for multiple patches is not yet implemented. Rerun the transcript patches step with patch_width=None."
        )

    settings.parallelization_backend = None
    proseg_run(patch_dir)

    sdata.attrs[SopaAttrs.BOUNDARIES] = key_added


class Proseg:
    def __init__(
        self, proseg_command: str, input_technology: str, capture_output: bool = True
    ):
        self.proseg_command = proseg_command
        self.input_technology = input_technology
        self.capture_output = capture_output

    def __call__(self, patch_dir: Path):
        import subprocess

        result = subprocess.run(
            self.proseg_command,
            cwd=patch_dir,
            shell=True,
            capture_output=self.capture_output,
        )

        if result.returncode != 0:
            raise subprocess.CalledProcessError(
                returncode=result.returncode,
                cmd=self.proseg_command,
                output=result.stdout,
                stderr=result.stderr,
            )


def _get_proseg_command(config: dict) -> str:
    prior_suffix = (
        f"--cell-id-column {prior_shapes_key}" if prior_shapes_key else ""
    )  # use a prior segmentation

    return f"proseg --{input_technology} transcripts.csv -x x -y y {prior_suffix}"


def _get_default_config(
    sdata: SpatialData,
    input_technology: str,
    prior_shapes_key: str | None,
) -> dict:
    points_key = sdata.attrs.get(SopaAttrs.TRANSCRIPTS)
    assert points_key, (
        f"Transcripts key not found in sdata.attrs['{SopaAttrs.TRANSCRIPTS}'], baysor config can't be inferred."
    )

    feature_key = get_feature_key(sdata[points_key], raise_error=True)

    config = {
        "data": {
            "input_technology": input_technology,
            "x": "x",
            "y": "y",
            "z": "z",
            "gene": str(feature_key),
        },
        "segmentation": {"prior_shapes_key": prior_shapes_key},
    }

    log.info(
        f"The Proseg config was not provided, using the following by default:\n{config}"
    )

    return config
