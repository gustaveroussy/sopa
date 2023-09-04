from pathlib import Path
from typing import List

import tifffile as tf


def write_ome_tif(
    output_path: Path,
    series: List[tf.TiffPageSeries],
    channel_names: List[str],
    pixelsize: float = 0.2125,
) -> None:
    with tf.TiffWriter(output_path, bigtiff=True) as tif:
        metadata = {
            "SignificantBits": 8,
            "PhysicalSizeX": pixelsize,
            "PhysicalSizeXUnit": "µm",
            "PhysicalSizeY": pixelsize,
            "PhysicalSizeYUnit": "µm",
            "Channel": {"Name": channel_names},
        }
        options = dict(
            photometric="minisblack",
            tile=(1024, 1024),
            compression="jpeg2000",
            resolutionunit="CENTIMETER",
        )
        tif.write(
            series[0].asarray(),
            subifds=len(series) - 1,
            resolution=(1e4 / pixelsize, 1e4 / pixelsize),
            metadata=metadata,
            **options,
        )

        for i in range(1, len(series)):
            print(f"Writing sub-resolution {i}...")
            tif.write(
                series[i].asarray(),
                subfiletype=1,
                resolution=(1e4 * 2**i / pixelsize, 1e4 * 2**i / pixelsize),
                **options,
            )
