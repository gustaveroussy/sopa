from pathlib import Path

import dask.array as da
import tiffslide
import xarray
from multiscale_spatial_image import MultiscaleSpatialImage
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Scale
from xarray import DataArray


def read_wsi(path: Path) -> SpatialData:
    image_name = Path(path).absolute().name.split(".")[0]
    tiff = tiffslide.open_slide(path)
    img = xarray.open_zarr(
        tiff.zarr_group.store,
        consolidated=False,
        mask_and_scale=False,
    )

    tiff_metadata = {
        "properties": tiff.properties,
        "dimensions": tiff.dimensions,
        "level_count": tiff.level_count,
        "level_dimensions": tiff.level_dimensions,
        "level_downsamples": tiff.level_downsamples,
    }

    images = {}
    for i, k in enumerate(list(img.keys())):
        ap = k if k != "0" else ""

        simg = SpatialImage(
            img[k].transpose("S", "Y" + ap, "X" + ap),
            dims=("c", "y", "x"),
            attrs={"metadata": tiff_metadata},
        ).chunk({"c": 3, "y": 256, "x": 256})

        sf = tiff.level_downsamples[i]
        images[f"scale{k}"] = Image2DModel.parse(
            simg,
            transformations={"pixels": Scale([sf, sf], axes=("x", "y"))},
            c_coords=("r", "g", "b"),
        )

    mimg = MultiscaleSpatialImage.from_dict(images)

    return SpatialData(images={image_name: mimg})


if __name__ == "__main__":
    # test file here: https://openslide.cs.cmu.edu/download/openslide-testdata/Hamamatsu/CMU-1.ndpi
    img = read_wsi("CMU-1.ndpi")
    print(img)
