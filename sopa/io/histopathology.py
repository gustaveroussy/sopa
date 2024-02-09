from pathlib import Path

import xarray
from multiscale_spatial_image import MultiscaleSpatialImage
from spatial_image import SpatialImage
from spatialdata import SpatialData
from spatialdata.models import Image2DModel
from spatialdata.transformations import Scale


def wsi(path: Path) -> SpatialData:
    import tiffslide

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
        scale_factor = tiff.level_downsamples[i]
        suffix = k if k != "0" else ""

        scale_image = SpatialImage(
            img[k].transpose("S", f"Y{suffix}", f"X{suffix}"),
            dims=("c", "y", "x"),
            attrs={"metadata": tiff_metadata},
        ).chunk((3, 256, 256))

        images[f"scale{k}"] = Image2DModel.parse(
            scale_image,
            transformations={"pixels": Scale([scale_factor, scale_factor], axes=("x", "y"))},
            c_coords=("r", "g", "b"),
        )

    multiscale_image = MultiscaleSpatialImage.from_dict(images)

    return SpatialData(images={image_name: multiscale_image})
