import xarray
import tiffslide
import dask.array as da
from pathlib import Path
from xarray import DataArray
from spatialdata import SpatialData
from spatial_image import SpatialImage
from spatialdata.models import Image2DModel
from spatialdata.transformations import Identity

def read_wsi(path: Path) -> SpatialData:
    image_name = Path(path).absolute().name.split(".")[0]
    tiff = tiffslide.open_slide(path)
    img = xarray.open_zarr(
        tiff.zarr_group.store,
        consolidated=False,
        mask_and_scale=False,
    )

    tiff_metadata = {
        'properties': tiff.properties,
        'dimensions': tiff.dimensions,
        'level_count': tiff.level_count,
        'level_dimensions': tiff.level_dimensions,
        'level_downsamples': tiff.level_downsamples
    }

    images = {}
    for k in list(img.keys()):
        ap = k if k!='0' else ''
        
        simg = SpatialImage(
            img[k].transpose("S","Y"+ap,"X"+ap),
            dims=('c','y','x'),
            attrs={'metadata': tiff_metadata}
        )

        images[k] = Image2DModel.parse(
            simg,
            transformations={"pixels": Identity()},
            c_coords=('r','g','b')
        )

    return SpatialData(images=images)

if __name__ == '__main__':
    # test file here: https://openslide.cs.cmu.edu/download/openslide-testdata/Hamamatsu/CMU-1.ndpi
    img = read_wsi('CMU-1.ndpi')
    print(img)