import xarray
import tiffslide
import dask.array as da
from pathlib import Path
from xarray import DataArray

def read_wsi(path: Path) -> DataArray:
    image_name = Path(path).absolute().name.split(".")[0]
    tiff = tiffslide.open_slide(path)
    image = xarray.open_zarr(
        tiff.zarr_group.store,
        consolidated=False,
        mask_and_scale=False,
    )

    image.attrs['properties'] = tiff.properties
    image.attrs['dimensions'] = tiff.dimensions
    image.attrs['level_count'] = tiff.level_count
    image.attrs['level_dimensions'] = tiff.level_dimensions
    image.attrs['level_downsamples'] = tiff.level_downsamples

    return image

if __name__ == '__main__':
    # test file here: https://openslide.cs.cmu.edu/download/openslide-testdata/Hamamatsu/CMU-1.ndpi
    img = read_wsi('CMU-1.ndpi')
    print(img)