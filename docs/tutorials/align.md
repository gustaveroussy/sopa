If you have an H&E image or immunofluorescence data that you want to align on the main image, this can be done with the Xenium Explorer, even if you don't work on Xenium data.

## Image conversion
Convert your image with QuPath as written in this [10x genomics webpage](https://www.10xgenomics.com/support/software/xenium-explorer/tutorials/xe-image-file-conversion).

If you are not familiar with QuPath, you can also use our API to write the image:
```python
import sopa.io

image = sopa.io.ome_tif("path/to/your/image.tif") # or use a different reader
sopa.io.write_image("where/to/save/image.ome.tif", image, is_dir=False)
```

!!! note "Xenium users"
    If using the Xenium machine, then you don't need conversion; the images provided by the Xenium machine already have the correct format.

## Open your data in the Xenium Explorer

If not done yet, convert your `SpatialData` object to the Xenium Explorer's inputs. This can be done as detailed in [this tutorial](../cli_usage/#visualization-xenium-explorer).

Double-click on the `experiment.xenium` file, or select it from the Xenium Explorer interface. It will display the data in the explorer.

## Keypoint placement

!!! warning
    Make sure your Xenium Explorer version is at least `1.3.0`

On the Xenium Explorer, under the "Images" panel, click "Add image" and follow the instructions on the screen.

<p align="center">
  <img src="../../assets/explorer/add_image.png" alt="add_image" width="300px"/>
</p>

Afterwards, the explorer will automatically align the images based on the key points you selected on both images.

## (Optional) Update the `SpatialData` object

After alignment, export the transformation matrix as a `.csv` file. For that, select your aligned image under the "Images" panel and click on "Download Alignment File":

<p align="center">
  <img src="../../assets/explorer/download_alignment.png" alt="add_image" width="300px"/>
</p>

Then, select only the "Transformation Matrix" box and download it:

<p align="center">
  <img src="../../assets/explorer/download_transformation_file.png" alt="add_image" width="400px"/>
</p>

Then, use the [CLI](../../cli/#sopa-explorer-add-aligned) to update your `SpatialData` object. You'll need the path to the `.zarr` directory corresponding to your `SpatialData` object (`SDATA_PATH`), the path to the `.ome.tif` image that you converted above (`IMAGE_PATH`), and the `.csv` transformation matrix that you exported from the Xenium Explorer (`TRANSFORMATION_MATRIX_PATH`):

```sh
sopa explorer add-aligned <SDATA_PATH> <IMAGE_PATH> <TRANSFORMATION_MATRIX_PATH>
```
