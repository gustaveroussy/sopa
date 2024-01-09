## Save the `SpatialData` object


For this tutorial, we use a generated dataset. The command below will generate it and save it on-disk (you can change the path `tuto.zarr` to save it somewhere else). See [here](`../../cli/#sopa-read`) for details to use your own technology.

```sh
# this generates a 'tuto.zarr' directory
sopa read . --sdata-path tuto.zarr --technology uniform
```

!!! Note
    This generates a `.zarr` directory corresponding to a [`SpatialData` object](https://github.com/scverse/spatialdata).

## (Optional) ROI selection

Sometimes, your slide may contain a region with low quality data, and we want to run the analysis only on the good quality region. For this, we can interactively select a region of interest (ROI), and Sopa will only run on the selected ROI.

=== "If working locally"
    Run the following command line, and follow the instructions displayed in the console:
    ```sh
    sopa crop --sdata-path tuto.zarr --channels DAPI
    ```

=== "If working on a machine without interative mode"
    When interactive mode is not available, the ROI selection will be performed in three steps.

    1. On the machine where the data is stored, save a light resized view of the original image (here, it will create a file called `image.zarr.zip`):
    ```sh
    sopa crop --sdata-path tuto.zarr --channels DAPI --intermediate-image image.zarr.zip
    ```

    2. Download the `image.zip` file locally (or on a machine with interactive mode), and select the ROI. Here, it will create a file called `roi.zarr.zip`:
    ```sh
    sopa crop --intermediate-image image.zarr.zip --intermediate-polygon roi.zarr.zip
    ```

    3. Upload the `roi.zarr.zip` file, and save it inside the `SpatialData` object:
    ```sh
    sopa crop --sdata-path tuto.zarr --intermediate-polygon roi.zarr.zip
    ```

## Run segmentation

### Option 1: Cellpose

Then, generate the bounding boxes of the patches on which Cellpose will be run. Here, the patches have a width and height of 1500 pixels, and an overlap of 50 pixels. We advise bigger sizes for real datasets (see our default parameters in one of our [config files](https://github.com/gustaveroussy/sopa/tree/master/workflow/config)). On the toy dataset, this will generate **4** patches.

```sh
sopa patchify image tuto.zarr --patch-width-pixel 1500 --patch-overlap-pixel 50
```

Now, we can run Cellpose on each individual patch, and for each "segmentation step" we want. On this toy example, we run 3 steps (don't forget to execute the three steps), with (i) DAPI + CK, (ii) DAPI + CD3, and (iii) DAPI + CD20.

!!! Advice
    Running manually the commands below can involve using many consecutive command, so we recommend automatizing it. For instance, this can be done using Snakemake or Nextflow. Mainly, this will help you parallelizing it, since you can run each task on seperate jobs, or using multithreading. You can also see how we do it in the [Sopa Snakemake pipeline](https://github.com/gustaveroussy/sopa/blob/master/workflow/Snakefile).

    To automatically get the number of patches, you can either open the `tuto.zarr/.sopa_cache/patches_file_image` file, or compute `len(sdata['sopa_patches'])` in Python.

Execute the following command line on all `patch-index` (i.e., `0`, `1`, `2`, and `3`) to run Cellpose using DAPI + CK

```sh
sopa segmentation cellpose tuto.zarr \
    --channels DAPI \
    --patch-dir tuto.zarr/.sopa_cache/cellpose \
    --diameter 35 \
    --min-area 2000 \
    --patch-index 0
```

!!! Note
    In the above commands, the `--diameter` and `--min-area` parameters are specific to the data type we work on. For your own data, consider using the default parameters from one of our [config files](https://github.com/gustaveroussy/sopa/tree/master/workflow/config). Here, `min-area` is in pixels^2.

At this stage, you executed 4 times Cellpose. Now, we need to resolve the conflict, i.e. where boundaries are overlapping due to segmentation on multiple patches.
```sh
sopa resolve cellpose tuto.zarr --patch-dir tuto.zarr/.sopa_cache/cellpose
```

### Option 2: Baysor

## Aggregation

To turn the data into an `AnnData` object, we can do count the transcript inside each cell, and/or average each channel intensity inside each cell boundary.

count the transcript inside each cell (by providing the name of the points dataframe, see `--gene-column genes` below), average the channels intensities inside each cell (using `average-intensities`).
```sh
sopa aggregate tuto.zarr --gene-column genes --average-intensities
```

## Annotation

Currently, we support Tangram for transcript-based annotation, and a simple scoring approach for channel-based annotation (called channel z-score).

=== "Tangram annotation"
    ...
=== "Channel Z-score annotation"
    ...   


## Pipeline report

You can create an HTML report of the pipeline run (on the example below, we save it under `report.html`). It contains some quality controls about your data.

```sh
sopa report tuto.zarr report.html
```

## Visualization (Xenium Explorer)
The Xenium Explorer is a software developed by 10X Genomics for visualizing spatial data, and it can be downloaded freely [here](https://www.10xgenomics.com/support/software/xenium-explorer/latest). Sopa allows the convertion to the Xenium Explorer, whatever the type of spatial data you worked on.

```sh
sopa explorer write tuto.zarr --gene-column genes
```

If you have downloaded the Xenium Explorer, you can now open the results in the explorer: `open tuto.explorer/experiment.xenium` (if using a Unix operating system), or double click on the latter file.

!!! note "Time efficiency"
    Creating the image needed by the Xenium Explorer can be time consuming. Therefore, we recommend to perform one run for the image generation (below) and another to save the transcripts/boundaries/observations.
    ```sh
    # this can be done directly after saving the raw data in a .zarr directory
    sopa explorer write tuto.zarr --mode '+i' --no-save-h5ad
    ```

    After running everything with Sopa, you can finally save all the other Xenium Explorer input (e.g. boundaries and cell categories):
    ```sh
    # this should be done after aggregation and an eventual annotation
    sopa explorer write tuto.zarr --mode '-i'
    ```
    For more details and customization, refer to the [command line helper](../../cli/#sopa-explorer-write).

## Geometric and spatial statistics

All functions to compute geometric and spatial statistics are detailed in the `sopa.stats` [API](../../api/stats). You can also read [this tutorial](../stats).

## Further analysis

- If you are familiar with the [`spatialdata` library](https://github.com/scverse/spatialdata), you can directly use the `tuto.zarr` directory, corresponding to a `SpatialData` object:
```python
import spatialdata

sdata = spatialdata.read_zarr("tuto.zarr")
```
- You can use [Squidpy](https://squidpy.readthedocs.io/en/latest/index.html) which operates on both the `SpatialData` object or the `AnnData` object, or use other tools of the `scverse` ecosystem such as [`scanpy`](https://scanpy.readthedocs.io/en/stable/index.html).
- You can also use the file `tuto.explorer/adata.h5ad` if you prefer the `AnnData` object instead of the full `SpatialData` object.