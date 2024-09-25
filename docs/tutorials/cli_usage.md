Here, we provide a minimal example of command line usage. For more details and to learn about other optional arguments, refer to the full [CLI documentation](../../cli).

!!! tip
    The [Snakemake pipeline](https://gustaveroussy.github.io/sopa/tutorials/snakemake/) is recommended to get started with Sopa. Using the CLI is advised if you want more flexibility, but you'll need to parallelize the segmentation yourself, as detailed below.

## Save the `SpatialData` object

For this tutorial, we use a generated dataset. You can expect a total runtime of a few minutes.

The command below will generate and save it on disk (you can change the path `tuto.zarr` to save it somewhere else). If you want to load your own data: choose the right panel below. For more information, refer to this [FAQ](../../faq/#what-kind-of-inputs-do-i-need-to-run-sopa) describing which data input you need, or see the [`sopa read`](`../../cli/#sopa-read`) documentation.

=== "Tutorial"
    ```sh
    # it will generate a 'tuto.zarr' directory
    sopa read . --sdata-path tuto.zarr --technology uniform
    ```
=== "Xenium"
    ```sh
    # it will generate a '/path/to/sample/directory.zarr' directory
    sopa read /path/to/sample/directory --technology xenium
    ```
=== "MERSCOPE"
    ```sh
    # it will generate a '/path/to/sample/directory.zarr' directory
    sopa read /path/to/sample/directory --technology merscope
    ```
=== "CosMX"
    ```sh
    # it will generate a '/path/to/sample/directory.zarr' directory
    sopa read /path/to/sample/directory --technology cosmx
    ```
=== "PhenoCycler"
    ```sh
    # it will generate a '/path/to/sample.zarr' directory
    sopa read /path/to/sample.qptiff --technology phenocycler
    ```
=== "MACSima"
    ```sh
    # it will generate a '/path/to/sample/directory.zarr' directory
    sopa read /path/to/sample/directory --technology macsima
    ```
=== "Hyperion"
    ```sh
    # it will generate a '/path/to/sample/directory.zarr' directory
    sopa read /path/to/sample/directory --technology hyperion
    ```


!!! info
    It has created a `.zarr` directory which stores a [`SpatialData` object](https://github.com/scverse/spatialdata) corresponding to your data sample. You can choose the location of the `.zarr` directory using the `--sdata-path` command line argument.

## (Optional) ROI selection

Sometimes, your slide may contain a region with low-quality data, and we want to run the analysis only on the good-quality region. For this, we can interactively select a region of interest (ROI), and Sopa will only run on the selected ROI.

=== "If working locally"
    Run the following command line and follow the instructions displayed in the console:
    ```sh
    sopa crop --sdata-path tuto.zarr --channels DAPI
    ```

=== "If working on a machine without interactive mode"
    When interactive mode is unavailable, the ROI selection will be performed in three steps.

    1. On the machine where the data is stored, save a light view of the original image (here, it will create a file called `image.zarr.zip`):
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

First, generate the bounding boxes of the patches on which Cellpose will be run. Here, the patches have a width and height of 1500 pixels and an overlap of 50 pixels. We advise bigger sizes for real datasets (see our default parameters in one of our [config files](https://github.com/gustaveroussy/sopa/tree/master/workflow/config)). On the toy dataset, this will generate **4** patches.

```sh
sopa patchify image tuto.zarr --patch-width-pixel 1500 --patch-overlap-pixel 50
```

Now, we can run Cellpose on each individual patch. Execute the following command line on all `patch-index` (i.e., `0`, `1`, `2`, and `3`) to run Cellpose using DAPI only (you can add an additional channel, for instance, `--channels DAPI --channels PolyT`):

!!! tip
    Manually running the commands below can involve using many consecutive commands, so we recommend automatizing it. For instance, this can be done using Snakemake or Nextflow. This will help you parallelize it since you can run each task on separate jobs or using multithreading. You can also see how we do it in our [Snakefile](https://github.com/gustaveroussy/sopa/blob/master/workflow/Snakefile). If you prefer using the already existing pipeline instead of the CLI, you can read our [Snakemake pipeline tutorial](https://gustaveroussy.github.io/sopa/tutorials/snakemake/).

    To automatically get the number of patches, you can either open the `tuto.zarr/.sopa_cache/patches_file_image` file, or compute `len(sdata['sopa_patches'])` in Python.

=== "Patch 0"
    ```sh
    sopa segmentation cellpose tuto.zarr \
        --channels DAPI \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 0
    ```
=== "Patch 1"
    ```sh
    sopa segmentation cellpose tuto.zarr \
        --channels DAPI \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 1
    ```
=== "Patch 2"
    ```sh
    sopa segmentation cellpose tuto.zarr \
        --channels DAPI \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 2
    ```
=== "Patch 3"
    ```sh
    sopa segmentation cellpose tuto.zarr \
        --channels DAPI \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 3
    ```

!!! Note
    In the above commands, the `--diameter` and `--min-area` parameters are specific to the data type we work on. For your own data, consider using the default parameters from one of our [config files](https://github.com/gustaveroussy/sopa/tree/master/workflow/config). Here, `min-area` is in pixels^2.

At this stage, you executed 4 times Cellpose (once per patch). Now, we need to resolve the conflict, i.e. where boundaries are overlapping due to segmentation on multiple patches:
```sh
sopa resolve cellpose tuto.zarr
```

### Option 2: Baysor

Baysor needs a config to be executed. You can find official config examples [here](https://github.com/kharchenkolab/Baysor/tree/master/configs).

!!! note
    You can also reuse the Baysor parameter we have defined for each machine, as in our [Snakemake config files](https://github.com/gustaveroussy/sopa/tree/master/workflow/config). Note that, our Snakemake config is a `.yaml` file, but the Baysor config should still be a `.toml` file.

For this tutorial, we will use the config below. Save this in a `config.toml` file.
```toml
[data]
force_2d = true
min_molecules_per_cell = 10
x = "x"
y = "y"
z = "z"
gene = "genes"
min_molecules_per_gene = 0
min_molecules_per_segment = 3
confidence_nn_id = 6

[segmentation]
scale = 30                          # typical cell radius
scale_std = "25%"                   # cell radius standard deviation
prior_segmentation_confidence = 0
estimate_scale_from_centers = false
n_clusters = 4
iters = 500
n_cells_init = 0
nuclei_genes = ""
cyto_genes = ""
```

Then, we generate the bounding boxes of the patches on which Baysor will be run. Here, the patches have a width and height of 1200 microns and an overlap of 50 microns. We advise bigger sizes for real datasets (see our default parameters in one of our [config files](https://github.com/gustaveroussy/sopa/tree/master/workflow/config)). On the toy dataset, this will generate **4** patches.

```sh
# config.toml is the Baysor config file you generated above
sopa patchify baysor tuto.zarr --config-path config.toml --patch-width-microns 1200 --patch-overlap-microns 50
```

Now, we can run Baysor on each individual patch. Execute the following command lines to run Baysor on each patch (i.e., `0`, `1`, `2`, and `3`).

!!! tip
    Manually running the commands below can involve using many consecutive commands, so we recommend automatizing it. For instance, this can be done using Snakemake or Nextflow. This will help you parallelize it since you can run each task on separate jobs or using multithreading. You can also see how we do it in the [Sopa Snakemake pipeline](https://github.com/gustaveroussy/sopa/blob/master/workflow/Snakefile).

    To automatically get the number of patches, you can open the `tuto.zarr/.sopa_cache/patches_file_baysor` file. This lists the names of the directories inside `tuto.zarr/.sopa_cache/baysor` related to each patch. If you selected an ROI, the excluded patches are effectively not in the `patches_file_baysor` file.

=== "Patch 0"
    ```sh
    cd tuto.zarr/.sopa_cache/baysor_boundaries/0

    # 'baysor' is the official baysor executable. If unavailable, replace it with your path to the executable
    baysor run --save-polygons GeoJSON -c config.toml transcripts.csv
    ```
=== "Patch 1"
    ```sh
    cd tuto.zarr/.sopa_cache/baysor_boundaries/1

    # 'baysor' is the official baysor executable. If unavailable, replace it with your path to the executable
    baysor run --save-polygons GeoJSON -c config.toml transcripts.csv
    ```
=== "Patch 2"
    ```sh
    cd tuto.zarr/.sopa_cache/baysor_boundaries/2

    # 'baysor' is the official baysor executable. If unavailable, replace it with your path to the executable
    baysor run --save-polygons GeoJSON -c config.toml transcripts.csv
    ```
=== "Patch 3"
    ```sh
    cd tuto.zarr/.sopa_cache/baysor_boundaries/3

    # 'baysor' is the official baysor executable. If unavailable, replace it with your path to the executable
    baysor run --save-polygons GeoJSON -c config.toml transcripts.csv
    ```

At this stage, you executed 4 times Baysor (once per patch). Now, we need to resolve the conflict, i.e. where boundaries are overlapping due to segmentation on multiple patches:
```sh
sopa resolve baysor tuto.zarr --gene-column genes
```

## Aggregation

This **mandatory** step turns the data into an `AnnData` object. We can count the transcript inside each cell, and/or average each channel intensity inside each cell boundary.

!!! info
    The `--gene-column` option below tells which column contains the gene names inside the transcript dataframe. If you don't know it, you can look to [our configs](https://github.com/gustaveroussy/sopa/tree/master/workflow/config) to find the right `gene-column` corresponding to your machine.

=== "Count transcripts + average intensities"
    ```sh
    sopa aggregate tuto.zarr --gene-column genes --average-intensities --min-transcripts 10
    ```
=== "Count transcripts"
    ```sh
    sopa aggregate tuto.zarr --gene-column genes --min-transcripts 10
    ```
=== "Average intensities"
    ```sh
    sopa aggregate tuto.zarr --average-intensities
    ```

!!! note "If using Baysor"
    Baysor already counts the transcripts inside each cell to create a cell-by-gene table. So you'll always have this table, and there is no need to use the `--gene-column` argument. If you don't want to average the intensities, you will still need to run `sopa aggregate tuto.zarr` before continuing.

## Annotation

If desired, cell-type annotation can be run. Currently, we support Tangram for transcript-based annotation and a simple scoring approach for channel-based annotation (called channel z-score).

=== "Channel Z-score annotation"
    For now, our fluorescence-based annotation is very simple. We provide a dictionary where a channel is associated with a population. Then, each cell is associated with the cell type whose corresponding channel is the brightest (according to a certain Z-score). In this tutorial example, we can annotate Tumoral cells, T cells, and B cells:
    ```sh
    sopa annotate fluorescence tuto.zarr --marker-cell-dict '{"CK": "Tumoral cell", "CD3": "T cell", "CD20": "B cell"}'
    ```
    !!! note "More complex annotation"
        If you have a large number of channels, it may be preferable to run clustering on your data, for instance, using [`Leiden clustering](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html). Then, you can annotate each cluster manually by plotting a heatmap of all channels expressions per cluster.
=== "Tangram annotation"
    [Tangram](https://github.com/broadinstitute/Tangram) is a transcript-based annotation that uses an annotated single-cell reference. Let's suppose your reference `AnnData` object is stored in a file called `adata_reference.h5ad` (preferably, keep raw counts), and the cell type is in `adata.obs["cell_type"]`. Then, you can annotate your spatial data as follows:
    ```sh
    sopa annotate tangram tuto.zarr --sc-reference-path adata_reference.h5ad --cell-type-key cell_type
    ```


## Pipeline report

You can optionally create an HTML report of the pipeline run (in the example below, we save it under `report.html`). It contains some quality controls for your data.

```sh
sopa report tuto.zarr report.html
```

## Visualization (Xenium Explorer)
The Xenium Explorer is a software developed by 10X Genomics for visualizing spatial data, and it can be downloaded freely [here](https://www.10xgenomics.com/support/software/xenium-explorer/latest). Sopa allows the conversion to the Xenium Explorer, whatever the type of spatial data you worked on. It will create some files under a new `tuto.explorer` directory:

```sh
sopa explorer write tuto.zarr --gene-column genes
```

If you have downloaded the Xenium Explorer, you can now open the results in the explorer: `open tuto.explorer/experiment.xenium` (if using a Unix operating system), or double-click on the latter file.

!!! info "Time efficiency"
    Creating the image needed by the Xenium Explorer can be time-consuming. Therefore, we recommend performing one run for the image generation (below) and another to save the transcripts/boundaries/observations.
    ```sh
    # this can be done directly after saving the raw data in a .zarr directory
    sopa explorer write tuto.zarr --mode "+i" --no-save-h5ad
    ```

    After running everything with Sopa, you can finally save all the other Xenium Explorer input (e.g. boundaries and cell categories):
    ```sh
    # this should be done after aggregation and an eventual annotation
    sopa explorer write tuto.zarr --mode "-i" --gene-column genes
    ```
    For more details and customization, refer to the [command line helper](../../cli/#sopa-explorer-write).

## Geometric and spatial statistics

All functions to compute geometric and spatial statistics are detailed in the `sopa.spatial` [API](../../api/spatial). You can also read [this tutorial](../spatial).

## Further analyses

- If you are familiar with the [`spatialdata` library](https://github.com/scverse/spatialdata), you can directly use the `tuto.zarr` directory, corresponding to a `SpatialData` object:
```python
import spatialdata

sdata = spatialdata.read_zarr("tuto.zarr")
```
- You can use [Squidpy](https://squidpy.readthedocs.io/en/latest/index.html) which operates on both the `SpatialData` object or the `AnnData` object, or use other tools of the `scverse` ecosystem such as [`Scanpy`](https://scanpy.readthedocs.io/en/stable/index.html).
- You can also use the file `tuto.explorer/adata.h5ad` if you prefer the `AnnData` object instead of the full `SpatialData` object.
