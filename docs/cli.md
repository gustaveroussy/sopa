# CLI (command-line-interface)

## Usage

When installing `sopa` as written in our [getting-started guidelines](../getting_started), a new command named `sopa` becomes available.

!!! note "CLI helper"
    Run `sopa --help` to get details about all the command line purposes. You can also use this helper on any subcommand, for instance, `sopa read --help`.

<div class="termy">
```console
// Run the Sopa CLI helper
$ sopa --help
 Usage: sopa [OPTIONS] COMMAND [ARGS]...
╭─ Commands ─────────────────────────────────────────────────────╮
│ aggregate     Aggregate transcripts/channels inside cells      │
│ annotate      Perform cell-type annotation                     │
│ check         Run some sanity checks                           │
│ crop          Crop an image based on a user-defined polygon    │
│ explorer      Conversion to the Xenium Explorer's inputs       │
│ patchify      Create patches with overlaps                     │
│ read          Read any technology + write a SpatialData object │
│ report        Create a web-report with figures/QCs             │
│ resolve       Resolve the segmentation conflicts over patches  │
│ segmentation  Perform cell segmentation on patches             │
╰────────────────────────────────────────────────────────────────╯
// Example: run cellpose segmentation
$ sopa segmentation cellpose sdata.zarr
... [Logs] ...
```
</div>

### Notes

If you don't know in which order to run these commands, refer to the image in the [homepage](..), or see our [CLI usage tutorial](../tutorials/cli_usage).

When running the `sopa` CLI, some arguments are required, while some are optional. For instance, for the `sopa read` command, `sdata_path` is an argument, and a path has to be given directly, while `technology` is an option (in this case, the `--technology` prefix has to be used). For instance, if you read MERSCOPE data, it will be:

```
sopa read /path/to/merscope/directory --technology merscope
```

Note that `/path/to/merscope/directory` refers to `sdata_path`, which is an argument. You **should not** add the suffix `--sdata_path`, as it is an argument.

!!! note "Required options"
    All the arguments are required, as shown by the `[required]` hint on the CLI helper. Note that some options may also be required too (in this case, the term `[required]` will appear on the CLI helper). But they still need to be called as a normal option.

## `sopa` commands
**Usage**:

```console
$ sopa [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `aggregate`: Create an `anndata` table containing the...
* `annotate`: Perform cell-type annotation (based on...
* `check`: Run some sanity checks (e.g., on the YAML...
* `crop`: Crop an image based on a user-defined...
* `explorer`: Convertion to the Xenium Explorer's...
* `patchify`: Create patches with overlaps.
* `read`: Read any technology data, and write a...
* `report`: Create a HTML report of the pipeline run...
* `resolve`: Resolve the segmentation conflicts over...
* `segmentation`: Perform cell segmentation on patches.

### `sopa aggregate`

Create an `anndata` table containing the transcript count and/or the channel intensities per cell

**Usage**:

```console
$ sopa aggregate [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--gene-column TEXT`: Column of the transcript dataframe representing the gene names. If not provided, it will not compute transcript count
* `--average-intensities / --no-average-intensities`: Whether to average the channel intensities inside each cell  [default: no-average-intensities]
* `--expand-radius-ratio FLOAT`: Cells polygons will be expanded by `expand_radius_ratio * mean_radius` for channels averaging **only**. This help better aggregate boundary stainings  [default: 0]
* `--min-transcripts INTEGER`: Cells with less transcript than this integer will be filtered  [default: 0]
* `--min-intensity-ratio FLOAT`: Cells whose mean channel intensity is less than `min_intensity_ratio * quantile_90` will be filtered  [default: 0]
* `--image-key TEXT`: Optional image key of the SpatialData object. By default, considers the only one image. It can be useful if another image is added later on
* `--method-name TEXT`: If segmentation was performed with a generic method, this is the name of the method used.
* `--help`: Show this message and exit.

### `sopa annotate`

Perform cell-type annotation (based on transcripts and/or channel intensities)

**Usage**:

```console
$ sopa annotate [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `fluorescence`: Simple annotation based on fluorescence,...
* `tangram`: Tangram segmentation (i.e., uses an...

#### `sopa annotate fluorescence`

Simple annotation based on fluorescence, where each provided channel corresponds to one cell type.

For each cell, one z-score statistic is computed and the population with the highest z-score is attributed.

**Usage**:

```console
$ sopa annotate fluorescence [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--marker-cell-dict TEXT`: [required]
* `--cell-type-key TEXT`: Key added in `adata.obs` corresponding to the cell type  [default: cell_type]
* `--help`: Show this message and exit.

#### `sopa annotate tangram`

Tangram segmentation (i.e., uses an annotated scRNAseq reference to transfer cell-types)

**Usage**:

```console
$ sopa annotate tangram [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--sc-reference-path TEXT`: Path to the scRNAseq annotated reference, as a `.h5ad` file  [required]
* `--cell-type-key TEXT`: Key of `adata_ref.obs` containing the cell-types  [required]
* `--reference-preprocessing TEXT`: Preprocessing method applied to the reference. Either None (raw counts), or `normalized` (sc.pp.normalize_total) or `log1p` (sc.pp.normalize_total and sc.pp.log1p)
* `--bag-size INTEGER`: Number of cells in each bag of the spatial table. Low values will decrease the memory usage  [default: 10000]
* `--max-obs-reference INTEGER`: Maximum samples to be considered in the reference for tangram. Low values will decrease the memory usage  [default: 10000]
* `--help`: Show this message and exit.

### `sopa check`

Run some sanity checks (e.g., on the YAML config, on the tangram reference, ...)

**Usage**:

```console
$ sopa check [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `config`: Perform sanity checks on a sopa yaml config
* `reference`: Perform sanity checks on a tangram...

#### `sopa check config`

Perform sanity checks on a sopa yaml config

**Usage**:

```console
$ sopa check config [OPTIONS] PATH
```

**Arguments**:

* `PATH`: Path to the YAML config  [required]

**Options**:

* `--help`: Show this message and exit.

#### `sopa check reference`

Perform sanity checks on a tangram scRNAseq reference

**Usage**:

```console
$ sopa check reference [OPTIONS] REFERENCE_PATH
```

**Arguments**:

* `REFERENCE_PATH`: Path to the scRNAseq reference as a `.h5ad` file  [required]

**Options**:

* `--cell-type-key TEXT`: Key of adata.obs containing the cell types  [required]
* `--help`: Show this message and exit.

### `sopa crop`

Crop an image based on a user-defined polygon (interactive mode).

!!! note "Usage"

    - [Locally] Only `--sdata-path` is required

    - [On a cluster] Run `sopa crop` with `--sdata-path` and `--intermediate-image` on the cluster. Then, download the image locally, and run `sopa crop` with `--intermediate-image` and `--intermediate-polygon`. Then, upload the polygon and run `sopa crop` on the cluster with `--sdata-path` and `--intermediate-polygon`.

**Usage**:

```console
$ sopa crop [OPTIONS]
```

**Options**:

* `--sdata-path TEXT`: Path to the SpatialData `.zarr` directory
* `--intermediate-image TEXT`: Path to the intermediate image, with a `.zip` extension. Use this only if the interactive mode is not available
* `--intermediate-polygon TEXT`: Path to the intermediate polygon, with a `.zip` extension. Use this locally, after downloading the intermediate_image
* `--channels TEXT`: List of channel names to be displayed. Optional if there are already only 1 or 3 channels
* `--scale-factor FLOAT`: Resize the image by this value (high value for a lower memory usage)  [default: 10]
* `--margin-ratio FLOAT`: Ratio of the image margin on the display (compared to the image size)  [default: 0.1]
* `--help`: Show this message and exit.

### `sopa explorer`

Convertion to the Xenium Explorer's inputs, and image alignment

**Usage**:

```console
$ sopa explorer [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `add-aligned`: After alignment on the Xenium Explorer,...
* `update-obs`: Update the cell categories for the Xenium...
* `write`: Convert a spatialdata object to Xenium...

#### `sopa explorer add-aligned`

After alignment on the Xenium Explorer, add an image to the SpatialData object

**Usage**:

```console
$ sopa explorer add-aligned [OPTIONS] SDATA_PATH IMAGE_PATH TRANSFORMATION_MATRIX_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]
* `IMAGE_PATH`: Path to the image file to be added (`.ome.tif` used in the explorer during alignment)  [required]
* `TRANSFORMATION_MATRIX_PATH`: Path to the `matrix.csv` file returned by the Explorer after alignment  [required]

**Options**:

* `--original-image-key TEXT`: Optional original-image key (of sdata.images) on which the new image will be aligned. This doesn't need to be provided if there is only one image
* `--overwrite / --no-overwrite`: Whether to overwrite the image if existing  [default: no-overwrite]
* `--help`: Show this message and exit.

#### `sopa explorer update-obs`

Update the cell categories for the Xenium Explorer's (i.e. what's in `adata.obs`). This is useful when you perform analysis and update your `AnnData` object

!!! note "Usage"
    Make sure you have already run `sopa explorer write` before. This command should only be used if you updated `adata.obs`

**Usage**:

```console
$ sopa explorer update-obs [OPTIONS] ADATA_PATH OUTPUT_PATH
```

**Arguments**:

* `ADATA_PATH`: Path to the anndata file (`zarr` or `h5ad`) containing the new observations  [required]
* `OUTPUT_PATH`: Path to the Xenium Explorer directory (it will update `analysis.zarr.zip`)  [required]

**Options**:

* `--help`: Show this message and exit.

#### `sopa explorer write`

Convert a spatialdata object to Xenium Explorer's inputs

**Usage**:

```console
$ sopa explorer write [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--output-path TEXT`: Path to a directory where Xenium Explorer's outputs will be saved. By default, writes to the same path as `sdata_path` but with the `.explorer` suffix
* `--gene-column TEXT`: Column name of the points dataframe containing the gene names
* `--shapes-key TEXT`: Sdata key for the boundaries. By default, uses the baysor boundaires, else the cellpose boundaries
* `--pixel-size FLOAT`: Number of microns in a pixel. Invalid value can lead to inconsistent scales in the Explorer.  [default: 0.2125]
* `--pixelsize FLOAT`: `pixelsize` is deprecated and will be removed in future versions. Use `pixel_size` instead.
* `--lazy / --no-lazy`: If `True`, will not load the full images in memory (except if the image memory is below `ram_threshold_gb`)  [default: lazy]
* `--ram-threshold-gb INTEGER`: Threshold (in gygabytes) from which image can be loaded in memory. If `None`, the image is never loaded in memory  [default: 4]
* `--mode TEXT`: String that indicated which files should be created. `'-ib'` means everything except images and boundaries, while `'+tocm'` means only transcripts/observations/counts/metadata (each letter corresponds to one explorer file). By default, keeps everything
* `--save-h5ad / --no-save-h5ad`: Whether to save the adata as h5ad in the explorer directory (for convenience only, since h5ad is faster to open than the original .zarr table)  [default: save-h5ad]
* `--help`: Show this message and exit.

### `sopa patchify`

Create patches with overlaps. Afterwards, segmentation will be run on each patch

**Usage**:

```console
$ sopa patchify [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `baysor`: Prepare patches for transcript-based...
* `comseg`: Prepare patches for transcript-based...
* `image`: Prepare patches for staining-based...

#### `sopa patchify baysor`

Prepare patches for transcript-based segmentation with Baysor

**Usage**:

```console
$ sopa patchify baysor [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--patch-width-microns FLOAT`: Width (and height) of each patch in microns  [required]
* `--patch-overlap-microns FLOAT`: Number of overlapping microns between the patches. We advise to choose approximately twice the diameter of a cell  [required]
* `--baysor-temp-dir TEXT`: Temporary directory where baysor inputs and outputs will be saved. By default, uses `.sopa_cache/baysor_boundaries`
* `--config-path TEXT`: Path to the baysor config (you can also directly provide the argument via the `config` option)
* `--config TEXT`: Dictionnary of baysor parameters, overwrite the config_path argument if provided  [default: {}]
* `--cell-key TEXT`: Optional column of the transcripts dataframe that indicates in which cell-id each transcript is, in order to use prior segmentation Default is 'cell' if cell_key=None
* `--unassigned-value INTEGER`: If --cell-key is provided, this is the value given to transcripts that are not inside any cell (if it's already 0, don't provide this argument)
* `--use-prior / --no-use-prior`: Whether to use cellpose segmentation as a prior for baysor (if True, make sure to first run Cellpose)  [default: no-use-prior]
* `--help`: Show this message and exit.

#### `sopa patchify comseg`

Prepare patches for transcript-based segmentation with ComSeg

**Usage**:

```console
$ sopa patchify comseg [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--patch-width-microns FLOAT`: Width (and height) of each patch in microns  [required]
* `--patch-overlap-microns FLOAT`: Number of overlapping microns between the patches. We advise to choose approximately twice the diameter of a cell  [required]
* `--comseg-temp-dir TEXT`: Temporary directory where baysor inputs and outputs will be saved. By default, uses `.sopa_cache/comseg_boundaries`
* `--config-path TEXT`: Path to the ComSeg json config file (you can also directly provide the argument via the `config` option)
* `--config TEXT`: Dictionnary of ComSeg parameters, overwrite the config_path argument if provided  [default: {}]
* `--cell-key TEXT`: Optional column of the transcripts dataframe that indicates in which cell-id each transcript is, in order to use prior segmentation. Default is cell if cell_key=None
* `--unassigned-value INTEGER`: If --cell-key is provided, this is the value given to transcripts that are not inside any cell (if it's already 0, don't provide this argument)
* `--help`: Show this message and exit.

#### `sopa patchify image`

Prepare patches for staining-based segmentation (including Cellpose)

**Usage**:

```console
$ sopa patchify image [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--patch-width-pixel FLOAT`: Width (and height) of each patch in pixels  [default: 5000]
* `--patch-overlap-pixel FLOAT`: Number of overlapping pixels between the patches. We advise to choose approximately twice the diameter of a cell  [default: 100]
* `--help`: Show this message and exit.

### `sopa read`

Read any technology data, and write a standardized SpatialData object.

Either `--technology` or `--config-path` has to be provided.

**Usage**:

```console
$ sopa read [OPTIONS] DATA_PATH
```

**Arguments**:

* `DATA_PATH`: Path to one data sample (most of the time, this corresponds to a directory with images files and eventually a transcript file)  [required]

**Options**:

* `--technology TEXT`: Name of the technology used to collected the data (`xenium`/`merfish`/`cosmx`/`phenocycler`/`macsima`/`hyperion`)
* `--sdata-path TEXT`: Optional path to write the SpatialData object. If not provided, will write to the `{data_path}.zarr` directory
* `--config-path TEXT`: Path to the snakemake config. This can be useful in order not to provide the `--technology` and the `--kwargs` arguments
* `--kwargs TEXT`: Dictionary provided to the reader function as kwargs  [default: {}]
* `--help`: Show this message and exit.

### `sopa report`

Create a HTML report of the pipeline run and some quality controls

**Usage**:

```console
$ sopa report [OPTIONS] SDATA_PATH PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]
* `PATH`: Path to the HTML report  [required]

**Options**:

* `--help`: Show this message and exit.

### `sopa resolve`

Resolve the segmentation conflicts over patches overlaps

**Usage**:

```console
$ sopa resolve [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `baysor`: Resolve patches conflicts after baysor...
* `cellpose`: Resolve patches conflicts after cellpose...
* `comseg`: Resolve patches conflicts after comseg...
* `generic`: Resolve patches conflicts after generic...

#### `sopa resolve baysor`

Resolve patches conflicts after baysor segmentation. Provide either `--baysor-temp-dir` or `--patches-dirs`

**Usage**:

```console
$ sopa resolve baysor [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--gene-column TEXT`: Column of the transcripts dataframe containing the genes names  [required]
* `--baysor-temp-dir TEXT`: Path to the directory containing all the baysor patches (see `sopa patchify`). By default, uses the `.sopa_cache/baysor_boundaries` directory
* `--min-area FLOAT`: Cells with an area less than this value (in microns^2) will be filtered  [default: 0]
* `--patches-dirs TEXT`: List of patches directories inside `baysor_temp_dir`
* `--help`: Show this message and exit.

#### `sopa resolve cellpose`

Resolve patches conflicts after cellpose segmentation

**Usage**:

```console
$ sopa resolve cellpose [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--patch-dir TEXT`: Directory containing the cellpose segmentation on patches (or multiple directories if using multi-step segmentation). By default, uses the `.sopa_cache/cellpose_boundaries` directory
* `--help`: Show this message and exit.

#### `sopa resolve comseg`

Resolve patches conflicts after comseg segmentation. Provide either `--comseg-temp-dir` or `--patches-dirs`

**Usage**:

```console
$ sopa resolve comseg [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--gene-column TEXT`: Column of the transcripts dataframe containing the genes names  [required]
* `--comseg-temp-dir TEXT`: Path to the directory containing all the comseg patches (see `sopa patchify`). By default, uses the `.sopa_cache/comseg_boundaries` directory
* `--min-area FLOAT`: Cells with an area less than this value (in microns^2) will be filtered  [default: 0]
* `--patches-dirs TEXT`: List of patches directories inside `comseg_temp_dir`
* `--help`: Show this message and exit.

#### `sopa resolve generic`

Resolve patches conflicts after generic segmentation

**Usage**:

```console
$ sopa resolve generic [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--method-name TEXT`: Name of the method used during segmentation. This is also the key correspnding to the boundaries in `sdata.shapes`  [required]
* `--patch-dir TEXT`: Directory containing the generic segmentation on patches (or multiple directories if using multi-step segmentation). By default, uses the `.sopa_cache/<method_name>` directory
* `--help`: Show this message and exit.

### `sopa segmentation`

Perform cell segmentation on patches. NB: for `baysor`, use directly the `baysor` command line.

**Usage**:

```console
$ sopa segmentation [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `cellpose`: Perform cellpose segmentation.
* `comseg`: Perform ComSeg segmentation.
* `generic-staining`: Perform generic staining-based segmentation.

#### `sopa segmentation cellpose`

Perform cellpose segmentation. This can be done on all patches directly, or on one individual patch.

!!! note "Usage"

    - [On one patch] Use this mode to run patches in parallel. Provide `--patch-index` to run one patch, and execute all patches in a parallel manner (you need to define your own parallelization, else, use the Snakemake pipeline).

    - [On all patches at once] For small images, you can run the segmentation method sequentially (`--patch-index` is not needed)

**Usage**:

```console
$ sopa segmentation cellpose [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--diameter FLOAT`: Cellpose diameter parameter  [required]
* `--channels TEXT`: Names of the channels used for Cellpose. If one channel, then provide just a nucleus channel. If two channels, this is the nucleus and then the cytoplasm channel  [required]
* `--flow-threshold FLOAT`: Cellpose `flow_threshold` parameter  [default: 2]
* `--cellprob-threshold FLOAT`: Cellpose `cellprob_threshold` parameter  [default: -6]
* `--model-type TEXT`: Name of the cellpose model  [default: cyto3]
* `--pretrained-model TEXT`: Path to the pretrained model to be loaded
* `--min-area INTEGER`: Minimum area (in pixels^2) for a cell to be considered as valid  [default: 0]
* `--clip-limit FLOAT`: Parameter for skimage.exposure.equalize_adapthist (applied before running cellpose)  [default: 0.2]
* `--gaussian-sigma FLOAT`: Parameter for scipy gaussian_filter (applied before running cellpose)  [default: 1]
* `--patch-index INTEGER`: Index of the patch on which cellpose should be run. NB: the number of patches is `len(sdata['sopa_patches'])`
* `--patch-dir TEXT`: Path to the temporary cellpose directory inside which we will store each individual patch segmentation. By default, saves into the `.sopa_cache/cellpose_boundaries` directory
* `--method-kwargs TEXT`: Kwargs for the cellpose method builder. This should be a dictionnary, in inline string format.  [default: {}]
* `--help`: Show this message and exit.

#### `sopa segmentation comseg`

Perform ComSeg segmentation. This can be done on all patches directly, or on one individual patch.

**Usage**:

```console
$ sopa segmentation comseg [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--patch-index INTEGER`: Index of the patch on which the segmentation method should be run.`
* `--patch-dir TEXT`: Path to the temporary the segmentation method directory inside which we will store each individual patch segmentation. By default, saves into the `.sopa_cache/comseg` directory
* `--help`: Show this message and exit.

#### `sopa segmentation generic-staining`

Perform generic staining-based segmentation. This can be done on all patches directly, or on one individual patch.

!!! note "Usage"
    First, define a new segmentation method, and write it under `sopa.segmentation.methods`. It should correspond to a function that is a "callable builder", i.e. kwargs will be provided to this function, and it will return a callable that will be applied on patches.

    As for Cellpose, two modes ara available:

    - [On one patch] Use this mode to run patches in parallel. Provide `--patch-index` to run one patch, and execute all patches in a parallel manner (you need to define your own parallelization, else, use the Snakemake pipeline).

    - [On all patches at once] For small images, you can run the segmentation method sequentially (`--patch-index` is not needed)

**Usage**:

```console
$ sopa segmentation generic-staining [OPTIONS] SDATA_PATH
```

**Arguments**:

* `SDATA_PATH`: Path to the SpatialData `.zarr` directory  [required]

**Options**:

* `--method-name TEXT`: Name of the segmentation method builder to use. The corresponding function (`sopa.segmentation.methods.<method_name>`) will be used, and the kwargs below will be used to instantiate the method.  [required]
* `--method-kwargs TEXT`: Kwargs for the method. This should be a dictionnary, in inline string format.  [default: {}]
* `--channels TEXT`: Names of the channels used for segmentation.  [required]
* `--min-area INTEGER`: Minimum area (in pixels^2) for a cell to be considered as valid  [default: 0]
* `--clip-limit FLOAT`: Parameter for skimage.exposure.equalize_adapthist (applied before running the segmentation method)  [default: 0.2]
* `--gaussian-sigma FLOAT`: Parameter for scipy gaussian_filter (applied before running the segmentation method)  [default: 1]
* `--patch-index INTEGER`: Index of the patch on which the segmentation method should be run. NB: the number of patches is `len(sdata['sopa_patches'])`
* `--patch-dir TEXT`: Path to the temporary the segmentation method directory inside which we will store each individual patch segmentation. By default, saves into the `.sopa_cache/<method_name>` directory
* `--help`: Show this message and exit.
