# Sopa command-line-interface (CLI)

## Usage

When installing `sopa` are written in our [getting-started guidelines](../getting_started), a new command named `sopa` becomes available.

!!! note "CLI helper"
    Run `sopa --help` to get details about all the command line purpose. You can also use this helper on any subcommand, for instance `sopa read --help`.

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
│ explorer      Convertion to the Xenium Explorer's inputs       │
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

When running the `sopa` CLI, some arguments are required while some are optional. For instance, for the `sopa read` command, `sdata_path` is an argument and a path has to be given directly, while `technology` is an option, and in this case the `--technology` prefix has to be used. For instance, if you read MERSCOPE data, it will be:

```
sopa read /path/to/merscope/directory --technology merscope
```

Note that `/path/to/merscope/directory` refers to `sdata_path`, which is an argument. You **should not** add the suffix `--sdata_path`, as it is an argument.

!!! note "Required options"
    All the arguments are required, as shown by the `[required]` hint on the CLI helper. Note that some options may also be required too (in this case, the term `[required]` will appear on the CLI helper). But they still need to be called as a normal option.

## Commands

#### `sopa read`
Read any technology data, and write a standardized SpatialData object. Either `--technology` or `--config-path` has to be provided.

> Usage: `sopa read [OPTIONS] DATA_PATH`

```txt
[Args]           
data_path: Path to one data sample (most of the time, this corresponds to a directory)
                                                                  
[Options]
technology: Name of the technology used to collected the data (e.g., 'xenium', 'merfish', ...)
sdata_path: Optional path to write the SpatialData object. If not provided, will write to the '{data_path}.zarr' directory
config_path: Path to the snakemake config. This can be useful in order not to provide the 'technology' and the 'kwargs' arguments
kwargs: Dictionary provided to the reader function.

╭─ Arguments ──────────────────────────────────────────────────────╮
│ *    data_path      TEXT  [default: None] [required]             │
╰──────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────╮
│ --technology         TEXT  [default: None]                       │
│ --sdata-path         TEXT  [default: None]                       │
│ --config-path        TEXT  [default: None]                       │
│ --kwargs             TEXT  [default: {}]                         │
│ --help                     Show this message and exit.           │
╰──────────────────────────────────────────────────────────────────╯
```

#### `sopa patchify cellpose`

Prepare patches for Cellpose segmentation

> Usage: `sopa patchify cellpose [OPTIONS] SDATA_PATH`

```
[Args]                                                             
sdata_path: Path to the SpatialData zarr directory

[Options]
patch_width_pixel: Width (and height) of each patch in pixels      
patch_overlap_pixel: Number of overlapping pixels between the patches. We advise to choose approximately twice the diameter of a cell

╭─ Arguments ──────────────────────────────────────────────────────╮
│ *    sdata_path      TEXT  [default: None] [required]            │
╰──────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────╮
│ --patch-width-pixel          FLOAT  [default: None]              │
│ --patch-overlap-pixel        FLOAT  [default: None]              │
│ --help                              Show this message and exit.  │
╰──────────────────────────────────────────────────────────────────╯
```

#### `sopa patchify baysor`

Prepare the patches for Baysor segmentation

> Usage: `sopa patchify baysor [OPTIONS] SDATA_PATH`

```txt
[Args]                                                                                    
sdata_path: Path to the SpatialData zarr directory                                        
                                                                                           
[Options]                                                                                 
patch_width_microns: Width (and height) of each patch in microns                          
patch_overlap_microns: Number of overlapping microns between the patches. We advise to choose approximately twice the diameter of a cell
baysor_temp_dir: Temporary directory where baysor inputs and outputs will be saved
config: Path to the snakemake config containing the baysor arguments
cell_key: Optional column of the transcripts dataframe that indicates in which cell-id each transcript is
unassigned_value: If 'cell_key' is provided, this is the value given to transcripts that are not inside any cell (if it's already 0, don't provide this argument)
use_prior: Whether to use cellpose segmentation as a prior for 

╭─ Arguments ─────────────────────────────────────────────────────────────────────────────╮
│ *    sdata_path      TEXT  [default: None] [required]                                   │
╰─────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ───────────────────────────────────────────────────────────────────────────────╮
│ *  --patch-width-microns                        FLOAT    [default: None] [required]     │
│ *  --patch-overlap-microns                      FLOAT    [default: None] [required]     │
│ *  --baysor-temp-dir                            TEXT     [default: None] [required]     │
│    --config                                     TEXT     [default: {}]                  │
│    --cell-key                                   TEXT     [default: None]                │
│    --unassigned-value                           INTEGER  [default: None]                │
│    --use-prior                --no-use-prior             [default: no-use-prior]        │
│    --help                                                Show this message and exit.    │
╰─────────────────────────────────────────────────────────────────────────────────────────╯
```

#### `sopa segmentation cellpose`

Perform cellpose segmentation. This can be done on all patches directly, or on one individual patch (provide `--patch-dir` and `--patch-index`).     

> Usage: `sopa segmentation cellpose [OPTIONS] SDATA_PATH`
                                 
```txt
[Args]
sdata_path: Path to the SpatialData zarr directory

[Options]
diameter: Cellpose diameter parameter
channels: Names of the channels used for Cellpose. If one channel, then provide just a
nucleus channel. If two channels, this is the nucleus and then the cytoplasm channel.
flow_threshold: Cellpose flow_threshold parameter
cellprob_threshold: Cellpose cellprob_threshold parameter
model_type: Name of the cellpose model
patch_width: Ignore this if you already run 'sopa patchify'. Patch width in pixels.
patch_overlap: Ignore this if you already run 'sopa patchify'. Patches overlaps in pixels.
expand_radius: Ignore this if you already run 'sopa patchify'. Cell boundaries radius expansion in pixels.
patch_index: Index of the patch on which cellpose should be run. NB: the number of patches is `len(sdata['sopa_patches'])`.
patch_dir: Path to the temporary cellpose directory inside which we will store each individual patch segmentation
       
╭─ Arguments ─────────────────────────────────────────────────────────╮
│ *    sdata_path      TEXT  [default: None] [required]               │
╰─────────────────────────────────────────────────────────────────────╯
╭─ Options ───────────────────────────────────────────────────────────╮
│ *  --diameter                  FLOAT    [default: None] [required]  │
│ *  --channels                  TEXT     [default: None] [required]  │
│ *  --flow-threshold            FLOAT    [default: None] [required]  │
│ *  --cellprob-threshold        FLOAT    [default: None] [required]  │
│    --model-type                TEXT     [default: cyto2]            │
│    --patch-width               INTEGER  [default: None]             │
│    --patch-overlap             INTEGER  [default: None]             │
│    --expand-radius             INTEGER  [default: 0]                │
│    --patch-index               INTEGER  [default: None]             │
│    --patch-dir                 TEXT     [default: None]             │
│    --help                               Show this message and exit. │
╰─────────────────────────────────────────────────────────────────────╯
```

#### `sopa resolve cellpose`

Resolve patches conflicts after cellpose segmentation

> Usage: `sopa resolve cellpose [OPTIONS] SDATA_PATH`

```txt
[Args]
sdata_path: Path to the SpatialData zarr directory

[Options]
patch_dir: Directory containing the cellpose segmentation on patches
expand_radius: Number of pixels for radius expansion of each cell boundary

╭─ Arguments ──────────────────────────────────────────────────╮
│ *    sdata_path      TEXT  [default: None] [required]        │
╰──────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────╮
│ *  --patch-dir            TEXT   [default: None] [required]  │
│    --expand-radius        FLOAT  [default: 0]                │
│    --help                        Show this message and exit. │
╰──────────────────────────────────────────────────────────────╯
```

#### `sopa resolve baysor`

Resolve patches conflicts after baysor segmentation. Provide either 'baysor_temp_dir' 'patches_dirs'

> Usage: sopa resolve baysor [OPTIONS] SDATA_PATH

```txt
[Args]
sdata_path: Path to the SpatialData zarr directory

[Options]
gene_column: Column of the transcripts dataframe containing the genes names
baysor_temp_dir: Path to the directory containing all the baysor patches (see 'sopa patchify')
min_area: Cells with an area less than this value (in microns^2) will be filtered
expand_radius: Number of microns for radius expansion of each cell boundary
patches_dirs: List of patches directories inside 'baysor_temp_dir'

╭─ Arguments ────────────────────────────────────────────────────╮
│ *    sdata_path      TEXT  [default: None] [required]          │
╰────────────────────────────────────────────────────────────────╯
╭─ Options ──────────────────────────────────────────────────────╮
│ *  --gene-column            TEXT   [default: None] [required]  │
│    --baysor-temp-dir        TEXT   [default: None]             │
│    --min-area               FLOAT  [default: 0]                │
│    --expand-radius          FLOAT  [default: 0]                │
│    --patches-dirs           TEXT   [default: None]             │
│    --help                          Show this message and exit. │
╰────────────────────────────────────────────────────────────────╯
```

#### `sopa aggregate`

Create an `anndata` table containing the transcript count and/or the channel intensities per cell

> Usage: `sopa aggregate [OPTIONS] SDATA_PATH`

```txt
[Args]
sdata_path: Path to the SpatialData zarr directory

[Options]
gene_column: Column of the transcript dataframe representing the gene names. If not provided, it will not compute transcript count
average_intensities: Whether to average the channel intensities inside each cell
min_transcripts: Cells with less transcript than this integer will be filtered
min_intensity_ratio: Cells whose mean channel intensity is less than min_intensity_ratio * quantile_90 will be filtered

╭─ Arguments ───────────────────────────────────────────────────────────────────────────────╮
│ *    sdata_path      TEXT  [default: None] [required]                                     │
╰───────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ─────────────────────────────────────────────────────────────────────────────────╮
│ --gene-column                                        TEXT     [default: None]             │
│ --average-intensities    --no-average-intensities             [default:                   │
│                                                               no-average-intensities]     │
│ --min-transcripts                                    INTEGER  [default: 0]                │
│ --min-intensity-ratio                                FLOAT    [default: 0]                │
│ --help                                                        Show this message and exit. │
╰───────────────────────────────────────────────────────────────────────────────────────────╯
```
#### `sopa annotate tangram`

Tangram segmentation (i.e., uses an annotated scRNAseq reference to transfer cell-types)

> Usage: `sopa annotate tangram [OPTIONS] SDATA_PATH`

[Args]
sdata_path: Path to the SpatialData zarr directory

[Options]
sc_reference_path: Path to the scRNAseq annotated reference
cell_type_key: Key of 'adata_ref.obs' containing the cell-types
reference_preprocessing: Preprocessing method applied to the reference. Either None (raw counts), or 'normalized' (sc.pp.normalize_total) or 'log1p' (sc.pp.normalize_total and sc.pp.log1p)
bag_size: Number of cells in each bag of the spatial table. Low values will decrease the memory usage
max_obs_reference: Maximum samples to be considered in the reference for tangram. Low values will decrease the memory usage

```txt
╭─ Arguments ──────────────────────────────────────────────────────────────╮
│ *    sdata_path      TEXT  [default: None] [required]                    │
╰──────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────╮
│ *  --sc-reference-path              TEXT     [default: None] [required]  │
│    --cell-type-key                  TEXT     [default: cell_type]        │
│    --reference-preprocessing        TEXT     [default: None]             │
│    --bag-size                       INTEGER  [default: 10000]            │
│    --max-obs-reference              INTEGER  [default: 10000]            │
│    --help                                    Show this message and exit. │
╰──────────────────────────────────────────────────────────────────────────╯
```

#### `sopa annotate fluorescence`

Simple annotation based on fluorescence, where each provided channel corresponds to one
cell type. For each cell, one z-score statistic is computed and the population with the
highest z-score is attributed.

> Usage: `sopa annotate fluorescence [OPTIONS] SDATA_PATH`

[Args]
sdata_path: Path to the SpatialData zarr directory

[Options]
marker_cell_dict: Dictionary chose keys are channel names, and values are the corresponding cell types
cell_type_key: Key added in 'adata.obs' corresponding to the cell type

```txt
╭─ Arguments ─────────────────────────────────────────────────╮
│ *    sdata_path      TEXT  [default: None] [required]       │
╰─────────────────────────────────────────────────────────────╯
╭─ Options ───────────────────────────────────────────────────╮
│ --marker-cell-dict        TEXT  [default: {}]               │
│ --cell-type-key           TEXT  [default: cell_type]        │
│ --help                          Show this message and exit. │
╰─────────────────────────────────────────────────────────────╯
```

#### `sopa report`

Create a HTML report of the pipeline run and some quality controls

> Usage: `sopa report [OPTIONS] SDATA_PATH PATH`

[Args]
sdata_path: Path to the SpatialData zarr directory
path: Path to the HTML report

```txt
╭─ Arguments ───────────────────────────────────────────╮
│ *    sdata_path      TEXT  [default: None] [required] │
│ *    path            TEXT  [default: None] [required] │
╰───────────────────────────────────────────────────────╯
╭─ Options ─────────────────────────────────────────────╮
│ --help          Show this message and exit.           │
╰───────────────────────────────────────────────────────╯
```

#### `sopa explorer`

Convert a spatialdata object to Xenium Explorer's inputs

> Usage: sopa explorer [OPTIONS] SDATA_PATH

[Args]
sdata_path: Path to the SpatialData zarr directory

[Options]
output_path: Path to a directory where Xenium Explorer's outputs will be saved. By default, writes to the same path as `sdata_path` but with the `.explorer` suffix
shapes_key: Key for the boundaries. By default, uses the baysor boundaires, else the cellpose boundaries.
gene_column: Column name of the points dataframe containing the gene names.
lazy: If `True`, will not load the full images in memory (except if the image memory is below `ram_threshold_gb`).
ram_threshold_gb: Threshold (in gygabytes) from which image can be loaded in memory. If `None`, the image is never loaded in memory.
save_image_mode: `1` is normal mode. `0` doesn't save the image. `2` saves **only** the image.

```txt
╭─ Arguments ─────────────────────────────────────────────────────────────╮
│ *    sdata_path      TEXT  [default: None] [required]                   │
╰─────────────────────────────────────────────────────────────────────────╯
╭─ Options ───────────────────────────────────────────────────────────────╮
│ --output-path                      TEXT     [default: None]             │
│ --gene-column                      TEXT     [default: None]             │
│ --shapes-key                       TEXT     [default: None]             │
│ --lazy                --no-lazy             [default: lazy]             │
│ --ram-threshold-gb                 INTEGER  [default: 4]                │
│ --save-image-mode                  INTEGER  [default: 1]                │
│ --help                                      Show this message and exit. │
╰─────────────────────────────────────────────────────────────────────────╯
```