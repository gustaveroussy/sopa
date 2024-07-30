## What kind of inputs do I need to run Sopa?

You need the raw inputs of your machine, that is:

- One or multiple image(s), usually corresponding to one or multiple `.tiff` file(s)

- Optionally, a file of transcript location, usually a `.csv` or `.parquet` file

In this documentation, `data_path` denotes the path to your raw data. Select the correct tab below to understand what is the right path to your raw data:

=== "Xenium"
    `data_path` is the path to the directory containing the following files: `morphology.ome.tif`, `experiment.xenium` and `transcripts.parquet`. In brief, you should have this file structure:

    ```txt
    .
    ├─ morphology_focus.ome.tif   # or a directory (for recent versions of the Xenium)
    ├─ experiment.xenium
    └─ transcripts.parquet
    ```

=== "MERSCOPE"
    `data_path` is the path to the "region" directory containing a `detected_transcripts.csv` file and an `images` directory. For instance, the directory may be called `region_0`. In brief, you should have this file structure:

    ```txt
    .
    ├─ detected_transcripts.csv
    └─ images
       ├─ mosaic_{stain}_z{z_layer}.tif
       └─ micron_to_mosaic_pixel_transform.csv
    ```
=== "CosMX"
    `data_path` is the path to the directory containing:

    - a transcript file `*_tx_file` (with columns `target`, `x_global_px`, `y_global_px`)
    - a FOV locations file `*_fov_positions_file` (with columns `FOV`, `X_mm`, `Y_mm`)
    - a `Morphology_ChannelID_Dictionary.txt` file containing channel names
    - a `Morphology2D` directory containing the images, end in `_F*.TIF`.

    These files must be exported as flat files in AtomX. That is: within a study, click on "Export" and then select files from the "Flat CSV Files" section (transcripts flat and FOV position flat). You should have this file structure:

    ```txt
    .
    ├─ <DATASET_ID>_tx_file.csv (or csv.gz)
    ├─ <DATASET_ID>_fov_positions_file.csv (or csv.gz)
    ├─ Morphology_ChannelID_Dictionary.txt
    └─ Morphology2D
       ├─ XXX_F001.TIF
       ├─ XXX_F002.TIF
       └─ ...
    ```
=== "MACSima"
    `data_path` is the path to the directory containing multiple `.ome.tif` files (one file per channel). In brief, you should have this file structure:

    ```txt
    .
    ├─ AAA.ome.tif
    ├─ BBB.ome.tif
    └─ CCC.ome.tif
    ```
=== "PhenoCycler"
    `data_path` is the path to one `.qptiff` file, or one `.tif` file (if exported from QuPath).
=== "Hyperion"
    `data_path` is path to the directory containing multiple `.ome.tiff` files (one file per channel). In brief, you should have this file structure:
    ```txt
    .
    ├─ AAA.ome.tiff
    ├─ BBB.ome.tiff
    └─ CCC.ome.tiff
    ```
=== "Others (CZI, ...)"
    Other file formats (ND2, CZI, LIF, or DV) are supported via the `aicsimageio` reader. In that case, you'll need to add new dependencies: `pip install aicsimageio` (and, for CZI data, also `pip install aicspylibczi`).

    This reader is called `aicsimageio`, i.e. you can use it via `sopa.io.aicsimageio(data_path)`, where `data_path` is the path to your data file containing your image(s). For the Snakemake pipeline, provide `aicsimageio` as a `technology` in the config file.

## I have small artifact cells, how do remove them?

You may have small cells that were segmented but that should be removed. For that, `Sopa` offers three filtering approaches: using their area, their transcript count, or their fluorescence intensity. Refer to the following config parameters from this [example config](https://github.com/gustaveroussy/sopa/blob/master/workflow/config/example_commented.yaml): `min_area`, `min_transcripts`, and `min_intensity_ratio`.

If using the CLI, `--min-area` can be provided to `sopa segmentation cellpose` or `sopa resolve baysor`, and `--min-transcripts`/`--min-intensity-ratio` can be provided to `sopa aggregate`.

## Cellpose is not segmenting enough cells; what should I do?

- The main Cellpose parameter to check is `diameter`, i.e. a typical cell diameter **in pixels**. Note that this is highly specific to the technology you're using since the micron-to-pixel ratio can differ. We advise you to start with the default parameter for your technology of interest (see the `diameter` parameter inside our config files [here](https://github.com/gustaveroussy/sopa/tree/master/workflow/config)).
- Maybe `min_area` is too high, and all the cells are filtered because they are smaller than this area. Remind that, when using Cellpose, the areas correspond to pixels^2.
- This can be due to a low image quality. If the image is too pixelated, consider increasing `gaussian_sigma` (e.g., `2`) under the cellpose parameters of our config. If the image has a low contrast, consider increasing `clip_limit` (e.g., `0.3`). These parameters are detailed in [this example config](https://github.com/gustaveroussy/sopa/blob/master/workflow/config/example_commented.yaml).
- Consider updating the official Cellpose parameters. In particular, try `cellprob_threshold=-6` and `flow_threshold=2`.

## How to use a custom Cellpose model?

You can use any existing [Cellpose model](https://cellpose.readthedocs.io/en/latest/models.html) with the `model_type` argument (via the API, CLI, or Snakemake pipeline). For the Snakemake pipeline, see [here](https://github.com/gustaveroussy/sopa/blob/master/workflow/config/example_commented.yaml) how to set this argument.
If you have a custom pretrained model, use the `pretrained_model` argument instead of `model_type`, and give the path to your cellpose model.

## How to provide other arguments to Cellpose?

When using the Snakemake pipeline, you can use `method_kwargs` to provide extra arguments to Cellpose. For instance, we use `resample=False` in the example below, which may significantly speed up the segmentation while not decreasing significantly the segmentation quality:

```yaml
segmentation:
  cellpose:
    diameter: 60
    channels: ["DAPI"]
    flow_threshold: 2
    cellprob_threshold: -6
    min_area: 2000
    method_kwargs:
      resample: False
```

## How to use a prior cell segmentation?

If you have MERSCOPE or Xenium data, you probably already have a cell segmentation. This can be used as a prior for Baysor, instead of running Cellpose with Sopa. For that, you have an existing config file for the Snakemake pipeline for both [MERSCOPE](https://github.com/gustaveroussy/sopa/blob/master/workflow/config/merscope/baysor_vizgen.yaml) and [Xenium](https://github.com/gustaveroussy/sopa/blob/master/workflow/config/xenium/baysor_multimodal.yaml) data. If using the API/CLI, consider using the `cell_key` and the `unassigned_value` arguments when creating the patches for the transcripts. For MERSCOPE data, `cell_key="cell_id"` and `unassigned_value=-1`. For Xenium data, `cell_key="cell_id"` and `unassigned_value="UNASSIGNED"`.

## How to provide dictionnaries to CLI arguments?

Some CLI arguments are optionnal dictionnaries. For instance, [`sopa read`](../cli/#sopa-read) has a `--kwargs` option. In that case, a dictionnary can be provided as an inline string, for instance:

`--kwargs "{'backend': 'rioxarray'}"`

## How to fix an "out-of-memory" issue on MERSCOPE data?

If using MERSCOPE data, images can be huge. To improve RAM efficiency, you can install `rioxarray` (`pip install rioxarray`). Then, the `rioxarray` will be used by default by the reader (no change needed, it will be detected automatically).

## Can I use Nextflow instead of Snakemake?

Nextflow is not supported yet, but we are working on it. You can also help re-write our Snakemake pipeline for Nextflow (see issue [#7](https://github.com/gustaveroussy/sopa/issues/7)).

## I have another issue; how do I fix it?

Don't hesitate to open an issue on [Sopa's Github repository](https://github.com/gustaveroussy/sopa/issues), and detail your issue with as much precision as possible for the maintainers to be able to reproduce it.
