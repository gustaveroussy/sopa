## What kind of inputs do I need to run Sopa?

You need the raw inputs of your machine, that is:

- One or multiple image(s), usually corresponding to one or multiple `.tiff` file(s)

- Optionally, a file of transcript location, usually a `.csv` or `.parquet` file

In this documentation, `data_path` denotes the path to your raw data. Select the correct tab below to understand what is the right path to your raw data:

=== "Xenium"
    `data_path` is the directory containing the following files: `morphology.ome.tif` and `transcripts.parquet`
=== "MERSCOPE"
    `data_path` is the "region" directory containing a `detected_transcripts.csv` file and an `image` directory. For instance, the directory can be called `region_0`.
=== "CosMX"
    (The CosMX data requires stitching the FOVs. It will be added soon, see [this issue](https://github.com/gustaveroussy/sopa/issues/5))
=== "MACSima"
    `data_path` is the directory containing multiple `.ome.tif` files (one file per channel)
=== "PhenoCycler"
    `data_path` corresponds to the path to one `.qptiff` file, or one `.tif` file (if exported from QuPath)
=== "Hyperion"
    `data_path` is the directory containing multiple `.ome.tiff` files (one file per channel)

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

## How to provide dictionnaries to CLI arguments?

Some CLI arguments are optionnal dictionnaries. For instance, [`sopa read`](../cli/#sopa-read) has a `--kwargs` option. In that case, a dictionnary can be provided as an inline string, for instance:

`--kwargs "{'backend': 'rioxarray'}"`

## Can I use Nextflow instead of Snakemake?

Nextflow is not supported yet, but we are working on it. You can also help re-write our Snakemake pipeline for Nextflow (see issue [#7](https://github.com/gustaveroussy/sopa/issues/7)).

## I have another issue; how do I fix it?

Don't hesitate to open an issue on [Sopa's Github repository](https://github.com/gustaveroussy/sopa/issues), and detail your issue with as much precision as possible for the maintainers to be able to reproduce it.
