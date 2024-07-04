# Snakemake pipeline

Sopa comes with an existing [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to get started quickly. This will not involve any coding but requires some setup specific to `snakemake`.

## Install `sopa`

Follow the ["Snakemake setup" instructions](../../getting_started/#snakemake-setup) of our installation page.

!!! warning "Baysor usage"
    Even though `pip install -e '.[baysor]'` will install some dependencies related to baysor, you still have to install the `baysor` command line (see the [official repository](https://github.com/kharchenkolab/Baysor)) if you want to use it.

    If your path to the baysor executable is not the default one (i.e. `~/.julia/bin/baysor`), you can do one of the following:

    - set a `BAYSOR_EXECUTABLE_PATH` environment variable with the path to your Baysor executable.
    - update the config described below to provide the right path to your Baysor executable.

## Choose a config

Our pipeline config is a YAML file that describes all the steps desired for the pipeline. It is flexible; for instance, if you remove the `baysor` arguments from the config, then it will not run baysor. Similarly, if you remove the `"annotation"` section, it will not run annotation.

You can choose a config among the existing ones [here](https://github.com/gustaveroussy/sopa/tree/master/workflow/config) or [create your own](./#create-your-own-config).

Note the relative path of your config since you'll need it later, e.g. `config/merscope/base.yaml`.

## Run the pipeline

1. First, locate the path to one sample's raw experiment file(s). This is usually a directory containing one or many image(s) and, eventually, a transcript file. If you don't know what data you need, see our [FAQ](../../faq/#what-kind-of-inputs-do-i-need-to-run-sopa).

2. Then, activate your environment that has the snakemake command, and go to the `workflow` directory inside the `sopa` directory that you cloned earlier:
```sh
conda activate sopa    # or an environment that has `snakemake`
cd workflow            # run this at the root of the 'sopa' directory
```

1. You can either execute the pipeline locally or on a high-performance-cluster (choose the right option below)

=== "Local execution (e.g., personal laptop)"

    You can execute the pipeline locally as below (in this example, we use only one core):

    ```sh
    # replace the configfile with yours
    # replace data_path with the path to your data directory

    snakemake --config data_path=/path/to/directory --configfile=config/merscope/base.yaml --cores 1 --use-conda
    ```

    !!! note "Faster pipeline"
        Even though Sopa can be run locally, we recommend to use it on high-performance-clusters to benefit from all the pipeline capabilities (see the second tab just above).

=== "High-performance-cluster (e.g., Slurm cluster)"

    To benefit from high-performance-cluster, you'll need a [Snakemake cluster profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles). By default, we use [this `slurm` profile](https://github.com/gustaveroussy/sopa/blob/master/workflow/slurm/config.yaml), but you can also update it or create your own profile under `sopa/workflow/<hpc-name>/config.yaml`.

    ```sh
    # replace the configfile with yours
    # replace data_path with the path to your data directory
    # replace profile with <hpc-name> as above, or keep 'slurm'

    snakemake --profile slurm --config data_path=/path/to/directory --configfile=config/merscope/base.yaml
    ```

    !!! note
        You may need to update the `partition` parameters inside the `sopa/workflow/Snakefile` file, according to the partition names of your cluster. You can also change `mem_mb`, depending on the RAM capabilities of your cluster.

For more customization, see the [snakemake CLI documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

## Toy example

In the example below, we run the pipeline on a generated toy dataset. Running it locally can help test a new pipeline or config.

Make sure you have installed everything as detailed in this tutorial, and then run the following command lines:

=== "Cellpose usage"
    Make sure you have installed sopa with the Cellpose extra
    ```sh
    conda activate sopa    # or an environment that has `snakemake`
    cd workflow            # run this at the root of the 'sopa' directory

    # you can replace tuto.zarr by another path where the data will be saved
    snakemake --config sdata_path=tuto.zarr --configfile=config/toy/uniform_cellpose.yaml --cores 1 --use-conda
    ```

=== "Baysor usage"
    Make sure you have installed sopa with the Baysor extra, and that you have installed the `baysor` command
    ```sh
    conda activate sopa    # or an environment that has `snakemake`
    cd workflow            # run this at the root of the 'sopa' directory

    # replace tuto.zarr by the path where you want the data to be saved
    snakemake --config sdata_path=tuto.zarr --configfile=config/toy/uniform_baysor.yaml --cores 1 --use-conda
    ```

!!! notes
    On the above example, it executes snakemake sequentially (one core), which is enough for debugging purposes.

You can then check `tuto.explorer` for output files. Notably, if you have installed the [Xenium Explorer](https://www.10xgenomics.com/support/software/xenium-explorer), double-click on `experiment.xenium` to visualize the results.

## Pipeline outputs

The pipeline outputs consists in two directories located next to your raw data directory. They have the same name as your raw directory, but with extension `.zarr` and `.explorer` respectively (see below for more details).

!!! info
    You can also change the path to the `.zarr` output, by providing `sdata_path=/your/path.zarr` just after `--config` on the snakemake execution line. This will also move the `.explorer` directory, that will be saved at `/your/path.explorer`

#### `SpatialData` directory

If you are familiar with the [`spatialdata` library](https://github.com/scverse/spatialdata), you can use the `.zarr` directory, corresponding to a `SpatialData` object:
```python
import spatialdata

sdata = spatialdata.read_zarr("/path/to/data.zarr")
```

#### Explorer directory

The `.explorer` directory contains the following files:

- `report.html` a short quality control of you data, as an HTML report

- `adata.h5ad` the AnnData object with spatial locations of the cells (see `adata.obsm['spatial']`), and also cell-by-gene table and/or the cell-by-channel table.

- `experiment.xenium` the Xenium Explorer file: double-click on it to open it on the Xenium Explorer ([download the software here](https://www.10xgenomics.com/support/software/xenium-explorer/downloads))

- The other files are data files related and required by the Xenium Explorer

## Create your own config

If the existing `config` files are not suited for your project, you can update an existing one or create a whole new one. For this, use [this commented config](https://github.com/gustaveroussy/sopa/blob/master/workflow/config/example_commented.yaml) to understand the purpose of each argument. Note that some sections are optional: in this case, remove the section or the argument, and Sopa will not run it.

When running snakemake, you will then need to provide the relative or absolute path to your `.yaml` config, for instance `--configfile=/path/to/your/config.yaml`.
