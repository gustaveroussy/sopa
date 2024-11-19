# Snakemake pipeline

Sopa comes with an existing [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow to get started quickly. This will not involve any coding but requires some setup specific to `snakemake`.

## Setup

### Installation

Follow our [installation instructions](../../getting_started/#snakemake-setup) until the end of the "Snakemake setup" section. At the end, you should have one `sopa` environment, one one environment with `snakemake` (it can be the same environment, if desired), and you should also have cloned the `sopa` repository.

### Choose a config file

Our pipeline config is a YAML file that describes all the steps desired for the pipeline. It is flexible; for instance, if you remove the `baysor` arguments from the config, then it will not run baysor. Similarly, if you remove the `"annotation"` section, it will not run annotation.

You can choose a config among the existing ones [here](https://github.com/gustaveroussy/sopa_workflow/tree/main/config) or [create your own](./#create-your-own-config).

Keep in mind the relative path of your config since you'll need it later, e.g. `config/merscope/base.yaml`.

## Run the pipeline

1. First, locate the path to one sample's raw experiment file(s). This is usually a directory containing one or many image(s) and, eventually, a transcript file. If you don't know what data you need, see our [FAQ](../../faq/#what-are-the-inputs-or-sopa).

2. Then, activate an environment that has the snakemake command:
```sh
conda activate snakemake    # or any environment that has `snakemake`
```

1. Go in the `workflow` directory
```sh
cd workflow   # move to the workflow directory inside the sopa repository
```

1. You can either execute the pipeline locally or on a high-performance-cluster (choose the right option below)

=== "Local execution (e.g., personal laptop)"

    You can execute the pipeline locally as below (in this example, we use only one core). Make sure to replace `data_path` with the path to your data directory, and `configfile` with the relative path to your config.

    ```sh
    snakemake \
        --config data_path=/path/to/directory \
        --configfile=config/merscope/base.yaml \
        --workflow-profile profile/local \
        --cores 1
    ```

    !!! note "Faster pipeline"
        Even though Sopa can be run locally, we recommend to use it on high-performance-clusters to benefit from all the pipeline capabilities (see the second tab just above).

=== "High-performance-cluster (e.g., Slurm cluster)"

    To benefit from high-performance-cluster, you'll need a [Snakemake cluster profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).

    If you have a Slurm cluster, you can use our default Slurm profile. For that, you'll need `snakemake>=8.0.0`, and you'll also need to install the [Slurm plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) via `pip install snakemake-executor-plugin-slurm`.

    Then, you can use the profile as below:

    ```sh
    snakemake \
        --config data_path=/path/to/directory \
        --configfile=config/merscope/base.yaml
        --workflow-profile profile/slurm  # or any profile you want
    ```

    !!! warning "Specify the slurm partition names"
        You may need to update the `slurm_partition` parameters inside the `workflow/profile/slurm/config.yaml` file according to the partition names of your cluster (else, it will always use the same partition). You can also change `mem_mb`, depending on the RAM capabilities of your cluster.

=== "Other"

    If you don't have a Slurm cluster, we recommend reading more about the [Snakemake profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles), and, especially, the different [executor plugins](https://snakemake.github.io/snakemake-plugin-catalog/index.html).

    Once you installed an executor plugin, you can use it via:

    ```sh
    snakemake \
        --config data_path=/path/to/directory \
        --configfile=config/merscope/base.yaml
        --executor my_executor  # your new executor
    ```

    Or, if you created a new config file under `workflow/profile/my_profile/config.yaml`, you can use it via:

    ```sh
    snakemake \
        --config data_path=/path/to/directory \
        --configfile=config/merscope/base.yaml
        --workflow-profile profile/my_profile  # your new profile
    ```

    !!! warning "RAM per rule"
        Some Snakemake rules may use more RAM, for instance `explorer`, or `to_spatialdata`. Consider adjusting the RAM and the walltime depending on the Snakemake rules.


For more customization, see the [snakemake CLI documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

## Toy example

In the example below, we run the pipeline on a generated toy dataset. Running it locally can help test a new pipeline or config.

Make sure you have installed everything as detailed in this tutorial, and then run the following command lines:

=== "Cellpose usage"
    Make sure you have installed sopa with the `cellpose` extra (for instance, this can be done via the following command: `pip install 'sopa[cellpose]'`).
    ```sh
    conda activate snakemake    # or any environment that has `snakemake`
    cd workflow   # move to the workflow directory inside the sopa repository

    # you can replace tuto.zarr by another path where the data will be saved
    snakemake \
        --config sdata_path=tuto.zarr \
        --configfile=config/toy/cellpose.yaml \
        --workflow-profile profile/local \
        --cores 1
    ```

=== "Baysor usage"
    Make sure you have installed sopa with the `baysor` extra, and that you have installed the `baysor` command (refer to our getting started).
    ```sh
    conda activate snakemake    # or any environment that has `snakemake`
    cd workflow   # move to the workflow directory inside the sopa repository

    # you can replace tuto.zarr by another path where the data will be saved
    snakemake \
        --config sdata_path=tuto.zarr \
        --configfile=config/toy/baysor.yaml \
        --workflow-profile profile/local \
        --cores 1
    ```

!!! notes
    On the above example, it executes snakemake sequentially (one core), which is enough for debugging purposes. You can remove the argument, or set a higher number of cores.

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

## Advanced usage

### Create your own config

If the existing `config` files are not suited for your project, you can update an existing one or create a whole new one. For this, use [this commented config](https://github.com/gustaveroussy/sopa/blob/master/workflow/config/example_commented.yaml) to understand the purpose of each argument. Note that some sections are optional: in this case, remove the section or the argument, and Sopa will not run it.

When running snakemake, you will then need to provide the relative or absolute path to your `.yaml` config, for instance `--configfile=/path/to/your/config.yaml`.

### Passing kwargs to the config

Internally, the Snakemake pipeline is calling [Sopa's CLI](../cli_usage). **Not all argument are used** in the default config files. You can pass additionnal kwargs to the config so that they are given to the CLI.

For instance, if you want to pass kwargs to the [MERSCOPE reader](../../api/readers/#sopa.io.merscope), you can update the config as below:
```
read:
  technology: merscope
  kwargs:
    z_layers: 2
```

### Create your own profile

As mentioned above, you can use a [Snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) to execute the pipeline. In Snakemake>=8.0.0, there are multiple [executor plugins](https://snakemake.github.io/snakemake-plugin-catalog/index.html) that can help you.

Save your new profile under `workflow/profile/my_profile/config.yaml`.

Then, to use the new profile, pass `--workflow-profile profile/my_profile` to your `snakemake` command.
