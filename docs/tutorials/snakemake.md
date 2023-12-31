# Snakemake pipeline

Sopa comes with an existing [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to get started quickly. This will not involve any coding, but requires some setup specific to `snakemake`.

## Setup

- Clone the `sopa` repository, i.e., `git clone https://github.com/gustaveroussy/sopa.git`
- Make sure you have a `conda` environment called `sopa` (see [getting started](../getting_started)), on which you installed `sopa` with the `snakemake` extra.

!!! Note
    You can also use a separate environment for `snakemake`. In this case, you don't need to install the `'snakemake'` extra when installing `sopa`. But you may still need to install other extras, for instance `'cellpose'` if you plan to run Cellpose.

## Choose a config

Our pipeline config is a YAML file that described all the steps desired for the pipeline. It is flexible, for instance if you remove the `baysor` arguments from the config, then it will not run baysor. Similarly, if you remove the `"annotation"` section, it will not run annotation.

You can choose a config among the existing one [here](https://github.com/gustaveroussy/sopa/tree/master/workflow/config), or create your own. Keep in mind the path of your config since you'll need it later, e.g. `config/merscope/base.yaml`.

!!! note "Baysor usage"
    Even though `pip install 'sopa[baysor]'` will install some dependencies related to baysor, you still have to install the `baysor` command line (see the [official repository](https://github.com/kharchenkolab/Baysor)) if you want to use it.

    If the path to the baysor executable is not the default one (i.e. `~/.julia/bin/baysor`), update the config to provide the right path to the executable

## Run the pipeline

1. First, locate the path to your experiment file(s) of one sample. Most of the time, this is a directory containing one or many images, and eventually a transcript file (except for PhenoCycler, which is only one `.qptiff` file). You don't need to change your raw data, `sopa` will take care of this.

2. Then, activate your environment that has the snakemake command, and go to the `workflow` directory inside the `sopa` directory that you cloned earlier:
```sh
conda activate sopa    # or an environment that has `snakemake`
cd sopa/workflow       # your own personal path to the workflow directory
```

3. You can either execute the pipeline locally, or on a high-performance-cluster (choose the right option below)

=== "Local execution (e.g., personal laptop)"

    Most of the time, when executing locally, you'll run the pipeline sequentially, i.e. with one core:

    ```sh
    # replace the configfile with yours
    # replace data_path with the path to your data directory

    snakemake --config data_path=/path/to/directory --configfile=config/merscope/base.yaml --cores 1 --use-conda
    ```

=== "High-performance-cluster (e.g., Slurm cluster)"

    You'll need a cluster profile. For instance, if on a Slurm cluster, it can look like [this file](https://github.com/gustaveroussy/sopa/blob/master/workflow/slurm/config.yaml). You can create your own in `sopa/workflow/<hpc-name>/config.yaml`, or simply re-use this file (as in the command below). For more details, see the [snakemake documentation on profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).

    ```sh
    # replace the configfile with yours
    # replace data_path with the path to your data directory
    # replace profile with <hpc-name> as above, or keep 'slurm'

    snakemake --profile slurm --config data_path=/path/to/directory --configfile=config/merscope/base.yaml
    ```

    !!! note
        You may need to update the params under `resources` inside the `sopa/workflow/Snakefile` file, according to the partition names of your cluster.

For more customization, see the [snakemake CLI documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

## Toy example

In the example below, we run the pipeline on a generated toy dataset. Running it locally can help testing a new pipeline or a new config.

Make sure you have setup everything as detailed in this tutorial, and then run the following command lines:

=== "Cellpose usage"
    Make sure you have installed sopa with the Cellpose extra
    ```sh
    conda activate sopa    # or an environment that has `snakemake`
    cd sopa/workflow       # your own personal path to the workflow directory

    # replace tuto.zarr by the path where you want the data to be saved
    snakemake --config data_path=. sdata_path=tuto.zarr --configfile=config/toy/uniform_cellpose.yaml --cores 1 --use-conda
    ```

=== "Baysor usage"
    Make sure you have installed sopa with the Baysor extra, and that you have installed the `baysor` command
    ```sh
    conda activate sopa    # or an environment that has `snakemake`
    cd sopa/workflow       # your own personal path to the workflow directory

    # replace tuto.zarr by the path where you want the data to be saved
    snakemake --config data_path=. sdata_path=tuto.zarr --configfile=config/toy/uniform_baysor.yaml --cores 1 --use-conda
    ```

!!! notes
    On the above example, it executes snakemake sequentially (one core), which is enough for debugging purposes

You can then check `toy.explorer` for output files. Notably, if you have installed the [Xenium Explorer](https://www.10xgenomics.com/support/software/xenium-explorer), double-click on `experiment.xenium` to visualize the results.

## Create your own config

If the existing `config` files are not suited for your project, you can update an existing one, or create a whole new one. For this, use [this commented config](https://github.com/gustaveroussy/sopa/blob/master/workflow/config/example_commented.yaml) to understand the purpose of each argument. Note that some sections are optional: in this case, just remove the section or the argument, and sopa will not run it.