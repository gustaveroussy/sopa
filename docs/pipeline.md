# Snakemake pipeline

If you don't want to dig into the CLI/API, you can directly use our existing [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. It will not involve any coding, but requires some setup for `snakemake`.

## Setup

- Clone the `sopa` repository, i.e., `git clone https://github.com/gustaveroussy/sopa.git`
- Make sure you have a `conda` environment called `sopa` (see [getting started](../getting_started)), on which you installed `sopa` with the `snakemake` extra.

!!! Note
    You can also use a separate environment for `snakemake`. In this case, you don't need to install the snakemake extra when installing `sopa`

## Choose a config

Our pipeline config is a YAML file that described all the steps desired for the pipeline. It is flexible, for instance if you remove the `baysor` arguments from the config, then it will not run baysor. Similarly, if you remove the `"annotation"` section, it will not run annotation.

You can choose a config among the existing one [here](https://github.com/gustaveroussy/sopa/tree/master/workflow/config), or create your own. Keep in mind the path of your config since you'll need it later, e.g. `config/merscope/base.yaml`.

## Run the pipeline

Activate your environment that has the snakemake command, and go to the `workflow` directory inside the `sopa` directory that you cloned earlier:
```sh
conda activate sopa    # or an environment that has `snakemake`
cd sopa/workflow       # your own personal path to the workflow directory
```

Then, you can either execute the pipeline locally, or on a high-performance-cluster (see below). 

### Local execution

First, locate the path to your experiment file(s) of one sample. Most of the time, this is a directory containing one or many images, and eventually a transcript file. You don't need to change your raw data, `sopa` will take care of this.

Most of the time, when executing locally, you'll run the pipeline sequentially with one core, i.e.:

```sh
# replace the configfile with yours
# replace data_path with the path to your data directory

snakemake --config data_path=/path/to/directory --configfile=config/merscope/base.yaml --cores 1
```

For more customization, see the [snakemake CLI documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

### High-performance-cluster execution

... TODO

## Configure your own config

... TODO