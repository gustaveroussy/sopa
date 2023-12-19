# Toy examples

Toy datasets can be useful for debugging purposes or to test some new configurations. In `sopa`, we provide a toy dataset composed of a multi-channel image, a transcript layer, and a cell boundaries layer (unused in sopa).

## Generate a toy SpatialData object

### Using the API
One solution to generate it is to use the API:

```
from sopa.utils.data import uniform

sdata = uniform()
```

For more details, see the [function documentation](../api/utils/data/#sopa.utils.data.uniform).

### Using the command line
You can also use the command line to generate a new toy dataset (it will be saved in the `/path/to/directory.zarr` provided below):
```sh
sopa read . --sdata-path /path/to/directory.zarr --technology uniform
```

For more details, [see the `sopa read` command](../cli/#sopa-read).

## Run the snakemake pipeline

Running the toy dataset locally can help testing and fixing a pipeline or a new config. Setup snakemake as detailed [here](../pipeline), and then run the following command lines:
```sh
conda activate sopa    # or an environment that has `snakemake`
cd sopa/workflow       # your own personal path to the workflow directory

# replace '/path/to/directory.zarr' by the path to the uniform dataset you generated above
snakemake --config sdata_path=/path/to/directory.zarr --configfile=config/toy/uniform_cellpose.yaml --cores 1 --use-conda
```

!!! notes
    It executes snakemake sequentially (one core), which is enough for debugging purposes

!!! Baysor usage
    If you want to test Baysor, use the `config/toy/uniform_baysor.yaml` config instead of `config/toy/uniform_cellpose.yaml`

You can then check `/path/to/directory.explorer` for output files. Notably, if you have installed the [Xenium Explorer](https://www.10xgenomics.com/support/software/xenium-explorer), double-click on `experiment.xenium` to visualize the results.