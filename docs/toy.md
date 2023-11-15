# Toy examples

Toy datasets can be useful for debugging purposes or to test some new configurations. In `sopa`, we provide a toy dataset composed of a multi-channel image, a transcript layer, and a cell boundaries layer (unused in sopa).

## Generate a toy SpatialData object

### Using the API
One solution to generate it is to use the API:

```
from sopa.utils.data import uniform

sdata = uniform()
```

For more details, see TODO.

### Using the command line
You can also use the command line to generate a new toy dataset: `sopa read . --sdata-path /path/to/directory.zarr --technology uniform` and play with it.

## Run the snakemake pipeline

Most of the time, you will run debugging locally. TODO

Clone the repository, and at the root of `/sopa`, run:
```
cd workflow
snakemake -c1 --configfile=config/toy/uniform.yaml --config data_path=/path/where/save/uniform
```