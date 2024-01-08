For staining-based segmentation, the Sopa CLI and pipeline are based on [Cellpose](https://github.com/MouseLand/cellpose). Yet, if desired, one can implement another staining-based segmentation algorithm, or use a multi-step segmentation process on multiple channels.

## Multi-step segmentation

### Using the CLI [WIP]

```sh
sopa read . --sdata-path tuto.zarr --technology uniform
```

```sh
sopa patchify image tuto.zarr --patch-width-pixel 1500 --patch-overlap-pixel 50
```

```sh
sopa segmentation cellpose tuto.zarr \
    --channels DAPI --channels CK \
    --patch-index 0 \
    --patch-dir tuto.zarr/.sopa_cache/cellpose_CK \
    --diameter 40 \
    --min-area 1000 --clip-limit 0.01 # these are optional parameters
```

Same for CD3 and CD20

```sh
sopa resolve cellpose tuto.zarr \
    --patch-dir tuto.zarr/.sopa_cache/cellpose_CK \
    --patch-dir tuto.zarr/.sopa_cache/cellpose_CD3 \
    --patch-dir tuto.zarr/.sopa_cache/cellpose_CD20
```

```sh
sopa aggregate tuto.zarr --gene-column genes --average-intensities --min-intensity-ratio 0.25
```

```sh
sopa explorer write tuto.zarr --gene-column genes
```

If you have downloaded the Xenium Explorer, you can now open the results in the explorer: `open tuto.explorer/experiment.xenium`

You can also use the file `tuto.explorer/adata.h5ad` if you prefer the `AnnData` object instead of the full `SpatialData` object.

This can be automatized in a snakemake pipeline, ...

## Custom staining-based segmentation

You can use your own segmentation model and plug it into Sopa to benefit from all the others functionnalities. Especially, it will scale the segmentation, since Sopa will be run on small patches.

For this, you need a python function as described below:

- The function input is an image of shape `(C, Y, X)` (`C` is the number of desired channels, it can be one if you want DAPI only)

- The function output is a mask of shape `(Y, X)`. This mask should contain positive values representing the segmented cells, and contain `0` outside of the cells. For instance, if 4 cells are segmented, the mask **should** contain the values 1, 2, 3, and eventually 0 (where there is no cell).

### Using the API

An example of custom segmentation using the API is detailed [here](../../api/segmentation/stainings/#sopa.segmentation.stainings.StainingSegmentation).

### Using the CLI

To use the CLI here, you'll need to clone the repository:
```sh
git clone https://github.com/gustaveroussy/sopa.git
cd sopa
```

Then, add your method inside the `sopa/segmentation/methods.py` file. An example function, called `dummy_method`, is given.

Now, install `sopa` to have your new method in the installation:
```sh
# install without extras
pip install -e .`

# install with extras
pip install -e '.[cellpose,baysor,...]'
```

Afterwards, simply call the CLI by providing the name of your function as the `<FUNCTION_NAME>` in the following commands:

- `sopa segmentation generic-staining <SDATA_PATH> --method-name <FUNCTION_NAME> ...` (see [here](../../cli/#sopa-segmentation-generic-staining) for CLI details)

- `sopa resolve generic <SDATA_PATH> --method-name <FUNCTION_NAME> --patch-dir <PATCH_DIR>` (see [here](../../cli/#sopa-resolve-generic) for CLI details)
