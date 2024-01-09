For staining-based segmentation, the Sopa CLI and pipeline are based on [Cellpose](https://github.com/MouseLand/cellpose). Yet, if desired, one can implement another staining-based segmentation algorithm, or use a multi-step segmentation process on multiple channels.

## Multi-step segmentation

Multi-step segmentation consist in running multiple times Cellpose over the whole slides with different parameters. For instance, we can first run a nucleus segmentation using DAPI, then another round using DAPI and a membrane staining, and finally using DAPI and a cell boundary staining. This can make the segmentation more robust. Note that the results of the multiple steps are combine into one segmentation result.

### 1. Save your data

For this tutorial, we use a generated dataset. The command below will generate it and save it on-disk (you can change the path `tuto.zarr` to save it somewhere else). See [here](`../../cli/#sopa-read`) for details to use your own technology.

```sh
sopa read . --sdata-path tuto.zarr --technology uniform
```

### 2. Run multi-step segmentation

Then, generate the bounding boxes of the patches on which Cellpose will be run. Here, the patches have a width and height of 1500 pixels, and an overlap of 50 pixels. We advise bigger sizes for real datasets (see our default parameters in one of our [config files](https://github.com/gustaveroussy/sopa/tree/master/workflow/config)). On the toy dataset, this will generate **4** patches.

```sh
sopa patchify image tuto.zarr --patch-width-pixel 1500 --patch-overlap-pixel 50
```

Now, we can run Cellpose on each individual patch, and for each "segmentation step" we want. On this toy example, we run 3 steps (don't forget to execute the three steps), with (i) DAPI + CK, (ii) DAPI + CD3, and (iii) DAPI + CD20.

!!! Advice
    Running the commands below manually can involve using many consecutive command, so we recommend automatizing it. For instance, this can be done using Snakemake or Nextflow. Mainly, this will help you parallelizing it, since you can run each task on seperate jobs, or using multithreading.

=== "Step 1"

    Execute the following command line on all `patch-index` (i.e., `0`, `1`, `2`, and `3`) to run Cellpose using DAPI + CK

    ```sh
    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CK \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CK \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 0
    ```

=== "Step 2"

    Execute the following command line on all `patch-index` (i.e., `0`, `1`, `2`, and `3`) to run Cellpose using DAPI + CD3

    ```sh
    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CD3 \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CD3 \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 0
    ```

=== "Step 3"

    Execute the following command line on all `patch-index` (i.e., `0`, `1`, `2`, and `3`) to run Cellpose using DAPI + CD20

    ```sh
    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CD20 \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CD20 \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 0
    ```

!!! Note
    In the above commands, the `--diameter` and `--min-area` parameters are specific to the data type we work on. For your own data, consider using the default parameters from one of our [config files](https://github.com/gustaveroussy/sopa/tree/master/workflow/config). Here, `min-area` is in pixels^2.

At this stage, you executed 12 times Cellpose (4 times on each of the three steps). Now, we need to resolve the conflict, i.e. merging the three segmentations into one. Note that we gave the paths to the temporary boundaries that we made above.
```sh
sopa resolve cellpose tuto.zarr \
    --patch-dir tuto.zarr/.sopa_cache/cellpose_CK \
    --patch-dir tuto.zarr/.sopa_cache/cellpose_CD3 \
    --patch-dir tuto.zarr/.sopa_cache/cellpose_CD20
```

### 3. Post-segmentation

Now, we can count the transcript inside each cell (by providing the name of the points dataframe, see `--gene-column genes` below), average the channels intensities inside each cell (using `average-intensities`). In the example below, we also filter cells whose average intensity if lower that `0.25 * Q90`, where Q90 is the 90th quantile (**Warning**: this may remove a lot of cells, we advise not to use this parameter when trying Sopa for the first time).
```sh
sopa aggregate tuto.zarr --gene-column genes --average-intensities --min-intensity-ratio 0.25
```

Other post-segmentations methods are available [here](../../cli). Among then, one can be used to convert the results for the Xenium Explorer:
```sh
sopa explorer write tuto.zarr --gene-column genes
```

If you have downloaded the Xenium Explorer, you can now open the results in the explorer: `open tuto.explorer/experiment.xenium` (if using a Unix operating system), or double click on the latter file.

!!! Note
    You can also use the file `tuto.explorer/adata.h5ad` if you prefer the `AnnData` object instead of the full `SpatialData` object.

## Custom staining-based segmentation

You can use your own segmentation model and plug it into Sopa to benefit from all the others functionnalities. Especially, it will scale the segmentation, since Sopa will be run on small patches.

### 1. Define your segmentation function

You need a python function as described below:

- The function input is an image of shape `(C, Y, X)` (`C` is the number of desired channels, it can be one if you want DAPI only)

- The function output is a mask of shape `(Y, X)`. This mask should contain positive values representing the segmented cells, and contain `0` outside of the cells. For instance, if 4 cells are segmented, the mask **should** contain the values 1, 2, 3, and eventually 0 (where there is no cell).

If you want to use our API, you can find a detailed example of custom segmentation [here](../../api/segmentation/stainings/#sopa.segmentation.stainings.StainingSegmentation). Else, if you want to use the CLI, continue below.

### 2. Setup

To use the CLI, you'll need to clone the repository:
```sh
git clone https://github.com/gustaveroussy/sopa.git
cd sopa
```

Then, add the function that you define above to the `sopa/segmentation/methods.py` file. An example function, called `dummy_method`, is given.

Now, install `sopa` to have your new method in the installation:
```sh
# install without extras
pip install -e .`

# install with extras
pip install -e '.[cellpose,baysor,...]'
```

### 3. Save your data

For this tutorial, we use a generated dataset. The command below will generate it and save it on-disk (you can change the path `tuto.zarr` to save it somewhere else). See [here](`../../cli/#sopa-read`) for details to use your own technology.

```sh
sopa read . --sdata-path tuto.zarr --technology uniform
```

### 4. Run your custom segmentation

Afterwards, simply call the CLI by providing the name of your function as the `<FUNCTION_NAME>` in the following commands:

- `sopa segmentation generic-staining <SDATA_PATH> --method-name <FUNCTION_NAME> ...` (see [here](../../cli/#sopa-segmentation-generic-staining) for CLI details)

- `sopa resolve generic <SDATA_PATH> --method-name <FUNCTION_NAME> --patch-dir <PATCH_DIR>` (see [here](../../cli/#sopa-resolve-generic) for CLI details)
