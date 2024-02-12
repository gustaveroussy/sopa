The Sopa CLI and pipeline are based on [Cellpose](https://github.com/MouseLand/cellpose) for staining-based segmentation. Yet, if desired, one can implement another staining-based segmentation algorithm or use a multi-step segmentation process on multiple channels.

## Multi-step segmentation

Multi-step segmentation consists of running multiple times Cellpose over the whole slides with different parameters. For instance, we can first run a nucleus segmentation using DAPI, then another round using DAPI and a membrane staining, and finally, DAPI and cell boundary staining. This can make the segmentation more robust. Note that the results of the multiple steps are combined into one final segmentation.

!!! warning
    Here, we only detail the multi-step segmentation. For the rest of the CLI usage, refer to our [CLI usage tutorial](../cli_usage), and only replace the "Run segmentation" section with the instructions below.

First, generate the bounding boxes of the patches on which Cellpose will be run. Here, the patches have a width and height of 1500 pixels, and an overlap of 50 pixels. We advise bigger sizes for real datasets (see our default parameters in one of our [config files](https://github.com/gustaveroussy/sopa/tree/master/workflow/config)). On the toy dataset, this will generate **4** patches.

```sh
sopa patchify image tuto.zarr --patch-width-pixel 1500 --patch-overlap-pixel 50
```

Now, we can run Cellpose on each of the four patches and for each "segmentation step" we want. In this toy example, we run 3 steps with (i) DAPI + CK, (ii) DAPI + CD3, and (iii) DAPI + CD20.

!!! tip
    Manually running the commands below can involve using many consecutive commands, so we recommend automatizing it. For instance, this can be done using Snakemake or Nextflow. Mainly, this will help you parallelize it since you can run each task on separate jobs or using multithreading. You can also see how we do it in the [Sopa Snakemake pipeline](https://github.com/gustaveroussy/sopa/blob/master/workflow/Snakefile).

    To automatically get the number of patches, you can either open the `tuto.zarr/.sopa_cache/patches_file_image` file or compute `len(sdata['sopa_patches'])` in Python.

=== "Patch 0"

    ```sh
    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CK \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CK \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 0

    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CD3 \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CD3 \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 0

    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CD20 \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CD20 \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 0
    ```
=== "Patch 1"

    ```sh
    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CK \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CK \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 1

    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CD3 \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CD3 \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 1

    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CD20 \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CD20 \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 1
    ```
=== "Patch 2"

    ```sh
    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CK \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CK \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 2

    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CD3 \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CD3 \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 2

    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CD20 \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CD20 \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 2
    ```
=== "Patch 3"

    ```sh
    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CK \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CK \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 3

    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CD3 \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CD3 \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 3

    sopa segmentation cellpose tuto.zarr \
        --channels DAPI --channels CD20 \
        --patch-dir tuto.zarr/.sopa_cache/cellpose_CD20 \
        --diameter 35 \
        --min-area 2000 \
        --patch-index 3
    ```

!!! Note
    In the above commands, the `--diameter` and `--min-area` parameters are specific to the data type we work on. For your own data, consider using the default parameters from one of our [config files](https://github.com/gustaveroussy/sopa/tree/master/workflow/config). Here, `min-area` is in pixels^2.

At this stage, you executed 12 times Cellpose (3 steps on each of the 4 patches). Now, we need to resolve the conflict, i.e., merge the three segmentations into one. Note that we gave the paths to the temporary boundaries we made above.
```sh
sopa resolve cellpose tuto.zarr \
    --patch-dir tuto.zarr/.sopa_cache/cellpose_CK \
    --patch-dir tuto.zarr/.sopa_cache/cellpose_CD3 \
    --patch-dir tuto.zarr/.sopa_cache/cellpose_CD20
```

Congrats, you have now merged the results of a three-step segmentation! You can now refer to our normal [CLI usage tutorial](../cli_usage) for all the other tasks.

## Custom staining-based segmentation

You can plug your segmentation model into Sopa to benefit from all the other functionalities. In particular, it will scale the segmentation since Sopa will run on small patches.

### 1. Define your segmentation function

You need a Python function as described below:

- The function input is an image of shape `(C, Y, X)` (`C` is the number of desired channels; it can be one if you want DAPI only)

- The function output is a mask of shape `(Y, X)`. This mask should contain positive values representing the segmented cells, and contain `0` outside of the cells. For instance, if 4 cells are segmented, the mask **should** contain the values 1, 2, 3, and eventually 0 (where there is no cell).

If you want to use our API, you can find a detailed example of custom segmentation [here](../../api/segmentation/stainings/#sopa.segmentation.stainings.StainingSegmentation). Else, if you want to use the CLI, continue below.

This function should be wrapped into a "method builder". The purpose is that you can provide arguments to the method builder that your desired function can then use. For instance, this is a dummy method builder whose segmentation function only creates one cell:

```python
def dummy_method(**method_kwargs):
    """A method builder (i.e. it returns a segmentation function).
    Kwargs can be provided and used in the below function"""

    def segmentation_function(image: np.ndarray) -> np.ndarray:
        """A dummy example of a custom segmentation method
        that creates one cell (with a padding of 10 pixels).

        Args:
            image: An image of shape `(C, Y, X)`

        Returns:
            A mask of shape `(Y, X)` containing one cell
        """
        mask = np.zeros(image.shape[1:], dtype=int)

        # one cell, corresponding to value 1
        mask[10:-10, 10:-10] = 1  # squared shaped

        return mask

    return segmentation_function
```

### 2. Setup

To use the CLI, you'll need to clone the repository. Also, we recommend installing Sopa in dev mode, allowing you to update your segmentation function without re-installing everything. For instance:

```sh
git clone https://github.com/gustaveroussy/sopa.git
cd sopa

# create an empty conda env
conda create --name sopa python=3.10
conda activate sopa

# install the extras you want
pip install -e ".[cellpose,baysor,tangram]"
```

### 3. Run your custom segmentation

!!! warning
    Here, we only detail the multi-step segmentation. For the rest of the CLI usage, refer to our [CLI usage tutorial](../cli_usage), and only replace the "Run segmentation" section with the instructions below.

Similarly to the tutorial for Cellpose, first start by generating patches:

```sh
sopa patchify image tuto.zarr --patch-width-pixel 1500 --patch-overlap-pixel 50
```

Then, call the CLI by providing the name of your method builder, i.e. replace `<BUILDER_NAME>` in the following commands:

!!! note
    The process below is similar to the one of the Cellpose tutorial in the [CLI usage tutorial](../cli_usage).

=== "Patch 0"
    ```sh
    sopa segmentation generic-staining tuto.zarr \
        --method-name <BUILDER_NAME> \
        --channels DAPI \
        --patch-index 0
    ```
=== "Patch 1"
    ```sh
    sopa segmentation generic-staining tuto.zarr \
        --method-name <BUILDER_NAME> \
        --channels DAPI \
        --patch-index 1
    ```
=== "Patch 2"
    ```sh
    sopa segmentation generic-staining tuto.zarr \
        --method-name <BUILDER_NAME> \
        --channels DAPI \
        --patch-index 2
    ```
=== "Patch 3"
    ```sh
    sopa segmentation generic-staining tuto.zarr \
        --method-name <BUILDER_NAME> \
        --channels DAPI \
        --patch-index 3
    ```

And, to resolve the conflicts:

```sh
sopa resolve generic tuto.zarr --method-name <BUILDER_NAME>
```

Finally, aggregate the results:

```sh
sopa aggregate tuto.zarr --method-name <BUILDER_NAME> --gene-column genes --average-intensities
```

Congrats, you have now run your custom segmentation method! You can now refer to our normal [CLI usage tutorial](../cli_usage) for all the other tasks.
