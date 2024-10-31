# Multi-step segmentation

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
