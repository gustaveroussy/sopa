
### Option 3: ComSeg


[ComSeg](https://github.com/fish-quant/ComSeg) is a transcript-based segmentation method. It uses a segmentation prior (here, Cellpose) and improves it using the transcripts information.

#### Run Cellpose to segment nuclei

```
sopa patchify image tuto.zarr --patch-width-pixel 1500 --patch-overlap-pixel 50
sopa segmentation cellpose tuto.zarr --channels DAPI --diameter 35 --min-area 2000
sopa resolve cellpose tuto.zarr
```

####  Save a ComSeg config file as config.jsons
More information on the parameters can be found in the [ComSeg documentation](https://comseg.readthedocs.io/en/latest/userguide/Minimal_example.html).
Below we display a minimal example of a ComSeg config file.


```json
{"dict_scale": {"x": 1, "y": 1, "z": 1},
"mean_cell_diameter": 15,
"max_cell_radius": 50,
"alpha": 0.5,
"min_rna_per_cell": 5,
"gene_column": "genes"}
```

####  Run ComSeg with the sopa command line tool

1) create the ComSeg patches
On the toy dataset, we will generate 4 patches.
```
sopa patchify comseg tuto.zarr --config-path config.json --patch-width-microns 200 --patch-overlap-microns 50
```

2) run ComSeg on all patches

!!! tip
    Manually running the commands below can involve using many consecutive commands, so we recommend automatizing it. For instance, this can be done using Snakemake or Nextflow. This will help you parallelize it since you can run each task on separate jobs or using multithreading. You can also see how we do it in the [Sopa Snakemake pipeline](https://github.com/gustaveroussy/sopa/blob/master/workflow/Snakefile).

    To automatically get the number of patches, you can open the `tuto.zarr/.sopa_cache/patches_file_comseg` file. This lists the names of the directories inside `tuto.zarr/.sopa_cache/comseg` related to each patch. If you selected an ROI, the excluded patches are effectively not in the `patches_file_comseg` file.

=== "Patch 0"
    ```sh
    cd tuto.zarr/.sopa_cache/comseg_boundaries/0

    # 'comseg' is the official comseg executable. If unavailable, replace it with your path to the executable
    comseg run --save-polygons GeoJSON -c config.toml transcripts.csv
    ```
=== "Patch 1"
    ```sh
    cd tuto.zarr/.sopa_cache/comseg_boundaries/1

    # 'comseg' is the official comseg executable. If unavailable, replace it with your path to the executable
    comseg run --save-polygons GeoJSON -c config.toml transcripts.csv
    ```
=== "Patch 2"
    ```sh
    cd tuto.zarr/.sopa_cache/comseg_boundaries/2

    # 'comseg' is the official comseg executable. If unavailable, replace it with your path to the executable
    comseg run --save-polygons GeoJSON -c config.toml transcripts.csv
    ```
=== "Patch 3"
    ```sh
    cd tuto.zarr/.sopa_cache/comseg_boundaries/3

    # 'comseg' is the official comseg executable. If unavailable, replace it with your path to the executable
    comseg run --save-polygons GeoJSON -c config.toml transcripts.csv
    ```

3) Merge the results
```sh
sopa resolve comseg tuto.zarr --gene-column genes
```