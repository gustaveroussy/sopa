# Misc

## Settings

### Disabling auto-save

By default, if your `SpatialData` object is stored on-disk, it will also store the new elements on disk.

You can disable this behavior as follow:

```python
sopa.settings.auto_save_on_disk = False
```

### Parallelization backends

Some methods (for instance `sopa.segmentation.cellpose`) may need a parallelization backend to run fast enough.

You can set it as below:

```python
sopa.settings.parallelization_backend = "dask" # using dask
sopa.settings.parallelization_backend = None # no backend (i.e., sequential)
```

!!! warning
    The `dask` backend is still experimental. You can add a comment to [this issue](https://github.com/gustaveroussy/sopa/issues/145) to help us improve it.

You can also pass some kwargs to the [dask Client](https://distributed.dask.org/en/stable/api.html#client):
```python
sopa.settings.dask_client_kwargs["n_workers"] = 4
```

Otherwise, if you don't use the API, you can also set the `SOPA_PARALLELIZATION_BACKEND` env variable, e.g.:
```sh
export SOPA_PARALLELIZATION_BACKEND=dask
```

### Gene filtering

Use `sopa.settings.gene_exclude_pattern` to filter out gene names during segmentation and aggregation. By default, we use the variable below:
```python
sopa.settings.gene_exclude_pattern: str | None = "negcontrol.*|blank.*|antisense.*|unassigned.*|deprecated.*|intergenic.*"
```
Use `sopa.settings.gene_exclude_pattern = None` to keep all genes.

## Shapes operations

::: sopa.shapes.expand_radius

::: sopa.shapes.remove_overlap

::: sopa.shapes.vectorize

## Xenium Explorer

::: sopa.io.explorer.write

::: sopa.io.explorer.align

::: sopa.io.explorer.add_explorer_selection

::: sopa.io.explorer.int_cell_id

::: sopa.io.explorer.str_cell_id

::: sopa.io.explorer.write_image

::: sopa.io.explorer.write_cell_categories

::: sopa.io.explorer.save_column_csv

## Report

::: sopa.io.write_report
