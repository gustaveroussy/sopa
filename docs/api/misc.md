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

## Xenium Explorer

::: sopa.io.explorer.write

::: sopa.io.explorer.align

::: sopa.io.explorer.add_explorer_selection

::: sopa.io.explorer.int_cell_id

::: sopa.io.explorer.str_cell_id

::: sopa.io.explorer.write_image

::: sopa.io.explorer.save_column_csv

## Report

::: sopa.io.write_report
