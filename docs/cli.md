# Sopa CLI
<div class="termy">
```console
// Run the Sopa CLI helper
$ sopa --help
 Usage: sopa [OPTIONS] COMMAND [ARGS]...                                                               
╭─ Commands ────────────────────────────────────────────────────────────╮
│ aggregate     Count transcripts and/or average channels inside cells  │
│ annotate      Perform cell-type annotation                            │
│ check         Run some sanity checks                                  │
│ crop          Crop an image based on a user-defined polygon           │
│ explorer      Convertion to the Xenium Explorer's inputs              │
│ patchify      Create patches with overlaps                            │
│ read          Read any technology data + write an SpatialData object  │
│ report        Create a HTML report of the pipeline run and some QCs   │
│ resolve       Resolve the segmentation conflicts over patches         │
│ segmentation  Perform cell segmentation on patches                    │
╰───────────────────────────────────────────────────────────────────────╯
// Example: run cellpose segmentation
$ sopa segmentation cellpose sdata.zarr
... [Logs] ...
```
</div>