# Sopa command-line-interface (CLI)

## Usage

When installing `sopa` are written in our [getting-started guidelines](../getting_started), a new command named `sopa` becomes available.

!!! note "CLI helper"
    Run `sopa --help` to get details about all the command line purpose. You can also use this helper on any subcommand, for instance `sopa read --help`.

<div class="termy">
```console
// Run the Sopa CLI helper
$ sopa --help
 Usage: sopa [OPTIONS] COMMAND [ARGS]...    
╭─ Commands ─────────────────────────────────────────────────────╮
│ aggregate     Aggregate transcripts/channels inside cells      │
│ annotate      Perform cell-type annotation                     │
│ check         Run some sanity checks                           │
│ crop          Crop an image based on a user-defined polygon    │
│ explorer      Convertion to the Xenium Explorer's inputs       │
│ patchify      Create patches with overlaps                     │
│ read          Read any technology + write a SpatialData object │
│ report        Create a web-report with figures/QCs             │
│ resolve       Resolve the segmentation conflicts over patches  │
│ segmentation  Perform cell segmentation on patches             │
╰────────────────────────────────────────────────────────────────╯
// Example: run cellpose segmentation
$ sopa segmentation cellpose sdata.zarr
... [Logs] ...
```
</div>

### Arguments vs options

When running the `sopa` CLI, some arguments are required while some are optional. For instance, for the `sopa read` command, `sdata_path` is an argument and a path has to be given directly, while `technology` is an option, and in this case the `--technology` prefix has to be used. For instance, if you read MERSCOPE data, it will be:

```
sopa read /path/to/merscope/directory --technology merscope
```

Note that `/path/to/merscope/directory` refers to `sdata_path`. Thus, you **should not** add the suffix `--sdata_path` for arguments.

!!! note "Required options"
    All the arguments are required, as shown by the `[required]` hint on the CLI helper. Note that some options may also be required too (in this case, the term `[required]` will appear on the CLI helper).


### Sub-commands

Some commands contains subcommands. This is the case e.g. for `sopa patchify`, which has two subcommands: `sopa patchify cellpose` and `sopa patchify baysor`.

## Commands

### `sopa read`

```sh
 Usage: sopa read [OPTIONS] DATA_PATH                                                                    
                                                                                                         
 Read any technology data, and write a standardized SpatialData object                                   
 Args:                                                                                                   
 data_path: Path to one data sample (most of the time, this corresponds to a directory)                  
 technology: Name of the technology used to collected the data (e.g., 'xenium', 'merfish', ...)          
 sdata_path: Optional path to write the SpatialData object. If not provided, will write to the '{data_path}.zarr' directory                                                                            
 config_path: Path to the snakemake config. This can be useful in order not to provide the 'technology' and the 'kwargs' arguments                                                                              
 kwargs: Dictionary provided to the reader function.                                                     
                                                                                                         
╭─ Arguments ───────────────────────────────────────────────────────────────────────────────────────────╮
│ *    data_path      TEXT  [default: None] [required]                                                  │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ─────────────────────────────────────────────────────────────────────────────────────────────╮
│ --technology         TEXT  [default: None]                                                            │
│ --sdata-path         TEXT  [default: None]                                                            │
│ --config-path        TEXT  [default: None]                                                            │
│ --kwargs             TEXT  [default: {}]                                                              │
│ --help                     Show this message and exit.                                                │
╰───────────────────────────────────────────────────────────────────────────────────────────────────────╯
```