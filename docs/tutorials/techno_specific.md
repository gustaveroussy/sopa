# Technology-specific advice

Sopa is designed to run on all spatial technologies at single-cell resolution, but some choices or parameters may be more or less suited to certain technologies.

!!! note
    All these advice are for the API usage, but they also apply for the CLI and the Snakemake pipeline (update the usage accordingly, or look at our [per-technology Snakemake configs](https://github.com/gustaveroussy/sopa/tree/master/workflow/config)).

## Xenium or MERSCOPE

For Xenium or MERSCOPE data, we recommend using a transcript based segmentation tool, e.g. Baysor or Proseg. With the Xenium 5k panel, Proseg is recommended due to fast results and high quality.

To use the prior 10X or Vizgen segmentation, you need to precise it in [`sopa.make_transcript_patches`](../../api/patches/#sopa.make_transcript_patches) - with `sopa>=2.0.4`, use `prior_shapes_key="auto"` to automatically detect the prior and use it for Baysor/Proseg.

```python
sopa.make_transcript_patches(sdata, prior_shapes_key="auto") # for proseg, also add patch_width=None
```

If using [Baysor](../../api/segmentation/#sopa.segmentation.baysor), you can also use `area=20` to filter cells with an area below 20 um².

To filter low quality cells, you can provide `min_transcripts=20` to the [aggregation](../../api/aggregation/#sopa.aggregate).

## CosMX

For CosMX data, the same advice as above are applicable, although you may experience some issue when reading the data with `sopa.io.cosmx` due to some frequent changes in the AtoMX exports. If so, please [open an issue](https://github.com/gustaveroussy/sopa/issues) to improve this.

## Visium HD

See [here](../visium_hd) the full Visium HD specific tutorial.

## MACSima or PhenoCycler

For spatial proteomics, [Cellpose](../../api/segmentation/#sopa.segmentation.cellpose) is needed. To avoid having artefacts outside of the tissue, we recommend running [`sopa.segmentation.tissue`](../../api/segmentation/#sopa.segmentation.tissue) before running Cellpose. This way, Cellpose will run only inside the tissue.

For MACSima, we recommend to use `diameter=35` pixels and `min_area=400` pixels² to [`sopa.segmentation.cellpose`](../../api/segmentation/#sopa.segmentation.cellpose). For PhenoCycler data, it will depend on the image resolution, so you'll need to choose a diameter (in pixels) that is biologically relevant.

During [aggregation](../../api/aggregation/#sopa.aggregate), you can for instance use `expand_radius_ratio=0.1` to expand the cells, `min_intensity_ratio=0.1` to filter cells with a too low intensity.

## Other

For other technologies not listed here, please [open a new GitHub issue](https://github.com/gustaveroussy/sopa/issues).
