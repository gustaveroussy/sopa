## Cellpose is not segmenting enough cells, what should I do?
- Maybe `min_area` is too high, and all the cells are filtered because they are smaller than this area. Remind that, when using Cellpose, the areas correspond to pixels^2.
- This can be due to a low image quality. If the image is too pixelated, consider increasing `gaussian_sigma` (e.g., `2`) under the cellpose parameters of our config. If the image has a low contrast, consider increasing `clip_limit` (e.g., `0.3`). These parameters are detailed in [this example config](https://github.com/gustaveroussy/sopa/blob/master/workflow/config/example_commented.yaml).
- Consider updating the official Cellpose parameters. In particular, try `cellprob_threshold=-6` and `flow_threshold=2`.

## How to align H&E or fluorescence data?

1. Convert your image with QuPath as written in this [10x genomics webpage](https://www.10xgenomics.com/support/software/xenium-explorer/tutorials/xe-image-file-conversion), or use our API to write the image:
```python
from sopa import io
from sopa.io.explorer import write_image

image = io.ome_tif("path/to/your/image.tif") # or use a different reader
write_image("where/to/save/image.ome.tif", image, is_dir=False)
```
1. On the Xenium Explorer, under the "Images" panel, click on "Add image", and follow the instructions.
2. [Optional] After alignment, you can export the transformation matrix as a `csv` file. Then, use [our cli](../cli/#sopa-explorer-add-aligned) to update the SpatialData object.

## Can I use Nextflow instead of Snakemake?

Nextflow is not supported yet, but we are working on it. If you want, you can also help re-writing our Snakemake pipeline for Nextflow.

## I have another issue, how to fix it?
- Make sure you have clearly read our CLI or API (depending on what you use).
- Don't hesitate to open an issue on [Sopa's Github repository](https://github.com/gustaveroussy/sopa/issues), and detail your issue with as much precision as possible, in order for the maintainers to be able to reproduce it.