# Aggregation

!!! tips "Recommendation"
    We recommend using the `sopa.aggregate` function below, which is a wrapper for all types of aggregation. Internally, it uses `aggregate_channels`, `count_transcripts`, and/or `aggregate_bins`, which are also documented below if needed. For marker expression prediction, see `sopa.aggregation.nimbus`.

::: sopa.aggregate

::: sopa.aggregation.aggregate_channels

::: sopa.aggregation.count_transcripts

::: sopa.aggregation.aggregate_bins

::: sopa.overlay_segmentation

## Nimbus aggregation

Use `sopa.aggregation.nimbus` to predict marker expression confidence scores.

!!! info "External dependency"
    Nimbus aggregation requires installing [`Nimbus-Inference`](https://nimbus-inference.readthedocs.io/en/latest/?badge=latest). You can install it via pip:

    ```bash
    pip install Nimbus-Inference
    ```

Example:

```python
anndata = sopa.aggregation.nimbus(
    sdata,
    image_key="raw_image",
    shapes_key="cellpose_boundaries",
)
```

::: sopa.aggregation.nimbus
