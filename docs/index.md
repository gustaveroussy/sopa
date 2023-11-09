# Spatial-omics pipeline and analysis

<p align="center">
  <img src="./assets/sopa.png" alt="sopa_logo" width="400px"/>
</p>

**S**patial-**o**mics **p**ipeline and **a**nalysis in Python. Built on top of [SpatialData](https://github.com/scverse/spatialdata), it enables processing and analyses of **any** imaging-based spatial-omics using a standard data structure and output. Sopa was designed for generability and low-memory consumption on large images (scales to `1TB+` images).

The pipeline outputs are composed of (i) Xenium Explorer files for interactive visualization, (ii) a HTML report for quick quality controls, and (iii) a SpatialData `.zarr` directory for further analyses.