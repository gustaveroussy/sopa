# Sopa docker containers

This directory contains the different Dockerfile for Sopa, depending on the segmentation tool to be used. The default `sopa` image contains no extra.

### Installation

The docker images can be found on [Docker Hub](https://hub.docker.com/repository/docker/quentinblampey/sopa/tags), and are available for **Linux** only (for now).

For instance, you can install the sopa images as follows:

```sh
docker pull quentinblampey/sopa:2.0.3 # image without any extra
docker pull quentinblampey/sopa:2.0.3-proseg # image with sopa+proseg
```

### Notes
- Numba issue (https://github.com/numba/numba/issues/4032) solved via `ENV NUMBA_CACHE_DIR=/tmp/numba_cache`
- `ffmpeg libsm6 libxext6` were installed because of the `cv2` dependency. This can be removed after [#239](https://github.com/gustaveroussy/sopa/pull/239) is merged, as it will remove the `cv2` dependency.
- `procps` was installed so that `nextflow` can run using these images
