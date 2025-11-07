FROM python:3.12-slim

ARG SOPA_VERSION

# command 'ps' for nextflow
RUN apt-get update && apt-get install -y procps && apt-get clean && rm -rf /var/lib/apt/lists/*

# to avoid permission issues with the numba cache
ENV NUMBA_CACHE_DIR=/tmp/numba_cache

RUN pip install --no-cache-dir sopa[stardist]==${SOPA_VERSION}

ENV PYTHONUNBUFFERED=1
