FROM python:3.10-slim

ARG SOPA_VERSION

# to avoid permission issues with the numba cache
ENV NUMBA_CACHE_DIR=/tmp/numba_cache

RUN pip install --no-cache-dir sopa==${SOPA_VERSION}

ENV PYTHONUNBUFFERED=1
