FROM rust:latest AS rust-builder

RUN cargo install proseg

FROM python:3.12-slim

ARG SOPA_VERSION

COPY --from=rust-builder /usr/local/cargo/bin/proseg /usr/local/bin/proseg

RUN apt-get update && apt-get install -y procps && apt-get clean && rm -rf /var/lib/apt/lists/*

# to avoid permission issues with the numba cache
ENV NUMBA_CACHE_DIR=/tmp/numba_cache

RUN pip install --no-cache-dir sopa==${SOPA_VERSION}

ENV PYTHONUNBUFFERED=1
