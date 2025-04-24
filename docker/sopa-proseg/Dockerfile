FROM rust:1-slim-bullseye AS rust-builder

RUN cargo install proseg

FROM python:3.10-slim

COPY --from=rust-builder /usr/local/cargo/bin/proseg /usr/local/bin/proseg

RUN apt-get update && apt-get install -y procps && apt-get clean && rm -rf /var/lib/apt/lists/*

# to avoid permission issues with the numba cache
ENV NUMBA_CACHE_DIR=/tmp/numba_cache

RUN pip install --no-cache-dir sopa==2.0.6

ENV PYTHONUNBUFFERED=1
