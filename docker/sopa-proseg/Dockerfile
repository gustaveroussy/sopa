FROM rust:1-slim-bullseye AS rust-builder

RUN cargo install proseg

FROM python:3.10-slim

COPY --from=rust-builder /usr/local/cargo/bin/proseg /usr/local/bin/proseg

RUN pip install --no-cache-dir sopa

ENV PYTHONUNBUFFERED=1
