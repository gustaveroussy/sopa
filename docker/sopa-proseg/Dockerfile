FROM rust:latest AS rust-builder

RUN cargo install proseg

FROM python:3.10

COPY --from=rust-builder /usr/local/cargo/bin/proseg /usr/local/bin/proseg

RUN pip install sopa

ENV PYTHONUNBUFFERED=1
