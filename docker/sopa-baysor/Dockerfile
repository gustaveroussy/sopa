FROM python:3.12-slim-bookworm

ARG SOPA_VERSION

RUN apt-get update && apt-get install -y wget unzip procps && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN cd /tmp && wget 'https://github.com/kharchenkolab/Baysor/releases/download/v0.7.1/baysor-x86_x64-linux-v0.7.1_build.zip' && unzip baysor-x86_x64-linux-v0.7.1_build.zip && mv ./bin/baysor /opt/ && rm ./baysor-x86_x64-linux-v0.7.1_build.zip
ENV PATH="/opt/baysor/bin:${PATH}"

# to avoid permission issues with the numba cache
ENV NUMBA_CACHE_DIR=/tmp/numba_cache

RUN pip install --no-cache-dir sopa[baysor]==${SOPA_VERSION}

ENV PYTHONUNBUFFERED=1
