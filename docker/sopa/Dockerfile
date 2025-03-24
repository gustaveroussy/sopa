FROM python:3.10-slim

# needed for the open-cv dependency
RUN apt-get update && apt-get install -y \
    ffmpeg \
    libsm6 \
    libxext6 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir sopa[cellpose]

ENV PYTHONUNBUFFERED=1
