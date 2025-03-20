FROM python:3.10

RUN apt-get update && apt-get install ffmpeg libsm6 libxext6  -y

RUN pip install sopa[stardist]

ENV PYTHONUNBUFFERED=1
