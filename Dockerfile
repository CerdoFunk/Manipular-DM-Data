FROM python:3
WORKDIR /workdir
COPY . .
RUN pip install \
    black