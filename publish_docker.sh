#!/bin/bash
VERSION=0.0.14

docker build --no-cache --tag genebe/pygenebe:$VERSION  . && \
docker push genebe/pygenebe:$VERSION
