#!/bin/bash
VERSION=`cat ./genebe/version.py | cut -f2 -d'"'`

echo "Sending version $VERSION"

docker build --no-cache --tag genebe/pygenebe:$VERSION  . && \
docker push genebe/pygenebe:$VERSION
