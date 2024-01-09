#
# build using
# docker build --no-cache --tag genebe/pygenebe:VERSION  .
# docker push genebe/pygenebe:VERSION
#

# Use an official Python runtime as a parent image
FROM python:3.11

# Set the working directory in the container
WORKDIR /app

# Install genebe package from PyPI
RUN pip install -U --no-cache-dir genebe[cpp]

# Run genebe script when the container launches
CMD ["genebe"]


