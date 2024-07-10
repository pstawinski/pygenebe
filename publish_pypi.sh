#!/bin/bash
source /ssd2/pio/safessd/workspace/workspace.python/genebe/.venv/bin/activate
python -m pip install build
python -m build --no-isolation
twine upload dist/*
rm -rf build
