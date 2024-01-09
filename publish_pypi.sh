#!/bin/bash
source /ssd2/pio/safessd/workspace/workspace.python/genebe/.venv/bin/activate
python -m build
twine upload dist/*
rm -rf build
