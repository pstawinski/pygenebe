#!/bin/bash

# Install pytest if not available
pip install pytest pytest-cov

# Run all tests with pytest
pytest tests/ -v

# Run with coverage if desired
echo ""
echo "Running tests with coverage:"
pytest tests/ --cov=genebe --cov-report=term-missing
