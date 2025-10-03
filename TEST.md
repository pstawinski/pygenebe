# Testing Guide

This document provides comprehensive instructions for running tests in the GeneBe Client project.

## Overview

The project uses [pytest](https://pytest.org/) as the testing framework, with tests organized in both the `tests/` directory and at the project root level.

## Test Structure

```
genebe-client/
├── tests/                    # Main test directory
│   ├── test_simple.py       # Basic functionality tests
│   ├── test_transcriptid.py # Transcript ID encoding/decoding tests
│   └── test_variants.py     # Variant annotation tests
├── test_genebe.py           # Legacy tests (still functional)
├── pytest.ini              # Pytest configuration
└── test.sh                 # Convenient test runner script
```

## Prerequisites

### Install Dependencies

1. **Install pytest and related packages:**
   ```bash
   pip install pytest pytest-cov pytest-mock
   ```

2. **Or install test dependencies via setup.py:**
   ```bash
   pip install -e .[test]
   ```

### Environment Setup

Since the project doesn't need to be installed to run tests, ensure the Python path includes the project directory:

```bash
export PYTHONPATH=/path/to/genebe-client:$PYTHONPATH
```

Or use the `-e` flag with pytest to automatically handle the path.

## Running Tests

### Quick Start

**Option 1: Use the provided test script**
```bash
./test.sh
```

**Option 2: Run pytest directly**
```bash
# Set Python path and run all tests
PYTHONPATH=. pytest -v
```

### Basic Test Commands

```bash
# Run all tests with verbose output
pytest -v

# Run tests in a specific directory
pytest tests/ -v

# Run a specific test file
pytest tests/test_simple.py -v

# Run a specific test function
pytest tests/test_simple.py::test_basic_encoding -v
```

### Test Categories and Markers

The project uses custom pytest markers to categorize tests:

- `@pytest.mark.unit` - Unit tests (fast, no external dependencies)
- `@pytest.mark.integration` - Integration tests
- `@pytest.mark.network` - Tests requiring network connectivity
- `@pytest.mark.slow` - Tests that may take longer to complete

#### Running Tests by Category

```bash
# Run only unit tests
pytest -m unit -v

# Run tests excluding network-dependent ones
pytest -m "not network" -v

# Run only fast tests (exclude slow ones)
pytest -m "not slow" -v

# Combine markers
pytest -m "unit or integration" -v
```

### Filtering Tests

```bash
# Run tests matching a pattern in the name
pytest -k "transcriptid" -v

# Exclude specific test patterns
pytest -k "not lift_over" -v

# Run tests excluding network-dependent functionality
pytest -k "not (lift_over or annotate_variants_list_liftover or parse_variants_multiple)" -v
```

## Test Coverage

### Generate Coverage Reports

```bash
# Run tests with coverage
pytest --cov=genebe --cov-report=term-missing

# Generate HTML coverage report
pytest --cov=genebe --cov-report=html

# Generate XML coverage report (useful for CI)
pytest --cov=genebe --cov-report=xml
```

The HTML report will be generated in `htmlcov/index.html`.

## Test Configuration

The project uses `pytest.ini` for configuration:

```ini
[tool:pytest]
testpaths = tests .                    # Look for tests in tests/ and root
python_files = test_*.py *_test.py     # Test file patterns
python_classes = Test*                 # Test class patterns
python_functions = test_*              # Test function patterns
addopts =
    -v                                 # Verbose output
    --tb=short                         # Short traceback format
    --strict-markers                   # Require marker registration
    --disable-warnings                 # Suppress warnings
    --color=yes                        # Colored output
```

## Network-Dependent Tests

Some tests require a GeneBe server running locally:

- **Server endpoint:** `http://localhost:7180`
- **Tests affected:** Liftover functionality, variant annotation with remote API

### Running Without Network Tests

If you don't have a local GeneBe server running:

```bash
# Method 1: Use keyword filtering
pytest -k "not (lift_over or annotate_variants_list_liftover)" -v

# Method 2: Use markers (when network tests are properly marked)
pytest -m "not network" -v
```

## Continuous Integration

### Example CI Commands

```bash
# Basic test run
pytest -v --tb=short

# With coverage for CI
pytest --cov=genebe --cov-report=xml --cov-report=term

# Exclude network tests in CI
pytest -m "not network" --cov=genebe --cov-report=xml
```

## Troubleshooting

### Common Issues

1. **ModuleNotFoundError: No module named 'genebe'**
   ```bash
   # Solution: Set PYTHONPATH
   PYTHONPATH=. pytest -v
   ```

2. **Connection refused errors**
   ```bash
   # These are expected if GeneBe server isn't running locally
   # Run without network tests:
   pytest -k "not lift_over" -v
   ```

3. **Import errors in tests**
   ```bash
   # Ensure you're in the project root directory
   cd /path/to/genebe-client
   pytest -v
   ```

### Debug Mode

```bash
# Run with debug information
pytest -v -s --tb=long

# Run with Python debugger on failures
pytest --pdb

# Stop on first failure
pytest -x
```

## Writing New Tests

### Test File Naming

- Place test files in the `tests/` directory
- Name files as `test_*.py` or `*_test.py`
- Name test functions as `test_*`

### Example Test Structure

```python
import pytest
from genebe import some_function

def test_basic_functionality():
    """Test basic functionality."""
    result = some_function("input")
    assert result == "expected_output"

@pytest.mark.network
def test_network_dependent():
    """Test requiring network access."""
    # This test will be skipped when running with -m "not network"
    pass

@pytest.fixture
def sample_data():
    """Provide test data."""
    return {"key": "value"}

def test_with_fixture(sample_data):
    """Test using fixture data."""
    assert sample_data["key"] == "value"
```

### Adding Test Markers

Register new markers in `pytest.ini`:

```ini
markers =
    slow: marks tests as slow (deselect with '-m "not slow"')
    integration: marks tests as integration tests
    unit: marks tests as unit tests
    network: marks tests that require network connectivity
    your_marker: description of your custom marker
```

## Performance

### Running Tests in Parallel

Install pytest-xdist for parallel execution:

```bash
pip install pytest-xdist

# Run tests in parallel
pytest -n auto  # Use all CPU cores
pytest -n 4     # Use 4 processes
```

### Test Selection for Development

```bash
# Quick feedback loop - run only fast unit tests
pytest -m "unit and not slow" -v

# Run specific test file you're working on
pytest tests/test_simple.py -v
```

## Summary of Key Commands

| Purpose | Command |
|---------|---------|
| Run all tests | `pytest -v` |
| Run quick tests | `pytest -m "not (network or slow)" -v` |
| Run with coverage | `pytest --cov=genebe --cov-report=html` |
| Run specific file | `pytest tests/test_simple.py -v` |
| Debug failing test | `pytest --pdb -x` |
| Use test script | `./test.sh` |

For more pytest options, run `pytest --help` or visit the [pytest documentation](https://docs.pytest.org/).
