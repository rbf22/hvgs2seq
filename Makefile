.PHONY: help install install-dev test lint format check-format typecheck clean clean-build clean-pyc clean-test

# Determine the Python interpreter to use
PYTHON = python
PIP = pip
POETRY = poetry

# Help message
help:
	@echo "Please use 'make <target>' where <target> is one of"
	@echo "  install           install the package in development mode"
	@echo "  install-dev       install development dependencies"
	@echo "  test              run tests quickly with the default Python"
	@echo "  test-all          run tests on every Python version with tox"
	@echo "  lint              check code style with ruff"
	@echo "  format            format code with black and isort"
	@echo "  typecheck         run type checking with mypy"
	@echo "  clean             remove all build, test, coverage and Python artifacts"

# Install the package in development mode
install:
	$(POETRY) install --with dev

# Install development dependencies
install-dev:
	$(POETRY) install --with dev
	$(POETRY) run pre-commit install

# Run tests
TEST_ARGS = -xvs --cov=hgvs2seq --cov-report=term-missing

test:
	$(POETRY) run pytest $(TEST_ARGS)

test-all:
	$(POETRY) run tox

# Lint the code
lint:
	$(POETRY) run ruff check .

# Format the code
format:
	$(POETRY) run black .
	$(POETRY) run isort .

# Type checking
typecheck:
	$(POETRY) run mypy hgvs2seq tests

# Clean up
clean: clean-build clean-pyc clean-test

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache/
	rm -fr .mypy_cache/
	rm -fr .ruff_cache/
