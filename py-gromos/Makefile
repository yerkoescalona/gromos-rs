.DEFAULT_GOAL := help

PYTHONPATH=
SHELL=bash
ifeq ($(VENV),)
VENV := ../.venv
endif

ifeq ($(OS),Windows_NT)
	VENV_BIN=$(VENV)/Scripts
else
	VENV_BIN=$(VENV)/bin
endif

##@ General

.PHONY: help
help: ## Display this help
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[36m<target>\033[0m\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)

##@ Setup

.PHONY: .venv
.venv: ## Set up virtual environment and install requirements
	@$(MAKE) -s -C .. $@

.PHONY: requirements
requirements: .venv ## Install/refresh Python project requirements
	@$(MAKE) -s -C .. $@

##@ Building

.PHONY: build
build: .venv ## Compile and install Python extension for development
	@$(MAKE) -s -C .. build-python

.PHONY: build-release
build-release: .venv ## Compile and install with optimizations
	@$(MAKE) -s -C .. build-release

##@ Testing

.PHONY: test
test: .venv build ## Run Python tests
	$(VENV_BIN)/pytest tests/ -v

.PHONY: test-verbose
test-verbose: .venv build ## Run tests with verbose output
	$(VENV_BIN)/pytest tests/ -vv -s

.PHONY: test-fast
test-fast: ## Run tests without rebuilding
	$(VENV_BIN)/pytest tests/ -v

.PHONY: test-coverage
test-coverage: .venv build ## Run tests with coverage
	$(VENV_BIN)/pytest tests/ --cov=gromos --cov-report=html --cov-report=term

##@ Code Quality

.PHONY: fmt
fmt: .venv ## Format Python code
	$(VENV_BIN)/ruff check --fix python/ tests/ examples/
	$(VENV_BIN)/ruff format python/ tests/ examples/

.PHONY: lint
lint: .venv ## Lint Python code
	$(VENV_BIN)/ruff check python/ tests/ examples/

.PHONY: pre-commit
pre-commit: fmt lint test ## Run pre-commit checks

##@ Notebooks

.PHONY: notebooks
notebooks: .venv build ## Run Jupyter notebooks
	$(VENV_BIN)/pip install jupyter notebook
	$(VENV_BIN)/jupyter notebook notebooks/

##@ Documentation

.PHONY: docs
docs: .venv ## Build static documentation site into site/
	$(VENV_BIN)/pip install --quiet mkdocs-material
	$(VENV_BIN)/mkdocs build --strict

.PHONY: docs-serve
docs-serve: .venv ## Serve documentation locally with live-reload (http://localhost:8000)
	$(VENV_BIN)/pip install --quiet mkdocs-material
	$(VENV_BIN)/mkdocs serve

##@ Cleaning

.PHONY: clean
clean: ## Clean build artifacts (preserves .venv)
	rm -rf .pytest_cache
	rm -rf .coverage htmlcov
	rm -rf .ruff_cache
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.so" -delete

.PHONY: clean-cache
clean-cache: ## Clean Python cache files only
	rm -rf .pytest_cache
	rm -rf .ruff_cache
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete

##@ Advanced

.PHONY: profile
profile: build-release ## Profile Python code
	$(VENV_BIN)/python -m cProfile -o profile.stats examples/basic_usage.py
	$(VENV_BIN)/python -m pstats profile.stats

.PHONY: type-check
type-check: venv ## Run type checking with mypy
	$(VENV_BIN)/pip install mypy
	$(VENV_BIN)/mypy python/ --ignore-missing-imports

.PHONY: publish-test
publish-test: build-wheel ## Publish to TestPyPI
	$(VENV_BIN)/pip install twine
	$(VENV_BIN)/twine upload --repository testpypi dist/*

.PHONY: publish
publish: build-wheel ## Publish to PyPI
	$(VENV_BIN)/pip install twine
	$(VENV_BIN)/twine upload dist/*

.PHONY: all
all: fmt clippy test build-release ## Run formatter, linter, tests and build release
