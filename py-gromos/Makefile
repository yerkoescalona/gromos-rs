.DEFAULT_GOAL := help

VENV := .venv
PYTHON := python3
VENV_BIN := $(VENV)/bin
MATURIN := $(VENV_BIN)/maturin

# Detect CPU features for optimization
RUST_FLAGS :=
ifeq ($(LTS_CPU),1)
    RUST_FLAGS := RUSTFLAGS='-C target-cpu=x86-64'
else
    RUST_FLAGS := RUSTFLAGS='-C target-cpu=native'
endif

##@ General

.PHONY: help
help: ## Display this help
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[36m<target>\033[0m\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)

##@ Setup

.PHONY: venv
venv: $(VENV) ## Create virtual environment

$(VENV):
	$(PYTHON) -m venv $(VENV)
	$(VENV_BIN)/pip install --upgrade pip
	$(VENV_BIN)/pip install maturin ruff pytest

.PHONY: install-dev
install-dev: venv ## Install development dependencies
	$(VENV_BIN)/pip install -e .[dev]

.PHONY: requirements
requirements: venv ## Install/refresh Python requirements
	$(VENV_BIN)/pip install --upgrade pip maturin
	$(VENV_BIN)/pip install -e .[dev]

##@ Building

.PHONY: build
build: venv ## Build Python extension (debug mode, fast compilation)
	$(RUST_FLAGS) $(MATURIN) develop

.PHONY: build-release
build-release: venv ## Build optimized release binary (slow compilation)
	$(RUST_FLAGS) $(MATURIN) develop --release

.PHONY: build-nodebug-release
build-nodebug-release: venv ## Build release without debug symbols
	$(RUST_FLAGS) $(MATURIN) develop --release --strip

.PHONY: build-debug-release
build-debug-release: venv ## Build release with debug symbols
	$(RUST_FLAGS) CARGO_PROFILE_RELEASE_DEBUG=true $(MATURIN) develop --release

.PHONY: build-dist-release
build-dist-release: venv ## Build for distribution with maximum optimization
	$(RUST_FLAGS) CARGO_PROFILE_RELEASE_LTO=fat CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1 $(MATURIN) build --release

.PHONY: build-wheel
build-wheel: venv ## Build Python wheel
	$(RUST_FLAGS) $(MATURIN) build --release

.PHONY: install
install: venv ## Install the package in the virtual environment
	$(RUST_FLAGS) $(MATURIN) develop --release

##@ Testing

.PHONY: test
test: build ## Run Python tests
	$(VENV_BIN)/pytest tests/ -v

.PHONY: test-verbose
test-verbose: build ## Run tests with verbose output
	$(VENV_BIN)/pytest tests/ -vv -s

.PHONY: test-fast
test-fast: ## Run tests without rebuilding
	$(VENV_BIN)/pytest tests/ -v

.PHONY: test-benchmark
test-benchmark: build-release ## Run benchmark tests
	$(VENV_BIN)/pytest tests/ -v --benchmark-only

.PHONY: test-coverage
test-coverage: build ## Run tests with coverage
	$(VENV_BIN)/pytest tests/ --cov=gromos --cov-report=html --cov-report=term

.PHONY: test-parallel
test-parallel: build ## Run tests in parallel
	$(VENV_BIN)/pytest tests/ -v -n auto

##@ Code Quality

.PHONY: fmt
fmt: venv ## Format Python code with ruff
	$(VENV_BIN)/ruff check --fix python/ tests/ examples/
	$(VENV_BIN)/ruff format python/ tests/ examples/

.PHONY: fmt-check
fmt-check: venv ## Check Python code formatting with ruff
	$(VENV_BIN)/ruff check python/ tests/ examples/
	$(VENV_BIN)/ruff format --check python/ tests/ examples/

.PHONY: lint
lint: venv ## Lint Python code
	$(VENV_BIN)/ruff check python/ tests/ examples/

.PHONY: clippy
clippy: ## Run Rust clippy on bindings
	cargo clippy --all-targets -- -D warnings

.PHONY: pre-commit
pre-commit: fmt clippy test ## Run pre-commit checks

##@ Documentation

.PHONY: doc
doc: venv ## Build Python documentation
	$(VENV_BIN)/pip install sphinx sphinx-rtd-theme
	cd docs && $(VENV_BIN)/sphinx-build -b html . _build/html

.PHONY: doc-open
doc-open: doc ## Build and open documentation
	$(PYTHON) -m webbrowser docs/_build/html/index.html

.PHONY: notebooks
notebooks: build ## Run Jupyter notebooks
	$(VENV_BIN)/pip install jupyter notebook
	$(VENV_BIN)/jupyter notebook notebooks/

##@ Examples

.PHONY: run-examples
run-examples: build ## Run all Python examples
	@for example in examples/*.py; do \
		echo "Running $$example..."; \
		$(VENV_BIN)/python $$example || exit 1; \
	done

##@ Cleaning

.PHONY: clean
clean: ## Clean build artifacts
	rm -rf $(VENV)
	rm -rf target/
	rm -rf dist/
	rm -rf build/
	rm -rf *.egg-info
	rm -rf .pytest_cache
	rm -rf .coverage
	rm -rf htmlcov
	rm -rf .ruff_cache
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
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
