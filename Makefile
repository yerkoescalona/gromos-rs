.DEFAULT_GOAL := help

PYTHONPATH=
SHELL=bash
ifeq ($(VENV),)
VENV := .venv
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

##@ Python Environment (Polars-style: .venv at repo root)

.venv: ## Set up Python virtual environment and install requirements
	@python3 -m venv $(VENV)
	@$(MAKE) requirements

.PHONY: requirements
requirements: ## Install/refresh Python project requirements
	@python3 -m venv $(VENV)
	@$(VENV_BIN)/pip install --upgrade pip
	@$(VENV_BIN)/pip install maturin numpy pytest matplotlib ruff

##@ Building

.PHONY: build
build: .venv build-rust build-python ## Build Rust binaries and Python extension

.PHONY: build-rust
build-rust: ## Build Rust workspace (all binaries)
	cargo build --release --bin md

.PHONY: build-python
build-python: .venv ## Compile and install Python bindings for development
	cd py-gromos && ../$(VENV_BIN)/maturin develop

.PHONY: build-release
build-release: .venv ## Build Python bindings with optimizations
	cd py-gromos && ../$(VENV_BIN)/maturin develop --release

##@ Testing

.PHONY: test
test: test-rust test-python ## Run all tests

.PHONY: test-rust
test-rust: ## Run Rust tests
	cargo test --workspace

.PHONY: test-references
test-references: ## Run gromosXX reference tests
	cargo test -p gromos-md --test test_gromosXX_references

.PHONY: test-python
test-python: build-python ## Run Python tests
	cd py-gromos && ../$(VENV_BIN)/pytest tests/ -v

.PHONY: test-python-fast
test-python-fast: ## Run Python tests without rebuilding
	cd py-gromos && ../$(VENV_BIN)/pytest tests/ -v

.PHONY: test-python-verbose
test-python-verbose: build-python ## Run Python tests with verbose output
	cd py-gromos && ../$(VENV_BIN)/pytest tests/ -vv -s

##@ Code Quality

.PHONY: fmt
fmt: ## Format Rust and Python code
	cargo fmt --all
	cd py-gromos && ../$(VENV_BIN)/ruff check --fix python/ tests/ examples/ || true
	cd py-gromos && ../$(VENV_BIN)/ruff format python/ tests/ examples/

.PHONY: lint
lint: ## Lint Rust and Python code
	cargo clippy --workspace
	cd py-gromos && ../$(VENV_BIN)/ruff check python/ tests/ examples/

.PHONY: check
check: ## Run cargo check on workspace
	cargo check --workspace

.PHONY: pre-commit
pre-commit: fmt lint test ## Run all pre-commit checks

##@ Documentation

.PHONY: doc
doc: ## Build Rust documentation
	cargo doc --workspace --no-deps --open

##@ Cleaning

.PHONY: clean
clean: ## Clean build artifacts (preserves .venv)
	cargo clean
	@rm -rf py-gromos/.pytest_cache
	@rm -rf py-gromos/.ruff_cache
	@find py-gromos -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true

.PHONY: clean-all
clean-all: clean ## Clean everything including .venv
	rm -rf $(VENV)

##@ Quick Commands

.PHONY: quick-test
quick-test: ## Quick test without full rebuild
	cargo test --lib
	cd py-gromos && ../$(VENV_BIN)/pytest tests/ -v

.PHONY: ci
ci: fmt lint test build-release ## Run full CI pipeline locally
