.DEFAULT_GOAL := help

##@ General

.PHONY: help
help: ## Display this help
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[36m<target>\033[0m\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)

##@ Complete Workflow

.PHONY: all
all: fmt lint test build ## Run complete workflow: format, lint, test, and build both projects

.PHONY: pre-commit
pre-commit: all ## Run all pre-commit checks

##@ Building

.PHONY: build
build: build-rust build-python ## Build both Rust and Python projects

.PHONY: build-rust
build-rust: ## Build gromos-rs library
	@echo "==> Building gromos-rs..."
	cd gromos-rs && $(MAKE) build

.PHONY: build-python
build-python: ## Build py-gromos Python bindings
	@echo "==> Building py-gromos..."
	cd py-gromos && $(MAKE) build

.PHONY: build-release
build-release: ## Build both projects in release mode
	@echo "==> Building gromos-rs (release)..."
	cd gromos-rs && $(MAKE) build-release
	@echo "==> Building py-gromos (release)..."
	cd py-gromos && $(MAKE) build-release

.PHONY: build-all-features
build-all-features: ## Build Rust with all features
	@echo "==> Building gromos-rs with all features..."
	cd gromos-rs && $(MAKE) build-all-features

##@ Testing

.PHONY: test
test: test-rust test-python ## Run all tests for both projects

.PHONY: test-rust
test-rust: ## Run gromos-rs tests
	@echo "==> Testing gromos-rs..."
	cd gromos-rs && $(MAKE) test

.PHONY: test-python
test-python: ## Run py-gromos tests
	@echo "==> Testing py-gromos..."
	cd py-gromos && $(MAKE) test

.PHONY: test-integration
test-integration: ## Run integration tests across both projects
	@echo "==> Running integration tests..."
	cd gromos-rs && $(MAKE) test-integration
	cd py-gromos && $(MAKE) test

.PHONY: test-all
test-all: ## Run all tests with all features
	@echo "==> Running all Rust tests with all features..."
	cd gromos-rs && $(MAKE) test-all-features
	@echo "==> Running all Python tests..."
	cd py-gromos && $(MAKE) test-verbose

.PHONY: test-coverage
test-coverage: ## Generate coverage reports for both projects
	@echo "==> Generating Rust coverage..."
	cd gromos-rs && $(MAKE) coverage || echo "Install cargo-llvm-cov for coverage"
	@echo "==> Generating Python coverage..."
	cd py-gromos && $(MAKE) test-coverage

##@ Code Quality

.PHONY: fmt
fmt: fmt-rust fmt-python ## Format code in both projects

.PHONY: fmt-rust
fmt-rust: ## Format Rust code
	@echo "==> Formatting gromos-rs..."
	cd gromos-rs && $(MAKE) fmt

.PHONY: fmt-python
fmt-python: ## Format Python code
	@echo "==> Formatting py-gromos..."
	cd py-gromos && $(MAKE) fmt

.PHONY: fmt-check
fmt-check: ## Check code formatting
	@echo "==> Checking gromos-rs formatting..."
	cd gromos-rs && $(MAKE) fmt-check
	@echo "==> Checking py-gromos formatting..."
	cd py-gromos && $(MAKE) fmt-check

.PHONY: lint
lint: lint-rust lint-python ## Lint code in both projects

.PHONY: lint-rust
lint-rust: ## Lint Rust code with clippy
	@echo "==> Linting gromos-rs..."
	cd gromos-rs && $(MAKE) clippy

.PHONY: lint-python
lint-python: ## Lint Python code
	@echo "==> Linting py-gromos..."
	cd py-gromos && $(MAKE) lint

.PHONY: check
check: ## Run cargo check on Rust code
	@echo "==> Checking gromos-rs..."
	cd gromos-rs && $(MAKE) check

##@ Documentation

.PHONY: doc
doc: doc-rust doc-python ## Build documentation for both projects

.PHONY: doc-rust
doc-rust: ## Build gromos-rs documentation
	@echo "==> Building gromos-rs docs..."
	cd gromos-rs && $(MAKE) doc

.PHONY: doc-python
doc-python: ## Build py-gromos documentation
	@echo "==> Building py-gromos docs..."
	cd py-gromos && $(MAKE) doc

.PHONY: doc-open
doc-open: ## Build and open Rust documentation
	cd gromos-rs && $(MAKE) doc-open

##@ Cleaning

.PHONY: clean
clean: clean-rust clean-python ## Clean all build artifacts

.PHONY: clean-rust
clean-rust: ## Clean gromos-rs artifacts
	@echo "==> Cleaning gromos-rs..."
	cd gromos-rs && $(MAKE) clean

.PHONY: clean-python
clean-python: ## Clean py-gromos artifacts
	@echo "==> Cleaning py-gromos..."
	cd py-gromos && $(MAKE) clean

.PHONY: clean-all
clean-all: clean ## Clean everything including root
	rm -rf target/
	rm -rf .pytest_cache/

##@ Advanced Features

.PHONY: bench
bench: ## Run benchmarks for both projects
	@echo "==> Running gromos-rs benchmarks..."
	cd gromos-rs && $(MAKE) bench

.PHONY: build-mpi
build-mpi: ## Build MPI-enabled binaries
	@echo "==> Building MPI binaries..."
	cd gromos-rs && $(MAKE) build-all-bins

.PHONY: test-mpi
test-mpi: ## Run MPI tests
	@echo "==> Testing MPI features..."
	cd gromos-rs && $(MAKE) test-mpi

.PHONY: audit
audit: ## Security audit of dependencies
	@echo "==> Auditing dependencies..."
	cd gromos-rs && $(MAKE) audit

##@ Quick Commands

.PHONY: quick-test
quick-test: ## Quick test without rebuilding everything
	cd gromos-rs && cargo test --lib
	cd py-gromos && $(MAKE) test-fast

.PHONY: quick-check
quick-check: fmt-check lint ## Quick format and lint check

##@ CI/CD

.PHONY: ci
ci: fmt-check lint test build-release ## Run full CI pipeline locally

.PHONY: ci-fast
ci-fast: fmt-check lint test ## Run CI without release build
