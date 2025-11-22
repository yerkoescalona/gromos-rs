.DEFAULT_GOAL := help

##@ General

.PHONY: help
help: ## Display this help
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[36m<target>\033[0m\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2 } /^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' $(MAKEFILE_LIST)

##@ Development

.PHONY: build
build: ## Build the Rust library (debug mode, fast compilation)
	cargo build

.PHONY: build-release
build-release: ## Build optimized release binary
	cargo build --release

.PHONY: build-all-features
build-all-features: ## Build with all features enabled
	cargo build --all-features

.PHONY: check
check: ## Run cargo check with all features
	cargo check --workspace --all-features --all-targets

.PHONY: clippy
clippy: ## Run clippy linter with all features
	cargo clippy --workspace --all-features --all-targets -- -D warnings -A clippy::upper_case_acronyms -A non_snake_case

.PHONY: clippy-default
clippy-default: ## Run clippy with default features only
	cargo clippy --workspace -- -D warnings

.PHONY: fmt
fmt: ## Format Rust code
	cargo fmt --all

.PHONY: fmt-check
fmt-check: ## Check Rust code formatting
	cargo fmt --all -- --check

##@ Testing

.PHONY: test
test: ## Run all tests
	cargo test --workspace

.PHONY: test-release
test-release: ## Run tests in release mode
	cargo test --workspace --release

.PHONY: test-all-features
test-all-features: ## Run tests with all features
	cargo test --workspace --all-features

.PHONY: test-mpi
test-mpi: ## Run MPI tests (requires MPI runtime)
	cargo test --workspace --features use-mpi

.PHONY: test-cuda
test-cuda: ## Run CUDA tests (requires CUDA runtime)
	cargo test --workspace --features use-cuda

.PHONY: test-integration
test-integration: ## Run integration tests only
	cargo test --test '*' --workspace

.PHONY: test-doc
test-doc: ## Run documentation tests
	cargo test --doc --workspace

.PHONY: test-verbose
test-verbose: ## Run tests with verbose output
	cargo test --workspace --verbose -- --nocapture

##@ Benchmarking

.PHONY: bench
bench: ## Run benchmarks
	cargo bench --workspace

.PHONY: bench-baseline
bench-baseline: ## Save benchmark baseline
	cargo bench --workspace -- --save-baseline main

.PHONY: bench-compare
bench-compare: ## Compare against baseline
	cargo bench --workspace -- --baseline main

##@ Binary Targets

.PHONY: build-md-mpi
build-md-mpi: ## Build MPI MD program
	cargo build --release --bin md_mpi --features use-mpi

.PHONY: build-repex-mpi
build-repex-mpi: ## Build MPI Replica Exchange MD
	cargo build --release --bin repex_mpi --features use-mpi

.PHONY: build-md-mpi-cuda
build-md-mpi-cuda: ## Build MPI+CUDA Hybrid MD
	cargo build --release --bin md_mpi_cuda --features use-mpi,use-cuda

.PHONY: build-all-bins
build-all-bins: build-md-mpi build-repex-mpi ## Build all binary targets

##@ Documentation

.PHONY: doc
doc: ## Build documentation
	cargo doc --workspace --all-features --no-deps

.PHONY: doc-open
doc-open: ## Build and open documentation
	cargo doc --workspace --all-features --no-deps --open

.PHONY: mkdocs-serve
mkdocs-serve: ## Serve MkDocs documentation locally
	mkdocs serve

.PHONY: mkdocs-build
mkdocs-build: ## Build MkDocs documentation
	mkdocs build

##@ Code Quality

.PHONY: pre-commit
pre-commit: fmt clippy test ## Run pre-commit checks (format, lint, test)

.PHONY: fix
fix: ## Auto-fix clippy warnings and format code
	cargo clippy --workspace --all-features --fix --allow-dirty --allow-staged
	cargo fmt --all

.PHONY: audit
audit: ## Security audit of dependencies
	cargo audit

.PHONY: outdated
outdated: ## Check for outdated dependencies
	cargo outdated

##@ Cleaning

.PHONY: clean
clean: ## Clean build artifacts
	cargo clean

.PHONY: clean-doc
clean-doc: ## Clean documentation
	rm -rf target/doc

##@ Advanced

.PHONY: bloat
bloat: ## Analyze binary size
	cargo bloat --release -n 20

.PHONY: expand
expand: ## Expand macros in main library
	cargo expand --lib

.PHONY: asm
asm: ## Show assembly for critical functions
	cargo asm --release --lib

.PHONY: coverage
coverage: ## Generate code coverage report (requires cargo-llvm-cov)
	cargo llvm-cov --workspace --html

.PHONY: all
all: fmt clippy test build-release ## Run formatter, linter, tests and build release
