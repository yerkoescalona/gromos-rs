# gromos-rs — local CI equivalent
# Usage: make check        (fast, before commit)
#        make ci           (full suite, matches GitHub CI)
#        make build-python (compile + install Python extension into .venv)
#        make test-python  (run py-gromos test suite)

VENV := .venv
VENV_BIN := $(VENV)/bin

.PHONY: check ci fmt fmt-check lint-rust lint test \
        .venv requirements build-python build-release test-python \
        docs docs-serve

# ── Rust ─────────────────────────────────────────────────────────────────────

# Run before every commit (matches .claude/overview.md rule)
check:
	cargo fmt
	cargo check --tests

# Full CI equivalent — run before pushing
ci: fmt-check lint-rust test

fmt:
	cargo fmt

fmt-check:
	cargo fmt --check

lint-rust:
	cargo clippy --workspace --all-targets -- \
		-W warnings \
		-A clippy::upper_case_acronyms \
		-A non_snake_case

test:
	cargo test --workspace

# ── Python ───────────────────────────────────────────────────────────────────

# Create virtual environment and install Python dev dependencies
.venv:
	python3 -m venv $(VENV)
	$(VENV_BIN)/pip install --quiet --upgrade pip
	$(VENV_BIN)/pip install --quiet maturin
	$(VENV_BIN)/pip install --quiet -e "py-gromos/[dev]" --no-build-isolation 2>/dev/null || true

# Reinstall Python dependencies without rebuilding the extension
requirements: .venv
	$(VENV_BIN)/pip install --quiet maturin
	$(VENV_BIN)/pip install --quiet "py-gromos/[dev]" --no-build-isolation 2>/dev/null || true

# Compile the Rust extension and install it into .venv (dev/debug build)
build-python: .venv
	cd py-gromos && $(CURDIR)/$(VENV_BIN)/maturin develop --extras dev

# Compile with optimisations
build-release: .venv
	cd py-gromos && $(CURDIR)/$(VENV_BIN)/maturin develop --release --extras dev

# Run the Python test suite (build first)
test-python: build-python
	$(VENV_BIN)/pytest py-gromos/tests/ -v

# ── Docs ─────────────────────────────────────────────────────────────────────

# Build static documentation site into py-gromos/site/
docs: .venv
	$(VENV_BIN)/pip install --quiet mkdocs-material
	cd py-gromos && $(CURDIR)/$(VENV_BIN)/mkdocs build --strict

# Serve documentation locally with live-reload (http://localhost:8000)
docs-serve: .venv
	$(VENV_BIN)/pip install --quiet mkdocs-material
	cd py-gromos && $(CURDIR)/$(VENV_BIN)/mkdocs serve
