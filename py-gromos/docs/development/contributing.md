# Contributing

## Setup

```bash
git clone https://github.com/yerkoescalona/gromos-rs.git
cd gromos-rs
make build-python
source .venv/bin/activate
```

## Workflow

Before committing, run the pre-commit checks defined in the root `Makefile`:

```bash
make check       # cargo fmt + cargo check --tests (fast, run every commit)
make test-python # build + full Python test suite
```

Before pushing:

```bash
make ci          # fmt-check + clippy + cargo test --workspace
```

## Where things live

| What you want to change | Where |
|------------------------|-------|
| Python API exposed to users | `crates/pyo3-gromos/src/` |
| Python module packaging | `py-gromos/python/gromos/__init__.py` |
| Type stubs | `py-gromos/python/gromos/gromos.pyi` |
| Python tests | `py-gromos/tests/` |
| Documentation | `py-gromos/docs/` |
| Rust MD engine | `crates/gromos-{core,forces,integrators,io}/` |
| Reference test data | `crates/gromos-md/tests/gromosXX_references/` |

## Adding a new Python-exposed feature

1. Implement in Rust under `crates/pyo3-gromos/src/`.
2. Register the class/function in `crates/pyo3-gromos/src/lib.rs`.
3. Re-export from `py-gromos/python/gromos/__init__.py`.
4. Add the type signature to `py-gromos/python/gromos/gromos.pyi`.
5. Write a test in `py-gromos/tests/`.
6. Document in `py-gromos/docs/api/reference.md`.

## Adding a reference test system

Reference data lives in `crates/gromos-md/tests/gromosXX_references/`. Each
system directory contains:

```
<system_name>/
    input.toml          ← paths + tolerances + reference values
    <system>.topo
    <system>.cnf
    <system>.in
    expected/
        energies.tre
        trajectory.trc
        forces.trf
        final.conf
```

See `crates/gromos-md/.claude/CONTEXT.md` for the full procedure.

## Code style

- **Rust**: `cargo fmt` (enforced in CI), `cargo clippy` with workspace-level lint flags.
- **Python**: `ruff check` + `ruff format` (`make fmt` inside `py-gromos/`).
- **Comments**: only when the *why* is non-obvious. No docstrings restating the signature.

## Documentation

```bash
cd py-gromos
make docs        # build static site → site/
make docs-serve  # live-reload at http://localhost:8000
```

Docs are written in Markdown under `py-gromos/docs/` and rendered with
[MkDocs Material](https://squidfunk.github.io/mkdocs-material/).
