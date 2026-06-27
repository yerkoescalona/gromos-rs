# gromos-rs task runner — install just: cargo install just

# Run all tests
test:
    cargo test --workspace

# Run only the GROMOS reference tests (ground truth)
ref:
    cargo test -p gromos-md --test test_gromosXX_references

# Build the md binary (release)
build:
    cargo build --release --bin md

# Open API docs in browser
doc:
    cargo doc -p gromos-core -p gromos-forces -p gromos-integrators -p gromos-io --no-deps --open

# Serve the scientific book at http://localhost:3000
book:
    mdbook serve docs/

# Clippy across the workspace
lint:
    cargo clippy --workspace --all-targets
