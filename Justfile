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

# Engine-vs-engine benchmark against gromosXX (see BENCHMARKING.md); needs `just build`
# and the gromosXX native build. Extra args are passed through, e.g.
#   just bench --system ch4_water_fep --repeats 3
#   just bench --ref-dir bench/systems --system ch4_water_7975 --steps 500 --threads 1,2,4,8 --engines rust,omp
bench *ARGS:
    python3 scripts/bench_engines.py --system water_216_box --system water_1000_spc_gridcell --repeats 5 {{ARGS}}
