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

# Clippy across the workspace, then the G6 drift gates (PLAN.md 3.9): one builder, one IMD
# reader, no process::exit in the library. Every grep must print nothing.
lint:
    #!/usr/bin/env bash
    set -uo pipefail
    cargo clippy --workspace --all-targets || exit 1
    fail=0
    gate() {  # gate <label> <grep args…> — a hit is a violation
        local label=$1; shift
        local out; out=$(grep -rn --include='*.rs' "$@" 2>/dev/null)
        if [ -n "$out" ]; then echo "G6 violation — $label:" >&2; echo "$out" >&2; fail=1; fi
    }
    # (1) no second builder: only gromos-run assembles an AlgorithmSequence (definition and test code excluded)
    gate "AlgorithmSequence::new() outside gromos-run" 'AlgorithmSequence::new()' crates --exclude-dir=gromos-run --exclude-dir=tests --exclude=algorithm.rs
    # (2) the binding pushes no algorithm of its own
    gate ".push(Box::new( in pyo3-gromos" '\.push(Box::new(' crates/pyo3-gromos
    # (3) one IMD reader: gromos_run::read_imd / Recipe (GAMD/EDS passthrough blocks parse their own text)
    gate "IMD read outside gromos-run/gromos-io" -e 'read_imd_file(' -e 'parse_imd_str(' crates --exclude-dir=gromos-run --exclude-dir=gromos-io --exclude-dir=tests
    # (4) gromos-run never exits the process (doc comments excluded)
    if grep -rn --include='*.rs' 'process::exit' crates/gromos-run/src | grep -v ':[[:space:]]*//' ; then
        echo "G6 violation — process::exit in gromos-run" >&2; fail=1
    fi
    [ "$fail" -eq 0 ] && echo "G6 gates: clean"
    exit $fail

# Engine-vs-engine benchmark against gromosXX (see BENCHMARKING.md); needs `just build`
# and the gromosXX native build. Extra args are passed through, e.g.
#   just bench --system ch4_water_fep --repeats 3
#   just bench --ref-dir bench/systems --system ch4_water_7975 --steps 500 --threads 1,2,4,8 --engines rust,omp
bench *ARGS:
    python3 scripts/bench_engines.py --system water_216_box --system water_1000_spc_gridcell --repeats 5 {{ARGS}}
