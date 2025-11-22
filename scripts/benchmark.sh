#!/bin/bash
# Run benchmarks and compare with baseline

set -e

BASELINE=${1:-main}

echo "==> Running benchmarks..."

# Run Rust benchmarks
echo "==> Running Rust benchmarks..."
cd gromos-rs
if [ "$BASELINE" != "main" ]; then
    echo "Comparing against baseline: $BASELINE"
    cargo bench --workspace -- --baseline "$BASELINE"
else
    echo "Saving baseline as: $BASELINE"
    cargo bench --workspace -- --save-baseline "$BASELINE"
fi
cd ..

# Run Python benchmarks
echo "==> Running Python benchmarks..."
cd py-gromos
if [ ! -d ".venv" ]; then
    echo "Virtual environment not found. Run 'make venv' first."
    exit 1
fi
.venv/bin/pip install pytest-benchmark 2>/dev/null || true
.venv/bin/pytest tests/ --benchmark-only || echo "No Python benchmarks available yet"
cd ..

echo ""
echo "==> Benchmarks complete!"
echo "Results saved in gromos-rs/target/criterion/"
