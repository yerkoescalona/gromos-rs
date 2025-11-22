#!/bin/bash
# Run CI checks locally before pushing

set -e

echo "==> Running CI checks locally..."
echo ""

# Format check
echo "==> 1. Checking code formatting..."
make fmt-check || {
    echo "❌ Format check failed. Run 'make fmt' to fix."
    exit 1
}
echo "✓ Format check passed"
echo ""

# Lint
echo "==> 2. Running linters..."
make lint || {
    echo "❌ Lint failed. Fix the issues or run 'make fix' for auto-fixes."
    exit 1
}
echo "✓ Lint passed"
echo ""

# Tests
echo "==> 3. Running tests..."
make test || {
    echo "❌ Tests failed. Fix the failing tests."
    exit 1
}
echo "✓ Tests passed"
echo ""

# Build release
echo "==> 4. Building release..."
make build-release || {
    echo "❌ Release build failed."
    exit 1
}
echo "✓ Release build passed"
echo ""

echo "==================================="
echo "✓ All CI checks passed!"
echo "==================================="
echo ""
echo "You can safely push your changes."
echo ""
