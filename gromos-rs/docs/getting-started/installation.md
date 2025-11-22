# Installation

This guide covers installing and building GROMOS-RS from source.

## Prerequisites

### Rust

GROMOS-RS requires Rust 1.70 or later. Install Rust using [rustup](https://rustup.rs):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Verify installation:

```bash
rustc --version  # Should be 1.70 or later
cargo --version
```

### Optional: C/C++ Compiler

For integration with md++ or gromos++, you'll need a C++ compiler:

=== "Linux"
    ```bash
    # Debian/Ubuntu
    sudo apt-get install build-essential cmake

    # Fedora/RHEL
    sudo dnf install gcc-c++ cmake

    # Arch
    sudo pacman -S base-devel cmake
    ```

=== "macOS"
    ```bash
    xcode-select --install
    brew install cmake
    ```

=== "Windows"
    Install [Visual Studio](https://visualstudio.microsoft.com/) with C++ support

## Building GROMOS-RS

### Clone Repository

```bash
git clone https://github.com/yerkoescalona/gromos-rs.git
cd gromosXX/gromos-rs
```

### Basic Build

```bash
# Debug build (faster compilation, slower runtime)
cargo build

# Release build (optimized, recommended for production)
cargo build --release
```

Binaries will be in:
- Debug: `target/debug/`
- Release: `target/release/`

### Build with Optimizations

For maximum performance:

```bash
# Native CPU optimizations (AVX2/AVX-512)
RUSTFLAGS="-C target-cpu=native" cargo build --release

# With mimalloc (faster allocator)
cargo build --release --features mimalloc

# Combined
RUSTFLAGS="-C target-cpu=native" cargo build --release --features mimalloc
```

!!! tip "Performance"
    Using `target-cpu=native` enables SIMD instructions for your specific CPU, providing 10-30% additional speedup.

### Profile-Guided Optimization (PGO)

For the absolute best performance:

```bash
# 1. Build with instrumentation
RUSTFLAGS="-Cprofile-generate=/tmp/pgo-data" cargo build --release

# 2. Run representative workload
cargo bench

# 3. Build optimized binary
RUSTFLAGS="-Cprofile-use=/tmp/pgo-data" cargo build --release
```

PGO can provide an additional 5-15% speedup.

## Verify Installation

### Run Tests

```bash
# All tests
cargo test --release

# Specific test
cargo test --release test_leap_frog

# With output
cargo test --release -- --nocapture
```

### Run Benchmarks

```bash
cargo bench
```

Results will be in `target/criterion/report/index.html`.

### Check Available Binaries

```bash
ls target/release/
```

You should see:
- `md` - Main MD simulation engine
- `remd` - Replica Exchange MD
- `eds` - Enveloping Distribution Sampling
- `gamd` - Gaussian Accelerated MD
- `make_pt_top` - FEP topology generator
- `ene_ana`, `rmsd`, `rmsf`, `rgyr`, etc. - Analysis tools

## System-Wide Installation

### Copy Binaries

```bash
# Install to ~/.local/bin (add to PATH if needed)
cargo install --path .

# Or manually copy
sudo cp target/release/md /usr/local/bin/
sudo cp target/release/remd /usr/local/bin/
# ... other binaries
```

### Add to PATH

Add to your `~/.bashrc` or `~/.zshrc`:

```bash
export PATH="$HOME/.cargo/bin:$PATH"
export PATH="$HOME/gromosXX/gromos-rs/target/release:$PATH"
```

## Integration with GROMOS

### Build GROMOS++ (Analysis Tools)

GROMOS++ provides 104 analysis tools. Build separately:

```bash
cd ../gromosPlusPlus/gromos++
./Config.sh
./configure --prefix=/opt/gromos
make -j$(nproc)
sudo make install
```

Add to PATH:

```bash
export PATH="/opt/gromos/bin:$PATH"
```

### Build md++ (Optional)

If you need the original md++ engine:

```bash
cd ../md++
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

Binary will be at `build/md++`.

## Docker Installation

A Dockerfile is provided for containerized use:

```bash
cd gromosXX/gromos-rs
docker build -t gromos-rs .
docker run -it gromos-rs md --help
```

## Troubleshooting

### Linking Errors

If you encounter linking errors:

```bash
# Linux: Install linker
sudo apt-get install lld  # Debian/Ubuntu
sudo dnf install lld       # Fedora/RHEL

# Use mold (faster linker)
cargo install mold
RUSTFLAGS="-C link-arg=-fuse-ld=mold" cargo build --release
```

### Out of Memory During Build

```bash
# Reduce parallel jobs
cargo build --release -j 2
```

### SIMD Not Working

Verify CPU features:

```bash
# Linux
cat /proc/cpuinfo | grep -E "avx2|avx512"

# macOS
sysctl -a | grep machdep.cpu.features
```

If AVX2/AVX-512 not available, SIMD will fall back to SSE4.2 or scalar code.

## Next Steps

- [Quick Start Guide](quick-start.md) - Run your first simulation
- [User Guide](../user-guide/introduction.md) - Detailed usage
- [Reference](../reference/integrators.md) - Algorithm documentation

## Getting Help

- **Issues**: [GitHub Issues](https://github.com/yerkoescalona/gromos-rs/issues)
- **GROMOS Forum**: [www.gromos.net/forum](https://www.gromos.net/forum)
- **Email**: biomos@igc.phys.chem.ethz.ch
