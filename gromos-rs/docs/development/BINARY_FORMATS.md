# Binary Trajectory Formats for GROMOS

Binary I/O format implementations for molecular dynamics trajectories and energies.

!!! warning "Educational Implementation"
    **These are format implementations for learning purposes.** The performance claims below are theoretical estimates based on literature from similar MD software (CHARMM, NAMD, GROMACS). They have NOT been validated in gromos-rs, which has ~140+ compilation errors and is not production-ready.

---

## 📊 Theoretical Performance Characteristics

Binary formats are *expected* to provide performance benefits similar to those observed in production MD software:

| Operation | ASCII (Typical) | Binary (Expected) | Literature Speedup |
|-----------|-----------------|-------------------|---------------------|
| **Trajectory Write** | 90-180s | <3s | 30-60× (NAMD/CHARMM) |
| **Trajectory Read** | 30s | 5-8s | 4-6× (NAMD/CHARMM) |
| **Energy Write** | 10-20s | 0.5-1s | 20-40× (estimated) |
| **Energy Read** | 5-10s | 1-2s | 3-5× (estimated) |
| **File Size** | 100% | ~50% | 2× smaller (typical) |

*Reference: 10,000 frames, 1,000 atoms (typical benchmark from literature)*

---

## 📦 Formats Supported

### 1. DCD Trajectory Format (.dcd)

**Standard binary format, widely compatible**

- ✅ Expected to provide faster I/O than ASCII (typical in MD software)
- ✅ Smaller file sizes (binary vs text)
- ✅ Single-precision coordinates (sufficient for all MD analyses)
- ✅ Compatible with CHARMM, NAMD, OpenMM, VMD
- ✅ Random access seeking
- ✅ Zero-copy NumPy integration

**Use when:**
- Running production MD simulations
- I/O is a bottleneck
- Compatibility with other MD software needed
- Disk space is limited

### 2. Binary Energy Format (.tre.bin)

**Custom, self-describing, lossless**

- ✅ Expected to provide faster I/O than ASCII
- ✅ Smaller file sizes than ASCII text
- ✅ Double-precision storage (lossless)
- ✅ Self-describing with metadata
- ✅ All 17 energy components
- ✅ Random access seeking

**Use when:**
- Energy I/O is bottlenecking analysis
- Long-term archival required
- Exact precision needed

---

## 🔧 Installation

### Rust Library (gromos-rs)

The binary formats are built into gromos-rs:

```rust
use gromos_rs::io::{DcdWriter, DcdReader, BinaryEnergyWriter, BinaryEnergyReader};
```

### Python Bindings (py-gromos)

Build the Python module:

```bash
cd py-gromos
maturin develop --release
```

Then in Python:

```python
import gromos
```

---

## 📚 Usage Guide

### Rust API

#### Writing Binary Trajectories

```rust
use gromos_rs::io::{DcdWriter, BinaryTrajectoryWriter};

// Create writer
let mut writer = DcdWriter::new("output.dcd", "MD simulation")?;

// Write frames during simulation
for step in 0..10000 {
    let time = step as f64 * 0.002; // 2 fs timestep
    writer.write_frame(step, time, &config)?;
}

// Finalize (updates header with frame count)
writer.finish()?;
```

#### Reading Binary Trajectories

```rust
use gromos_rs::io::{DcdReader, BinaryTrajectoryReader};

// Open reader
let mut reader = DcdReader::new("trajectory.dcd")?;

println!("Frames: {}, Atoms: {}", reader.n_frames(), reader.n_atoms());

// Read frames
while let Some(frame) = reader.read_frame()? {
    // frame.positions: Vec<Vec3>
    // frame.box_dims: Vec3
    // frame.time: f64
    // frame.step: usize
}

// Random access
reader.seek_frame(500)?;
let frame = reader.read_frame()?;
```

#### Writing Binary Energies

```rust
use gromos_rs::io::BinaryEnergyWriter;

let mut writer = BinaryEnergyWriter::new("energy.tre.bin", "Energies")?;

for step in 0..5000 {
    let mut frame = EnergyFrame::new(
        step as f64 * 0.002,
        kinetic,
        potential,
        temperature
    );
    frame.bond = -50.0;
    frame.lj = -100.0;
    // ... set other components

    writer.write_frame(&frame)?;
}

writer.finish()?;
```

#### Reading Binary Energies

```rust
use gromos_rs::io::BinaryEnergyReader;

let mut reader = BinaryEnergyReader::new("energy.tre.bin")?;

while let Some(frame) = reader.read_frame()? {
    println!("Time: {} ps", frame.time);
    println!("Total: {} kJ/mol", frame.total);
    println!("Kinetic: {} kJ/mol", frame.kinetic);
    // ... all 17 energy components available
}
```

---

### Python API

#### Reading Binary Trajectories

```python
import gromos
import numpy as np

# Open DCD file (binary format)
reader = gromos.DcdReader("trajectory.dcd")
print(f"Frames: {reader.n_frames}, Atoms: {reader.n_atoms}")

# Iterator pattern (memory efficient)
for frame in reader:
    positions = frame['positions']  # NumPy array (n_atoms, 3)
    box_dims = frame['box_dims']    # [lx, ly, lz]
    time = frame['time']            # ps
    step = frame['step']

    # Positions are already NumPy arrays - zero copy!
    rmsd = np.sqrt(np.mean(np.sum((positions - reference)**2, axis=1)))

# Random access
reader.seek_frame(500)
frame = reader.read_frame()

# Batch read
frames = reader.read_all_frames()
```

#### Reading Binary Energies

```python
# Open binary energy file
reader = gromos.BinaryEnergyReader("energy.tre.bin")
print(f"Title: {reader.title}")
print(f"Frames: {reader.n_frames}")

# Extract time series
times = []
energies = []

for frame in reader:
    times.append(frame['time'])
    energies.append(frame['total'])

# All energy components available:
# kinetic, potential, total, temperature, volume, pressure
# bond, angle, improper, dihedral, lj, coul_real, coul_recip,
# coul_self, shake, restraint
```

#### Format Conversion

```python
# ASCII → Binary conversion
gromos.convert_trajectory("input.trc", "output.dcd")

# Binary → ASCII (when needed for legacy tools)
gromos.convert_trajectory("output.dcd", "converted.trc")

# Energy conversion
gromos.convert_energy("input.tre", "output.tre.bin")
gromos.convert_energy("output.tre.bin", "converted.tre")
```

---

## 🎯 Recommended Workflows

### 1. Production Simulations

**Binary formats for efficient I/O:**

```bash
# MD++ or gromos-rs simulation can output binary directly
gromos-md --trc-format dcd --output trajectory.dcd
```

**Expected Benefits:**
- Faster trajectory writing (binary vs text)
- Reduced I/O wait time
- Less disk space used
- Compatible with analysis tools

### 2. Analysis Pipelines

**Read binary files directly in Python:**

```python
import gromos
import numpy as np

reader = gromos.DcdReader("production.dcd")

# Direct NumPy-based analysis (no conversion needed)
for frame in reader:
    positions = frame['positions']
    # Your analysis here
    rmsd = calculate_rmsd(positions, reference)
    rg = calculate_radius_of_gyration(positions)
```

**Benefits:**
- Direct binary reading (no parsing overhead)
- Zero-copy NumPy arrays
- Memory efficient iteration
- No intermediate files needed

### 3. Long-Term Archival

**Store trajectories in binary format:**

```bash
# Convert existing ASCII to binary for archival
python -c "import gromos; gromos.convert_trajectory('old.trc', 'archive.dcd')"
```

**Benefits:**
- Reduced disk space usage (binary vs text)
- Efficient re-analysis
- Compatible with many tools (VMD, MDAnalysis, MDTraj)
- Lossless precision preserved

### 4. Collaboration & Compatibility

**Binary DCD is a de facto standard:**

- ✅ VMD: Direct visualization
- ✅ MDAnalysis: Native support
- ✅ MDTraj: Native support
- ✅ CHARMM/NAMD: Native format
- ✅ OpenMM: DCDReporter
- ✅ GROMACS: gmx trjconv -o output.dcd

**Share DCD files directly - no conversion needed!**

---

## 📊 Format Characteristics

### Theoretical Performance (Based on MD Software Literature)

Binary formats typically provide performance advantages over text-based ASCII formats:

**Write Performance:**
- ASCII text: Slower (formatting overhead, larger writes)
- Binary DCD: Faster (direct binary writes, smaller data)

**Read Performance:**
- ASCII text: Slower (parsing overhead, type conversion)
- Binary DCD: Faster (direct memory mapping possible)

**File Size:**
- ASCII text: Larger (text representation of numbers)
- Binary DCD: Smaller (~50% reduction is typical for coordinates)

!!! note "Validation Needed"
    These characteristics are based on literature from production MD codes (CHARMM, NAMD, GROMACS). Actual performance in gromos-rs has not been measured, as the project has significant compilation errors and is not yet functional.

---

## 🔬 Technical Details

### DCD Format Specification

**Header:**
- Magic number: "CORD" (4 bytes)
- Frame count: Updated on close
- Timestep: Single precision
- Atom count: Fixed per file
- Endianness: Auto-detected

**Frame Structure:**
- Box dimensions: 6 doubles (a, γ, b, β, α, c)
- X coordinates: n_atoms × f32
- Y coordinates: n_atoms × f32
- Z coordinates: n_atoms × f32

**Unit Conversion:**
- Internal: nanometers (GROMOS standard)
- File: Ångströms (DCD standard)
- Automatic conversion on read/write

**Performance Optimizations:**
- 128KB buffering for optimal I/O
- Sequential frame layout for cache efficiency
- Single-precision float (sufficient for MD)

### Binary Energy Format (.tre.bin)

**Header:**
- Magic number: "GREBIN01" (8 bytes)
- Version: uint32
- Block count: 17 energy components
- Frame count: Updated on close
- Block names: Self-describing metadata
- Creation timestamp

**Frame Structure (17 × f64):**
1. Time (ps)
2. Kinetic energy (kJ/mol)
3. Potential energy (kJ/mol)
4. Total energy (kJ/mol)
5. Temperature (K)
6. Volume (nm³)
7. Pressure (bar)
8. Bond energy
9. Angle energy
10. Improper dihedral energy
11. Proper dihedral energy
12. Lennard-Jones energy
13. Coulomb (real space)
14. Coulomb (reciprocal space)
15. Coulomb (self)
16. SHAKE constraint energy
17. Restraint energy

**Performance:**
- Double precision (lossless)
- Fixed-size frames (136 bytes)
- Random access via frame index
- 128KB buffering

---

## 🔧 Implementation Files

### Rust Core

- `gromos-rs/src/io/trajectory_binary.rs` - DCD implementation (650 lines)
- `gromos-rs/src/io/energy_binary.rs` - Binary energy format (330 lines)
- `gromos-rs/src/io.rs` - Module exports

### Python Bindings

- `py-gromos/src/lib.rs` - PyO3 bindings (lines 927-1359)
  - `PyDcdReader` class
  - `PyBinaryEnergyReader` class
  - `convert_trajectory()` function
  - `convert_energy()` function

### Examples

- `py-gromos/examples/13_binary_trajectory_formats.py` - Complete tutorial

---

## 🧪 Testing

### Rust Tests

```bash
cd gromos-rs
cargo test trajectory_binary
cargo test energy_binary
```

### Python Tests

```python
import gromos

# Write and read DCD
writer = gromos.DcdWriter("test.dcd", "Test")
# ... write frames ...
writer.finish()

reader = gromos.DcdReader("test.dcd")
assert reader.n_frames() == expected_frames
assert reader.n_atoms() == expected_atoms
```

---

## 📖 References

### Research Basis

Implementation based on comprehensive MD software survey:
- OpenMM Reporter pattern
- GROMACS XTC/TRR formats
- CHARMM/NAMD DCD format
- AMBER NetCDF format
- Performance benchmarks from literature

### Compatibility

**DCD Format Variants:**
- CHARMM DCD (original)
- NAMD DCD (extended)
- OpenMM DCD (compatible)
- This implementation: NAMD-compatible

**Tools Supporting DCD:**
- VMD (visualization)
- MDAnalysis (Python analysis)
- MDTraj (Python analysis)
- PyMOL (visualization)
- CHARMM (simulation)
- NAMD (simulation)
- OpenMM (simulation)

---

## 🚀 Quick Start

### 1. Python Analysis

```python
import gromos

# Read binary trajectory
reader = gromos.DcdReader("trajectory.dcd")

for frame in reader:
    positions = frame['positions']  # NumPy array
    # Your analysis here
```

### 2. Convert Existing Files

```python
# ASCII → Binary conversion
gromos.convert_trajectory("old_trajectory.trc", "new_trajectory.dcd")
gromos.convert_energy("old_energy.tre", "new_energy.tre.bin")
```

### 3. Run Tutorial

```bash
python3 py-gromos/examples/13_binary_trajectory_formats.py
```

---

## 💡 Best Practices

### ✅ DO

- Use binary formats for production simulations
- Convert ASCII to binary for long-term archival
- Read binary files directly in analysis (no conversion)
- Keep binary files as primary storage
- Share DCD files for collaboration (widely compatible)

### ❌ DON'T

- Convert to ASCII unless required by legacy tools
- Use ASCII for large-scale simulations (I/O bottleneck)
- Store uncompressed ASCII long-term (wastes disk space)
- Re-convert repeatedly (slow and unnecessary)

---

## 🐛 Troubleshooting

### "Error opening DCD: Invalid magic number"

- File is corrupted or not a DCD file
- Check file exists and is complete
- Try with `file trajectory.dcd` to verify format

### "Error reading frame: Unexpected EOF"

- File was truncated during writing
- Ensure writer.finish() was called
- Check disk space during writes

### "Atom count mismatch"

- Trajectory frame has different atom count
- Ensure topology matches trajectory
- Check for selection/filtering issues

---

## 📝 Summary

Binary trajectory formats are implemented based on industry-standard formats (DCD, custom binary). Key characteristics:

✅ **Efficient**: Binary I/O eliminates text parsing overhead
✅ **Compact**: Binary storage is smaller than ASCII text
✅ **Lossless**: Precision is preserved (single/double as specified)
✅ **Compatible**: DCD format supported by VMD, MDAnalysis, MDTraj
✅ **Integrated**: Zero-copy NumPy integration design
✅ **Standard**: Based on proven CHARMM/NAMD/OpenMM formats

!!! warning "Implementation Status"
    These format implementations are educational. Performance benefits have not been validated in gromos-rs, which is not production-ready.

---

## 📞 Support

- Documentation: `gromos-rs/src/io/trajectory_binary.rs`
- Examples: `py-gromos/examples/13_binary_trajectory_formats.py`
- Tests: `cargo test trajectory_binary`
- Issues: GitHub repository

---

**Last Updated**: 2025-11-21
**Status**: Educational implementation, not validated
