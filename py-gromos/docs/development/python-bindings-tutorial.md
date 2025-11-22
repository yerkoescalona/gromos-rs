# Python Bindings Architecture Tutorial

A comprehensive guide to understanding and extending the py-gromos Python bindings.

## Table of Contents

1. [Introduction](#introduction)
2. [PyO3 Fundamentals](#pyo3-fundamentals)
3. [Architecture Overview](#architecture-overview)
4. [Binding Construction](#binding-construction)
5. [Zero-Copy NumPy Integration](#zero-copy-numpy-integration)
6. [Extending the Bindings](#extending-the-bindings)
7. [Advanced Patterns](#advanced-patterns)
8. [Best Practices](#best-practices)

---

## Introduction

py-gromos uses [PyO3](https://pyo3.rs) to create Python bindings for the Rust-based gromos-rs library. This tutorial explains how these bindings work internally and how you can extend them.

### Why PyO3?

PyO3 provides:
- **Zero-overhead**: Direct FFI calls with minimal abstraction
- **Memory safety**: Rust's ownership system prevents memory leaks
- **Performance**: Native code execution without Python interpreter overhead
- **Ergonomic**: Clean APIs for both Rust and Python sides

### The Polars Inspiration

Our architecture is inspired by [Polars](https://github.com/pola-rs/polars):

```
Python Layer (py-gromos)
    ↓ PyO3 bindings (zero-copy when possible)
Rust Core (gromos-rs)
    ↓ Native computation (SIMD, parallelism)
Hardware
```

This separation provides:
- **High-level Python API**: Pythonic, familiar interface
- **Low-level Rust core**: Maximum performance, memory safety
- **Seamless integration**: Zero-copy data sharing via NumPy

---

## PyO3 Fundamentals

### The Module Definition

Every PyO3 module starts with a module function:

```rust
use pyo3::prelude::*;

#[pymodule]
fn gromos(_py: Python, m: &PyModule) -> PyResult<()> {
    // Register classes
    m.add_class::<PyVec3>()?;
    m.add_class::<PyState>()?;

    // Register functions
    m.add_function(wrap_pyfunction!(calculate_rmsd, m)?)?;

    Ok(())
}
```

**What happens here:**
1. `#[pymodule]` macro generates Python module initialization code
2. `m.add_class` registers Rust structs as Python classes
3. `m.add_function` exposes Rust functions to Python
4. The module name (`gromos`) becomes the Python import name

### Basic Type Wrapping

The simplest binding wraps a Rust type:

```rust
// Rust side
#[pyclass(name = "Vec3")]
#[derive(Clone)]
pub struct PyVec3 {
    inner: RustVec3,  // Wrap the native Rust type
}

#[pymethods]
impl PyVec3 {
    #[new]
    fn new(x: f32, y: f32, z: f32) -> Self {
        Self {
            inner: RustVec3::new(x, y, z),
        }
    }
}
```

**Python usage:**
```python
import gromos
v = gromos.Vec3(1.0, 2.0, 3.0)  # Calls Rust constructor
```

### Key PyO3 Attributes

| Attribute | Purpose | Example |
|-----------|---------|---------|
| `#[pyclass]` | Define Python class | `#[pyclass(name = "State")]` |
| `#[pymethods]` | Define methods block | `impl PyState { ... }` |
| `#[new]` | Constructor | `fn new() -> Self` |
| `#[getter]` | Property getter | `fn x(&self) -> f32` |
| `#[setter]` | Property setter | `fn set_x(&mut self, val: f32)` |
| `#[staticmethod]` | Static method | `fn from_numpy(...)` |
| `#[pyo3(name = "...")]` | Rename in Python | `#[pyo3(name = "from_array")]` |

---

## Architecture Overview

### Project Structure

```
py-gromos/
├── src/
│   └── lib.rs           # All bindings (PyO3 wrappers)
├── Cargo.toml           # Rust dependencies + PyO3
├── pyproject.toml       # Python package metadata
└── examples/            # Python usage examples
```

### Layered Design

```
┌─────────────────────────────────────────────────┐
│  Python API (gromos.*)                          │
│  - Pythonic interface                           │
│  - Type hints, docstrings                       │
│  - NumPy integration                            │
└─────────────────┬───────────────────────────────┘
                  │ PyO3 FFI layer
┌─────────────────┴───────────────────────────────┐
│  PyO3 Wrappers (PyVec3, PyState, etc.)          │
│  - Python class definitions                     │
│  - Type conversions                             │
│  - Error handling                               │
└─────────────────┬───────────────────────────────┘
                  │ Direct calls (zero overhead)
┌─────────────────┴───────────────────────────────┐
│  Rust Core (gromos-rs)                          │
│  - Native data structures (Vec3, State)         │
│  - Algorithms (MD integrators)                  │
│  - SIMD operations (via glam)                   │
└─────────────────────────────────────────────────┘
```

### Data Flow Example

**Python call:**
```python
v1 = gromos.Vec3(1.0, 2.0, 3.0)
v2 = gromos.Vec3(4.0, 5.0, 6.0)
distance = v1.distance(v2)  # ← This is the call
```

**What happens:**
1. Python invokes `v1.distance(v2)`
2. PyO3 FFI translates to Rust: `PyVec3::distance(&self, &PyVec3)`
3. Rust extracts inner values: `self.inner.distance(other.inner)`
4. Native Rust computation (SIMD-accelerated)
5. Result returned through FFI as Python `float`

**Performance:** Near-native speed! No Python interpreter in the hot path.

---

## Binding Construction

### Pattern 1: Simple Wrapper

For value types like vectors:

```rust
/// 3D vector with SIMD acceleration
#[pyclass(name = "Vec3")]
#[derive(Clone)]
pub struct PyVec3 {
    inner: RustVec3,  // Wraps gromos_rs::Vec3
}

#[pymethods]
impl PyVec3 {
    /// Constructor
    #[new]
    fn new(x: f32, y: f32, z: f32) -> Self {
        Self {
            inner: RustVec3::new(x, y, z),
        }
    }

    /// Method: compute length
    fn length(&self) -> f32 {
        self.inner.length()  // Delegate to Rust
    }

    /// Method: dot product
    fn dot(&self, other: &Self) -> f32 {
        self.inner.dot(other.inner)
    }

    /// Operator overload: addition
    fn __add__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner + other.inner,
        }
    }

    /// String representation
    fn __repr__(&self) -> String {
        format!("Vec3({:.4}, {:.4}, {:.4})",
                self.inner.x, self.inner.y, self.inner.z)
    }
}
```

**Python usage:**
```python
v1 = gromos.Vec3(1.0, 2.0, 3.0)
v2 = gromos.Vec3(4.0, 5.0, 6.0)
v3 = v1 + v2              # __add__
print(v1.length())        # length()
print(v1.dot(v2))         # dot()
print(v1)                 # __repr__
```

### Pattern 2: Properties (Getters/Setters)

For mutable state:

```rust
#[pymethods]
impl PyVec3 {
    /// Getter: read x component
    #[getter]
    fn x(&self) -> f32 {
        self.inner.x
    }

    /// Setter: write x component
    #[setter]
    fn set_x(&mut self, value: f32) {
        self.inner.x = value;
    }
}
```

**Python usage:**
```python
v = gromos.Vec3(1.0, 2.0, 3.0)
print(v.x)      # Calls getter → 1.0
v.x = 5.0       # Calls setter
print(v.x)      # → 5.0
```

### Pattern 3: Static Methods

For constructors and utilities:

```rust
#[pymethods]
impl PyVec3 {
    /// Create from NumPy array
    #[staticmethod]
    fn from_numpy(arr: PyReadonlyArray1<f32>) -> PyResult<Self> {
        let slice = arr.as_slice()?;
        if slice.len() != 3 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Array must have exactly 3 elements"
            ));
        }
        Ok(Self {
            inner: RustVec3::new(slice[0], slice[1], slice[2]),
        })
    }
}
```

**Python usage:**
```python
import numpy as np
arr = np.array([1.0, 2.0, 3.0], dtype=np.float32)
v = gromos.Vec3.from_numpy(arr)  # Static method
```

### Pattern 4: Complex Types with Lifetimes

For types that hold large amounts of data:

```rust
/// Molecular system state
#[pyclass(name = "State")]
pub struct PyState {
    inner: RustState,  // Contains Vec<Vec3> for positions, etc.
}

#[pymethods]
impl PyState {
    /// Constructor
    #[new]
    fn new(num_atoms: usize, num_temp_groups: usize,
           num_energy_groups: usize) -> Self {
        Self {
            inner: RustState::new(
                num_atoms,
                num_temp_groups,
                num_energy_groups
            ),
        }
    }

    /// Get number of atoms
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }

    /// Get positions as NumPy array (zero-copy)
    fn positions<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
        // IMPORTANT: Zero-copy view into Rust memory
        let data = self.inner.positions_slice();  // &[f32]
        let shape = (self.inner.num_atoms(), 3);

        // Create NumPy array that views Rust memory
        unsafe {
            PyArray2::from_ptr(
                py,
                data.as_ptr(),
                shape,
                [3 * std::mem::size_of::<f32>(), std::mem::size_of::<f32>()],
            )
        }
    }

    /// Set positions from NumPy array
    fn set_positions(&mut self, arr: PyReadonlyArray2<f32>) -> PyResult<()> {
        let arr = arr.as_array();
        if arr.shape() != [self.inner.num_atoms(), 3] {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Expected shape ({}, 3)", self.inner.num_atoms())
            ));
        }

        // Copy from NumPy to Rust
        self.inner.set_positions_from_slice(arr.as_slice().unwrap());
        Ok(())
    }
}
```

**Python usage:**
```python
state = gromos.State(num_atoms=1000, num_temp_groups=1, num_energy_groups=1)

# Zero-copy access (fast!)
positions = state.positions()  # NumPy array viewing Rust memory
print(positions.shape)  # (1000, 3)

# Modify in-place (modifies Rust memory directly)
positions[0, :] = [1.0, 2.0, 3.0]

# Or set new data (copies from NumPy to Rust)
new_pos = np.random.rand(1000, 3).astype(np.float32)
state.set_positions(new_pos)
```

---

## Zero-Copy NumPy Integration

### The Problem

Copying large arrays between Python and Rust is expensive:

```python
# BAD: Copies data twice
positions = state.get_positions()  # Copy Rust → NumPy
positions[0, 0] = 5.0
state.set_positions(positions)     # Copy NumPy → Rust
```

For 1M atoms × 3 coords × 4 bytes = 12 MB copied **twice** = 24 MB!

### The Solution: Zero-Copy Views

Share memory between Rust and NumPy:

```rust
fn positions<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
    let data = self.inner.positions_slice();  // &[f32] in Rust
    let shape = (self.inner.num_atoms(), 3);

    // Create NumPy array pointing to Rust memory
    unsafe {
        PyArray2::from_ptr(
            py,
            data.as_ptr(),      // Pointer to Rust data
            shape,              // Shape tuple
            strides,            // Memory layout
        )
    }
}
```

**Result:**
```python
positions = state.positions()     # Zero-copy! Just returns a view
positions[0, 0] = 5.0             # Modifies Rust memory directly
# No need to call set_positions!
```

### When to Use Zero-Copy

✅ **Use zero-copy when:**
- Data is large (> 1 KB)
- Data is frequently accessed
- Python needs read/write access
- Lifetime is clear (data won't be freed while Python uses it)

❌ **Don't use zero-copy when:**
- Data is small (< 100 bytes) → copying is faster
- Lifetime is complex → safety concerns
- Data needs transformation → copy anyway

### Memory Layout Considerations

For efficient zero-copy:

```rust
// GOOD: Contiguous memory
struct State {
    positions: Vec<f32>,  // [x0, y0, z0, x1, y1, z1, ...]
}

// BAD: Non-contiguous (Vec of Vecs)
struct State {
    positions: Vec<Vec3>,  // Each Vec3 might not be contiguous
}
```

**Why?** NumPy expects contiguous memory for optimal performance.

---

## Extending the Bindings

### Example: Adding a New Type

Let's add bindings for a `Force` calculator:

**Step 1: Define the Rust wrapper**

```rust
/// Force calculation for molecular systems
#[pyclass(name = "ForceCalculator")]
pub struct PyForceCalculator {
    inner: RustForceCalculator,
}

#[pymethods]
impl PyForceCalculator {
    /// Create a new force calculator
    #[new]
    fn new(cutoff: f32, use_pme: bool) -> Self {
        Self {
            inner: RustForceCalculator::new(cutoff, use_pme),
        }
    }

    /// Calculate forces for a state
    fn calculate(&self, state: &mut PyState, topology: &PyTopology) -> PyResult<()> {
        self.inner
            .calculate(&mut state.inner, &topology.inner)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                format!("Force calculation failed: {}", e)
            ))
    }

    /// Get cutoff distance
    #[getter]
    fn cutoff(&self) -> f32 {
        self.inner.cutoff()
    }

    /// Set cutoff distance
    #[setter]
    fn set_cutoff(&mut self, value: f32) {
        self.inner.set_cutoff(value);
    }
}
```

**Step 2: Register in the module**

```rust
#[pymodule]
fn gromos(_py: Python, m: &PyModule) -> PyResult<()> {
    // ... existing classes ...
    m.add_class::<PyForceCalculator>()?;  // ← Add this line
    Ok(())
}
```

**Step 3: Use from Python**

```python
import gromos

# Create calculator
calc = gromos.ForceCalculator(cutoff=1.2, use_pme=True)

# Calculate forces
state = gromos.State(num_atoms=1000, ...)
topology = gromos.Topology()
calc.calculate(state, topology)

# Access/modify properties
print(calc.cutoff)  # 1.2
calc.cutoff = 1.5
```

### Example: Adding a Module Function

For functions that don't belong to a class:

**Step 1: Define the function**

```rust
/// Calculate RMSD between two structures
#[pyfunction]
fn calculate_rmsd(
    positions1: PyReadonlyArray2<f32>,
    positions2: PyReadonlyArray2<f32>,
    masses: Option<PyReadonlyArray1<f32>>,
) -> PyResult<f32> {
    let pos1 = positions1.as_array();
    let pos2 = positions2.as_array();

    if pos1.shape() != pos2.shape() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Position arrays must have same shape"
        ));
    }

    let rmsd = if let Some(m) = masses {
        gromos_rs::analysis::rmsd_weighted(
            pos1.as_slice().unwrap(),
            pos2.as_slice().unwrap(),
            m.as_slice()?,
        )
    } else {
        gromos_rs::analysis::rmsd(
            pos1.as_slice().unwrap(),
            pos2.as_slice().unwrap(),
        )
    };

    Ok(rmsd)
}
```

**Step 2: Register in module**

```rust
#[pymodule]
fn gromos(_py: Python, m: &PyModule) -> PyResult<()> {
    // ... existing code ...
    m.add_function(wrap_pyfunction!(calculate_rmsd, m)?)?;
    Ok(())
}
```

**Step 3: Use from Python**

```python
import gromos
import numpy as np

pos1 = np.random.rand(100, 3).astype(np.float32)
pos2 = np.random.rand(100, 3).astype(np.float32)

# Without masses (unweighted)
rmsd = gromos.calculate_rmsd(pos1, pos2)

# With masses (weighted)
masses = np.ones(100, dtype=np.float32)
rmsd = gromos.calculate_rmsd(pos1, pos2, masses)
```

---

## Advanced Patterns

### Error Handling

Convert Rust errors to Python exceptions:

```rust
fn load_topology(&self, filename: &str) -> PyResult<PyTopology> {
    let topology = gromos_rs::io::read_topology(filename)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
            format!("Failed to load topology: {}", e)
        ))?;

    Ok(PyTopology { inner: topology })
}
```

**Python side:**
```python
try:
    topo = gromos.load_topology("missing.top")
except IOError as e:
    print(f"Error: {e}")
```

### Iterators

Expose Rust iterators to Python:

```rust
#[pymethods]
impl PyTrajectory {
    fn __iter__(slf: PyRef<Self>) -> PyResult<PyTrajectoryIterator> {
        Ok(PyTrajectoryIterator {
            inner: slf.inner.iter(),
        })
    }
}

#[pyclass]
struct PyTrajectoryIterator {
    inner: gromos_rs::TrajectoryIter,
}

#[pymethods]
impl PyTrajectoryIterator {
    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<Self>) -> Option<PyFrame> {
        slf.inner.next().map(|f| PyFrame { inner: f })
    }
}
```

**Python usage:**
```python
traj = gromos.Trajectory.load("traj.dcd")
for frame in traj:
    print(frame.time)
```

### Context Managers

For resource cleanup:

```rust
#[pymethods]
impl PyTrajectoryWriter {
    fn __enter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }

    fn __exit__(
        &mut self,
        _exc_type: Option<&PyAny>,
        _exc_value: Option<&PyAny>,
        _traceback: Option<&PyAny>,
    ) -> PyResult<bool> {
        self.inner.flush()?;
        Ok(false)
    }
}
```

**Python usage:**
```python
with gromos.TrajectoryWriter("output.dcd") as writer:
    writer.write_frame(frame)
# Automatically flushes on exit
```

### Threading and the GIL

Release the GIL for expensive operations:

```rust
use pyo3::Python;

#[pymethods]
impl PySimulation {
    fn run(&mut self, py: Python, num_steps: usize) -> PyResult<()> {
        // Release GIL during expensive computation
        py.allow_threads(|| {
            self.inner.run(num_steps)
        }).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
            format!("Simulation failed: {}", e)
        ))
    }
}
```

**Benefit:** Other Python threads can run while simulation is computing!

---

## Best Practices

### 1. Type Annotations

Always add Python type hints:

```python
from typing import Optional
import numpy as np
import numpy.typing as npt

def calculate_rmsd(
    positions1: npt.NDArray[np.float32],
    positions2: npt.NDArray[np.float32],
    masses: Optional[npt.NDArray[np.float32]] = None,
) -> float:
    """Calculate RMSD between two structures.

    Args:
        positions1: First structure (N, 3) array
        positions2: Second structure (N, 3) array
        masses: Optional atomic masses for weighting

    Returns:
        RMSD value in nm
    """
    ...
```

### 2. Error Messages

Make them actionable:

```rust
// BAD
return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Invalid"));

// GOOD
return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
    format!("Invalid array shape: expected ({}, 3), got {:?}",
            expected_atoms, actual_shape)
));
```

### 3. Docstrings

Document in Rust (appears in Python):

```rust
#[pymethods]
impl PyVec3 {
    /// Calculate dot product with another vector.
    ///
    /// The dot product is the sum of element-wise products:
    /// dot(v1, v2) = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
    ///
    /// Args:
    ///     other: Another Vec3 instance
    ///
    /// Returns:
    ///     Dot product as float
    ///
    /// Example:
    ///     >>> v1 = gromos.Vec3(1.0, 0.0, 0.0)
    ///     >>> v2 = gromos.Vec3(0.0, 1.0, 0.0)
    ///     >>> v1.dot(v2)
    ///     0.0
    fn dot(&self, other: &Self) -> f32 {
        self.inner.dot(other.inner)
    }
}
```

### 4. NumPy Dtype Validation

Always check dtypes:

```rust
fn set_positions(&mut self, arr: PyReadonlyArray2<f32>) -> PyResult<()> {
    // PyReadonlyArray2<f32> automatically validates dtype = float32
    let arr = arr.as_array();

    // But still check shape!
    if arr.shape()[0] != self.inner.num_atoms() || arr.shape()[1] != 3 {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("Expected shape ({}, 3), got {:?}",
                    self.inner.num_atoms(), arr.shape())
        ));
    }

    Ok(())
}
```

### 5. Avoid Unnecessary Copies

```rust
// BAD: Allocates temporary Vec
fn get_positions(&self) -> Vec<f32> {
    self.inner.positions().to_vec()
}

// GOOD: Returns view (zero-copy)
fn positions<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
    // Return view into existing memory
    ...
}
```

### 6. Version Compatibility

Specify minimum Python version:

```toml
# pyproject.toml
[project]
requires-python = ">=3.12"
```

### 7. Testing

Write tests for both Rust and Python:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vec3_creation() {
        let v = PyVec3::new(1.0, 2.0, 3.0);
        assert_eq!(v.length(), (1.0f32 + 4.0 + 9.0).sqrt());
    }
}
```

```python
# tests/test_vec3.py
import gromos
import numpy as np
import pytest

def test_vec3_creation():
    v = gromos.Vec3(1.0, 2.0, 3.0)
    assert abs(v.length() - np.sqrt(14)) < 1e-6

def test_vec3_numpy_conversion():
    arr = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    v = gromos.Vec3.from_numpy(arr)
    assert np.allclose(v.to_numpy(), arr)
```

---

## Summary

### Key Takeaways

1. **PyO3 provides zero-overhead bindings** between Rust and Python
2. **Wrapper pattern**: PyO3 types wrap Rust types (`PyVec3` wraps `Vec3`)
3. **Zero-copy is crucial** for performance with large arrays
4. **Follow Polars' architecture**: Rust core + thin Python wrapper
5. **Extension is straightforward**: Add Rust wrapper → Register → Use from Python

### Checklist for Adding Bindings

- [ ] Define Rust wrapper struct with `#[pyclass]`
- [ ] Implement methods with `#[pymethods]`
- [ ] Add constructor with `#[new]`
- [ ] Add getters/setters with `#[getter]`/`#[setter]`
- [ ] Implement `__repr__`, `__str__` for debugging
- [ ] Add NumPy conversions if applicable
- [ ] Register in `#[pymodule]` function
- [ ] Write docstrings
- [ ] Add Python type hints
- [ ] Write tests (Rust + Python)
- [ ] Update API documentation

### Resources

- **PyO3 Guide**: https://pyo3.rs
- **Polars Source**: https://github.com/pola-rs/polars
- **NumPy C API**: https://numpy.org/doc/stable/reference/c-api/
- **Maturin**: https://www.maturin.rs/

### Next Steps

1. Study the existing bindings in `src/lib.rs`
2. Try adding a simple new type (e.g., `Angle` or `Distance`)
3. Experiment with NumPy zero-copy patterns
4. Contribute improvements back to the project!

---

*This tutorial is part of the py-gromos educational project. For production-quality Python bindings, study the Polars and PyO3 documentation.*
