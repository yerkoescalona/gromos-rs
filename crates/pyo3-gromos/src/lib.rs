//! PyO3 bindings for GROMOS-RS
//!
//! This crate provides the low-level PyO3 bindings that wrap the Rust
//! GROMOS-RS library for use from Python. The `py-gromos` crate builds
//! on top of this to provide the final Python extension module.
//!
//! # Architecture
//!
//! ```text
//! Python (gromos package)
//!     ↓
//! py-gromos (cdylib, extension module)
//!     ↓
//! pyo3-gromos (rlib, binding types)
//!     ↓
//! gromos-* crates (pure Rust)
//! ```

use pyo3::prelude::*;
use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use glam::Vec3A as Vec3;

// Re-export core types
pub use gromos_core::{Configuration, Topology, Energy, State};
pub use gromos_forces::{ForceEnergy, LJParameters};
pub use gromos_integrators::{LeapFrog, VelocityVerlet, StochasticDynamics};
pub use gromos_io::{EnergyFrame, EnergyWriter};

/// Python-wrapped 3D vector
#[pyclass(name = "Vec3")]
#[derive(Clone, Debug)]
pub struct PyVec3 {
    pub inner: Vec3,
}

#[pymethods]
impl PyVec3 {
    #[new]
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self {
            inner: Vec3::new(x, y, z),
        }
    }

    #[getter]
    pub fn x(&self) -> f32 {
        self.inner.x
    }

    #[getter]
    pub fn y(&self) -> f32 {
        self.inner.y
    }

    #[getter]
    pub fn z(&self) -> f32 {
        self.inner.z
    }

    pub fn length(&self) -> f32 {
        self.inner.length()
    }

    pub fn dot(&self, other: &PyVec3) -> f32 {
        self.inner.dot(other.inner)
    }

    pub fn cross(&self, other: &PyVec3) -> PyVec3 {
        PyVec3 {
            inner: self.inner.cross(other.inner),
        }
    }

    pub fn normalize(&self) -> PyVec3 {
        PyVec3 {
            inner: self.inner.normalize(),
        }
    }

    pub fn __repr__(&self) -> String {
        format!("Vec3({}, {}, {})", self.inner.x, self.inner.y, self.inner.z)
    }

    pub fn __add__(&self, other: &PyVec3) -> PyVec3 {
        PyVec3 {
            inner: self.inner + other.inner,
        }
    }

    pub fn __sub__(&self, other: &PyVec3) -> PyVec3 {
        PyVec3 {
            inner: self.inner - other.inner,
        }
    }

    pub fn __mul__(&self, scalar: f32) -> PyVec3 {
        PyVec3 {
            inner: self.inner * scalar,
        }
    }
}

/// Python-wrapped Energy container
#[pyclass(name = "Energy")]
#[derive(Clone, Debug)]
pub struct PyEnergy {
    pub total: f64,
    pub kinetic: f64,
    pub potential: f64,
    pub bond: f64,
    pub angle: f64,
    pub dihedral: f64,
    pub improper: f64,
    pub lj: f64,
    pub coulomb: f64,
}

#[pymethods]
impl PyEnergy {
    #[new]
    pub fn new() -> Self {
        Self {
            total: 0.0,
            kinetic: 0.0,
            potential: 0.0,
            bond: 0.0,
            angle: 0.0,
            dihedral: 0.0,
            improper: 0.0,
            lj: 0.0,
            coulomb: 0.0,
        }
    }

    #[getter]
    pub fn total(&self) -> f64 {
        self.total
    }

    #[getter]
    pub fn kinetic(&self) -> f64 {
        self.kinetic
    }

    #[getter]
    pub fn potential(&self) -> f64 {
        self.potential
    }

    pub fn __repr__(&self) -> String {
        format!(
            "Energy(total={:.4}, kinetic={:.4}, potential={:.4})",
            self.total, self.kinetic, self.potential
        )
    }
}

/// Python-wrapped simulation frame
#[pyclass(name = "Frame")]
#[derive(Clone)]
pub struct PyFrame {
    pub time: f64,
    pub step: u64,
    pub positions: Vec<Vec3>,
    pub velocities: Option<Vec<Vec3>>,
    pub box_vectors: Option<[Vec3; 3]>,
}

#[pymethods]
impl PyFrame {
    #[new]
    pub fn new(time: f64, step: u64) -> Self {
        Self {
            time,
            step,
            positions: Vec::new(),
            velocities: None,
            box_vectors: None,
        }
    }

    #[getter]
    pub fn time(&self) -> f64 {
        self.time
    }

    #[getter]
    pub fn step(&self) -> u64 {
        self.step
    }

    #[getter]
    pub fn n_atoms(&self) -> usize {
        self.positions.len()
    }

    /// Get positions as numpy array (N x 3)
    pub fn positions_array<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let n_atoms = self.positions.len();
        let data: Vec<f32> = self.positions.iter()
            .flat_map(|v| [v.x, v.y, v.z])
            .collect();
        
        Ok(PyArray2::from_vec2_bound(py, &data.chunks(3).map(|c| c.to_vec()).collect::<Vec<_>>())?)
    }

    pub fn __repr__(&self) -> String {
        format!("Frame(time={:.3}, step={}, n_atoms={})", self.time, self.step, self.positions.len())
    }
}

/// Calculate RMSD between two sets of positions
#[pyfunction]
pub fn rmsd<'py>(
    positions: PyReadonlyArray2<'py, f32>,
    reference: PyReadonlyArray2<'py, f32>,
) -> PyResult<f64> {
    let pos = positions.as_array();
    let ref_arr = reference.as_array();
    
    if pos.shape() != ref_arr.shape() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Position and reference arrays must have the same shape"
        ));
    }
    
    let n_atoms = pos.shape()[0];
    let mut sum_sq = 0.0f64;
    
    for i in 0..n_atoms {
        let dx = (pos[[i, 0]] - ref_arr[[i, 0]]) as f64;
        let dy = (pos[[i, 1]] - ref_arr[[i, 1]]) as f64;
        let dz = (pos[[i, 2]] - ref_arr[[i, 2]]) as f64;
        sum_sq += dx * dx + dy * dy + dz * dz;
    }
    
    Ok((sum_sq / n_atoms as f64).sqrt())
}

/// Calculate radial distribution function
#[pyfunction]
pub fn rdf<'py>(
    py: Python<'py>,
    positions: PyReadonlyArray2<'py, f32>,
    group1: Vec<usize>,
    group2: Vec<usize>,
    n_bins: usize,
    r_max: f64,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
    let pos = positions.as_array();
    let dr = r_max / n_bins as f64;
    
    let mut hist = vec![0u64; n_bins];
    
    // Calculate distances
    for &i in &group1 {
        for &j in &group2 {
            if i != j {
                let dx = (pos[[i, 0]] - pos[[j, 0]]) as f64;
                let dy = (pos[[i, 1]] - pos[[j, 1]]) as f64;
                let dz = (pos[[i, 2]] - pos[[j, 2]]) as f64;
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                
                let bin = (r / dr) as usize;
                if bin < n_bins {
                    hist[bin] += 1;
                }
            }
        }
    }
    
    // Normalize to g(r)
    let n1 = group1.len() as f64;
    let n2 = group2.len() as f64;
    let volume = (4.0 / 3.0) * std::f64::consts::PI * r_max.powi(3);
    let rho = n2 / volume;
    
    let mut r_values = Vec::with_capacity(n_bins);
    let mut g_values = Vec::with_capacity(n_bins);
    
    for bin in 0..n_bins {
        let r_inner = bin as f64 * dr;
        let r_outer = (bin + 1) as f64 * dr;
        let shell_volume = (4.0 / 3.0) * std::f64::consts::PI * (r_outer.powi(3) - r_inner.powi(3));
        
        let r = (r_inner + r_outer) / 2.0;
        let g = (hist[bin] as f64) / (n1 * rho * shell_volume);
        
        r_values.push(r);
        g_values.push(g);
    }
    
    Ok((
        PyArray1::from_vec_bound(py, r_values),
        PyArray1::from_vec_bound(py, g_values),
    ))
}

/// Register all bindings in a Python module
pub fn register_bindings(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyVec3>()?;
    m.add_class::<PyEnergy>()?;
    m.add_class::<PyFrame>()?;
    m.add_function(wrap_pyfunction!(rmsd, m)?)?;
    m.add_function(wrap_pyfunction!(rdf, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pyvec3() {
        let v = PyVec3::new(1.0, 2.0, 3.0);
        assert!((v.length() - 3.7416575).abs() < 1e-5);
    }
}
