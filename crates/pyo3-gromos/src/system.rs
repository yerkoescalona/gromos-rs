//! PySystem — paired Topology + Configuration with atom-count validation.

use numpy::PyArray2;
use pyo3::prelude::*;

use gromos_io::g96::write_g96;

use crate::py_conf::PyConfiguration;
use crate::topology::PyTopology;

/// A molecular system: topology paired with a matching configuration.
///
/// Validates at construction that both objects describe the same number of
/// atoms, mirroring the gromosXX invariant that `.top` and `.cnf` must agree.
///
/// # Example (Python)
///
/// ```python
/// system = System.from_files("system.topo", "initial.cnf")
/// print(system.n_atoms)          # 648
/// print(system.charge)           # 0
/// print(system.positions.shape)  # (648, 3)
/// ```
#[pyclass(name = "System")]
#[derive(Debug)]
pub struct PySystem {
    pub(crate) topology: PyTopology,
    pub(crate) configuration: PyConfiguration,
}

#[pymethods]
impl PySystem {
    /// Create a System from pre-loaded Topology and Configuration objects.
    ///
    /// Raises ValueError if atom counts do not match.
    #[new]
    fn new(topology: &Bound<'_, PyTopology>, configuration: &Bound<'_, PyConfiguration>) -> PyResult<Self> {
        let topo_ref = topology.borrow();
        let conf_ref = configuration.borrow();
        validate_atom_count_match(topo_ref.inner.num_atoms(), conf_ref.pos_data.len())
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))?;
        Ok(Self {
            topology: PyTopology {
                inner: topo_ref.inner.clone(),
                physical_constants: topo_ref.physical_constants,
            },
            configuration: PyConfiguration {
                pos_data: conf_ref.pos_data.clone(),
                vel_data: conf_ref.vel_data.clone(),
                box_dims: conf_ref.box_dims,
            },
        })
    }

    /// Load a System directly from topology and coordinate file paths.
    #[staticmethod]
    fn from_files(topo_file: &str, conf_file: &str) -> PyResult<Self> {
        let topo = PyTopology::from_file(topo_file)?;
        let conf = PyConfiguration::from_file(conf_file)?;
        validate_atom_count_match(topo.inner.num_atoms(), conf.pos_data.len())
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))?;
        Ok(Self { topology: topo, configuration: conf })
    }

    // ── Properties ────────────────────────────────────────────────────────────

    /// Total number of atoms.
    #[getter]
    fn n_atoms(&self) -> usize {
        self.topology.inner.num_atoms()
    }

    /// Total integer charge of the system (e).
    #[getter]
    fn charge(&self) -> i32 {
        self.topology.charge()
    }

    /// Positions as an Nx3 numpy array (nm).
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        self.configuration.positions(py)
    }

    /// Velocities as an Nx3 numpy array (nm/ps).
    #[getter]
    fn velocities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        self.configuration.velocities(py)
    }

    /// Simulation box dimensions (Lx, Ly, Lz) in nm.
    #[getter]
    fn r#box(&self) -> (f64, f64, f64) {
        self.configuration.box_dimensions()
    }

    /// The underlying Topology object.
    #[getter]
    fn topology(&self) -> PyTopology {
        PyTopology {
            inner: self.topology.inner.clone(),
            physical_constants: self.topology.physical_constants,
        }
    }

    /// The underlying Configuration object.
    #[getter]
    fn configuration(&self) -> PyConfiguration {
        PyConfiguration {
            pos_data: self.configuration.pos_data.clone(),
            vel_data: self.configuration.vel_data.clone(),
            box_dims: self.configuration.box_dims,
        }
    }

    /// Write the configuration to a GROMOS coordinate file (.cnf).
    fn write(&self, path: &str) -> PyResult<()> {
        let pos = &self.configuration.pos_data;
        let vel = &self.configuration.vel_data;
        let box_v = self.configuration.box_dims;

        let vels = if vel.iter().any(|v| *v != gromos_core::math::Vec3::ZERO) {
            Some(vel.as_slice())
        } else {
            None
        };
        let box_opt = if box_v.x > 0.0 { Some(box_v) } else { None };

        write_g96(path, "gromos-rs", pos, vels, box_opt, Some(&self.topology.inner))
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e))
    }

    fn __repr__(&self) -> String {
        let (lx, ly, lz) = self.configuration.box_dimensions();
        format!(
            "System(n_atoms={}, charge={}, box=({:.3}, {:.3}, {:.3}))",
            self.n_atoms(),
            self.charge(),
            lx,
            ly,
            lz,
        )
    }
}

pub fn register_system(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySystem>()?;
    Ok(())
}

/// Rust-only validation logic extracted for unit testing without pyo3 linker symbols.
pub(crate) fn validate_atom_count_match(n_topo: usize, n_conf: usize) -> Result<(), String> {
    if n_topo != n_conf {
        Err(format!(
            "Topology has {} atoms but configuration has {} — they must match",
            n_topo, n_conf
        ))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atom_count_match_ok() {
        assert!(validate_atom_count_match(648, 648).is_ok());
    }

    #[test]
    fn atom_count_mismatch_error() {
        let err = validate_atom_count_match(648, 3).unwrap_err();
        assert!(err.contains("must match"), "message: {err}");
    }
}
