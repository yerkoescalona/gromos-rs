//! Python wrapper for GROMOS configuration (system state).

use numpy::PyArray2;
use pyo3::prelude::*;

use gromos_core::math::Vec3;
use gromos_io::coordinate::read_coordinates;

/// System configuration: positions, velocities, and simulation box.
///
/// Loaded from a GROMOS coordinate file (.cnf/.g96). Represents the
/// initial state of the system before simulation.
///
/// # Example (Python)
///
/// ```python
/// conf = Configuration("initial.cnf")
/// print(conf.n_atoms)
/// print(conf.positions.shape)      # (N, 3) numpy array
/// print(conf.box_dimensions)       # (Lx, Ly, Lz) in nm
/// ```
#[pyclass(name = "Configuration")]
#[derive(Debug)]
pub struct PyConfiguration {
    pub(crate) pos_data: Vec<Vec3>,
    pub(crate) vel_data: Vec<Vec3>,
    pub(crate) box_dims: Vec3,
}

#[pymethods]
impl PyConfiguration {
    /// Load a configuration from a GROMOS coordinate file.
    ///
    /// Args:
    ///     conf_file: Path to coordinate file (.cnf, .g96)
    #[new]
    fn new(conf_file: &str) -> PyResult<Self> {
        Self::from_file(conf_file)
    }

    /// Load a configuration from a GROMOS coordinate file.
    #[staticmethod]
    pub fn from_file(conf_file: &str) -> PyResult<Self> {
        let coord_data = read_coordinates(conf_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read coordinates '{}': {}",
                conf_file, e
            ))
        })?;
        Ok(Self {
            pos_data: coord_data.positions,
            vel_data: coord_data.velocities,
            box_dims: coord_data.box_dims,
        })
    }

    /// Number of atoms with coordinates.
    #[getter]
    fn n_atoms(&self) -> usize {
        self.pos_data.len()
    }

    /// Positions as an Nx3 numpy array (nm).
    #[getter]
    pub fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let data: Vec<Vec<f64>> = self.pos_data.iter().map(|v| vec![v.x, v.y, v.z]).collect();
        Ok(PyArray2::from_vec2_bound(py, &data)?)
    }

    /// Velocities as an Nx3 numpy array (nm/ps).
    #[getter]
    pub fn velocities<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let data: Vec<Vec<f64>> = self.vel_data.iter().map(|v| vec![v.x, v.y, v.z]).collect();
        Ok(PyArray2::from_vec2_bound(py, &data)?)
    }

    /// Simulation box dimensions (Lx, Ly, Lz) in nm.
    #[getter]
    pub fn box_dimensions(&self) -> (f64, f64, f64) {
        (self.box_dims.x, self.box_dims.y, self.box_dims.z)
    }

    fn __repr__(&self) -> String {
        format!(
            "Configuration(n_atoms={}, box=({:.3}, {:.3}, {:.3}))",
            self.pos_data.len(),
            self.box_dims.x,
            self.box_dims.y,
            self.box_dims.z,
        )
    }
}

pub fn register_configuration(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyConfiguration>()?;
    Ok(())
}
