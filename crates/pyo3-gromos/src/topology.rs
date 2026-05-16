//! Python wrapper for GROMOS molecular topology.

use numpy::PyArray1;
use pyo3::prelude::*;

use gromos_core::Topology;
use gromos_io::topology::{build_topology, read_topology_file};

/// Molecular topology loaded from a GROMOS topology file.
///
/// Contains atom types, masses, charges, bonded parameters, LJ parameters,
/// exclusions, and molecule definitions.
///
/// # Example (Python)
///
/// ```python
/// topo = Topology("system.topo")
/// print(topo.n_atoms)          # solute atoms (before solvation)
/// print(topo.n_solute_atoms)
/// print(topo.masses[:5])       # first 5 masses
/// topo.solvate(216)            # add 216 solvent molecules
/// print(topo.n_atoms)          # now includes solvent
/// ```
#[pyclass(name = "Topology")]
pub struct PyTopology {
    pub(crate) inner: Topology,
}

#[pymethods]
impl PyTopology {
    /// Load a topology from a GROMOS topology file.
    ///
    /// Args:
    ///     topo_file: Path to topology file (.topo, .top)
    #[new]
    fn new(topo_file: &str) -> PyResult<Self> {
        let topo_data = read_topology_file(topo_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read topology '{}': {}",
                topo_file, e
            ))
        })?;
        let topo = build_topology(topo_data);
        Ok(Self { inner: topo })
    }

    /// Total number of atoms (solute + solvent).
    #[getter]
    fn n_atoms(&self) -> usize {
        self.inner.num_atoms()
    }

    /// Number of solute atoms.
    #[getter]
    fn n_solute_atoms(&self) -> usize {
        self.inner.num_solute_atoms()
    }

    /// Number of solvent atoms.
    #[getter]
    fn n_solvent_atoms(&self) -> usize {
        self.inner.num_atoms() - self.inner.num_solute_atoms()
    }

    /// Atomic masses as a numpy array.
    #[getter]
    fn masses<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec_bound(py, self.inner.mass.clone())
    }

    /// Atomic charges as a numpy array.
    #[getter]
    fn charges<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec_bound(py, self.inner.charge.clone())
    }

    /// Add solvent molecules to the topology.
    ///
    /// Args:
    ///     nsm: Number of solvent molecules to add
    fn solvate(&mut self, nsm: usize) {
        self.inner.solvate(nsm);
    }

    fn __repr__(&self) -> String {
        format!(
            "Topology(n_atoms={}, n_solute={}, n_solvent={})",
            self.inner.num_atoms(),
            self.inner.num_solute_atoms(),
            self.inner.num_atoms() - self.inner.num_solute_atoms(),
        )
    }

    fn __len__(&self) -> usize {
        self.inner.num_atoms()
    }
}

pub fn register_topology(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTopology>()?;
    Ok(())
}
