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
    fn new(
        topology: &Bound<'_, PyTopology>,
        configuration: &Bound<'_, PyConfiguration>,
    ) -> PyResult<Self> {
        let topo_ref = topology.borrow();
        let conf_ref = configuration.borrow();
        validate_atom_count_match(&topo_ref.inner, conf_ref.pos_data.len())
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
        validate_atom_count_match(&topo.inner, conf.pos_data.len())
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))?;
        Ok(Self {
            topology: topo,
            configuration: conf,
        })
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

        write_g96(
            path,
            "gromos-rs",
            pos,
            vels,
            box_opt,
            Some(&self.topology.inner),
        )
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
/// The configuration must hold the topology's atoms — or, for a topology that is not yet
/// solvated, the solute plus a whole number of solvent molecules: `prepare_system` solvates
/// from the coordinate count (the `md` binary's rule; the recipe's NSM is only a hint), so a
/// `System` of a solute topology and a solvated coordinate file is a valid one.
pub(crate) fn validate_atom_count_match(
    topo: &gromos_core::Topology,
    n_conf: usize,
) -> Result<(), String> {
    atom_count_ok(
        topo.num_atoms(),
        topo.num_solute_atoms(),
        topo.solvent_atom_template.len(),
        n_conf,
    )
}

fn atom_count_ok(
    n_topo: usize,
    n_solute: usize,
    per_solvent: usize,
    n_conf: usize,
) -> Result<(), String> {
    if n_topo == n_conf {
        return Ok(());
    }
    let unsolvated = n_topo == n_solute && per_solvent > 0;
    if unsolvated && n_conf >= n_solute && (n_conf - n_solute) % per_solvent == 0 {
        return Ok(());
    }
    Err(if unsolvated {
        format!(
            "Topology has {n_solute} solute atoms and {per_solvent}-atom solvent molecules, but \
             the configuration has {n_conf} atoms — not solute plus whole solvent molecules"
        )
    } else {
        format!("Topology has {n_topo} atoms but configuration has {n_conf} — they must match")
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atom_count_match_ok() {
        assert!(atom_count_ok(648, 0, 3, 648).is_ok());
    }

    #[test]
    fn atom_count_mismatch_error() {
        let err = atom_count_ok(648, 0, 3, 3).unwrap_err();
        assert!(err.contains("must match"), "message: {err}");
    }

    #[test]
    fn unsolvated_topology_accepts_solute_plus_whole_solvent_molecules() {
        // nacl_1water_box: 2 solute atoms, 3-atom SPC solvent, 5 atoms in the coordinates.
        assert!(atom_count_ok(2, 2, 3, 5).is_ok());
        assert!(atom_count_ok(2, 2, 3, 2).is_ok());
        let err = atom_count_ok(2, 2, 3, 6).unwrap_err();
        assert!(err.contains("whole solvent molecules"), "message: {err}");
        // A solvated topology (all 648 atoms present) is still an exact match only.
        assert!(atom_count_ok(648, 0, 3, 651).is_err());
    }
}
