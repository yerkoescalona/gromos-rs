//! Python wrapper for GROMOS simulation input parameters.

use pyo3::prelude::*;

use gromos_io::imd::{read_imd_file, ImdParameters};

/// Simulation parameters loaded from a GROMOS input file (.imd/.in).
///
/// Contains all MD control parameters: timestep, number of steps,
/// cutoffs, thermostat/barostat settings, constraints, etc.
///
/// # Example (Python)
///
/// ```python
/// params = InputParameters("run.imd")
/// print(f"dt={params.dt} ps, steps={params.nstlim}")
/// print(f"T={params.temperature} K, cutoff={params.cutoff} nm")
/// print(f"solvent molecules: {params.nsm}")
/// ```
#[pyclass(name = "InputParameters")]
pub struct PyInputParameters {
    pub(crate) inner: ImdParameters,
}

#[pymethods]
impl PyInputParameters {
    /// Load parameters from a GROMOS input file.
    ///
    /// Args:
    ///     input_file: Path to input file (.imd, .in)
    #[new]
    fn new(input_file: &str) -> PyResult<Self> {
        let imd = read_imd_file(input_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Failed to read input file '{}': {}",
                input_file, e
            ))
        })?;
        Ok(Self { inner: imd })
    }

    /// Time step in picoseconds.
    #[getter]
    fn dt(&self) -> f64 {
        self.inner.dt
    }

    /// Number of simulation steps.
    #[getter]
    fn nstlim(&self) -> usize {
        self.inner.nstlim
    }

    /// Number of solvent molecules.
    #[getter]
    fn nsm(&self) -> usize {
        self.inner.nsm
    }

    /// Long-range nonbonded cutoff (nm).
    #[getter]
    fn cutoff(&self) -> f64 {
        self.inner.rcutl
    }

    /// Short-range pairlist cutoff (nm).
    #[getter]
    fn rcutp(&self) -> f64 {
        self.inner.rcutp
    }

    /// Target temperature (K) of the first temperature bath.
    #[getter]
    fn temperature(&self) -> f64 {
        if !self.inner.temp_bath.is_empty() && !self.inner.temp_bath[0].temp0.is_empty() {
            self.inner.temp_bath[0].temp0[0]
        } else {
            0.0
        }
    }

    /// SHAKE constraint mode (1=none, 2=H-bonds, 3=all bonds).
    #[getter]
    fn ntc(&self) -> i32 {
        self.inner.ntc
    }

    /// Boundary condition type (0=vacuum, 1=rectangular).
    #[getter]
    fn ntb(&self) -> i32 {
        self.inner.ntb
    }

    /// Pairlist update frequency (steps).
    #[getter]
    fn nsnb(&self) -> usize {
        self.inner.nsnb
    }

    /// Trajectory write frequency (steps).
    #[getter]
    fn ntwx(&self) -> usize {
        self.inner.ntwx
    }

    /// Energy write frequency (steps).
    #[getter]
    fn ntwe(&self) -> usize {
        self.inner.ntwe
    }

    fn __repr__(&self) -> String {
        format!(
            "InputParameters(dt={}, nstlim={}, nsm={}, cutoff={:.3}, T={:.1})",
            self.inner.dt,
            self.inner.nstlim,
            self.inner.nsm,
            self.inner.rcutl,
            self.temperature(),
        )
    }
}

pub fn register_parameters(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyInputParameters>()?;
    Ok(())
}
