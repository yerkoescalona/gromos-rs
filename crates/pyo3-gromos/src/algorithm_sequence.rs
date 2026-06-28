//! Compositional AlgorithmSequence API for Python.
//!
//! Allows users to build custom MD algorithm sequences with full control
//! over ordering, parameters, and composition. Algorithms are described as
//! lightweight Python objects and resolved into the real Rust algorithms
//! when a Simulation is created.
//!
//! # Example (Python)
//!
//! ```python
//! from gromos import (
//!     Topology, Configuration, InputParameters, Simulation,
//!     AlgorithmSequence, Forcefield, LeapFrogIntegrator,
//!     BerendsenThermostat, BerendsenBarostat, ShakeConstraints,
//!     TemperatureCalculation, PressureCalculation, EnergyCalculation,
//!     RemoveCOMMotion,
//! )
//!
//! # Build from parameters (standard preset)
//! topo = Topology("system.topo")
//! params = InputParameters("run.imd")
//! seq = AlgorithmSequence.from_parameters(topo, params)
//!
//! # Build custom sequence
//! seq = AlgorithmSequence()
//! seq.add(Forcefield(cutoff=0.8, epsilon_rf=61.0))
//! seq.add(LeapFrogIntegrator())
//! seq.add(BerendsenThermostat(temperature=300.0, tau=0.1))
//! seq.add(ShakeConstraints(tolerance=1e-4))
//! seq.add(TemperatureCalculation())
//! seq.add(EnergyCalculation())
//!
//! # Presets
//! seq = AlgorithmSequence.nve(topo, params)
//! seq = AlgorithmSequence.nvt(topo, params)
//! seq = AlgorithmSequence.npt(topo, params)
//!
//! # Modify
//! seq.insert_after("LeapFrogVelocity", BerendsenThermostat(temperature=300.0, tau=0.1))
//! seq.remove("PressureCalculation")
//!
//! # Use with Simulation
//! sim = Simulation(topo, conf, sequence=seq)
//! ```

use pyo3::prelude::*;

use gromos_core::algorithm::AlgorithmSequence;
use gromos_core::configuration::BoxType;
use gromos_core::configuration::Configuration;
use gromos_core::math::{Periodicity, Rectangular, Vacuum, Vec3};
use gromos_core::pairlist::{PairlistAlgorithm, PairlistContainer};
use gromos_core::topology::Topology;
use gromos_forces::nonbonded::CRFParameters;
use gromos_integrators::algorithms::{
    BerendsenBarostat as RustBarostat, BerendsenBarostatParams,
    BerendsenThermostat as RustThermostat, EnergyCalculation as RustEnergy,
    Forcefield as RustForcefield, LeapFrogPosition, LeapFrogVelocity,
    PressureCalculation as RustPressure, RemoveCOMMotion as RustRemoveCOM,
    ShakeAlgorithm as RustShake, TemperatureCalculation as RustTemperature, VirialType,
};
use gromos_integrators::constraints::{NtcMode, ShakeParameters};
use gromos_io::imd::ImdParameters;

use super::parameters::PyInputParameters;
use super::topology::PyTopology;

// ============================================================================
// Algorithm descriptor types (lightweight Python objects holding config)
// ============================================================================

/// Forcefield algorithm: computes bonded and nonbonded forces.
///
/// This is the main force calculation step in the MD sequence.
/// It handles pairlist construction, LJ interactions, Coulomb reaction field,
/// bond/angle/dihedral forces, and virial computation.
///
/// # Example (Python)
///
/// ```python
/// ff = Forcefield(cutoff=0.8, epsilon_rf=61.0)
/// ff = Forcefield(cutoff=1.4, rcutp=0.8, epsilon_rf=61.0, kappa=0.0)
/// ff = Forcefield()  # use defaults from InputParameters
/// ```
#[pyclass(name = "Forcefield")]
#[derive(Clone, Debug)]
pub struct PyForcefield {
    /// Long-range cutoff (nm). None = use from InputParameters.
    #[pyo3(get, set)]
    pub cutoff: Option<f64>,
    /// Short-range (pairlist) cutoff (nm). None = same as cutoff.
    #[pyo3(get, set)]
    pub rcutp: Option<f64>,
    /// Reaction field dielectric constant. None = use from InputParameters.
    #[pyo3(get, set)]
    pub epsilon_rf: Option<f64>,
    /// Inverse Debye screening length. None = use from InputParameters.
    #[pyo3(get, set)]
    pub kappa: Option<f64>,
    /// Pairlist update frequency (steps). None = use from InputParameters.
    #[pyo3(get, set)]
    pub pairlist_update: Option<usize>,
    /// Virial type: "none", "atomic", "molecular". None = from InputParameters.
    #[pyo3(get, set)]
    pub virial: Option<String>,
    /// NTF bond flag (compute bond forces).
    #[pyo3(get, set)]
    pub ntf_bond: Option<bool>,
    /// NTF angle flag (compute angle forces).
    #[pyo3(get, set)]
    pub ntf_angle: Option<bool>,
    /// NTF improper flag (compute improper dihedral forces).
    #[pyo3(get, set)]
    pub ntf_improper: Option<bool>,
    /// NTF dihedral flag (compute proper dihedral forces).
    #[pyo3(get, set)]
    pub ntf_dihedral: Option<bool>,
}

#[pymethods]
impl PyForcefield {
    #[new]
    #[pyo3(signature = (
        cutoff=None, rcutp=None, epsilon_rf=None, kappa=None,
        pairlist_update=None, virial=None,
        ntf_bond=None, ntf_angle=None, ntf_improper=None, ntf_dihedral=None,
    ))]
    fn new(
        cutoff: Option<f64>,
        rcutp: Option<f64>,
        epsilon_rf: Option<f64>,
        kappa: Option<f64>,
        pairlist_update: Option<usize>,
        virial: Option<String>,
        ntf_bond: Option<bool>,
        ntf_angle: Option<bool>,
        ntf_improper: Option<bool>,
        ntf_dihedral: Option<bool>,
    ) -> Self {
        Self {
            cutoff,
            rcutp,
            epsilon_rf,
            kappa,
            pairlist_update,
            virial,
            ntf_bond,
            ntf_angle,
            ntf_improper,
            ntf_dihedral,
        }
    }

    fn __repr__(&self) -> String {
        let mut parts = Vec::new();
        if let Some(c) = self.cutoff {
            parts.push(format!("cutoff={:.3}", c));
        }
        if let Some(e) = self.epsilon_rf {
            parts.push(format!("epsilon_rf={:.1}", e));
        }
        if let Some(v) = &self.virial {
            parts.push(format!("virial='{}'", v));
        }
        format!("Forcefield({})", parts.join(", "))
    }
}

/// Leap-Frog integrator (velocity + position update).
///
/// Standard GROMOS Leap-Frog scheme:
///   - velocity half-step: v(t+dt/2) = v(t-dt/2) + f(t)/m * dt
///   - position full step: r(t+dt) = r(t) + v(t+dt/2) * dt
///
/// **Important:** In the algorithm sequence, this expands to TWO internal steps:
/// `LeapFrogVelocity` (which also exchanges state) and `LeapFrogPosition`.
/// They are placed adjacently. If you need a thermostat BETWEEN them
/// (GROMOS convention), use the individual `LeapFrogVelocity` and
/// `LeapFrogPosition` classes instead.
///
/// # Example (Python)
///
/// ```python
/// # Simple (both steps together):
/// seq.add(LeapFrogIntegrator())
///
/// # Advanced (thermostat between velocity and position update):
/// seq.add(LeapFrogVelocity())
/// seq.add(BerendsenThermostat(temperature=300.0, tau=0.1))
/// seq.add(LeapFrogPosition())
/// ```
#[pyclass(name = "LeapFrogIntegrator")]
#[derive(Clone, Debug)]
pub struct PyLeapFrogIntegrator {}

#[pymethods]
impl PyLeapFrogIntegrator {
    #[new]
    fn new() -> Self {
        Self {}
    }

    fn __repr__(&self) -> String {
        "LeapFrogIntegrator()".to_string()
    }
}

/// Leap-Frog velocity half-step (advanced).
///
/// Performs state exchange and velocity update: v(t+dt/2) = v(t-dt/2) + f(t)/m * dt.
/// Use this with `LeapFrogPosition` when you need to place algorithms (e.g. thermostat)
/// between the velocity and position updates.
///
/// # Example (Python)
///
/// ```python
/// seq.add(LeapFrogVelocity())
/// seq.add(BerendsenThermostat(temperature=300.0, tau=0.1))
/// seq.add(LeapFrogPosition())
/// ```
#[pyclass(name = "LeapFrogVelocity")]
#[derive(Clone, Debug)]
pub struct PyLeapFrogVelocity {}

#[pymethods]
impl PyLeapFrogVelocity {
    #[new]
    fn new() -> Self {
        Self {}
    }

    fn __repr__(&self) -> String {
        "LeapFrogVelocity()".to_string()
    }
}

/// Leap-Frog position full step (advanced).
///
/// Performs position update: r(t+dt) = r(t) + v(t+dt/2) * dt.
/// Use together with `LeapFrogVelocity` for full control over algorithm placement.
///
/// # Example (Python)
///
/// ```python
/// seq.add(LeapFrogVelocity())
/// seq.add(BerendsenThermostat(temperature=300.0, tau=0.1))
/// seq.add(LeapFrogPosition())
/// ```
#[pyclass(name = "LeapFrogPosition")]
#[derive(Clone, Debug)]
pub struct PyLeapFrogPosition {}

#[pymethods]
impl PyLeapFrogPosition {
    #[new]
    fn new() -> Self {
        Self {}
    }

    fn __repr__(&self) -> String {
        "LeapFrogPosition()".to_string()
    }
}

/// Berendsen weak-coupling thermostat.
///
/// Rescales velocities to couple the system to a heat bath at target temperature.
/// Placed between LeapFrogVelocity and LeapFrogPosition in the sequence.
///
/// Scaling factor: λ = sqrt(1 + dt/τ · (T₀/T - 1))
///
/// # Example (Python)
///
/// ```python
/// thermostat = BerendsenThermostat(temperature=300.0, tau=0.1)
/// ```
#[pyclass(name = "BerendsenThermostat")]
#[derive(Clone, Debug)]
pub struct PyBerendsenThermostat {
    /// Target temperature T₀ (K).
    #[pyo3(get, set)]
    pub temperature: f64,
    /// Coupling time constant τ (ps). Smaller = stronger coupling.
    #[pyo3(get, set)]
    pub tau: f64,
}

#[pymethods]
impl PyBerendsenThermostat {
    #[new]
    #[pyo3(signature = (temperature=300.0, tau=0.1))]
    fn new(temperature: f64, tau: f64) -> Self {
        Self { temperature, tau }
    }

    fn __repr__(&self) -> String {
        format!(
            "BerendsenThermostat(temperature={:.1}, tau={:.3})",
            self.temperature, self.tau
        )
    }
}

/// Berendsen weak-coupling barostat (isotropic).
///
/// Scales box and positions to couple the system to a pressure bath.
/// Placed after PressureCalculation in the sequence.
///
/// Scaling: μ = (1 - κ·dt/τ·(P₀ - P))^(1/3)
///
/// # Example (Python)
///
/// ```python
/// barostat = BerendsenBarostat(pressure=1.0, tau=0.5, compressibility=4.575e-4)
/// ```
#[pyclass(name = "BerendsenBarostat")]
#[derive(Clone, Debug)]
pub struct PyBerendsenBarostat {
    /// Target pressure P₀ (kJ/(mol·nm³)).
    #[pyo3(get, set)]
    pub pressure: f64,
    /// Coupling time constant τ_P (ps).
    #[pyo3(get, set)]
    pub tau: f64,
    /// Isothermal compressibility κ ((kJ/(mol·nm³))⁻¹).
    #[pyo3(get, set)]
    pub compressibility: f64,
    /// Virial type: "none", "atomic", "molecular".
    #[pyo3(get, set)]
    pub virial: Option<String>,
}

#[pymethods]
impl PyBerendsenBarostat {
    #[new]
    #[pyo3(signature = (pressure=1.0, tau=0.5, compressibility=4.575e-4, virial=None))]
    fn new(pressure: f64, tau: f64, compressibility: f64, virial: Option<String>) -> Self {
        Self {
            pressure,
            tau,
            compressibility,
            virial,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "BerendsenBarostat(pressure={:.2}, tau={:.3}, compressibility={:.3e})",
            self.pressure, self.tau, self.compressibility
        )
    }
}

/// SHAKE bond constraint algorithm.
///
/// Enforces rigid bond lengths via iterative constraint satisfaction.
/// Applied after the position update step.
///
/// # Example (Python)
///
/// ```python
/// shake = ShakeConstraints(tolerance=1e-4, mode="hydrogen")
/// shake = ShakeConstraints(tolerance=1e-4, mode="all")
/// ```
#[pyclass(name = "ShakeConstraints")]
#[derive(Clone, Debug)]
pub struct PyShakeConstraints {
    /// Convergence tolerance for SHAKE iterations.
    #[pyo3(get, set)]
    pub tolerance: f64,
    /// Maximum number of iterations.
    #[pyo3(get, set)]
    pub max_iterations: usize,
    /// Constraint mode: "solvent" (NTC=1), "hydrogen" (NTC=2), "all" (NTC=3).
    #[pyo3(get, set)]
    pub mode: String,
}

#[pymethods]
impl PyShakeConstraints {
    #[new]
    #[pyo3(signature = (tolerance=1e-4, max_iterations=1000, mode="hydrogen".to_string()))]
    fn new(tolerance: f64, max_iterations: usize, mode: String) -> Self {
        Self {
            tolerance,
            max_iterations,
            mode,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "ShakeConstraints(tolerance={:.0e}, mode='{}')",
            self.tolerance, self.mode
        )
    }
}

/// Temperature calculation (kinetic energy → temperature).
///
/// Computes kinetic energy from velocities and stores temperature.
/// Required for thermostat operation and energy output.
///
/// # Example (Python)
///
/// ```python
/// seq.add(TemperatureCalculation())
/// ```
#[pyclass(name = "TemperatureCalculation")]
#[derive(Clone, Debug)]
pub struct PyTemperatureCalculation {}

#[pymethods]
impl PyTemperatureCalculation {
    #[new]
    fn new() -> Self {
        Self {}
    }

    fn __repr__(&self) -> String {
        "TemperatureCalculation()".to_string()
    }
}

/// Pressure calculation (virial → pressure tensor).
///
/// Computes the pressure tensor from kinetic energy tensor and virial.
/// Required for barostat operation.
///
/// # Example (Python)
///
/// ```python
/// seq.add(PressureCalculation(virial="molecular"))
/// ```
#[pyclass(name = "PressureCalculation")]
#[derive(Clone, Debug)]
pub struct PyPressureCalculation {
    /// Virial type: "none", "atomic", "molecular".
    #[pyo3(get, set)]
    pub virial: String,
}

#[pymethods]
impl PyPressureCalculation {
    #[new]
    #[pyo3(signature = (virial="molecular".to_string()))]
    fn new(virial: String) -> Self {
        Self { virial }
    }

    fn __repr__(&self) -> String {
        format!("PressureCalculation(virial='{}')", self.virial)
    }
}

/// Total energy calculation.
///
/// Finalizes potential + kinetic = total energy.
/// Typically the last algorithm in the sequence.
///
/// # Example (Python)
///
/// ```python
/// seq.add(EnergyCalculation())
/// ```
#[pyclass(name = "EnergyCalculation")]
#[derive(Clone, Debug)]
pub struct PyEnergyCalculation {}

#[pymethods]
impl PyEnergyCalculation {
    #[new]
    fn new() -> Self {
        Self {}
    }

    fn __repr__(&self) -> String {
        "EnergyCalculation()".to_string()
    }
}

/// Centre of mass motion removal.
///
/// Removes translational COM velocity from all atoms.
/// Placed first in the sequence (before Forcefield, GROMOS convention).
///
/// # Example (Python)
///
/// ```python
/// com = RemoveCOMMotion(nscm=10)           # remove every 10 steps
/// com = RemoveCOMMotion(initial=True)       # also remove at step 0
/// ```
#[pyclass(name = "RemoveCOMMotion")]
#[derive(Clone, Debug)]
pub struct PyRemoveCOMMotion {
    /// Remove COM motion at step 0 (NTICOM≥1).
    #[pyo3(get, set)]
    pub initial: bool,
    /// Signed periodic removal control: >0 translation every `nscm` steps,
    /// <0 translation+rotation every `|nscm|` steps, 0 = off.
    #[pyo3(get, set)]
    pub nscm: i32,
}

#[pymethods]
impl PyRemoveCOMMotion {
    #[new]
    #[pyo3(signature = (initial=true, nscm=10))]
    fn new(initial: bool, nscm: i32) -> Self {
        Self { initial, nscm }
    }

    fn __repr__(&self) -> String {
        format!(
            "RemoveCOMMotion(initial={}, nscm={})",
            self.initial, self.nscm
        )
    }
}

// ============================================================================
// Algorithm descriptor enum (internal, not exposed to Python)
// ============================================================================

/// Internal representation of an algorithm in the sequence.
#[derive(Clone, Debug)]
pub(crate) enum AlgorithmDescriptor {
    Forcefield(PyForcefield),
    LeapFrogIntegrator(PyLeapFrogIntegrator),
    LeapFrogVelocity(PyLeapFrogVelocity),
    LeapFrogPosition(PyLeapFrogPosition),
    BerendsenThermostat(PyBerendsenThermostat),
    BerendsenBarostat(PyBerendsenBarostat),
    ShakeConstraints(PyShakeConstraints),
    TemperatureCalculation(PyTemperatureCalculation),
    PressureCalculation(PyPressureCalculation),
    EnergyCalculation(PyEnergyCalculation),
    RemoveCOMMotion(PyRemoveCOMMotion),
}

impl AlgorithmDescriptor {
    fn name(&self) -> &str {
        match self {
            Self::Forcefield(_) => "Forcefield",
            Self::LeapFrogIntegrator(_) => "LeapFrogIntegrator",
            Self::LeapFrogVelocity(_) => "LeapFrogVelocity",
            Self::LeapFrogPosition(_) => "LeapFrogPosition",
            Self::BerendsenThermostat(_) => "BerendsenThermostat",
            Self::BerendsenBarostat(_) => "BerendsenBarostat",
            Self::ShakeConstraints(_) => "ShakeConstraints",
            Self::TemperatureCalculation(_) => "TemperatureCalculation",
            Self::PressureCalculation(_) => "PressureCalculation",
            Self::EnergyCalculation(_) => "EnergyCalculation",
            Self::RemoveCOMMotion(_) => "RemoveCOMMotion",
        }
    }
}

// ============================================================================
// PyAlgorithmSequence — the main compositional type
// ============================================================================

/// An ordered sequence of MD algorithms defining a single simulation step.
///
/// Build custom sequences with full control over algorithm ordering and
/// parameters, or use preset constructors for standard ensembles.
///
/// The sequence determines EXACTLY what happens each MD step and in what order.
/// This is the most powerful way to configure a GROMOS simulation.
///
/// # Standard MD sequence (NVT)
///
/// ```python
/// seq = AlgorithmSequence()
/// seq.add(RemoveCOMMotion(nscm=10))
/// seq.add(Forcefield(cutoff=0.8, epsilon_rf=61.0))
/// seq.add(LeapFrogIntegrator())
/// seq.add(BerendsenThermostat(temperature=300.0, tau=0.1))
/// seq.add(ShakeConstraints(tolerance=1e-4))
/// seq.add(TemperatureCalculation())
/// seq.add(EnergyCalculation())
/// ```
///
/// # Presets
///
/// ```python
/// seq = AlgorithmSequence.nve(topo, params)  # microcanonical
/// seq = AlgorithmSequence.nvt(topo, params)  # canonical (Berendsen)
/// seq = AlgorithmSequence.npt(topo, params)  # isothermal-isobaric
/// ```
#[pyclass(name = "AlgorithmSequence")]
#[derive(Clone, Debug)]
pub struct PyAlgorithmSequence {
    pub(crate) algorithms: Vec<AlgorithmDescriptor>,
}

#[pymethods]
impl PyAlgorithmSequence {
    /// Create an empty algorithm sequence.
    #[new]
    fn new() -> Self {
        Self {
            algorithms: Vec::new(),
        }
    }

    /// Add an algorithm to the end of the sequence.
    ///
    /// Args:
    ///     algorithm: Any algorithm descriptor object (Forcefield, LeapFrogIntegrator, etc.)
    ///
    /// Example:
    ///     seq.add(Forcefield(cutoff=0.8))
    ///     seq.add(LeapFrogIntegrator())
    fn add(&mut self, algorithm: &Bound<'_, PyAny>) -> PyResult<()> {
        let desc = extract_descriptor(algorithm)?;
        self.algorithms.push(desc);
        Ok(())
    }

    /// Insert an algorithm after a named algorithm in the sequence.
    ///
    /// Args:
    ///     after: Name of the algorithm to insert after.
    ///     algorithm: The algorithm to insert.
    ///
    /// Raises:
    ///     ValueError: If the named algorithm is not found.
    ///
    /// Example:
    ///     seq.insert_after("LeapFrogIntegrator", BerendsenThermostat(temperature=300.0))
    fn insert_after(&mut self, after: &str, algorithm: &Bound<'_, PyAny>) -> PyResult<()> {
        let desc = extract_descriptor(algorithm)?;
        // Find the position (searching from the start, matching by name)
        let pos = self
            .algorithms
            .iter()
            .position(|a| a.name() == after)
            .ok_or_else(|| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Algorithm '{}' not found in sequence. Available: {:?}",
                    after,
                    self.algorithm_names()
                ))
            })?;
        self.algorithms.insert(pos + 1, desc);
        Ok(())
    }

    /// Insert an algorithm before a named algorithm in the sequence.
    ///
    /// Args:
    ///     before: Name of the algorithm to insert before.
    ///     algorithm: The algorithm to insert.
    ///
    /// Raises:
    ///     ValueError: If the named algorithm is not found.
    ///
    /// Example:
    ///     seq.insert_before("Forcefield", RemoveCOMMotion(nscm=10))
    fn insert_before(&mut self, before: &str, algorithm: &Bound<'_, PyAny>) -> PyResult<()> {
        let desc = extract_descriptor(algorithm)?;
        let pos = self
            .algorithms
            .iter()
            .position(|a| a.name() == before)
            .ok_or_else(|| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Algorithm '{}' not found in sequence. Available: {:?}",
                    before,
                    self.algorithm_names()
                ))
            })?;
        self.algorithms.insert(pos, desc);
        Ok(())
    }

    /// Remove an algorithm by name.
    ///
    /// Args:
    ///     name: Name of the algorithm to remove.
    ///
    /// Raises:
    ///     ValueError: If the named algorithm is not found.
    ///
    /// Example:
    ///     seq.remove("BerendsenThermostat")
    fn remove(&mut self, name: &str) -> PyResult<()> {
        let pos = self
            .algorithms
            .iter()
            .position(|a| a.name() == name)
            .ok_or_else(|| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Algorithm '{}' not found in sequence. Available: {:?}",
                    name,
                    self.algorithm_names()
                ))
            })?;
        self.algorithms.remove(pos);
        Ok(())
    }

    /// Replace an algorithm by name with a new one.
    ///
    /// Args:
    ///     name: Name of the algorithm to replace.
    ///     algorithm: The new algorithm to use in its place.
    ///
    /// Raises:
    ///     ValueError: If the named algorithm is not found.
    ///
    /// Example:
    ///     seq.replace("BerendsenThermostat", BerendsenThermostat(temperature=350.0, tau=0.05))
    fn replace(&mut self, name: &str, algorithm: &Bound<'_, PyAny>) -> PyResult<()> {
        let desc = extract_descriptor(algorithm)?;
        let pos = self
            .algorithms
            .iter()
            .position(|a| a.name() == name)
            .ok_or_else(|| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Algorithm '{}' not found in sequence. Available: {:?}",
                    name,
                    self.algorithm_names()
                ))
            })?;
        self.algorithms[pos] = desc;
        Ok(())
    }

    /// Number of algorithms in the sequence.
    fn __len__(&self) -> usize {
        self.algorithms.len()
    }

    /// List of algorithm names in order.
    #[getter]
    fn names(&self) -> Vec<String> {
        self.algorithm_names()
    }

    /// Check if a named algorithm is in the sequence.
    fn __contains__(&self, name: &str) -> bool {
        self.algorithms.iter().any(|a| a.name() == name)
    }

    fn __repr__(&self) -> String {
        let names: Vec<&str> = self.algorithms.iter().map(|a| a.name()).collect();
        format!("AlgorithmSequence([\n  {}\n])", names.join(",\n  "))
    }

    // ========================================================================
    // Preset constructors
    // ========================================================================

    /// Create a standard NVE (microcanonical) algorithm sequence.
    ///
    /// Uses parameters from Topology and InputParameters for forcefield configuration.
    /// No thermostat or barostat.
    ///
    /// Args:
    ///     topo: Topology object
    ///     params: InputParameters object
    ///
    /// Example:
    ///     seq = AlgorithmSequence.nve(topo, params)
    #[staticmethod]
    fn nve(_topo: &PyTopology, params: &PyInputParameters) -> PyResult<Self> {
        let imd = &params.inner;
        let mut seq = Self::new();

        // COM removal
        if imd.nticom >= 1 || imd.nscm != 0 {
            seq.algorithms
                .push(AlgorithmDescriptor::RemoveCOMMotion(PyRemoveCOMMotion {
                    initial: imd.nticom >= 1,
                    nscm: imd.nscm,
                }));
        }

        // Forcefield
        seq.algorithms
            .push(AlgorithmDescriptor::Forcefield(forcefield_from_imd(imd)));

        // Leap-Frog integrator
        seq.algorithms.push(AlgorithmDescriptor::LeapFrogIntegrator(
            PyLeapFrogIntegrator {},
        ));

        // SHAKE if needed
        if imd.ntc > 1 || (imd.ntcs > 0 && imd.nsm > 0) {
            seq.algorithms
                .push(AlgorithmDescriptor::ShakeConstraints(shake_from_imd(imd)));
        }

        // Temperature + Energy calculation
        seq.algorithms
            .push(AlgorithmDescriptor::TemperatureCalculation(
                PyTemperatureCalculation {},
            ));
        seq.algorithms.push(AlgorithmDescriptor::EnergyCalculation(
            PyEnergyCalculation {},
        ));

        Ok(seq)
    }

    /// Create a standard NVT (canonical) algorithm sequence with Berendsen thermostat.
    ///
    /// Args:
    ///     topo: Topology object
    ///     params: InputParameters object
    ///
    /// Example:
    ///     seq = AlgorithmSequence.nvt(topo, params)
    #[staticmethod]
    fn nvt(_topo: &PyTopology, params: &PyInputParameters) -> PyResult<Self> {
        let imd = &params.inner;
        let mut seq = Self::new();

        // COM removal
        if imd.nticom >= 1 || imd.nscm != 0 {
            seq.algorithms
                .push(AlgorithmDescriptor::RemoveCOMMotion(PyRemoveCOMMotion {
                    initial: imd.nticom >= 1,
                    nscm: imd.nscm,
                }));
        }

        // Forcefield
        seq.algorithms
            .push(AlgorithmDescriptor::Forcefield(forcefield_from_imd(imd)));

        // Leap-Frog velocity
        seq.algorithms
            .push(AlgorithmDescriptor::LeapFrogVelocity(PyLeapFrogVelocity {}));

        // Berendsen thermostat (between velocity and position, GROMOS convention)
        let temperature = if !imd.temp_bath.is_empty() && !imd.temp_bath[0].temp0.is_empty() {
            imd.temp_bath[0].temp0[0]
        } else {
            300.0
        };
        let tau = if !imd.temp_bath.is_empty() && !imd.temp_bath[0].tau.is_empty() {
            imd.temp_bath[0].tau[0]
        } else {
            0.1
        };
        seq.algorithms
            .push(AlgorithmDescriptor::BerendsenThermostat(
                PyBerendsenThermostat { temperature, tau },
            ));

        // Leap-Frog position
        seq.algorithms
            .push(AlgorithmDescriptor::LeapFrogPosition(PyLeapFrogPosition {}));

        // SHAKE if needed
        if imd.ntc > 1 || (imd.ntcs > 0 && imd.nsm > 0) {
            seq.algorithms
                .push(AlgorithmDescriptor::ShakeConstraints(shake_from_imd(imd)));
        }

        // Temperature + Energy
        seq.algorithms
            .push(AlgorithmDescriptor::TemperatureCalculation(
                PyTemperatureCalculation {},
            ));
        seq.algorithms.push(AlgorithmDescriptor::EnergyCalculation(
            PyEnergyCalculation {},
        ));

        Ok(seq)
    }

    /// Create a standard NPT (isothermal-isobaric) algorithm sequence.
    ///
    /// Uses Berendsen thermostat + Berendsen barostat.
    ///
    /// Args:
    ///     topo: Topology object
    ///     params: InputParameters object
    ///
    /// Example:
    ///     seq = AlgorithmSequence.npt(topo, params)
    #[staticmethod]
    fn npt(_topo: &PyTopology, params: &PyInputParameters) -> PyResult<Self> {
        let imd = &params.inner;
        let mut seq = Self::new();

        // COM removal
        if imd.nticom >= 1 || imd.nscm != 0 {
            seq.algorithms
                .push(AlgorithmDescriptor::RemoveCOMMotion(PyRemoveCOMMotion {
                    initial: imd.nticom >= 1,
                    nscm: imd.nscm,
                }));
        }

        // Forcefield with virial
        let mut ff = forcefield_from_imd(imd);
        let virial_str = if let Some(ref pp) = imd.pressure_parameters {
            match pp.virial {
                2 => "molecular",
                1 => "atomic",
                _ => "none",
            }
        } else {
            "molecular"
        };
        ff.virial = Some(virial_str.to_string());
        seq.algorithms.push(AlgorithmDescriptor::Forcefield(ff));

        // Leap-Frog velocity
        seq.algorithms
            .push(AlgorithmDescriptor::LeapFrogVelocity(PyLeapFrogVelocity {}));

        // Berendsen thermostat (between velocity and position)
        let temperature = if !imd.temp_bath.is_empty() && !imd.temp_bath[0].temp0.is_empty() {
            imd.temp_bath[0].temp0[0]
        } else {
            300.0
        };
        let tau = if !imd.temp_bath.is_empty() && !imd.temp_bath[0].tau.is_empty() {
            imd.temp_bath[0].tau[0]
        } else {
            0.1
        };
        seq.algorithms
            .push(AlgorithmDescriptor::BerendsenThermostat(
                PyBerendsenThermostat { temperature, tau },
            ));

        // Leap-Frog position
        seq.algorithms
            .push(AlgorithmDescriptor::LeapFrogPosition(PyLeapFrogPosition {}));

        // SHAKE if needed
        if imd.ntc > 1 || (imd.ntcs > 0 && imd.nsm > 0) {
            seq.algorithms
                .push(AlgorithmDescriptor::ShakeConstraints(shake_from_imd(imd)));
        }

        // Temperature + Pressure + Barostat + Energy
        seq.algorithms
            .push(AlgorithmDescriptor::TemperatureCalculation(
                PyTemperatureCalculation {},
            ));

        let pp = imd.pressure_parameters.as_ref();
        seq.algorithms
            .push(AlgorithmDescriptor::PressureCalculation(
                PyPressureCalculation {
                    virial: virial_str.to_string(),
                },
            ));
        seq.algorithms.push(AlgorithmDescriptor::BerendsenBarostat(
            PyBerendsenBarostat {
                pressure: pp.map(|p| p.pressure0[0][0]).unwrap_or(1.0),
                tau: pp.map(|p| p.tau_p).unwrap_or(0.5),
                compressibility: pp.map(|p| p.compressibility[0][0]).unwrap_or(4.575e-4),
                virial: Some(virial_str.to_string()),
            },
        ));

        seq.algorithms.push(AlgorithmDescriptor::EnergyCalculation(
            PyEnergyCalculation {},
        ));

        Ok(seq)
    }

    /// Create an algorithm sequence from InputParameters (auto-detect ensemble).
    ///
    /// Inspects parameters to determine NVE/NVT/NPT and builds the appropriate
    /// sequence. Equivalent to what the md binary does internally.
    ///
    /// Args:
    ///     topo: Topology object
    ///     params: InputParameters object
    ///
    /// Example:
    ///     seq = AlgorithmSequence.from_parameters(topo, params)
    #[staticmethod]
    fn from_parameters(_topo: &PyTopology, params: &PyInputParameters) -> PyResult<Self> {
        let imd = &params.inner;
        if imd.couple_pressure {
            Self::npt(_topo, params)
        } else {
            let thermostat_on = if !imd.temp_bath.is_empty() && !imd.temp_bath[0].tau.is_empty() {
                imd.temp_bath[0].tau[0] > 0.0
            } else {
                false
            };
            if thermostat_on {
                Self::nvt(_topo, params)
            } else {
                Self::nve(_topo, params)
            }
        }
    }
}

// Private helpers
impl PyAlgorithmSequence {
    fn algorithm_names(&self) -> Vec<String> {
        self.algorithms
            .iter()
            .map(|a| a.name().to_string())
            .collect()
    }
}

// ============================================================================
// Resolution: descriptors → real Rust AlgorithmSequence
// ============================================================================

/// Resolve a PyAlgorithmSequence into a real Rust AlgorithmSequence.
///
/// This is called by the Simulation constructor. It uses the topology and
/// configuration to build the actual algorithm objects with full state.
pub(crate) fn resolve_algorithm_sequence(
    descriptors: &PyAlgorithmSequence,
    topo: &Topology,
    conf: &Configuration,
    imd: &ImdParameters,
    box_dims: Vec3,
) -> Result<AlgorithmSequence, String> {
    let n_atoms = topo.num_atoms();
    let mut md_sequence = AlgorithmSequence::new();

    for desc in &descriptors.algorithms {
        match desc {
            AlgorithmDescriptor::RemoveCOMMotion(d) => {
                let nticom = if d.initial { 1 } else { 0 };
                md_sequence.push(Box::new(RustRemoveCOM::new(nticom, d.nscm)));
            },

            AlgorithmDescriptor::Forcefield(d) => {
                let cutoff = d.cutoff.unwrap_or(imd.rcutl);
                let rcutp = d.rcutp.unwrap_or(imd.rcutp);
                let epsilon_rf = d.epsilon_rf.unwrap_or(imd.epsrf);
                let kappa = d.kappa.unwrap_or(imd.appak);
                let pairlist_update = d.pairlist_update.unwrap_or(imd.nsnb);

                let crf_params = CRFParameters::new(cutoff, 1.0, epsilon_rf, kappa);
                let lj_params = RustForcefield::convert_lj_parameters(topo);

                let mut pairlist = PairlistContainer::new(rcutp, cutoff, 0.0);
                pairlist.update_frequency = pairlist_update;

                let periodicity = if box_dims.x == 0.0 && box_dims.y == 0.0 && box_dims.z == 0.0 {
                    Periodicity::Vacuum(Vacuum)
                } else {
                    Periodicity::Rectangular(Rectangular::new(box_dims))
                };
                let box_type = match &periodicity {
                    Periodicity::Rectangular(_) => BoxType::Rectangular,
                    Periodicity::Triclinic(_) => BoxType::Triclinic,
                    Periodicity::Vacuum(_) => BoxType::Vacuum,
                };
                let pairlist_algorithm = PairlistAlgorithm::from_imd(
                    imd.algorithm,
                    topo.num_atoms(),
                    box_type,
                    !topo.chargegroups.is_empty(),
                );

                pairlist_algorithm.update(topo, conf, &mut pairlist, &periodicity);

                let mut forcefield = RustForcefield::new(
                    lj_params,
                    crf_params,
                    periodicity,
                    pairlist,
                    pairlist_algorithm,
                );
                forcefield.ntf_bond = d.ntf_bond.unwrap_or(imd.ntf[0] != 0);
                forcefield.ntf_angle = d.ntf_angle.unwrap_or(imd.ntf[1] != 0);
                forcefield.ntf_improper = d.ntf_improper.unwrap_or(imd.ntf[2] != 0);
                forcefield.ntf_dihedral = d.ntf_dihedral.unwrap_or(imd.ntf[3] != 0);
                if !topo.solvent_atom_template.is_empty() {
                    forcefield.atoms_per_solvent = topo.solvent_atom_template.len();
                }

                // Virial type
                let virial_type = match d.virial.as_deref() {
                    Some("molecular") => VirialType::Molecular,
                    Some("atomic") => VirialType::Atomic,
                    Some("none") => VirialType::None,
                    _ => {
                        if imd.couple_pressure {
                            match imd
                                .pressure_parameters
                                .as_ref()
                                .map(|p| p.virial)
                                .unwrap_or(0)
                            {
                                2 => VirialType::Molecular,
                                1 => VirialType::Atomic,
                                _ => VirialType::None,
                            }
                        } else {
                            VirialType::None
                        }
                    },
                };
                forcefield.virial_type = virial_type;

                md_sequence.push(Box::new(forcefield));
            },

            AlgorithmDescriptor::LeapFrogIntegrator(_) => {
                md_sequence.push(Box::new(LeapFrogVelocity::new()));
                md_sequence.push(Box::new(LeapFrogPosition::new()));
            },

            AlgorithmDescriptor::LeapFrogVelocity(_) => {
                md_sequence.push(Box::new(LeapFrogVelocity::new()));
            },

            AlgorithmDescriptor::LeapFrogPosition(_) => {
                md_sequence.push(Box::new(LeapFrogPosition::new()));
            },

            AlgorithmDescriptor::BerendsenThermostat(d) => {
                let n_solute = topo.num_solute_atoms();
                let atoms_per_solvent = if !topo.solvent_atom_template.is_empty() {
                    topo.solvent_atom_template.len()
                } else {
                    1
                };
                let n_solvent_molecules = if atoms_per_solvent > 0 && n_atoms > n_solute {
                    (n_atoms - n_solute) / atoms_per_solvent
                } else {
                    0
                };
                let shake_enabled = imd.ntc > 1 || (imd.ntcs > 0 && imd.nsm > 0);
                let solvent_constraint_dof = if shake_enabled {
                    n_solvent_molecules * topo.solvent_constraint_template.len()
                } else {
                    0
                };
                let total_dof = (3 * n_atoms - solvent_constraint_dof) as f64 - imd.ndfmin as f64;

                md_sequence.push(Box::new(RustThermostat::new_single_bath(
                    d.temperature,
                    d.tau,
                    total_dof,
                    n_atoms,
                )));
            },

            AlgorithmDescriptor::ShakeConstraints(d) => {
                let ntc_mode = match d.mode.as_str() {
                    "all" => NtcMode::AllBonds,
                    "hydrogen" => NtcMode::HydrogenBonds,
                    _ => NtcMode::SolventOnly,
                };
                let mut shake_alg = RustShake::new(ShakeParameters {
                    tolerance: d.tolerance,
                    max_iterations: d.max_iterations,
                    ntc: ntc_mode,
                });
                if imd.ntishk >= 1 {
                    shake_alg.shake_initial_positions = true;
                }
                if imd.ntishk >= 2 {
                    shake_alg.shake_initial_velocities = true;
                }
                md_sequence.push(Box::new(shake_alg));
            },

            AlgorithmDescriptor::TemperatureCalculation(_) => {
                md_sequence.push(Box::new(RustTemperature::new()));
            },

            AlgorithmDescriptor::PressureCalculation(d) => {
                let virial_type = match d.virial.as_str() {
                    "molecular" => VirialType::Molecular,
                    "atomic" => VirialType::Atomic,
                    _ => VirialType::None,
                };
                md_sequence.push(Box::new(RustPressure::new(virial_type)));
            },

            AlgorithmDescriptor::BerendsenBarostat(d) => {
                md_sequence.push(Box::new(RustBarostat::new(BerendsenBarostatParams {
                    pressure0: d.pressure,
                    compressibility: d.compressibility,
                    tau: d.tau,
                })));
            },

            AlgorithmDescriptor::EnergyCalculation(_) => {
                md_sequence.push(Box::new(RustEnergy::new()));
            },
        }
    }

    Ok(md_sequence)
}

// ============================================================================
// Helper functions
// ============================================================================

/// Extract an AlgorithmDescriptor from a Python object.
fn extract_descriptor(obj: &Bound<'_, PyAny>) -> PyResult<AlgorithmDescriptor> {
    if let Ok(d) = obj.extract::<PyForcefield>() {
        return Ok(AlgorithmDescriptor::Forcefield(d));
    }
    if let Ok(d) = obj.extract::<PyLeapFrogIntegrator>() {
        return Ok(AlgorithmDescriptor::LeapFrogIntegrator(d));
    }
    if let Ok(d) = obj.extract::<PyLeapFrogVelocity>() {
        return Ok(AlgorithmDescriptor::LeapFrogVelocity(d));
    }
    if let Ok(d) = obj.extract::<PyLeapFrogPosition>() {
        return Ok(AlgorithmDescriptor::LeapFrogPosition(d));
    }
    if let Ok(d) = obj.extract::<PyBerendsenThermostat>() {
        return Ok(AlgorithmDescriptor::BerendsenThermostat(d));
    }
    if let Ok(d) = obj.extract::<PyBerendsenBarostat>() {
        return Ok(AlgorithmDescriptor::BerendsenBarostat(d));
    }
    if let Ok(d) = obj.extract::<PyShakeConstraints>() {
        return Ok(AlgorithmDescriptor::ShakeConstraints(d));
    }
    if let Ok(d) = obj.extract::<PyTemperatureCalculation>() {
        return Ok(AlgorithmDescriptor::TemperatureCalculation(d));
    }
    if let Ok(d) = obj.extract::<PyPressureCalculation>() {
        return Ok(AlgorithmDescriptor::PressureCalculation(d));
    }
    if let Ok(d) = obj.extract::<PyEnergyCalculation>() {
        return Ok(AlgorithmDescriptor::EnergyCalculation(d));
    }
    if let Ok(d) = obj.extract::<PyRemoveCOMMotion>() {
        return Ok(AlgorithmDescriptor::RemoveCOMMotion(d));
    }

    Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(format!(
        "Expected an algorithm descriptor (Forcefield, LeapFrogIntegrator, \
         LeapFrogVelocity, LeapFrogPosition, BerendsenThermostat, \
         BerendsenBarostat, ShakeConstraints, TemperatureCalculation, \
         PressureCalculation, EnergyCalculation, RemoveCOMMotion), got: {}",
        obj.get_type().name().unwrap_or_default()
    )))
}

/// Build forcefield descriptor from ImdParameters.
fn forcefield_from_imd(imd: &ImdParameters) -> PyForcefield {
    PyForcefield {
        cutoff: Some(imd.rcutl),
        rcutp: Some(imd.rcutp),
        epsilon_rf: Some(imd.epsrf),
        kappa: Some(imd.appak),
        pairlist_update: Some(imd.nsnb),
        virial: None,
        ntf_bond: Some(imd.ntf[0] != 0),
        ntf_angle: Some(imd.ntf[1] != 0),
        ntf_improper: Some(imd.ntf[2] != 0),
        ntf_dihedral: Some(imd.ntf[3] != 0),
    }
}

/// Build SHAKE descriptor from ImdParameters.
fn shake_from_imd(imd: &ImdParameters) -> PyShakeConstraints {
    let mode = match imd.ntc {
        3 => "all",
        2 => "hydrogen",
        _ => "solvent",
    };
    PyShakeConstraints {
        tolerance: imd.shake_tol,
        max_iterations: 1000,
        mode: mode.to_string(),
    }
}

/// Register all algorithm sequence classes in the Python module.
pub fn register_algorithm_sequence(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyAlgorithmSequence>()?;
    m.add_class::<PyForcefield>()?;
    m.add_class::<PyLeapFrogIntegrator>()?;
    m.add_class::<PyLeapFrogVelocity>()?;
    m.add_class::<PyLeapFrogPosition>()?;
    m.add_class::<PyBerendsenThermostat>()?;
    m.add_class::<PyBerendsenBarostat>()?;
    m.add_class::<PyShakeConstraints>()?;
    m.add_class::<PyTemperatureCalculation>()?;
    m.add_class::<PyPressureCalculation>()?;
    m.add_class::<PyEnergyCalculation>()?;
    m.add_class::<PyRemoveCOMMotion>()?;
    Ok(())
}
