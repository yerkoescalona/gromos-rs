//! Python bindings for GROMOS-RS molecular dynamics engine
//!
//! Inspired by Polars' architecture:
//! - Zero-copy data sharing where possible
//! - Rust core with Python wrapper
//! - High-performance parallel execution
//! - Memory-safe operations

use numpy::{
    PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods, ToPyArray,
};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

// Re-export from gromos-rs
use gromos_rs::{
    // Virtual atoms
    algorithm::virtual_atoms::{
        VirtualAtom as RustVirtualAtom, VirtualAtomsManager as RustVirtualAtomsManager,
    },
    // FEP
    fep::{
        InteractionLambdas as RustInteractionLambdas, LambdaController as RustLambdaController,
        PerturbedAtom as RustPerturbedAtom,
    },
    // Local elevation
    interaction::local_elevation::{
        CoordinateType as RustCoordinateType, LECoordinate as RustLECoordinate,
        Umbrella as RustUmbrella,
    },
    // Perturbed nonbonded
    interaction::nonbonded::{
        perturbed_lj_crf_interaction, PerturbedLambdaParams as RustPerturbedLambdaParams,
    },
    // Polarization
    interaction::polarization::{
        PolarizabilityParameters as RustPolarizabilityParameters,
        PolarizationCalculator as RustPolarizationCalculator,
        PolarizationModel as RustPolarizationModel,
    },
    // QM/MM
    interaction::qmmm::{
        CouplingScheme as RustCouplingScheme, QMMMCalculator as RustQMMMCalculator,
        QMMethod as RustQMMethod, QMRegion as RustQMRegion,
    },
    // NMR restraints
    interaction::restraints::{
        JValueRestraint as RustJValueRestraint, RDCRestraint as RustRDCRestraint,
    },
    io::energy_binary::BinaryEnergyReader,
    io::topology::{build_topology, read_topology_file},
    io::trajectory::TrajectoryReader,
    io::trajectory_binary::{BinaryTrajectoryReader, DcdReader},
    selection::AtomSelection,
    Configuration as RustConfiguration,
    Energy as RustEnergy,
    Mat3 as RustMat3,
    State as RustState,
    Topology as RustTopology,
    Vec3 as RustVec3,
};

/// Python module for GROMOS molecular dynamics
#[pymodule]
fn gromos(_py: Python, m: &PyModule) -> PyResult<()> {
    // Math types
    m.add_class::<PyVec3>()?;
    m.add_class::<PyMat3>()?;

    // Core data structures
    m.add_class::<PyBox>()?;
    m.add_class::<PyEnergy>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyConfiguration>()?;
    m.add_class::<PyTopology>()?;

    // Binary I/O readers
    m.add_class::<PyDcdReader>()?;
    m.add_class::<PyBinaryEnergyReader>()?;

    // NMR Restraints
    m.add_class::<PyJValueRestraint>()?;
    m.add_class::<PyRDCRestraint>()?;

    // Virtual Atoms
    m.add_class::<PyVirtualAtom>()?;
    m.add_class::<PyVirtualAtomsManager>()?;

    // Local Elevation / Metadynamics
    m.add_class::<PyUmbrella>()?;
    m.add_class::<PyLECoordinate>()?;
    m.add_class::<PyCoordinateType>()?;

    // Polarization
    m.add_class::<PyPolarizationCalculator>()?;
    m.add_class::<PyPolarizationModel>()?;
    m.add_class::<PyPolarizabilityParameters>()?;

    // QM/MM
    m.add_class::<PyQMMMCalculator>()?;
    m.add_class::<PyQMRegion>()?;
    m.add_class::<PyQMMethod>()?;
    m.add_class::<PyCouplingScheme>()?;

    // Free Energy Perturbation
    m.add_class::<PyLambdaController>()?;
    m.add_class::<PyPerturbedAtom>()?;
    m.add_class::<PyInteractionLambdas>()?;
    m.add_class::<PyPerturbedLambdaParams>()?;

    // Analysis functions
    m.add_function(wrap_pyfunction!(calculate_rmsd, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_rmsf, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_rgyr, m)?)?;
    m.add_function(wrap_pyfunction!(analyze_trajectory, m)?)?;

    // Conversion functions
    m.add_function(wrap_pyfunction!(convert_trajectory, m)?)?;
    m.add_function(wrap_pyfunction!(convert_energy, m)?)?;

    // Perturbed nonbonded interaction
    m.add_function(wrap_pyfunction!(py_perturbed_lj_crf_interaction, m)?)?;

    Ok(())
}

//==============================================================================
// MATH TYPES
//==============================================================================

/// 3D vector with SIMD acceleration
///
/// Wraps Rust's Vec3A (SIMD-accelerated vector from glam crate).
/// Provides zero-copy conversion to/from NumPy arrays where possible.
#[pyclass(name = "Vec3")]
#[derive(Clone)]
pub struct PyVec3 {
    inner: RustVec3,
}

#[pymethods]
impl PyVec3 {
    /// Create a new Vec3
    #[new]
    fn new(x: f32, y: f32, z: f32) -> Self {
        Self {
            inner: RustVec3::new(x, y, z),
        }
    }

    /// Get x component
    #[getter]
    fn x(&self) -> f32 {
        self.inner.x
    }

    /// Get y component
    #[getter]
    fn y(&self) -> f32 {
        self.inner.y
    }

    /// Get z component
    #[getter]
    fn z(&self) -> f32 {
        self.inner.z
    }

    /// Set x component
    #[setter]
    fn set_x(&mut self, value: f32) {
        self.inner.x = value;
    }

    /// Set y component
    #[setter]
    fn set_y(&mut self, value: f32) {
        self.inner.y = value;
    }

    /// Set z component
    #[setter]
    fn set_z(&mut self, value: f32) {
        self.inner.z = value;
    }

    /// Vector length (magnitude)
    fn length(&self) -> f32 {
        self.inner.length()
    }

    /// Squared length (faster than length)
    fn length_squared(&self) -> f32 {
        self.inner.length_squared()
    }

    /// Normalize the vector (unit length)
    fn normalize(&self) -> Self {
        Self {
            inner: self.inner.normalize(),
        }
    }

    /// Dot product
    fn dot(&self, other: &Self) -> f32 {
        self.inner.dot(other.inner)
    }

    /// Cross product
    fn cross(&self, other: &Self) -> Self {
        Self {
            inner: self.inner.cross(other.inner),
        }
    }

    /// Distance to another vector
    fn distance(&self, other: &Self) -> f32 {
        self.inner.distance(other.inner)
    }

    /// Add two vectors
    fn __add__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner + other.inner,
        }
    }

    /// Subtract two vectors
    fn __sub__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner - other.inner,
        }
    }

    /// Multiply by scalar
    fn __mul__(&self, scalar: f32) -> Self {
        Self {
            inner: self.inner * scalar,
        }
    }

    /// String representation
    fn __repr__(&self) -> String {
        format!(
            "Vec3({:.4}, {:.4}, {:.4})",
            self.inner.x, self.inner.y, self.inner.z
        )
    }

    /// Convert to NumPy array
    fn to_numpy<'py>(&self, py: Python<'py>) -> &'py PyArray1<f32> {
        let arr = [self.inner.x, self.inner.y, self.inner.z];
        arr.to_pyarray(py)
    }

    /// Create from NumPy array
    #[staticmethod]
    fn from_numpy(arr: PyReadonlyArray1<f32>) -> PyResult<Self> {
        let slice = arr.as_slice()?;
        if slice.len() != 3 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Array must have exactly 3 elements",
            ));
        }
        Ok(Self {
            inner: RustVec3::new(slice[0], slice[1], slice[2]),
        })
    }
}

/// 3x3 matrix with SIMD acceleration
#[pyclass(name = "Mat3")]
#[derive(Clone)]
pub struct PyMat3 {
    inner: RustMat3,
}

#[pymethods]
impl PyMat3 {
    /// Create identity matrix
    #[staticmethod]
    fn identity() -> Self {
        Self {
            inner: RustMat3::IDENTITY,
        }
    }

    /// Create from column vectors
    #[staticmethod]
    fn from_cols(x: &PyVec3, y: &PyVec3, z: &PyVec3) -> Self {
        Self {
            inner: RustMat3::from_cols(x.inner, y.inner, z.inner),
        }
    }

    /// Matrix determinant
    fn determinant(&self) -> f32 {
        self.inner.determinant()
    }

    /// Matrix inverse
    fn inverse(&self) -> Self {
        Self {
            inner: self.inner.inverse(),
        }
    }

    /// Matrix transpose
    fn transpose(&self) -> Self {
        Self {
            inner: self.inner.transpose(),
        }
    }

    /// Matrix-vector multiplication
    fn mul_vec3(&self, v: &PyVec3) -> PyVec3 {
        PyVec3 {
            inner: self.inner * v.inner,
        }
    }

    /// Convert to NumPy array (3x3)
    fn to_numpy<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
        let flat = vec![
            self.inner.x_axis.x,
            self.inner.y_axis.x,
            self.inner.z_axis.x,
            self.inner.x_axis.y,
            self.inner.y_axis.y,
            self.inner.z_axis.y,
            self.inner.x_axis.z,
            self.inner.y_axis.z,
            self.inner.z_axis.z,
        ];
        PyArray2::from_vec2(py, &[flat])
            .unwrap()
            .reshape([3, 3])
            .unwrap()
    }

    fn __repr__(&self) -> String {
        format!("Mat3([[{:.4}, {:.4}, {:.4}],\n     [{:.4}, {:.4}, {:.4}],\n     [{:.4}, {:.4}, {:.4}]])",
            self.inner.x_axis.x, self.inner.y_axis.x, self.inner.z_axis.x,
            self.inner.x_axis.y, self.inner.y_axis.y, self.inner.z_axis.y,
            self.inner.x_axis.z, self.inner.y_axis.z, self.inner.z_axis.z)
    }
}

//==============================================================================
// CORE DATA STRUCTURES
//==============================================================================

/// Simulation box representation
#[pyclass(name = "Box")]
#[derive(Clone)]
pub struct PyBox {
    inner: gromos_rs::configuration::Box,
}

#[pymethods]
impl PyBox {
    /// Create vacuum box (no periodicity)
    #[staticmethod]
    fn vacuum() -> Self {
        Self {
            inner: gromos_rs::configuration::Box::vacuum(),
        }
    }

    /// Create rectangular box
    #[staticmethod]
    fn rectangular(lx: f32, ly: f32, lz: f32) -> Self {
        Self {
            inner: gromos_rs::configuration::Box::rectangular(lx, ly, lz),
        }
    }

    /// Create triclinic box from vectors
    #[staticmethod]
    fn triclinic(mat: &PyMat3) -> Self {
        Self {
            inner: gromos_rs::configuration::Box::triclinic(mat.inner),
        }
    }

    /// Get box volume
    fn volume(&self) -> f64 {
        self.inner.volume()
    }

    /// Get box dimensions (for rectangular box)
    fn dimensions(&self) -> PyVec3 {
        PyVec3 {
            inner: self.inner.dimensions(),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Box(type={:?}, volume={:.2})",
            self.inner.box_type,
            self.volume()
        )
    }
}

/// Energy storage for molecular system
///
/// Stores total energies and per-group energies for detailed accounting.
/// All energies are in kJ/mol.
#[pyclass(name = "Energy")]
#[derive(Clone)]
pub struct PyEnergy {
    inner: RustEnergy,
}

#[pymethods]
impl PyEnergy {
    /// Create new energy object
    #[new]
    fn new(num_temperature_groups: usize, num_energy_groups: usize) -> Self {
        Self {
            inner: RustEnergy::new(num_temperature_groups, num_energy_groups),
        }
    }

    /// Total energy (kinetic + potential)
    fn total(&self) -> f64 {
        self.inner.total()
    }

    /// Kinetic energy
    #[getter]
    fn kinetic(&self) -> f64 {
        self.inner.kinetic_total
    }

    /// Potential energy
    #[getter]
    fn potential(&self) -> f64 {
        self.inner.potential_total
    }

    /// Bond energy
    #[getter]
    fn bond(&self) -> f64 {
        self.inner.bond_total
    }

    /// Angle energy
    #[getter]
    fn angle(&self) -> f64 {
        self.inner.angle_total
    }

    /// Dihedral energy
    #[getter]
    fn dihedral(&self) -> f64 {
        self.inner.dihedral_total
    }

    /// Lennard-Jones energy
    #[getter]
    fn lj(&self) -> f64 {
        self.inner.lj_total
    }

    /// Coulomb reaction field energy
    #[getter]
    fn coulomb(&self) -> f64 {
        self.inner.crf_total
    }

    /// Clear all energies to zero
    fn clear(&mut self) {
        self.inner.clear();
    }

    /// Get energy as dictionary
    fn to_dict(&self, py: Python) -> PyResult<PyObject> {
        let dict = PyDict::new(py);
        dict.set_item("total", self.total())?;
        dict.set_item("kinetic", self.inner.kinetic_total)?;
        dict.set_item("potential", self.inner.potential_total)?;
        dict.set_item("bond", self.inner.bond_total)?;
        dict.set_item("angle", self.inner.angle_total)?;
        dict.set_item("dihedral", self.inner.dihedral_total)?;
        dict.set_item("lj", self.inner.lj_total)?;
        dict.set_item("coulomb", self.inner.crf_total)?;
        Ok(dict.into())
    }

    fn __repr__(&self) -> String {
        format!(
            "Energy(total={:.2} kJ/mol, kinetic={:.2}, potential={:.2})",
            self.total(),
            self.inner.kinetic_total,
            self.inner.potential_total
        )
    }
}

/// System state (positions, velocities, forces)
///
/// Uses zero-copy sharing with NumPy arrays for efficient data access.
#[pyclass(name = "State")]
pub struct PyState {
    inner: RustState,
}

#[pymethods]
impl PyState {
    /// Create new state for N atoms
    #[new]
    fn new(num_atoms: usize, num_temp_groups: usize, num_energy_groups: usize) -> Self {
        Self {
            inner: RustState::new(num_atoms, num_temp_groups, num_energy_groups),
        }
    }

    /// Number of atoms
    fn num_atoms(&self) -> usize {
        self.inner.pos.len()
    }

    /// Get positions as NumPy array (N x 3)
    fn positions<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
        let n = self.inner.pos.len();
        let mut arr = Vec::with_capacity(n * 3);
        for v in &self.inner.pos {
            arr.push(v.x);
            arr.push(v.y);
            arr.push(v.z);
        }
        PyArray2::from_vec2(py, &[arr])
            .unwrap()
            .reshape([n, 3])
            .unwrap()
    }

    /// Get velocities as NumPy array (N x 3)
    fn velocities<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
        let n = self.inner.vel.len();
        let mut arr = Vec::with_capacity(n * 3);
        for v in &self.inner.vel {
            arr.push(v.x);
            arr.push(v.y);
            arr.push(v.z);
        }
        PyArray2::from_vec2(py, &[arr])
            .unwrap()
            .reshape([n, 3])
            .unwrap()
    }

    /// Get forces as NumPy array (N x 3)
    fn forces<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
        let n = self.inner.force.len();
        let mut arr = Vec::with_capacity(n * 3);
        for v in &self.inner.force {
            arr.push(v.x);
            arr.push(v.y);
            arr.push(v.z);
        }
        PyArray2::from_vec2(py, &[arr])
            .unwrap()
            .reshape([n, 3])
            .unwrap()
    }

    /// Set positions from NumPy array (N x 3)
    fn set_positions(&mut self, arr: PyReadonlyArray2<f32>) -> PyResult<()> {
        let shape = arr.shape();
        if shape[1] != 3 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Array must have shape (N, 3)",
            ));
        }

        let slice = arr.as_slice()?;
        self.inner.pos.clear();
        for i in 0..shape[0] {
            let idx = i * 3;
            self.inner
                .pos
                .push(RustVec3::new(slice[idx], slice[idx + 1], slice[idx + 2]));
        }
        Ok(())
    }

    /// Set velocities from NumPy array (N x 3)
    fn set_velocities(&mut self, arr: PyReadonlyArray2<f32>) -> PyResult<()> {
        let shape = arr.shape();
        if shape[1] != 3 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Array must have shape (N, 3)",
            ));
        }

        let slice = arr.as_slice()?;
        self.inner.vel.clear();
        for i in 0..shape[0] {
            let idx = i * 3;
            self.inner
                .vel
                .push(RustVec3::new(slice[idx], slice[idx + 1], slice[idx + 2]));
        }
        Ok(())
    }

    fn __repr__(&self) -> String {
        format!("State({} atoms)", self.num_atoms())
    }
}

/// System configuration (combines state, energy, topology)
#[pyclass(name = "Configuration")]
pub struct PyConfiguration {
    inner: RustConfiguration,
}

#[pymethods]
impl PyConfiguration {
    /// Create new configuration
    #[new]
    fn new(num_atoms: usize, num_temp_groups: usize, num_energy_groups: usize) -> Self {
        Self {
            inner: RustConfiguration::new(num_atoms, num_temp_groups, num_energy_groups),
        }
    }

    /// Get current state
    fn current_state(&self) -> PyState {
        PyState {
            inner: self.inner.current().clone(),
        }
    }

    /// Get current energy
    fn current_energy(&self) -> PyEnergy {
        PyEnergy {
            inner: self.inner.current().energies.clone(),
        }
    }

    fn __repr__(&self) -> String {
        format!("Configuration({} atoms)", self.inner.current().pos.len())
    }
}

/// Molecular topology (atoms, bonds, parameters)
#[pyclass(name = "Topology")]
pub struct PyTopology {
    inner: RustTopology,
}

#[pymethods]
impl PyTopology {
    /// Create new empty topology
    #[new]
    fn new() -> Self {
        Self {
            inner: RustTopology::new(),
        }
    }

    /// Number of atoms
    fn num_atoms(&self) -> usize {
        self.inner.solute.atoms.len()
    }

    /// Number of bonds
    fn num_bonds(&self) -> usize {
        self.inner.solute.bonds.len()
    }

    /// Number of angles
    fn num_angles(&self) -> usize {
        self.inner.solute.angles.len()
    }

    /// Number of dihedrals (proper)
    fn num_dihedrals(&self) -> usize {
        self.inner.solute.proper_dihedrals.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "Topology({} atoms, {} bonds)",
            self.num_atoms(),
            self.num_bonds()
        )
    }
}

//==============================================================================
// ANALYSIS FUNCTIONS
//==============================================================================

/// Calculate RMSD (Root Mean Square Deviation) for trajectory
///
/// Parameters
/// ----------
/// topology_file : str
///     Path to topology file
/// trajectory_file : str
///     Path to trajectory file
/// reference_frame : int
///     Reference frame index (0-based)
/// atom_selection : str, optional
///     Atom selection string (default: "all")
/// do_fit : bool
///     Perform rotational fit before RMSD calculation
///
/// Returns
/// -------
/// dict
///     Dictionary with 'times' and 'rmsd' arrays
#[pyfunction]
#[pyo3(signature = (topology_file, trajectory_file, reference_frame=0, atom_selection="all"))]
fn calculate_rmsd<'py>(
    py: Python<'py>,
    topology_file: &str,
    trajectory_file: &str,
    reference_frame: usize,
    atom_selection: &str,
) -> PyResult<&'py PyDict> {
    // Read topology
    let blocks = read_topology_file(topology_file).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading topology: {}", e))
    })?;
    let topo = build_topology(blocks);

    // Parse atom selection
    let selection = AtomSelection::from_string(atom_selection, &topo).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Error parsing selection: {}", e))
    })?;

    // Read trajectory
    let mut reader = TrajectoryReader::new(trajectory_file).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error opening trajectory: {}", e))
    })?;

    let frames = reader.read_all_frames().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading frames: {}", e))
    })?;

    if frames.is_empty() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "No frames in trajectory",
        ));
    }

    if reference_frame >= frames.len() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Reference frame {} out of range (0-{})",
            reference_frame,
            frames.len() - 1
        )));
    }

    // Get reference positions
    let reference_positions = &frames[reference_frame].positions;

    // Calculate RMSD for each frame
    let mut times = Vec::with_capacity(frames.len());
    let mut rmsds = Vec::with_capacity(frames.len());

    for frame in &frames {
        times.push(frame.time);

        // Calculate RMSD (simplified - no fitting for now)
        let rmsd =
            calculate_rmsd_simple(&frame.positions, reference_positions, selection.indices());
        rmsds.push(rmsd);
    }

    // Return as dictionary
    let result = PyDict::new(py);
    result.set_item("times", times.to_pyarray(py))?;
    result.set_item("rmsd", rmsds.to_pyarray(py))?;
    Ok(result)
}

/// Calculate RMSF (Root Mean Square Fluctuation) for trajectory
///
/// Parameters
/// ----------
/// topology_file : str
///     Path to topology file
/// trajectory_file : str
///     Path to trajectory file
/// atom_selection : str, optional
///     Atom selection string (default: "all")
/// skip_frames : int
///     Number of initial frames to skip
///
/// Returns
/// -------
/// dict
///     Dictionary with 'atom_indices' and 'rmsf' arrays
#[pyfunction]
#[pyo3(signature = (topology_file, trajectory_file, atom_selection="all", skip_frames=0))]
fn calculate_rmsf<'py>(
    py: Python<'py>,
    topology_file: &str,
    trajectory_file: &str,
    atom_selection: &str,
    skip_frames: usize,
) -> PyResult<&'py PyDict> {
    // Read topology
    let blocks = read_topology_file(topology_file).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading topology: {}", e))
    })?;
    let topo = build_topology(blocks);

    // Parse atom selection
    let selection = AtomSelection::from_string(atom_selection, &topo).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Error parsing selection: {}", e))
    })?;

    // Read trajectory
    let mut reader = TrajectoryReader::new(trajectory_file).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error opening trajectory: {}", e))
    })?;

    let frames = reader.read_all_frames().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading frames: {}", e))
    })?;

    if frames.len() <= skip_frames {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Not enough frames after skipping",
        ));
    }

    let analysis_frames = &frames[skip_frames..];

    // Calculate mean position for each atom
    let n_atoms = topo.num_atoms();
    let n_frames = analysis_frames.len() as f32;
    let mut mean_positions = vec![RustVec3::ZERO; n_atoms];

    for frame in analysis_frames {
        for (i, pos) in frame.positions.iter().enumerate() {
            mean_positions[i] = mean_positions[i] + *pos;
        }
    }

    for pos in &mut mean_positions {
        *pos = *pos / n_frames;
    }

    // Calculate RMSF for selected atoms
    let mut rmsf_values = Vec::new();
    let atom_indices: Vec<usize> = selection.indices().to_vec();

    for &atom_idx in &atom_indices {
        let mut sum_sq = 0.0f32;

        for frame in analysis_frames {
            if atom_idx < frame.positions.len() {
                let diff = frame.positions[atom_idx] - mean_positions[atom_idx];
                sum_sq += diff.length_squared();
            }
        }

        let rmsf = (sum_sq / n_frames).sqrt();
        rmsf_values.push(rmsf);
    }

    // Return as dictionary
    let result = PyDict::new(py);
    result.set_item("atom_indices", atom_indices.to_pyarray(py))?;
    result.set_item("rmsf", rmsf_values.to_pyarray(py))?;
    Ok(result)
}

/// Calculate radius of gyration for trajectory
///
/// Parameters
/// ----------
/// topology_file : str
///     Path to topology file
/// trajectory_file : str
///     Path to trajectory file
/// atom_selection : str, optional
///     Atom selection string (default: "all")
///
/// Returns
/// -------
/// dict
///     Dictionary with 'times' and 'rgyr' arrays
#[pyfunction]
#[pyo3(signature = (topology_file, trajectory_file, atom_selection="all"))]
fn calculate_rgyr<'py>(
    py: Python<'py>,
    topology_file: &str,
    trajectory_file: &str,
    atom_selection: &str,
) -> PyResult<&'py PyDict> {
    // Read topology
    let blocks = read_topology_file(topology_file).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading topology: {}", e))
    })?;
    let topo = build_topology(blocks);

    // Parse atom selection
    let selection = AtomSelection::from_string(atom_selection, &topo).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Error parsing selection: {}", e))
    })?;

    // Read trajectory
    let mut reader = TrajectoryReader::new(trajectory_file).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error opening trajectory: {}", e))
    })?;

    let frames = reader.read_all_frames().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading frames: {}", e))
    })?;

    // Calculate Rg for each frame
    let mut times = Vec::with_capacity(frames.len());
    let mut rgyr_values = Vec::with_capacity(frames.len());

    let atom_indices = selection.indices();

    for frame in &frames {
        times.push(frame.time);

        // Calculate center of mass
        let mut com = RustVec3::ZERO;
        let n_atoms = atom_indices.len() as f32;

        for &idx in atom_indices {
            if idx < frame.positions.len() {
                com = com + frame.positions[idx];
            }
        }
        com = com / n_atoms;

        // Calculate sum of squared distances from COM
        let mut sum_sq = 0.0f32;
        for &idx in atom_indices {
            if idx < frame.positions.len() {
                let diff = frame.positions[idx] - com;
                sum_sq += diff.length_squared();
            }
        }

        let rgyr = (sum_sq / n_atoms).sqrt();
        rgyr_values.push(rgyr);
    }

    // Return as dictionary
    let result = PyDict::new(py);
    result.set_item("times", times.to_pyarray(py))?;
    result.set_item("rgyr", rgyr_values.to_pyarray(py))?;
    Ok(result)
}

/// Analyze trajectory and return basic statistics
///
/// Parameters
/// ----------
/// topology_file : str
///     Path to topology file
/// trajectory_file : str
///     Path to trajectory file
///
/// Returns
/// -------
/// dict
///     Dictionary with trajectory statistics
#[pyfunction]
fn analyze_trajectory<'py>(
    py: Python<'py>,
    topology_file: &str,
    trajectory_file: &str,
) -> PyResult<&'py PyDict> {
    // Read topology
    let blocks = read_topology_file(topology_file).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading topology: {}", e))
    })?;
    let topo = build_topology(blocks);

    // Read trajectory
    let mut reader = TrajectoryReader::new(trajectory_file).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error opening trajectory: {}", e))
    })?;

    let frames = reader.read_all_frames().map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading frames: {}", e))
    })?;

    // Calculate statistics
    let n_frames = frames.len();
    let n_atoms = if n_frames > 0 {
        frames[0].positions.len()
    } else {
        0
    };
    let start_time = if n_frames > 0 {
        frames[0].time as f64
    } else {
        0.0
    };
    let end_time = if n_frames > 0 {
        frames[n_frames - 1].time as f64
    } else {
        0.0
    };
    let time_step = if n_frames > 1 {
        (end_time - start_time) / (n_frames as f64 - 1.0)
    } else {
        0.0
    };

    // Return as dictionary
    let result = PyDict::new(py);
    result.set_item("n_frames", n_frames)?;
    result.set_item("n_atoms", n_atoms)?;
    result.set_item("start_time", start_time)?;
    result.set_item("end_time", end_time)?;
    result.set_item("time_step", time_step)?;
    result.set_item("title", reader.title())?;
    Ok(result)
}

// Helper function for RMSD calculation
fn calculate_rmsd_simple(pos1: &[RustVec3], pos2: &[RustVec3], atom_indices: &[usize]) -> f32 {
    let indices = if atom_indices.is_empty() {
        (0..pos1.len()).collect::<Vec<_>>()
    } else {
        atom_indices.to_vec()
    };

    let n = indices.len();
    if n == 0 {
        return 0.0;
    }

    let mut sum_sq = 0.0f32;
    for &i in &indices {
        if i >= pos1.len() || i >= pos2.len() {
            continue;
        }
        let diff = pos1[i] - pos2[i];
        sum_sq += diff.length_squared();
    }

    (sum_sq / n as f32).sqrt()
}

//==============================================================================
// BINARY I/O READERS
//==============================================================================

/// Binary DCD trajectory reader (fast, lossless)
///
/// Reads DCD format trajectories with 30-60× faster performance than ASCII.
/// Compatible with CHARMM, NAMD, and OpenMM DCD files.
///
/// Examples
/// --------
/// >>> reader = gromos.DcdReader("trajectory.dcd")
/// >>> print(f"Frames: {reader.n_frames}, Atoms: {reader.n_atoms}")
/// >>> frame = reader.read_frame()
/// >>> print(frame['positions'].shape)  # (n_atoms, 3)
#[pyclass(name = "DcdReader")]
pub struct PyDcdReader {
    reader: DcdReader,
}

#[pymethods]
impl PyDcdReader {
    /// Open a DCD trajectory file
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     Path to DCD file
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let reader = DcdReader::new(path).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error opening DCD: {}", e))
        })?;
        Ok(Self { reader })
    }

    /// Total number of frames in trajectory
    #[getter]
    fn n_frames(&self) -> usize {
        self.reader.n_frames()
    }

    /// Number of atoms per frame
    #[getter]
    fn n_atoms(&self) -> usize {
        self.reader.n_atoms()
    }

    /// Read the next frame
    ///
    /// Returns
    /// -------
    /// dict or None
    ///     Dictionary with 'step', 'time', 'positions', 'box_dims' or None if EOF
    fn read_frame<'py>(&mut self, py: Python<'py>) -> PyResult<Option<&'py PyDict>> {
        match self.reader.read_frame().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading frame: {}", e))
        })? {
            Some(frame) => {
                let dict = PyDict::new(py);
                dict.set_item("step", frame.step)?;
                dict.set_item("time", frame.time)?;

                // Convert positions to NumPy array (n_atoms × 3)
                let n_atoms = frame.positions.len();
                let mut pos_data = Vec::with_capacity(n_atoms * 3);
                for v in &frame.positions {
                    pos_data.push(v.x);
                    pos_data.push(v.y);
                    pos_data.push(v.z);
                }
                let positions = PyArray2::from_vec2(py, &[pos_data])
                    .unwrap()
                    .reshape([n_atoms, 3])
                    .unwrap();
                dict.set_item("positions", positions)?;

                // Box dimensions
                let box_dims = [frame.box_dims.x, frame.box_dims.y, frame.box_dims.z];
                dict.set_item("box_dims", box_dims.to_pyarray(py))?;

                Ok(Some(dict))
            },
            None => Ok(None),
        }
    }

    /// Read all frames
    ///
    /// Returns
    /// -------
    /// list
    ///     List of frame dictionaries
    fn read_all_frames<'py>(&mut self, py: Python<'py>) -> PyResult<Vec<&'py PyDict>> {
        let mut frames = Vec::new();
        while let Some(frame) = self.read_frame(py)? {
            frames.push(frame);
        }
        Ok(frames)
    }

    /// Seek to a specific frame
    ///
    /// Parameters
    /// ----------
    /// frame : int
    ///     Frame index (0-based)
    fn seek_frame(&mut self, frame: usize) -> PyResult<()> {
        self.reader.seek_frame(frame).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error seeking: {}", e))
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "DcdReader({} frames, {} atoms)",
            self.n_frames(),
            self.n_atoms()
        )
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(
        mut slf: PyRefMut<'_, Self>,
        py: Python<'py>,
    ) -> PyResult<Option<&'py PyDict>> {
        slf.read_frame(py)
    }
}

/// Binary energy trajectory reader (fast, lossless)
///
/// Reads binary energy files with 20-40× faster performance than ASCII.
///
/// Examples
/// --------
/// >>> reader = gromos.BinaryEnergyReader("energy.tre.bin")
/// >>> print(f"Frames: {reader.n_frames}")
/// >>> frame = reader.read_frame()
/// >>> print(f"Total energy: {frame['total']:.2f} kJ/mol")
#[pyclass(name = "BinaryEnergyReader")]
pub struct PyBinaryEnergyReader {
    reader: BinaryEnergyReader,
}

#[pymethods]
impl PyBinaryEnergyReader {
    /// Open a binary energy file
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     Path to .tre.bin file
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let reader = BinaryEnergyReader::new(path).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "Error opening energy file: {}",
                e
            ))
        })?;
        Ok(Self { reader })
    }

    /// Total number of frames
    #[getter]
    fn n_frames(&self) -> usize {
        self.reader.n_frames()
    }

    /// Trajectory title
    #[getter]
    fn title(&self) -> &str {
        self.reader.title()
    }

    /// Read the next energy frame
    ///
    /// Returns
    /// -------
    /// dict or None
    ///     Dictionary with all energy components or None if EOF
    fn read_frame<'py>(&mut self, py: Python<'py>) -> PyResult<Option<&'py PyDict>> {
        match self.reader.read_frame().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading frame: {}", e))
        })? {
            Some(frame) => {
                let dict = PyDict::new(py);
                dict.set_item("time", frame.time)?;
                dict.set_item("kinetic", frame.kinetic)?;
                dict.set_item("potential", frame.potential)?;
                dict.set_item("total", frame.total)?;
                dict.set_item("temperature", frame.temperature)?;
                dict.set_item("volume", frame.volume)?;
                dict.set_item("pressure", frame.pressure)?;
                dict.set_item("bond", frame.bond)?;
                dict.set_item("angle", frame.angle)?;
                dict.set_item("improper", frame.improper)?;
                dict.set_item("dihedral", frame.dihedral)?;
                dict.set_item("lj", frame.lj)?;
                dict.set_item("coul_real", frame.coul_real)?;
                dict.set_item("coul_recip", frame.coul_recip)?;
                dict.set_item("coul_self", frame.coul_self)?;
                dict.set_item("shake", frame.shake)?;
                dict.set_item("restraint", frame.restraint)?;
                Ok(Some(dict))
            },
            None => Ok(None),
        }
    }

    /// Read all energy frames
    ///
    /// Returns
    /// -------
    /// list
    ///     List of energy frame dictionaries
    fn read_all_frames<'py>(&mut self, py: Python<'py>) -> PyResult<Vec<&'py PyDict>> {
        let mut frames = Vec::new();
        while let Some(frame) = self.read_frame(py)? {
            frames.push(frame);
        }
        Ok(frames)
    }

    /// Seek to a specific frame
    ///
    /// Parameters
    /// ----------
    /// frame : int
    ///     Frame index (0-based)
    fn seek_frame(&mut self, frame: usize) -> PyResult<()> {
        self.reader.seek_frame(frame).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error seeking: {}", e))
        })
    }

    fn __repr__(&self) -> String {
        format!("BinaryEnergyReader({} frames)", self.n_frames())
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(
        mut slf: PyRefMut<'_, Self>,
        py: Python<'py>,
    ) -> PyResult<Option<&'py PyDict>> {
        slf.read_frame(py)
    }
}

//==============================================================================
// CONVERSION FUNCTIONS
//==============================================================================

/// Convert trajectory between ASCII and binary formats
///
/// Parameters
/// ----------
/// input_file : str
///     Input trajectory file (.trc or .dcd)
/// output_file : str
///     Output trajectory file (.trc or .dcd)
/// topology_file : str, optional
///     Topology file for metadata
///
/// Examples
/// --------
/// >>> # ASCII to binary
/// >>> gromos.convert_trajectory("input.trc", "output.dcd")
///
/// >>> # Binary to ASCII
/// >>> gromos.convert_trajectory("input.dcd", "output.trc")
#[pyfunction]
#[pyo3(signature = (input_file, output_file, topology_file=None))]
fn convert_trajectory(
    input_file: &str,
    output_file: &str,
    topology_file: Option<&str>,
) -> PyResult<()> {
    use gromos_rs::io::{BinaryTrajectoryWriter, DcdWriter};
    use std::path::Path;

    let input_ext = Path::new(input_file)
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    let output_ext = Path::new(output_file)
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("");

    // Determine conversion direction
    match (input_ext, output_ext) {
        // ASCII to DCD
        ("trc" | "trj" | "g96", "dcd") => {
            let mut reader = TrajectoryReader::new(input_file).map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error opening input: {}", e))
            })?;

            let mut writer = DcdWriter::new(output_file, "Converted from ASCII").map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                    "Error creating output: {}",
                    e
                ))
            })?;

            let frames = reader.read_all_frames().map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading frames: {}", e))
            })?;

            // Convert frames
            for frame in &frames {
                // Create a temporary configuration for DCD writer
                let mut config = RustConfiguration::new(frame.positions.len(), 1, 1);
                config.current_mut().pos = frame.positions.clone();
                config.current_mut().box_config = gromos_rs::configuration::Box::rectangular(
                    frame.box_dims.x,
                    frame.box_dims.y,
                    frame.box_dims.z,
                );

                writer
                    .write_frame(frame.step, frame.time as f64, &config)
                    .map_err(|e| {
                        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                            "Error writing frame: {}",
                            e
                        ))
                    })?;
            }

            writer.finish().map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error finalizing: {}", e))
            })?;

            println!(
                "✓ Converted {} frames: {} → {}",
                frames.len(),
                input_file,
                output_file
            );
        },

        // DCD to ASCII
        ("dcd", "trc" | "trj") => {
            use gromos_rs::io::TrajectoryWriter;

            let mut reader = DcdReader::new(input_file).map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error opening input: {}", e))
            })?;

            let mut writer = TrajectoryWriter::new(output_file, "Converted from DCD", false, false)
                .map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                        "Error creating output: {}",
                        e
                    ))
                })?;

            let mut frame_count = 0;
            while let Some(frame) = reader.read_frame().map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading frame: {}", e))
            })? {
                // Create temporary configuration
                let mut config = RustConfiguration::new(frame.positions.len(), 1, 1);
                config.current_mut().pos = frame.positions;
                config.current_mut().box_config = gromos_rs::configuration::Box::rectangular(
                    frame.box_dims.x,
                    frame.box_dims.y,
                    frame.box_dims.z,
                );

                writer
                    .write_frame(frame.step, frame.time as f64, &config)
                    .map_err(|e| {
                        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                            "Error writing frame: {}",
                            e
                        ))
                    })?;

                frame_count += 1;
            }

            writer.flush().map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error flushing: {}", e))
            })?;

            println!(
                "✓ Converted {} frames: {} → {}",
                frame_count, input_file, output_file
            );
        },

        _ => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Unsupported conversion: {} → {}. Supported: .trc↔.dcd",
                input_ext, output_ext
            )));
        },
    }

    Ok(())
}

/// Convert energy trajectory between ASCII and binary formats
///
/// Parameters
/// ----------
/// input_file : str
///     Input energy file (.tre or .tre.bin)
/// output_file : str
///     Output energy file (.tre or .tre.bin)
///
/// Examples
/// --------
/// >>> # ASCII to binary
/// >>> gromos.convert_energy("input.tre", "output.tre.bin")
///
/// >>> # Binary to ASCII
/// >>> gromos.convert_energy("input.tre.bin", "output.tre")
#[pyfunction]
fn convert_energy(input_file: &str, output_file: &str) -> PyResult<()> {
    use gromos_rs::io::{BinaryEnergyWriter, EnergyReader, EnergyWriter};
    use std::path::Path;

    let input_ext = Path::new(input_file)
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    let output_name = Path::new(output_file)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    let is_output_binary = output_name.contains(".bin");

    // Determine format
    if input_ext == "tre" && !input_file.contains(".bin") && is_output_binary {
        // ASCII to binary
        let mut reader = EnergyReader::new(input_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error opening input: {}", e))
        })?;

        let mut writer = BinaryEnergyWriter::new(output_file, reader.title()).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error creating output: {}", e))
        })?;

        let frames = reader.read_all_frames().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading frames: {}", e))
        })?;

        for frame in &frames {
            writer.write_frame(frame).map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error writing frame: {}", e))
            })?;
        }

        writer.finish().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error finalizing: {}", e))
        })?;

        println!(
            "✓ Converted {} energy frames: {} → {}",
            frames.len(),
            input_file,
            output_file
        );
    } else if is_output_binary == false && input_file.contains(".bin") {
        // Binary to ASCII
        let mut reader = BinaryEnergyReader::new(input_file).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error opening input: {}", e))
        })?;

        let mut writer = EnergyWriter::new(output_file, reader.title()).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error creating output: {}", e))
        })?;

        let mut frame_count = 0;
        while let Some(frame) = reader.read_frame().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error reading frame: {}", e))
        })? {
            writer.write_frame(&frame).map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error writing frame: {}", e))
            })?;
            frame_count += 1;
        }

        writer.finalize().map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error finalizing: {}", e))
        })?;

        println!(
            "✓ Converted {} energy frames: {} → {}",
            frame_count, input_file, output_file
        );
    } else {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Unsupported conversion. Supported: .tre ↔ .tre.bin",
        ));
    }

    Ok(())
}

//==============================================================================
// NMR RESTRAINTS
//==============================================================================

/// J-value restraint for NMR refinement
///
/// Restrains J-coupling constants based on dihedral angles using the Karplus relation.
///
/// Examples
/// --------
/// >>> jval = gromos.JValueRestraint(
/// ...     i=0, j=1, k=2, l=3,
/// ...     karplus_a=6.4, karplus_b=-1.4, karplus_c=1.9,
/// ...     target_j=7.5, force_constant=10.0
/// ... )
/// >>> print(jval.j_current)
#[pyclass(name = "JValueRestraint")]
#[derive(Clone)]
pub struct PyJValueRestraint {
    inner: RustJValueRestraint,
}

#[pymethods]
impl PyJValueRestraint {
    #[new]
    #[pyo3(signature = (i, j, k_atom, l_atom, karplus_a, karplus_b, karplus_c, target_j, force_constant, delta=0.0, flat_bottom=0.0, half_harmonic=0, tau=0.0))]
    fn new(
        i: usize,
        j: usize,
        k_atom: usize,
        l_atom: usize,
        karplus_a: f64,
        karplus_b: f64,
        karplus_c: f64,
        target_j: f64,
        force_constant: f64,
        delta: f64,
        flat_bottom: f64,
        half_harmonic: i32,
        tau: f64,
    ) -> Self {
        Self {
            inner: RustJValueRestraint {
                i,
                j,
                k_atom,
                l_atom,
                a: karplus_a,
                b: karplus_b,
                c: karplus_c,
                delta,
                j0: target_j,
                k: force_constant,
                flat_bottom_width: flat_bottom,
                half_harmonic,
                tau,
                j_average: 0.0,
                j_current: 0.0,
            },
        }
    }

    #[getter]
    fn j_current(&self) -> f64 {
        self.inner.j_current
    }

    #[getter]
    fn j_average(&self) -> f64 {
        self.inner.j_average
    }

    fn __repr__(&self) -> String {
        format!(
            "JValueRestraint(atoms=[{},{},{},{}], J={:.2} Hz, J0={:.2} Hz)",
            self.inner.i,
            self.inner.j,
            self.inner.k_atom,
            self.inner.l_atom,
            self.inner.j_current,
            self.inner.j0
        )
    }
}

/// RDC (Residual Dipolar Coupling) restraint for NMR refinement
///
/// Restrains based on alignment tensor and internuclear vectors.
///
/// Examples
/// --------
/// >>> rdc = gromos.RDCRestraint(
/// ...     i=0, j=1,
/// ...     d_max=10000.0,
/// ...     target_rdc=5.2,
/// ...     force_constant=10.0,
/// ...     saupe_matrix=[0.1, 0.2, -0.3, 0.0, 0.0]
/// ... )
#[pyclass(name = "RDCRestraint")]
#[derive(Clone)]
pub struct PyRDCRestraint {
    inner: RustRDCRestraint,
}

#[pymethods]
impl PyRDCRestraint {
    #[new]
    #[pyo3(signature = (i, j, d_max, target_rdc, force_constant, saupe_matrix, tau=0.0, flat_bottom_width=0.0, rdc_type=0))]
    fn new(
        i: usize,
        j: usize,
        d_max: f64,
        target_rdc: f64,
        force_constant: f64,
        saupe_matrix: [f64; 5],
        tau: f64,
        flat_bottom_width: f64,
        rdc_type: i32,
    ) -> Self {
        Self {
            inner: RustRDCRestraint {
                i,
                j,
                d_max,
                d0: target_rdc,
                k: force_constant,
                flat_bottom_width,
                rdc_type,
                saupe_matrix,
                tau,
                rdc_average: 0.0,
                rdc_current: 0.0,
            },
        }
    }

    #[getter]
    fn rdc_current(&self) -> f64 {
        self.inner.rdc_current
    }

    #[getter]
    fn rdc_average(&self) -> f64 {
        self.inner.rdc_average
    }

    fn __repr__(&self) -> String {
        format!(
            "RDCRestraint(atoms=[{},{}], RDC={:.2}, RDC0={:.2})",
            self.inner.i, self.inner.j, self.inner.rdc_current, self.inner.d0
        )
    }
}

//==============================================================================
// VIRTUAL ATOMS
//==============================================================================

/// Virtual atom for coarse-grained models
///
/// Virtual atoms have no mass but are constructed from real atoms.
/// Used in united-atom models, TIP4P/TIP5P water, Martini force field, etc.
///
/// Examples
/// --------
/// >>> # Create TIP4P water virtual site (type 3 - planar)
/// >>> virt = gromos.VirtualAtom(
/// ...     atom_index=3,
/// ...     virt_type=3,
/// ...     parent_atoms=[0, 1, 2],  # O, H1, H2
/// ...     parameters=[0.15, 0.5],   # distance, weight
/// ...     masses=[16.0, 1.0, 1.0]
/// ... )
#[pyclass(name = "VirtualAtom")]
#[derive(Clone)]
pub struct PyVirtualAtom {
    inner: RustVirtualAtom,
}

#[pymethods]
impl PyVirtualAtom {
    #[new]
    fn new(
        atom_index: usize,
        virt_type: i32,
        parent_atoms: Vec<usize>,
        parameters: Vec<f64>,
        masses: Vec<f64>,
    ) -> Self {
        Self {
            inner: RustVirtualAtom {
                atom_index,
                virt_type,
                parent_atoms,
                parameters,
                masses,
            },
        }
    }

    #[getter]
    fn atom_index(&self) -> usize {
        self.inner.atom_index
    }

    #[getter]
    fn virt_type(&self) -> i32 {
        self.inner.virt_type
    }

    fn __repr__(&self) -> String {
        format!(
            "VirtualAtom(index={}, type={}, parents={:?})",
            self.inner.atom_index, self.inner.virt_type, self.inner.parent_atoms
        )
    }
}

/// Manager for all virtual atoms in a system
#[pyclass(name = "VirtualAtomManager")]
pub struct PyVirtualAtomsManager {
    inner: RustVirtualAtomsManager,
}

#[pymethods]
impl PyVirtualAtomsManager {
    #[new]
    fn new() -> Self {
        Self {
            inner: RustVirtualAtomsManager::new(),
        }
    }

    fn add_virtual_atom(&mut self, virt_atom: PyVirtualAtom) {
        self.inner.add_virtual_atom(virt_atom.inner.clone());
    }

    fn num_virtual_atoms(&self) -> usize {
        self.inner.virtual_atoms.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "VirtualAtomsManager({} virtual atoms)",
            self.num_virtual_atoms()
        )
    }
}

//==============================================================================
// LOCAL ELEVATION / METADYNAMICS
//==============================================================================

/// Coordinate type for local elevation
#[pyclass(name = "CoordinateType")]
#[derive(Clone)]
pub struct PyCoordinateType {
    inner: RustCoordinateType,
}

#[pymethods]
impl PyCoordinateType {
    #[new]
    fn new(coord_type: i32) -> Self {
        let inner = match coord_type {
            1 => RustCoordinateType::Distance,
            2 => RustCoordinateType::Angle,
            3 => RustCoordinateType::Dihedral,
            4 => RustCoordinateType::RMSD,
            5 => RustCoordinateType::RadiusOfGyration,
            _ => RustCoordinateType::Distance,
        };
        Self { inner }
    }

    /// Distance between two atoms
    #[staticmethod]
    fn distance() -> Self {
        Self {
            inner: RustCoordinateType::Distance,
        }
    }

    /// Angle between three atoms
    #[staticmethod]
    fn angle() -> Self {
        Self {
            inner: RustCoordinateType::Angle,
        }
    }

    /// Dihedral angle between four atoms
    #[staticmethod]
    fn dihedral() -> Self {
        Self {
            inner: RustCoordinateType::Dihedral,
        }
    }

    /// RMSD from reference structure
    #[staticmethod]
    fn rmsd() -> Self {
        Self {
            inner: RustCoordinateType::RMSD,
        }
    }

    /// Radius of gyration
    #[staticmethod]
    fn radius_of_gyration() -> Self {
        Self {
            inner: RustCoordinateType::RadiusOfGyration,
        }
    }

    fn __repr__(&self) -> String {
        format!("CoordinateType::{:?}", self.inner)
    }
}

/// Local elevation coordinate definition
#[pyclass(name = "LECoordinate")]
#[derive(Clone)]
pub struct PyLECoordinate {
    inner: RustLECoordinate,
}

#[pymethods]
impl PyLECoordinate {
    #[new]
    #[pyo3(signature = (umbrella_id, coord_type, atoms, reference_positions=None))]
    fn new(
        umbrella_id: usize,
        coord_type: PyCoordinateType,
        atoms: Vec<usize>,
        reference_positions: Option<Vec<Vec<f64>>>,
    ) -> Self {
        let ref_pos = match reference_positions {
            Some(positions) => positions
                .iter()
                .map(|v| RustVec3::new(v[0] as f32, v[1] as f32, v[2] as f32))
                .collect(),
            None => Vec::new(),
        };
        Self {
            inner: RustLECoordinate {
                umbrella_id,
                coord_type: coord_type.inner,
                atoms,
                reference_positions: ref_pos,
                value: 0.0,
                force: 0.0,
            },
        }
    }

    #[getter]
    fn value(&self) -> f64 {
        self.inner.value
    }

    #[getter]
    fn force(&self) -> f64 {
        self.inner.force
    }

    fn __repr__(&self) -> String {
        format!(
            "LECoordinate(umbrella={}, type={:?}, atoms={:?}, value={:.3})",
            self.inner.umbrella_id, self.inner.coord_type, self.inner.atoms, self.inner.value
        )
    }
}

/// Local elevation umbrella for metadynamics
///
/// Implements adaptive biasing force method with Gaussian hills.
///
/// Examples
/// --------
/// >>> # Create 1D umbrella for dihedral angle
/// >>> coord = gromos.LECoordinate(
/// ...     umbrella_id=1,
/// ...     coord_type=gromos.CoordinateType.dihedral(),
/// ...     atoms=[0, 1, 2, 3]
/// ... )
/// >>> umbrella = gromos.Umbrella(
/// ...     umbrella_id=1,
/// ...     coordinates=[coord],
/// ...     grid_sizes=[360],
/// ...     grid_mins=[-180.0],
/// ...     grid_maxs=[180.0],
/// ...     grid_spacings=[1.0],
/// ...     hill_height=0.5,
/// ...     gaussian_widths=[5.0],
/// ...     deposition_frequency=100
/// ... )
#[pyclass(name = "Umbrella")]
pub struct PyUmbrella {
    inner: RustUmbrella,
}

#[pymethods]
impl PyUmbrella {
    #[new]
    #[pyo3(signature = (umbrella_id, coordinates, grid_sizes, grid_mins, grid_maxs, grid_spacings, hill_height, gaussian_widths, deposition_frequency, enabled=true, building=true, periodic=None))]
    fn new(
        umbrella_id: usize,
        coordinates: Vec<PyLECoordinate>,
        grid_sizes: Vec<usize>,
        grid_mins: Vec<f64>,
        grid_maxs: Vec<f64>,
        grid_spacings: Vec<f64>,
        hill_height: f64,
        gaussian_widths: Vec<f64>,
        deposition_frequency: usize,
        enabled: bool,
        building: bool,
        periodic: Option<Vec<bool>>,
    ) -> Self {
        use gromos_rs::interaction::local_elevation::UmbrellaWeightMethod;

        let rust_coords: Vec<RustLECoordinate> =
            coordinates.iter().map(|c| c.inner.clone()).collect();
        let dimensionality = rust_coords.len();
        let total_grid_size: usize = grid_sizes.iter().product();

        let periodic_vec = periodic.unwrap_or_else(|| vec![false; dimensionality]);

        Self {
            inner: RustUmbrella {
                id: umbrella_id,
                dimensionality,
                coordinates: rust_coords,
                grid_sizes,
                grid_mins,
                grid_maxs,
                grid_spacings,
                potential_grid: vec![0.0; total_grid_size],
                visit_counts: vec![0; total_grid_size],
                gaussian_widths,
                hill_height,
                deposition_frequency,
                enabled,
                building,
                periodic: periodic_vec,
                weight_method: UmbrellaWeightMethod::NumberOfVisits,
            },
        }
    }

    #[getter]
    fn dimensionality(&self) -> usize {
        self.inner.dimensionality
    }

    #[getter]
    fn building(&self) -> bool {
        self.inner.building
    }

    #[setter]
    fn set_building(&mut self, value: bool) {
        self.inner.building = value;
    }

    fn __repr__(&self) -> String {
        format!(
            "Umbrella(id={}, dims={}, building={})",
            self.inner.id, self.inner.dimensionality, self.inner.building
        )
    }
}

//==============================================================================
// POLARIZATION
//==============================================================================

/// Polarization model type
#[pyclass(name = "PolarizationModel")]
#[derive(Clone)]
pub struct PyPolarizationModel {
    inner: RustPolarizationModel,
}

#[pymethods]
impl PyPolarizationModel {
    /// No polarization (fixed charges)
    #[staticmethod]
    fn none() -> Self {
        Self {
            inner: RustPolarizationModel::None,
        }
    }

    /// Point induced dipoles (classical Drude)
    #[staticmethod]
    fn point_dipole() -> Self {
        Self {
            inner: RustPolarizationModel::PointDipole,
        }
    }

    /// Shell model (Dick-Overhauser)
    #[staticmethod]
    fn shell_model() -> Self {
        Self {
            inner: RustPolarizationModel::ShellModel,
        }
    }

    /// Drude oscillators with auxiliary particles
    #[staticmethod]
    fn drude_oscillator() -> Self {
        Self {
            inner: RustPolarizationModel::DrudeOscillator,
        }
    }

    /// Fluctuating charges (QEq/EEM)
    #[staticmethod]
    fn fluctuating_charge() -> Self {
        Self {
            inner: RustPolarizationModel::FluctuatingCharge,
        }
    }

    fn __repr__(&self) -> String {
        format!("PolarizationModel::{:?}", self.inner)
    }
}

/// Polarizability parameters for an atom
#[pyclass(name = "PolarizabilityParameters")]
#[derive(Clone)]
pub struct PyPolarizabilityParameters {
    inner: RustPolarizabilityParameters,
}

#[pymethods]
impl PyPolarizabilityParameters {
    #[new]
    #[pyo3(signature = (alpha, thole_a=2.6, drude_charge=0.0, drude_k=0.0, drude_mass=0.4))]
    fn new(alpha: f64, thole_a: f64, drude_charge: f64, drude_k: f64, drude_mass: f64) -> Self {
        Self {
            inner: RustPolarizabilityParameters {
                alpha,
                thole_a,
                drude_charge,
                drude_k,
                drude_mass,
            },
        }
    }

    #[getter]
    fn alpha(&self) -> f64 {
        self.inner.alpha
    }

    #[getter]
    fn thole_a(&self) -> f64 {
        self.inner.thole_a
    }

    fn __repr__(&self) -> String {
        format!(
            "PolarizabilityParameters(α={:.4}, thole={:.2})",
            self.inner.alpha, self.inner.thole_a
        )
    }
}

/// Polarization force field calculator
///
/// Implements polarizable force fields with induced dipoles or Drude oscillators.
///
/// Examples
/// --------
/// >>> # Create polarization calculator
/// >>> calc = gromos.PolarizationCalculator(
/// ...     n_atoms=100,
/// ...     model=gromos.PolarizationModel.point_dipole()
/// ... )
/// >>> # Add polarizable atoms
/// >>> params = gromos.PolarizabilityParameters(alpha=1.0, thole_a=2.6)
/// >>> calc.add_atom(params)
#[pyclass(name = "PolarizationCalculator")]
pub struct PyPolarizationCalculator {
    inner: RustPolarizationCalculator,
}

#[pymethods]
impl PyPolarizationCalculator {
    #[new]
    fn new(n_atoms: usize, model: PyPolarizationModel) -> Self {
        Self {
            inner: RustPolarizationCalculator::new(n_atoms, model.inner),
        }
    }

    fn add_atom(&mut self, params: PyPolarizabilityParameters) {
        self.inner.polarizabilities.push(params.inner);
    }

    fn num_atoms(&self) -> usize {
        self.inner.polarizabilities.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "PolarizationCalculator({:?}, {} atoms, tol={:.1e})",
            self.inner.model,
            self.num_atoms(),
            self.inner.scf_tolerance
        )
    }
}

//==============================================================================
// QM/MM
//==============================================================================

/// QM method selection
#[pyclass(name = "QMMethod")]
#[derive(Clone)]
pub struct PyQMMethod {
    inner: RustQMMethod,
}

#[pymethods]
impl PyQMMethod {
    #[staticmethod]
    fn gfn1_xtb() -> Self {
        Self {
            inner: RustQMMethod::GFN1XTB,
        }
    }

    #[staticmethod]
    fn gfn2_xtb() -> Self {
        Self {
            inner: RustQMMethod::GFN2XTB,
        }
    }

    fn __repr__(&self) -> String {
        format!("QMMethod::{:?}", self.inner)
    }
}

/// QM/MM coupling scheme
#[pyclass(name = "CouplingScheme")]
#[derive(Clone)]
pub struct PyCouplingScheme {
    inner: RustCouplingScheme,
}

#[pymethods]
impl PyCouplingScheme {
    #[staticmethod]
    fn mechanical() -> Self {
        Self {
            inner: RustCouplingScheme::Mechanical,
        }
    }

    #[staticmethod]
    fn electrostatic() -> Self {
        Self {
            inner: RustCouplingScheme::Electrostatic,
        }
    }

    #[staticmethod]
    fn polarized() -> Self {
        Self {
            inner: RustCouplingScheme::Polarized,
        }
    }

    fn __repr__(&self) -> String {
        format!("CouplingScheme::{:?}", self.inner)
    }
}

/// QM region definition
#[pyclass(name = "QMRegion")]
#[derive(Clone)]
pub struct PyQMRegion {
    inner: RustQMRegion,
}

#[pymethods]
impl PyQMRegion {
    #[new]
    #[pyo3(signature = (atoms, charge=0, multiplicity=1))]
    fn new(atoms: Vec<usize>, charge: i32, multiplicity: i32) -> Self {
        Self {
            inner: RustQMRegion {
                atoms,
                charge,
                multiplicity,
                link_atoms: Vec::new(),
            },
        }
    }

    #[getter]
    fn atoms(&self) -> Vec<usize> {
        self.inner.atoms.clone()
    }

    #[getter]
    fn charge(&self) -> i32 {
        self.inner.charge
    }

    #[getter]
    fn multiplicity(&self) -> i32 {
        self.inner.multiplicity
    }

    fn __repr__(&self) -> String {
        format!(
            "QMRegion({} atoms, charge={}, mult={})",
            self.inner.atoms.len(),
            self.inner.charge,
            self.inner.multiplicity
        )
    }
}

/// QM/MM hybrid calculator
///
/// Combines quantum mechanics for reactive region with classical MD for environment.
///
/// Examples
/// --------
/// >>> # Define QM region (e.g., enzyme active site)
/// >>> qm_region = gromos.QMRegion(
/// ...     atoms=[10, 11, 12, 13, 14],  # Active site residue
/// ...     charge=0,
/// ...     multiplicity=1
/// ... )
/// >>> # Create QM/MM calculator with xTB
/// >>> qmmm = gromos.QMMMCalculator(
/// ...     qm_region=qm_region,
/// ...     qm_method=gromos.QMMethod.gfn2_xtb(),
/// ...     coupling=gromos.CouplingScheme.electrostatic(),
/// ...     mm_cutoff=1.2
/// ... )
#[pyclass(name = "QMMMCalculator")]
pub struct PyQMMMCalculator {
    atoms: Vec<usize>,
    charge: i32,
    multiplicity: i32,
    qm_method: String,
    coupling: String,
    mm_cutoff: f64,
}

#[pymethods]
impl PyQMMMCalculator {
    #[new]
    fn new(
        qm_region: PyQMRegion,
        qm_method: PyQMMethod,
        coupling: PyCouplingScheme,
        mm_cutoff: f64,
    ) -> Self {
        Self {
            atoms: qm_region.inner.atoms.clone(),
            charge: qm_region.inner.charge,
            multiplicity: qm_region.inner.multiplicity,
            qm_method: format!("{:?}", qm_method.inner),
            coupling: format!("{:?}", coupling.inner),
            mm_cutoff,
        }
    }

    #[getter]
    fn qm_atoms(&self) -> Vec<usize> {
        self.atoms.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "QMMMCalculator({} QM atoms, method={}, coupling={})",
            self.atoms.len(),
            self.qm_method,
            self.coupling
        )
    }
}

//==============================================================================
// FREE ENERGY PERTURBATION
//==============================================================================

/// Lambda controller for FEP simulations
///
/// Manages lambda coupling parameter(s) for free energy calculations.
///
/// Examples
/// --------
/// >>> # Create lambda controller for TI
/// >>> lambda_ctrl = gromos.LambdaController(
/// ...     lambda_value=0.5,
/// ...     lambda_exponent=1
/// ... )
/// >>> print(lambda_ctrl.lambda_value)
/// 0.5
#[pyclass(name = "LambdaController")]
#[derive(Clone)]
pub struct PyLambdaController {
    inner: RustLambdaController,
}

#[pymethods]
impl PyLambdaController {
    #[new]
    #[pyo3(signature = (lambda_value, lambda_exponent=1, dlambda_dt=0.0))]
    fn new(lambda_value: f64, lambda_exponent: i32, dlambda_dt: f64) -> Self {
        Self {
            inner: RustLambdaController {
                lambda: lambda_value,
                dlambda_dt,
                lambda_exponent,
                interaction_lambdas: RustInteractionLambdas::default(),
            },
        }
    }

    #[getter]
    fn lambda_value(&self) -> f64 {
        self.inner.lambda
    }

    #[setter]
    fn set_lambda_value(&mut self, value: f64) {
        self.inner.lambda = value;
    }

    #[getter]
    fn lambda_exponent(&self) -> i32 {
        self.inner.lambda_exponent
    }

    #[getter]
    fn dlambda_dt(&self) -> f64 {
        self.inner.dlambda_dt
    }

    #[setter]
    fn set_dlambda_dt(&mut self, value: f64) {
        self.inner.dlambda_dt = value;
    }

    fn lambda_derivative(&self) -> f64 {
        self.inner.lambda_derivative()
    }

    fn update(&mut self, dt: f64) {
        self.inner.update(dt);
    }

    fn __repr__(&self) -> String {
        format!(
            "LambdaController(λ={:.4}, n={}, dλ/dt={:.4})",
            self.inner.lambda, self.inner.lambda_exponent, self.inner.dlambda_dt
        )
    }
}

/// Perturbed atom for dual-topology FEP
///
/// Stores parameters for both state A and state B.
///
/// Examples
/// --------
/// >>> # Mutation: ASP- (charged) -> ASN (neutral)
/// >>> perturbed = gromos.PerturbedAtom(
/// ...     atom_index=10,
/// ...     a_charge=-1.0,  # ASP negative
/// ...     b_charge=0.0,   # ASN neutral
/// ...     a_iac=5,
/// ...     b_iac=6,
/// ...     lj_softcore=0.5,
/// ...     crf_softcore=0.5
/// ... )
/// >>> print(perturbed.charge_at_lambda(0.5))  # Interpolated charge
/// -0.5
#[pyclass(name = "PerturbedAtom")]
#[derive(Clone)]
pub struct PyPerturbedAtom {
    inner: RustPerturbedAtom,
}

#[pymethods]
impl PyPerturbedAtom {
    #[new]
    #[pyo3(signature = (atom_index, a_charge, b_charge, a_iac, b_iac, a_mass=0.0, b_mass=0.0, lj_softcore=0.0, crf_softcore=0.0))]
    fn new(
        atom_index: usize,
        a_charge: f64,
        b_charge: f64,
        a_iac: usize,
        b_iac: usize,
        a_mass: f64,
        b_mass: f64,
        lj_softcore: f64,
        crf_softcore: f64,
    ) -> Self {
        Self {
            inner: RustPerturbedAtom {
                atom_index,
                a_iac,
                a_mass,
                a_charge,
                b_iac,
                b_mass,
                b_charge,
                lj_softcore,
                crf_softcore,
            },
        }
    }

    #[getter]
    fn atom_index(&self) -> usize {
        self.inner.atom_index
    }

    #[getter]
    fn a_charge(&self) -> f64 {
        self.inner.a_charge
    }

    #[getter]
    fn b_charge(&self) -> f64 {
        self.inner.b_charge
    }

    #[getter]
    fn lj_softcore(&self) -> f64 {
        self.inner.lj_softcore
    }

    #[getter]
    fn crf_softcore(&self) -> f64 {
        self.inner.crf_softcore
    }

    /// Get interpolated charge at given lambda
    fn charge_at_lambda(&self, lambda: f64) -> f64 {
        self.inner.charge_at_lambda(lambda)
    }

    /// Get interpolated mass at given lambda
    fn mass_at_lambda(&self, lambda: f64) -> f64 {
        self.inner.mass_at_lambda(lambda)
    }

    fn __repr__(&self) -> String {
        format!(
            "PerturbedAtom(index={}, qA={:.2}, qB={:.2}, α_LJ={:.2}, α_CRF={:.2})",
            self.inner.atom_index,
            self.inner.a_charge,
            self.inner.b_charge,
            self.inner.lj_softcore,
            self.inner.crf_softcore
        )
    }
}

/// Interaction-specific lambda values
#[pyclass(name = "InteractionLambdas")]
#[derive(Clone)]
pub struct PyInteractionLambdas {
    inner: RustInteractionLambdas,
}

#[pymethods]
impl PyInteractionLambdas {
    #[new]
    #[pyo3(signature = (lj=1.0, lj_softness=1.0, crf=1.0, crf_softness=1.0, bond=1.0, angle=1.0, dihedral=1.0, improper=1.0, mass=1.0))]
    fn new(
        lj: f64,
        lj_softness: f64,
        crf: f64,
        crf_softness: f64,
        bond: f64,
        angle: f64,
        dihedral: f64,
        improper: f64,
        mass: f64,
    ) -> Self {
        Self {
            inner: RustInteractionLambdas {
                lj,
                lj_softness,
                crf,
                crf_softness,
                bond,
                angle,
                dihedral,
                improper,
                mass,
            },
        }
    }

    #[getter]
    fn lj(&self) -> f64 {
        self.inner.lj
    }

    #[getter]
    fn crf(&self) -> f64 {
        self.inner.crf
    }

    #[getter]
    fn bond(&self) -> f64 {
        self.inner.bond
    }

    #[getter]
    fn improper(&self) -> f64 {
        self.inner.improper
    }

    #[getter]
    fn mass(&self) -> f64 {
        self.inner.mass
    }

    fn __repr__(&self) -> String {
        format!(
            "InteractionLambdas(LJ={:.2}, CRF={:.2}, bond={:.2}, improper={:.2}, mass={:.2})",
            self.inner.lj, self.inner.crf, self.inner.bond, self.inner.improper, self.inner.mass
        )
    }
}

/// Perturbed nonbonded lambda parameters
#[pyclass(name = "PerturbedLambdaParams")]
#[derive(Clone)]
pub struct PyPerturbedLambdaParams {
    inner: RustPerturbedLambdaParams,
}

#[pymethods]
impl PyPerturbedLambdaParams {
    #[new]
    #[pyo3(signature = (lj_lambda, ljs_lambda, crf_lambda, crfs_lambda, lj_deriv=1.0, ljs_deriv=1.0, crf_deriv=1.0, crfs_deriv=1.0, lambda_exp=1))]
    fn new(
        lj_lambda: f64,
        ljs_lambda: f64,
        crf_lambda: f64,
        crfs_lambda: f64,
        lj_deriv: f64,
        ljs_deriv: f64,
        crf_deriv: f64,
        crfs_deriv: f64,
        lambda_exp: i32,
    ) -> Self {
        Self {
            inner: RustPerturbedLambdaParams::from_lambda(
                lj_lambda,
                ljs_lambda,
                crf_lambda,
                crfs_lambda,
                lj_deriv,
                ljs_deriv,
                crf_deriv,
                crfs_deriv,
                lambda_exp,
            ),
        }
    }

    #[getter]
    fn a_lj_lambda_n(&self) -> f64 {
        self.inner.a_lj_lambda_n
    }

    #[getter]
    fn b_lj_lambda_n(&self) -> f64 {
        self.inner.b_lj_lambda_n
    }

    fn __repr__(&self) -> String {
        format!(
            "PerturbedLambdaParams(λ_LJ^n=[{:.4}, {:.4}], λ_CRF^n=[{:.4}, {:.4}])",
            self.inner.a_lj_lambda_n,
            self.inner.b_lj_lambda_n,
            self.inner.a_crf_lambda_n,
            self.inner.b_crf_lambda_n
        )
    }
}

/// Calculate perturbed LJ+CRF interaction with soft-core
///
/// Parameters
/// ----------
/// r : Vec3
///     Distance vector between atoms i and j
/// a_c6, a_c12 : float
///     State A Lennard-Jones parameters
/// b_c6, b_c12 : float
///     State B Lennard-Jones parameters
/// a_q, b_q : float
///     State A and B charge products (qi*qj)
/// alpha_lj, alpha_crf : float
///     Soft-core parameters for LJ and electrostatics
/// lambda_params : PerturbedLambdaParams
///     Lambda values and derivatives
/// crf_cut, crf_2cut3i, crf_cut3i : float
///     CRF parameters
///
/// Returns
/// -------
/// tuple
///     (force_magnitude, e_lj, e_crf, de_lj, de_crf)
///     Force magnitude, LJ energy, CRF energy, and lambda derivatives
///
/// Examples
/// --------
/// >>> # Calculate interaction for charge appearing
/// >>> r = gromos.Vec3(0.3, 0.0, 0.0)  # 0.3 nm separation
/// >>> a_c6, a_c12 = 0.0, 0.0  # Dummy state
/// >>> b_c6, b_c12 = 0.001, 0.0001  # Real particle
/// >>> a_q, b_q = 0.0, 0.25  # No charge -> +0.5e charge product
/// >>> alpha_lj, alpha_crf = 0.5, 0.5  # Soft-core parameters
/// >>> lambda_params = gromos.PerturbedLambdaParams(0.5, 0.5, 0.5, 0.5, lambda_exp=1)
/// >>> force, e_lj, e_crf, de_lj, de_crf = gromos.perturbed_lj_crf_interaction(
/// ...     r, a_c6, a_c12, b_c6, b_c12, a_q, b_q,
/// ...     alpha_lj, alpha_crf, lambda_params,
/// ...     crf_cut=1.4, crf_2cut3i=0.182, crf_cut3i=0.13
/// ... )
/// >>> print(f"Energy: LJ={e_lj:.2f}, CRF={e_crf:.2f} kJ/mol")
/// >>> print(f"dH/dλ: LJ={de_lj:.2f}, CRF={de_crf:.2f} kJ/mol")
#[pyfunction]
fn py_perturbed_lj_crf_interaction(
    r: PyVec3,
    a_c6: f64,
    a_c12: f64,
    b_c6: f64,
    b_c12: f64,
    a_q: f64,
    b_q: f64,
    alpha_lj: f64,
    alpha_crf: f64,
    lambda_params: PyPerturbedLambdaParams,
    crf_cut: f64,
    crf_2cut3i: f64,
    crf_cut3i: f64,
) -> (f64, f64, f64, f64, f64) {
    use gromos_rs::interaction::nonbonded::CRFParameters;

    let crf = CRFParameters {
        crf_cut,
        crf_2cut3i,
        crf_cut3i,
    };

    perturbed_lj_crf_interaction(
        r.inner,
        a_c6,
        a_c12,
        b_c6,
        b_c12,
        a_q,
        b_q,
        alpha_lj,
        alpha_crf,
        &lambda_params.inner,
        &crf,
    )
}

//==============================================================================
// NOTE: Integrators and advanced sampling methods (GaMD, EDS, REMD) are
// available through the GROMOS-RS command-line binaries.
// Python bindings focus on data structures, I/O, and analysis tools.
//==============================================================================
