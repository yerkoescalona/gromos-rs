//! The one error type of the run assembly.
//!
//! Every failure the `md` binary used to report with `eprintln!` + `process::exit(1)` is a
//! variant here, so the binary and the Python binding render the same condition the same way
//! (PLAN.md 3.9 C10). Python maps variants onto builtin exception types in `pyo3-gromos`.

use std::fmt;

use gromos_core::validation::ValidationReport;
use gromos_io::IoError;

pub enum RunError {
    /// Reading an input file failed.
    Io { what: String, source: IoError },
    /// The parameters ask for an input file that was not given
    /// (e.g. `NTPOR=1` without a position-restraint specification).
    MissingInput { flag: &'static str, reason: String },
    /// Topology and coordinate file disagree on the number of atoms.
    AtomCountMismatch { topology: usize, coordinates: usize },
    /// The coordinate file's non-solute atoms are not a whole number of solvent molecules.
    SolventCount {
        coordinates: usize,
        solute: usize,
        atoms_per_solvent: usize,
    },
    /// A validation stage found fatal problems; the report carries them.
    Validation {
        stage: &'static str,
        report: ValidationReport,
    },
    /// Inputs contradict each other (e.g. a perturbed topology with `NTG=0`).
    Inconsistent(String),
    /// `AlgorithmSequence::init` failed.
    Init(String),
    /// `AlgorithmSequence::run_step` failed at `step`.
    Step { step: usize, message: String },
}

impl fmt::Display for RunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RunError::Io { what, source } => write!(f, "error reading {what}: {source}"),
            RunError::MissingInput { flag, reason } => {
                write!(f, "{reason} (no @{flag} file specified)")
            },
            RunError::AtomCountMismatch {
                topology,
                coordinates,
            } => write!(
                f,
                "number of atoms in topology ({topology}) != coordinates ({coordinates})"
            ),
            RunError::SolventCount {
                coordinates,
                solute,
                atoms_per_solvent,
            } => write!(
                f,
                "({coordinates} coords - {solute} solute) is not divisible by \
                 {atoms_per_solvent} atoms/solvent"
            ),
            RunError::Validation { stage, report } => write!(
                f,
                "fatal errors in {stage} validation ({} finding(s)) - cannot continue",
                report.errors.len()
            ),
            RunError::Inconsistent(msg) => write!(f, "inconsistent inputs: {msg}"),
            RunError::Init(msg) => write!(f, "error initializing algorithm sequence: {msg}"),
            RunError::Step { step, message } => write!(f, "error at step {step}: {message}"),
        }
    }
}

impl fmt::Debug for RunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "RunError({self})")
    }
}

impl std::error::Error for RunError {}
