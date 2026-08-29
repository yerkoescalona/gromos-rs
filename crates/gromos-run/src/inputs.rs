//! Auxiliary run inputs and execution options.
//!
//! An `.imd` file describes a run only together with its auxiliary files (PLAN.md 3.9 A5);
//! [`RunInputs`] carries those. [`RunOptions`] carries what is *not* physics and therefore has
//! no IMD field (A15, `execution` group of the future recipe): the parallel policy.

use std::path::PathBuf;

/// Paths of the optional GROMOS input files — the `md` binary's `@pttopo`, `@posresspec`,
/// `@refpos` and `@distrest` flags.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RunInputs {
    /// Perturbation topology (`.ptp`); only read when `NTG != 0`.
    pub pttopo: Option<PathBuf>,
    /// Position-restraint atom specification (`POSRESSPEC`); required when `NTPOR > 0`.
    pub posresspec: Option<PathBuf>,
    /// Reference positions for position restraints; required when `NTPORB >= 1`.
    pub refpos: Option<PathBuf>,
    /// Distance-restraint specification; required when `NTDIR != 0`.
    pub distrest: Option<PathBuf>,
}

/// Whether the nonbonded kernels run on the rayon pool.
///
/// A platform concern, not physics: the parallel kernels are reference-verified at the suite
/// tolerance (1e-8) and bit-identical run-to-run at a fixed thread count, but differ from the
/// serial kernels in the last digit (BENCHMARKING.md, "Kernel determinism").
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ParallelPolicy {
    /// What the `md` binary has always done: parallel above 100 atoms.
    #[default]
    Auto,
    /// Never parallel — the setting for exact (`==`) cross-front-end comparisons.
    Serial,
    /// Always parallel.
    Parallel,
}

impl ParallelPolicy {
    /// Resolve the policy for a system of `n_atoms` atoms.
    pub fn resolve(self, n_atoms: usize) -> bool {
        match self {
            ParallelPolicy::Auto => n_atoms > 100,
            ParallelPolicy::Serial => false,
            ParallelPolicy::Parallel => true,
        }
    }
}

/// Execution options that are not part of the physics description.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RunOptions {
    pub parallel: ParallelPolicy,
}
