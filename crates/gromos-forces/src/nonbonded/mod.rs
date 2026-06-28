//! Nonbonded interactions: Lennard-Jones and Coulomb (CRF)

pub mod excluded;
pub mod innerloops;
pub mod params;
pub mod perturbed;

pub use excluded::*;
pub use innerloops::*;
pub use params::*;
pub use perturbed::*;

/// MPI-parallel nonbonded force calculation
///
/// # Status
/// Stub — pending a dedicated `gromos-mpi` crate that holds `MpiControl` /
/// `MpiNonbondedMaster`. `gromos-forces` cannot import from `gromos-integrators`
/// (would be circular). Wire up when the MPI support crate is extracted.
///
/// Ported from `md++/src/interaction/nonbonded/interaction/mpi_nonbonded_master.cc`
#[cfg(feature = "use-mpi")]
pub mod mpi_nonbonded {
    use super::*;

    /// Calculate nonbonded interactions with MPI parallelisation (stub).
    ///
    /// Signature mirrors `MPI_Nonbonded_Master::calculate_interactions`.
    /// Returns a zeroed `ForceStorage` until the MPI layer is wired.
    pub fn calculate_mpi<BC: gromos_core::math::BoundaryCondition>(
        _positions: &[gromos_core::math::Vec3],
        _charges: &[f64],
        _iac: &[u32],
        _pairlist: &[(u32, u32)],
        _lj_params: &LJParamMatrix,
        _crf: &CRFParameters,
        _periodicity: &BC,
    ) -> ForceStorage {
        // TODO: rewire — needs gromos-mpi crate with MpiControl + MpiNonbondedMaster
        ForceStorage::new(_positions.len())
    }
}
