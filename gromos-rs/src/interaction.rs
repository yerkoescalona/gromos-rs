//! Interaction calculations (nonbonded, bonded, electrostatics, restraints)

pub mod bonded;
pub mod electrostatics;
pub mod local_elevation;
pub mod nonbonded;
pub mod pme_mpi;
pub mod polarization;
pub mod qmmm;
pub mod restraints;

// Re-export commonly used types and functions
pub use bonded::{
    calculate_angle_forces, calculate_bond_forces_harmonic, calculate_bond_forces_quartic,
    calculate_bonded_forces, calculate_dihedral_forces, calculate_improper_dihedral_forces,
    calculate_perturbed_angle_forces, calculate_perturbed_bond_forces,
    calculate_perturbed_dihedral_forces, ForceEnergy, ForceEnergyLambda,
};

pub use electrostatics::{
    pme_real_space_interaction, pme_reciprocal_space, pme_self_energy, reaction_field_interaction,
    ElectrostaticsMethod, PMEParameters, ReactionFieldParameters,
};

pub use restraints::{
    AngleRestraint, AngleRestraints, DihedralRestraint, DihedralRestraints, DistanceRestraint,
    DistanceRestraints, JValueRestraint, JValueRestraints, PositionRestraint, PositionRestraints,
    RDCRestraint, RDCRestraints,
};

pub use local_elevation::{
    CoordinateType, LECoordinate, LocalElevation, Umbrella, UmbrellaWeightMethod,
};

pub use polarization::{
    PolarizabilityParameters, PolarizationCalculator, PolarizationModel, PolarizationState,
    PolarizationTopology,
};

pub use qmmm::{
    CouplingScheme, DFTBPlusEngine, LinkAtom, QMEngine, QMEngineInterface, QMMMCalculator,
    QMMethod, QMRegion, QMResult, XTBEngine,
};
