//! Physical constants and unit system for GROMOS-RS.
//!
//! # Unit system
//!
//! All quantities are expressed in the GROMOS internal unit system:
//!
//! | Quantity    | Unit              |
//! |-------------|-------------------|
//! | Length      | nm                |
//! | Time        | ps                |
//! | Mass        | u (amu)           |
//! | Energy      | kJ mol⁻¹          |
//! | Charge      | e                 |
//! | Temperature | K                 |
//! | Force       | kJ mol⁻¹ nm⁻¹    |
//!
//! # Constants: gromosXX defaults
//!
//! The constants here are the **fallback defaults** used when no topology
//! PHYSICALCONSTANTS block is present (e.g. in analysis tools or unit tests).
//! They match the values hardcoded in gromosXX `math/math.cc`.
//!
//! For simulations the values come from the topology file — `gromos-io` parses
//! the PHYSICALCONSTANTS block and sets `PhysicalConstants` accordingly, which
//! `Forcefield` then uses.  The constants here are never the authoritative
//! source during an MD run.
//!
//! CODATA 2018 exact values for reference:
//!   four_pi_eps_i = 138.9354859  kJ mol⁻¹ nm e⁻²
//!   kB            = 0.008314462618  kJ mol⁻¹ K⁻¹
//!
//! # Naming convention
//!
//! Established physics symbols (`kB`, `hBar`) use scientific notation
//! rather than Rust's `SCREAMING_SNAKE_CASE`. The lint is suppressed at
//! file scope only — it does not propagate to importers.

// Physics symbols follow scientific notation, not SCREAMING_SNAKE_CASE.
#![allow(non_upper_case_globals)]

// ── Physical constants (gromosXX defaults) ───────────────────────────────────

/// Boltzmann constant (kJ mol⁻¹ K⁻¹).
///
/// gromosXX `math/math.cc`: `k_Boltzmann = 0.00831441`.
/// CODATA 2018 exact: 0.008314462618.
pub const kB: f64 = 0.00831441;

/// Coulomb prefactor 1/(4πε₀) (kJ mol⁻¹ nm e⁻²).
///
/// gromosXX `math/math.cc`: `four_pi_eps_i = 138.9354`.
/// CODATA 2018 exact: 138.9354859.
pub const four_pi_eps_i: f64 = 138.9354;

/// ε₀⁻¹ in GROMOS units (kJ mol⁻¹ nm e⁻²).
///
/// Derived as 4π × `four_pi_eps_i`.
pub const eps0_i: f64 = 4.0 * std::f64::consts::PI * four_pi_eps_i;

/// Reduced Planck constant ħ (kJ ps mol⁻¹).
///
/// gromosXX `math/math.cc`: `h_bar = 0.0635078`.
pub const hBar: f64 = 0.0635078;

/// Speed of light (nm ps⁻¹).
///
/// gromosXX `math/math.cc`: `spd_l = 299792.458`.
pub const spd_l: f64 = 299792.458;

// ── Unit conversions ─────────────────────────────────────────────────────────

/// e·nm → Debye.
pub const ENM_TO_DEBYE: f64 = 48.0321;

/// kJ mol⁻¹ → kcal mol⁻¹.
pub const KJ_TO_KCAL: f64 = 0.239006;

/// kcal mol⁻¹ → kJ mol⁻¹.
pub const KCAL_TO_KJ: f64 = 4.184;

/// nm → Å.
pub const NM_TO_ANGSTROM: f64 = 10.0;

/// Å → nm.
pub const ANGSTROM_TO_NM: f64 = 0.1;

/// ps → fs.
pub const PS_TO_FS: f64 = 1000.0;

/// fs → ps.
pub const FS_TO_PS: f64 = 0.001;

// ── Physical constants struct ─────────────────────────────────────────────────

/// Bundle of physical constants for explicit dependency injection.
///
/// The primary source is the topology PHYSICALCONSTANTS block, parsed by
/// `gromos-io` and stored here.  These values flow into `Forcefield` and
/// every force/thermostat function that needs them.
///
/// The module-level consts (`four_pi_eps_i`, `kB`) are used as fallback
/// defaults when constructing this struct without a topology.
#[derive(Debug, Clone, Copy)]
#[allow(non_snake_case)]
pub struct PhysicalConstants {
    /// Coulomb prefactor 1/(4πε₀) in GROMOS units (kJ mol⁻¹ nm e⁻²).
    pub four_pi_eps_i: f64,
    /// Boltzmann constant (kJ mol⁻¹ K⁻¹).
    pub kB: f64,
    /// Reduced Planck constant ħ (kJ ps mol⁻¹).
    pub hBar: f64,
    /// Speed of light (nm ps⁻¹).
    pub spd_l: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Guard that the default constants are exactly the gromosXX legacy values.
    ///
    /// CODATA 2018 references (for human audit):
    ///   four_pi_eps_i  138.9354859   (we use 138.9354,   Δ ~ 6e-7)
    ///   kB             0.008314462618 (we use 0.00831441,  Δ ~ 6e-6)
    ///   hBar           0.06350779927  (we use 0.0635078,   Δ ~ 1e-8)
    ///   spd_l          299792.458     (exact, same)
    #[test]
    fn constants_are_gromosxx_legacy_values() {
        assert_eq!(
            four_pi_eps_i, 138.9354,
            "four_pi_eps_i must be the gromosXX legacy value"
        );
        assert_eq!(kB, 0.00831441, "kB must be the gromosXX legacy value");
        assert_eq!(hBar, 0.0635078, "hBar must be the gromosXX legacy value");
        assert_eq!(spd_l, 299792.458, "spd_l must be the gromosXX legacy value");
    }
}
