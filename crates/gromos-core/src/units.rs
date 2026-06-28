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
//! # Source
//!
//! gromosXX defines these in `math/math.cc` as **mutable globals** with
//! truncated values. gromos-rs fixes both problems: `const` items with
//! CODATA 2018 precision. Changed values are noted per constant.
//!
//! # Naming convention
//!
//! Established physics symbols (`kB`, `hBar`) use scientific notation
//! rather than Rust's `SCREAMING_SNAKE_CASE`. The lint is suppressed at
//! file scope only — it does not propagate to importers.

// Physics symbols follow scientific notation, not SCREAMING_SNAKE_CASE.
#![allow(non_upper_case_globals)]

// ── Physical constants (CODATA 2018) ────────────────────────────────────────

/// Boltzmann constant (kJ mol⁻¹ K⁻¹).
///
/// gromosXX `k_Boltzmann = 0.00831441` — off by ~5 ppm; corrected here.
pub const kB: f64 = 0.008314462618;

/// Coulomb prefactor 1/(4πε₀) (kJ mol⁻¹ nm e⁻²).
///
/// gromosXX `four_pi_eps_i = 138.9354` — truncated; corrected here.
pub const four_pi_eps_i: f64 = 138.9354859;

/// ε₀⁻¹ in GROMOS units (kJ mol⁻¹ nm e⁻²).
///
/// Defined as 4π × `four_pi_eps_i`. gromosXX `eps0_i = 1745.9137`.
pub const eps0_i: f64 = 4.0 * std::f64::consts::PI * four_pi_eps_i;

/// Reduced Planck constant ħ (kJ ps mol⁻¹).
///
/// gromosXX `h_bar = 0.0635078`.
pub const hBar: f64 = 0.0635078;

/// Speed of light (nm ps⁻¹).
///
/// gromosXX `spd_l = 299792.458`.
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
