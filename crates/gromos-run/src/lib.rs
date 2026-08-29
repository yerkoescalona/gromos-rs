//! gromos-run — the one place that turns "what to run" into a running `AlgorithmSequence`.
//!
//! L3 orchestration (see `.claude/architecture.md`). Before this crate existed, the `md`
//! binary and the Python binding each carried their own copy of the IMD → algorithm-sequence
//! assembly and had drifted (PLAN.md 3.8). Both now call the functions here, so the gromosXX
//! reference suite (which drives the binary) and the Python suite (which drives the binding)
//! exercise the same code from two sides — drift guard G1 of PLAN.md 3.9.
//!
//! Stage of the plan: **step 1** — the assembly is lifted verbatim from `md.rs` and still
//! IMD-driven (`ImdParameters` in, `AlgorithmSequence` out). Step 2 introduces the plain-data
//! `RunRecipe` → `Vec<AlgorithmSpec>` plan in front of it; nothing in the public surface here
//! is meant to survive that unchanged except [`prepare_system`], [`start`] and [`RunError`].
//!
//! Pipeline:
//!
//! ```text
//! ImdParameters + Topology + Coordinates + RunInputs
//!        │ prepare_system      (pttopo merge, truncated-octahedron transform, NSM from the
//!        ▼                      coordinate file, initial velocities, validation)
//!     Prepared { topology, configuration, physical_constants, … }
//!        │ build_sequence_from_imd  (GROMOS step order, once)
//!        ▼
//!     Built { sequence: AlgorithmSequence, summary }
//!        │ start               (init + step 0 — NTISHK moves positions here, so it is physics)
//!        ▼
//!     caller loops run_step / advance
//! ```
//!
//! Rules (crate contract):
//! - No `println!`, no `process::exit`: everything the binary prints comes back as data
//!   ([`Prepared::notes`], [`BuildSummary`]) or as a [`RunError`].
//! - No second builder: if a caller needs a different sequence, it edits the plan (step 2),
//!   it does not assemble algorithms itself.

pub mod build;
pub mod bundle;
pub mod constraints;
pub mod dof;
pub mod error;
pub mod inputs;
#[cfg(test)]
mod legacy_builder;
pub mod plan;
pub mod prepare;
pub mod recipe;

pub use build::{
    build_sequence_from_imd, build_sequence_from_plan, build_sequence_from_recipe, instantiate,
    periodicity_of, start, BarostatSummary, BuildSummary, Built, Instantiated, ThermostatSummary,
};
pub use bundle::{load_bundle, read_imd, write_bundle, RunBundle};
pub use constraints::ConstraintSelection;
pub use dof::total_dof;
pub use error::RunError;
pub use inputs::{ParallelPolicy, RunInputs, RunOptions};
pub use plan::{build_plan, plan_from_json, plan_to_json, validate_plan, AlgorithmSpec};
pub use prepare::{prepare_system, Coordinates, PrepareNote, Prepared};
pub use recipe::{Diagnostics, PassthroughPolicy, RunRecipe, TermSpec, RECIPE_VERSION};
