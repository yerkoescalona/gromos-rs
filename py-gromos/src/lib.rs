//! Python extension module for GROMOS-RS
//!
//! This is the main entry point for the Python `gromos` package.
//! It uses the pyo3-gromos crate for the low-level bindings.

use pyo3::prelude::*;

// Re-export bindings from pyo3-gromos
use pyo3_gromos::register_bindings;

/// GROMOS-RS Python Module
///
/// High-performance molecular dynamics engine written in Rust.
#[pymodule]
fn gromos(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Register all bindings from pyo3-gromos
    register_bindings(m)?;
    
    // Add module metadata
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__author__", "GROMOS Developers")?;
    
    Ok(())
}
