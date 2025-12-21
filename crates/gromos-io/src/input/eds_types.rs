//! EDS type definitions for I/O
//! 
//! These types are used for parsing EDS input blocks.
//! The actual simulation types are in gromos-integrators.

/// EDS functional form
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EDSForm {
    /// Single smoothing parameter s
    SingleS,
    /// Multiple smoothing parameters (multi-s EDS)
    MultiS,
    /// Pair-based s values
    PairS,
    /// Accelerated EDS with adaptive parameters
    AEDS,
}

impl Default for EDSForm {
    fn default() -> Self {
        EDSForm::SingleS
    }
}

/// EDS parameters from input file
#[derive(Debug, Clone, Default)]
pub struct EDSParameters {
    pub num_states: usize,
    pub form: EDSForm,
    pub s_values: Vec<f64>,
    pub e_offsets: Vec<f64>,
    pub temperature: f64,
}

/// Accelerated EDS parameters
#[derive(Debug, Clone, Default)]
pub struct AEDSParameters {
    pub search_enabled: bool,
    pub e_max: f64,
    pub e_min: f64,
    pub search_interval: u64,
}
