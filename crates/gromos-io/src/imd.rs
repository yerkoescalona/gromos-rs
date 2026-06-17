//! Parser for GROMOS input parameter files (.imd)
//!
//! The IMD format uses block-based keywords to specify simulation parameters.
//!
//! # Format Example
//! ```text
//! TITLE
//!   My simulation
//! END
//! SYSTEM
//!   NPM     1
//!   NSM     0
//! END
//! STEP
//!   NSTLIM  1000
//!   T       0.0
//!   DT      0.002
//! END
//! ```
//!
//! # Major Blocks
//! - TITLE: Simulation description
//! - SYSTEM: System composition (NPM, NSM)
//! - STEP: Integration parameters (NSTLIM, DT, T)
//! - BOUNDCOND: Boundary conditions (NTB, NDFMIN)
//! - MULTIBATH: Temperature coupling groups
//! - PRESSURESCALE: Pressure coupling
//! - FORCE: Force field terms to compute
//! - CONSTRAINT: Constraint algorithm settings
//! - PAIRLIST: Neighbor list parameters
//! - NONBONDED: Cutoffs and long-range electrostatics
//! - INITIALISE: Initial velocities
//! - WRITETRAJ: Trajectory output settings
//! - PRINTOUT: Energy output settings

use crate::IoError;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// GROMOS simulation parameters from .imd file
#[derive(Debug, Clone)]
pub struct ImdParameters {
    /// Simulation title
    pub title: String,

    // SYSTEM block
    pub npm: usize, // Number of (identical) protein molecules
    pub nsm: usize, // Number of (identical) solvent molecules

    // STEP block
    pub nstlim: usize, // Number of MD steps
    pub t0: f64,       // Initial time (ps)
    pub dt: f64,       // Time step (ps)

    // BOUNDCOND block
    pub ntb: i32, // Boundary type (0=vacuum, 1=rectangular, 2=triclinic, -1=truncated octahedron)
    pub ndfmin: i32, // Minimum number of degrees of freedom

    // MULTIBATH block
    pub num_temp_baths: usize,
    pub temp_bath: Vec<TempBathParameters>,

    // PRESSURESCALE block
    pub couple_pressure: bool,
    pub pressure_parameters: Option<PressureParameters>,

    // FORCE block
    pub force_groups: Vec<Vec<(usize, usize)>>, // Energy group pairs

    // CONSTRAINT block
    pub ntc: i32,       // SHAKE constraints (1=none, 2=H-bonds, 3=all bonds, 4=all)
    pub ntcp: i32,      // P-SHAKE (pressure-SHAKE)
    pub ntcs: i32,      // Solvent SHAKE/SETTLE
    pub shake_tol: f64, // SHAKE tolerance
    pub lincs_order_solute: usize, // LINCS expansion order for solute (NTCP0 when NTCP=lincs)
    pub lincs_order_solvent: usize, // LINCS expansion order for solvent (NTCS0 when NTCS=lincs)

    // PAIRLIST block
    pub algorithm: i32, // Pairlist algorithm
    pub nsnb: usize,    // Update frequency
    pub rcutp: f64,     // Short-range cutoff (nm)
    pub rcutl: f64,     // Long-range cutoff (nm)
    pub size: f64,      // Grid cell size (nm)
    pub type_: i32,     // Pairlist type

    // NONBONDED block
    pub nlrele: i32,   // Long-range electrostatics (0=cutoff, 1=RF, 2=PME, 3=P3M)
    pub appak: f64,    // Reaction field κ (nm⁻¹)
    pub rcrf: f64,     // Reaction field cutoff (nm)
    pub epsrf: f64,    // Reaction field permittivity
    pub nslfexcl: i32, // Exclusions

    // FORCE block
    /// NTF force flags: [bonds, angles, improper, dihedral, charge, nonbonded]
    /// 0=off, 1=on for each term
    pub ntf: [i32; 6],
    /// Number of energy groups
    pub negr: usize,
    /// Last atom of each energy group (1-indexed)
    pub nre: Vec<usize>,

    // COMTRANSROT
    pub nscm: i32, // COM motion removal: >0 translation every NSCM steps, <0 translation+rotation every |NSCM| steps, 0=off

    // PME-specific parameters
    pub grid_x: usize,    // PME grid size X
    pub grid_y: usize,    // PME grid size Y
    pub grid_z: usize,    // PME grid size Z
    pub pme_order: usize, // PME spline order (typically 4)
    pub pme_alpha: f64,   // PME Ewald parameter

    // INITIALISE block
    pub ntivel: i32, // Initial velocities (0=read, 1=generate)
    pub ntishk: i32, // SHAKE initial configuration
    pub ntinht: i32, // Initial temperature assignment
    pub ntinhb: i32, // Initial bond constraints
    pub ntishi: i32, // Initial SHAKE iterations
    pub nticom: i32, // Initial COM motion removal (0=off, 1=translation, 2=trans+rot)
    pub ig: i64,     // Random seed
    pub tempi: f64,  // Initial temperature (K)

    // WRITETRAJ block
    pub ntwx: usize,        // Trajectory write frequency
    pub ntwe: usize,        // Energy write frequency
    pub ntwv: bool,         // Write velocities
    pub ntwf: bool,         // Write forces
    pub ntwe_special: bool, // Special energy format

    // PRINTOUT block
    pub ntpr: usize, // Print frequency

    // POSITIONRES block
    /// Position restraints (0=off, 1=restrain with CPOR, 2=use B-factors, 3=constrain)
    pub ntpor: i32,
    /// Position restraint read B-factors (0=no, 1=yes)
    pub ntporb: i32,
    /// Position restraint scaling (0=off, 1=on)
    pub ntpors: i32,
    /// Position restraint force constant (kJ mol⁻¹ nm⁻²)
    pub cpor: f64,

    // DISTANCERES block
    /// NTDIR: -2=time-avg+force, -1=time-avg, 0=off, 1=instantaneous, 2=instantaneous+w0
    pub ntdir: i32,
    /// NTDIRA: time-averaging flag
    pub ntdira: i32,
    /// CDIR: force constant for distance restraints (kJ mol⁻¹ nm⁻²)
    pub cdir: f64,
    /// DIR0: linear region threshold (nm)
    pub dir0: f64,
    /// TAUDIR: time constant for time averaging (ps)
    pub taudir: f64,
    /// FORCESCALE: force scaling mode for time-averaged restraints (0,1,2)
    pub forcescale: i32,
    /// VDIR: averaging exponent (unused in plain instantaneous mode)
    pub vdir: i32,
    /// NTWDIR: trajectory write frequency for restraints (0=off)
    pub ntwdir: i32,

    // PERTURBATION block (for FEP/lambda calculations)
    /// NTG: perturbation flag (0=off, 1=on)
    pub ntg: i32,
    /// NRDGL: read lambda from startup (0=no)
    pub nrdgl: i32,
    /// RLAM: starting lambda value
    pub rlam: f64,
    /// DLAMT: lambda change per step
    pub dlamt: f64,
    /// ALPHLJ: LJ soft-core alpha
    pub alphlj: f64,
    /// ALPHC: Coulomb soft-core alpha
    pub alphc: f64,
    /// NLAM: lambda exponent
    pub nlam: i32,
    /// NSCALE: energy group scaling (0=off)
    pub nscale: i32,

    // ENERGYMIN block
    /// Energy minimization method (0=off, 1=steepest descent, 2=Fletcher-Reeves CG)
    pub ntem: i32,
    /// Minimum steps before convergence check
    pub nmin: usize,
    /// Energy convergence tolerance (kJ/mol)
    pub dele: f64,
    /// Initial step size (nm)
    pub dx0: f64,
    /// Maximum step size (nm)
    pub dxm: f64,
    /// Maximum force magnitude limit (0=unlimited)
    pub flim: f64,

    /// Raw blocks for custom parsing
    pub raw_blocks: HashMap<String, Vec<String>>,
}

/// Temperature bath parameters (MULTIBATH block)
#[derive(Debug, Clone)]
pub struct TempBathParameters {
    pub algorithm: i32, // Coupling algorithm: 0=Berendsen, 1=NHC single, N>=2=NHC chain length N
    pub nhc_chain: usize, // NHC chain length (only used when algorithm >= 2)
    pub num_bath_groups: usize,
    pub temp0: Vec<f64>, // Reference temperatures (K)
    pub tau: Vec<f64>,   // Coupling times (ps)
    pub dof: Vec<usize>, // Degrees of freedom per group
}

/// Pressure bath parameters (PRESSURESCALE block)
#[derive(Debug, Clone)]
pub struct PressureParameters {
    pub algorithm: i32, // Coupling algorithm (1=Berendsen, 2=Parrinello-Rahman)
    pub pressure0: [[f64; 3]; 3], // Reference pressure tensor (GPa/bar)
    pub compressibility: [[f64; 3]; 3], // Isothermal compressibility
    pub tau_p: f64,     // Coupling time (ps)
    pub virial: i32,    // Virial calculation method
}

impl Default for ImdParameters {
    fn default() -> Self {
        Self {
            title: String::from("GROMOS simulation"),
            npm: 1,
            nsm: 0,
            nstlim: 1000,
            t0: 0.0,
            dt: 0.002,
            ntb: 1,
            ndfmin: 0,
            num_temp_baths: 1,
            temp_bath: vec![TempBathParameters::default()],
            couple_pressure: false,
            pressure_parameters: None,
            force_groups: Vec::new(),
            ntc: 1,
            ntcp: 0,
            ntcs: 1,
            shake_tol: 1e-4,
            lincs_order_solute: 4,
            lincs_order_solvent: 4,
            algorithm: 0,
            nsnb: 5,
            rcutp: 0.8,
            rcutl: 1.4,
            size: 0.4,
            type_: 0,
            nlrele: 1,
            appak: 0.0,
            rcrf: 1.4,
            epsrf: 0.0,
            nslfexcl: 0,
            ntf: [1, 1, 1, 1, 1, 1], // All force terms on by default
            negr: 1,
            nre: Vec::new(),
            nscm: 100000,
            grid_x: 64,
            grid_y: 64,
            grid_z: 64,
            pme_order: 4,
            pme_alpha: 0.0,
            ntivel: 0,
            ntishk: 0,
            ntinht: 0,
            ntinhb: 0,
            ntishi: 1000,
            nticom: 0,
            ig: 12345,
            tempi: 300.0,
            ntwx: 100,
            ntwe: 100,
            ntwv: false,
            ntwf: false,
            ntwe_special: false,
            ntpr: 100,
            ntpor: 0,
            ntporb: 0,
            ntpors: 0,
            cpor: 0.0,
            ntdir: 0,
            ntdira: 0,
            cdir: 0.0,
            dir0: 0.0,
            taudir: 0.0,
            forcescale: 0,
            vdir: 0,
            ntwdir: 0,
            ntg: 0,
            nrdgl: 0,
            rlam: 0.0,
            dlamt: 0.0,
            alphlj: 0.0,
            alphc: 0.0,
            nlam: 1,
            nscale: 0,
            ntem: 0,
            nmin: 1,
            dele: 0.1,
            dx0: 0.01,
            dxm: 0.05,
            flim: 0.0,
            raw_blocks: HashMap::new(),
        }
    }
}

impl Default for TempBathParameters {
    fn default() -> Self {
        Self {
            algorithm: 0,
            nhc_chain: 1,
            num_bath_groups: 1,
            temp0: vec![300.0],
            tau: vec![0.1],
            dof: vec![0],
        }
    }
}

/// Parse GROMOS .imd/.in parameter file.
///
/// Handles both key-value format and gromosXX positional format.
/// In positional format, comment lines starting with `#` describe fields,
/// and the following data line contains values in that order.
pub fn read_imd_file<P: AsRef<Path>>(path: P) -> Result<ImdParameters, IoError> {
    let file = File::open(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    let reader = BufReader::new(file);

    let mut params = ImdParameters::default();
    let mut current_block = String::new();
    let mut block_lines: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim().to_string();

        if trimmed.is_empty() {
            continue;
        }

        // Check for block start/end
        if trimmed == "END" {
            if !current_block.is_empty() {
                // Filter out comment lines for positional parsing
                let data_lines: Vec<String> = block_lines
                    .iter()
                    .filter(|l| !l.starts_with('#'))
                    .cloned()
                    .collect();
                parse_block(&mut params, &current_block, &data_lines)?;
                params
                    .raw_blocks
                    .insert(current_block.clone(), block_lines.clone());
                block_lines.clear();
            }
            current_block.clear();
        } else if !current_block.is_empty() {
            // Inside a block - collect ALL lines (including comments for raw_blocks)
            block_lines.push(trimmed);
        } else if !trimmed.starts_with('#') {
            // New block starting (not a comment)
            current_block = trimmed;
        }
    }

    Ok(params)
}

/// Parse a specific IMD block using gromosXX positional format.
///
/// Data lines are the non-comment lines within the block, in order.
/// The format follows gromosXX conventions where each data line
/// corresponds to a group of parameters described by preceding comments.
fn parse_block(
    params: &mut ImdParameters,
    block_name: &str,
    data_lines: &[String],
) -> Result<(), IoError> {
    match block_name {
        "TITLE" => {
            params.title = data_lines.join(" ");
        },
        "SYSTEM" => {
            // Line 0: NPM NSM
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if v.len() >= 1 {
                    params.npm = parse_usize(&v[0]);
                }
                if v.len() >= 2 {
                    params.nsm = parse_usize(&v[1]);
                }
            }
        },
        "STEP" => {
            // Line 0: NSTLIM T DT
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if v.len() >= 1 {
                    params.nstlim = parse_usize(&v[0]);
                }
                if v.len() >= 2 {
                    params.t0 = parse_f64(&v[1]);
                }
                if v.len() >= 3 {
                    params.dt = parse_f64(&v[2]);
                }
            }
        },
        "BOUNDCOND" => {
            // Line 0: NTB NDFMIN
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if v.len() >= 1 {
                    params.ntb = parse_i32(&v[0]);
                }
                if v.len() >= 2 {
                    params.ndfmin = parse_i32(&v[1]);
                }
            }
        },
        "MULTIBATH" => {
            // gromosXX format:
            //   Line 0: algorithm name (e.g. "weak-coupling")
            //   Line 1: NBATHS
            //   Line 2: TEMP0 TAU (per bath, may have multiple values)
            //   Line 3: DOFSET (number of DOF sets)
            //   Line 4+: LAST COM-BATH IR-BATH (per DOF set)
            let mut bath = TempBathParameters::default();
            if data_lines.is_empty() {
                return Ok(());
            }

            let mut line_idx = 0;

            // Line 0: algorithm (string like "weak-coupling" or number)
            // Matches gromosXX in_parameter.cc logic:
            //   "weak-coupling" / 0  → algorithm = 0  (Berendsen)
            //   "nose-hoover"   / 1  → algorithm = 1  (NHC single)
            //   "nose-hoover-chains" / 2 → reads next token as chain length N, algorithm = N
            if line_idx < data_lines.len() {
                let v0 = parse_values(&data_lines[line_idx]);
                if let Some(first) = v0.first() {
                    match first.as_str() {
                        "weak-coupling" => {
                            bath.algorithm = 0;
                        },
                        "nose-hoover" => {
                            bath.algorithm = 1;
                        },
                        "nose-hoover-chains" => {
                            // The chain length N follows on the same line, e.g. "nose-hoover-chains  3"
                            let n = if v0.len() > 1 { parse_usize(&v0[1]) } else { 3 };
                            let n = n.max(2); // chain length must be >= 2
                            bath.algorithm = n as i32;
                            bath.nhc_chain = n;
                        },
                        _ => {
                            bath.algorithm = parse_i32(first);
                        },
                    }
                }
                line_idx += 1;
            }

            // Line 1: NBATHS
            if line_idx < data_lines.len() {
                let v = parse_values(&data_lines[line_idx]);
                if let Some(nb) = v.first() {
                    bath.num_bath_groups = parse_usize(nb);
                }
                line_idx += 1;
            }

            // Line 2: TEMP0 TAU (repeated per bath)
            if line_idx < data_lines.len() {
                let v = parse_values(&data_lines[line_idx]);
                bath.temp0.clear();
                bath.tau.clear();
                // Format: temp0_1 tau_1 [temp0_2 tau_2 ...]
                let mut i = 0;
                while i + 1 < v.len() {
                    bath.temp0.push(parse_f64(&v[i]));
                    bath.tau.push(parse_f64(&v[i + 1]));
                    i += 2;
                }
                line_idx += 1;
            }

            // Line 3: DOFSET (number of DOF groups)
            if line_idx < data_lines.len() {
                // Just skip, we don't use it directly
                line_idx += 1;
            }

            // Lines 4+: LAST COM-BATH IR-BATH
            // (skip, used internally by gromosXX)
            let _ = line_idx;

            params.temp_bath = vec![bath];
        },
        "PRESSURESCALE" => {
            // gromosXX PRESSURESCALE block format:
            //   Line 0: COUPLE SCALE COMP TAUP VIRIAL
            //     COUPLE: off/calc/scale (or 0/1/2)
            //     SCALE: off/iso/aniso/full/semianiso (or 0/1/2/3/4)
            //     COMP: compressibility value
            //     TAUP: pressure coupling time
            //     VIRIAL: none/atomic/molecular (or 0/1/2)
            //   Line 1: SEMI (3 values, for semianisotropic)
            //   Lines 2-4: reference pressure (3x3)
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if v.len() >= 1 {
                    // Parse COUPLE keyword/number
                    let couple = match v[0].to_lowercase().as_str() {
                        "off" => 0,
                        "calc" => 1,
                        "scale" => 2,
                        other => parse_i32(other),
                    };
                    params.couple_pressure = couple > 0;
                    if couple > 0 {
                        // Parse SCALE keyword
                        let scale = if v.len() >= 2 {
                            match v[1].to_lowercase().as_str() {
                                "off" => 0,
                                "iso" | "isotropic" => 1,
                                "aniso" | "anisotropic" => 2,
                                "full" | "full_anisotropic" => 3,
                                "semianiso" | "semi_anisotropic" => 4,
                                other => parse_i32(other),
                            }
                        } else {
                            0
                        };
                        // Parse COMP (compressibility)
                        let comp = if v.len() >= 3 {
                            parse_f64(&v[2])
                        } else {
                            4.575e-4
                        };
                        // Parse TAUP
                        let tau_p = if v.len() >= 4 { parse_f64(&v[3]) } else { 0.5 };
                        // Parse VIRIAL keyword
                        let virial = if v.len() >= 5 {
                            match v[4].to_lowercase().as_str() {
                                "none" => 0,
                                "atomic" => 1,
                                "molecular" => 2,
                                other => parse_i32(other),
                            }
                        } else {
                            0
                        };

                        let mut pp = PressureParameters {
                            algorithm: scale, // store SCALE mode (iso/aniso/etc.)
                            pressure0: [[0.0; 3]; 3],
                            compressibility: [[comp; 3]; 3],
                            tau_p,
                            virial,
                        };

                        // Line 1: SEMI (semianisotropic couplings), skip
                        let mut dl = 1;
                        if dl < data_lines.len() {
                            dl += 1;
                        }

                        // Lines 2-4: reference pressure (3x3)
                        for row in 0..3 {
                            if dl < data_lines.len() {
                                let pv = parse_values(&data_lines[dl]);
                                for col in 0..3.min(pv.len()) {
                                    pp.pressure0[row][col] = parse_f64(&pv[col]);
                                }
                                dl += 1;
                            }
                        }
                        params.pressure_parameters = Some(pp);
                    }
                }
            }
        },
        "CONSTRAINT" => {
            // gromosXX format:
            //   Line 0: NTC
            //   Line 1: NTCP (string: "shake" or number)
            //   Line 2: NTCP0 (tolerance)
            //   Line 3: NTCS (string: "shake" or number)
            //   Line 4: NTCS0 (tolerance)
            if data_lines.is_empty() {
                return Ok(());
            }
            let mut idx = 0;
            if idx < data_lines.len() {
                params.ntc = parse_i32(&parse_values(&data_lines[idx])[0]);
                idx += 1;
            }
            if idx < data_lines.len() {
                let v = parse_values(&data_lines[idx]);
                // NTCP can be "shake" (string) or a number
                params.ntcp = match v[0].as_str() {
                    "shake" => 1,
                    "lincs" => 2,
                    "settle" => 3,
                    _ => parse_i32(&v[0]),
                };
                idx += 1;
            }
            if idx < data_lines.len() {
                let v = parse_values(&data_lines[idx]);
                // NTCP0: SHAKE tolerance, or LINCS expansion order when NTCP=lincs
                if params.ntcp == 2 {
                    params.lincs_order_solute = parse_i32(&v[0]) as usize;
                } else {
                    params.shake_tol = parse_f64(&v[0]);
                }
                idx += 1;
            }
            if idx < data_lines.len() {
                let v = parse_values(&data_lines[idx]);
                params.ntcs = match v[0].as_str() {
                    "shake" => 1,
                    "lincs" => 2,
                    "settle" => 3,
                    _ => parse_i32(&v[0]),
                };
                idx += 1;
            }
            if idx < data_lines.len() {
                let v = parse_values(&data_lines[idx]);
                // NTCS0: SHAKE tolerance, or LINCS expansion order when NTCS=lincs
                // (settle reads no NTCS0 parameter — line absent from .in for NTCS=settle)
                if params.ntcs == 2 {
                    params.lincs_order_solvent = parse_i32(&v[0]) as usize;
                }
            }
        },
        "FORCE" => {
            // gromosXX format:
            //   Line 0: bonds angles imp dih charge nonbonded (6 NTF values)
            //   Line 1: NEGR NRE(1) NRE(2) ... NRE(NEGR)
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                for i in 0..6.min(v.len()) {
                    params.ntf[i] = parse_i32(&v[i]);
                }
            }
            if data_lines.len() >= 2 {
                let v = parse_values(&data_lines[1]);
                if !v.is_empty() {
                    params.negr = parse_usize(&v[0]);
                    params.nre = v[1..].iter().map(|s| parse_usize(s)).collect();
                }
            }
        },
        "PAIRLIST" => {
            // gromosXX format:
            //   Line 0: ALGORITHM NSNB RCUTP RCUTL SIZE TYPE
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if v.len() >= 1 {
                    // ALGORITHM can be "standard" or a number
                    params.algorithm = match v[0].as_str() {
                        "standard" => 0,
                        "grid" => 1,
                        _ => parse_i32(&v[0]),
                    };
                }
                if v.len() >= 2 {
                    params.nsnb = parse_usize(&v[1]);
                }
                if v.len() >= 3 {
                    params.rcutp = parse_f64(&v[2]);
                }
                if v.len() >= 4 {
                    params.rcutl = parse_f64(&v[3]);
                }
                if v.len() >= 5 {
                    // SIZE can be "auto" or a number
                    params.size = match v[4].as_str() {
                        "auto" => 0.0,
                        _ => parse_f64(&v[4]),
                    };
                }
                if v.len() >= 6 {
                    params.type_ = match v[5].as_str() {
                        "chargegroup" => 0,
                        "atomic" => 1,
                        _ => parse_i32(&v[5]),
                    };
                }
            }
        },
        "NONBONDED" => {
            // gromosXX format:
            //   Line 0: NLRELE
            //   Line 1: APPAK RCRF EPSRF NSLFEXCL
            //   Line 2: NSHAPE ASHAPE NA2CLC TOLA2 EPSLS
            //   Line 3: NKX NKY NKZ NK2 (optional, for PME)
            //   Line 4: NGX NGY NGZ NASORD NFDORD NALIAS NSPORD (optional)
            //   ...
            if let Some(line) = data_lines.first() {
                params.nlrele = parse_i32(&parse_values(line)[0]);
            }
            if data_lines.len() >= 2 {
                let v = parse_values(&data_lines[1]);
                if v.len() >= 1 {
                    params.appak = parse_f64(&v[0]);
                }
                if v.len() >= 2 {
                    params.rcrf = parse_f64(&v[1]);
                }
                if v.len() >= 3 {
                    params.epsrf = parse_f64(&v[2]);
                }
                if v.len() >= 4 {
                    params.nslfexcl = parse_i32(&v[3]);
                }
            }
            // Line 3+ for PME parameters
            if data_lines.len() >= 4 {
                let v = parse_values(&data_lines[3]);
                if v.len() >= 3 {
                    params.grid_x = parse_usize(&v[0]);
                    params.grid_y = parse_usize(&v[1]);
                    params.grid_z = parse_usize(&v[2]);
                }
            }
        },
        "INITIALISE" => {
            // gromosXX format:
            //   Line 0: NTIVEL NTISHK NTINHT NTINHB
            //   Line 1: NTISHI NTIRTC NTICOM
            //   Line 2: NTISTI
            //   Line 3: IG TEMPI
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if v.len() >= 1 {
                    params.ntivel = parse_i32(&v[0]);
                }
                if v.len() >= 2 {
                    params.ntishk = parse_i32(&v[1]);
                }
                if v.len() >= 3 {
                    params.ntinht = parse_i32(&v[2]);
                }
                if v.len() >= 4 {
                    params.ntinhb = parse_i32(&v[3]);
                }
            }
            if data_lines.len() >= 2 {
                let v = parse_values(&data_lines[1]);
                if v.len() >= 1 {
                    params.ntishi = parse_i32(&v[0]);
                }
                if v.len() >= 3 {
                    params.nticom = parse_i32(&v[2]);
                }
            }
            // Line 2: NTISTI (skip)
            if data_lines.len() >= 4 {
                let v = parse_values(&data_lines[3]);
                if v.len() >= 1 {
                    params.ig = v[0].parse::<i64>().unwrap_or(12345);
                }
                if v.len() >= 2 {
                    params.tempi = parse_f64(&v[1]);
                }
            }
        },
        "WRITETRAJ" => {
            // gromosXX format:
            //   Line 0: NTWX NTWSE NTWV NTWF NTWE NTWG NTWB
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if v.len() >= 1 {
                    params.ntwx = parse_usize(&v[0]);
                }
                // NTWSE at v[1] - skip
                if v.len() >= 3 {
                    params.ntwv = parse_i32(&v[2]) != 0;
                }
                if v.len() >= 4 {
                    params.ntwf = parse_i32(&v[3]) != 0;
                }
                if v.len() >= 5 {
                    params.ntwe = parse_usize(&v[4]);
                }
            }
        },
        "PRINTOUT" => {
            // gromosXX format:
            //   Line 0: NTPR NTPP
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if v.len() >= 1 {
                    params.ntpr = parse_usize(&v[0]);
                }
            }
        },
        "COMTRANSROT" => {
            // Line 0: NSCM
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if v.len() >= 1 {
                    params.nscm = parse_i32(&v[0]);
                }
            }
        },
        "POSITIONRES" => {
            // gromosXX format:
            //   Line 0: NTPOR NTPORB NTPORS CPOR
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if v.len() >= 1 {
                    params.ntpor = parse_i32(&v[0]);
                }
                if v.len() >= 2 {
                    params.ntporb = parse_i32(&v[1]);
                }
                if v.len() >= 3 {
                    params.ntpors = parse_i32(&v[2]);
                }
                if v.len() >= 4 {
                    params.cpor = parse_f64(&v[3]);
                }
            }
        },
        "DISTANCERES" => {
            // gromosXX format (one data line):
            //   NTDIR NTDIRA CDIR DIR0 TAUDIR FORCESCALE VDIR NTWDIR
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if v.len() >= 1 { params.ntdir      = parse_i32(&v[0]); }
                if v.len() >= 2 { params.ntdira     = parse_i32(&v[1]); }
                if v.len() >= 3 { params.cdir       = parse_f64(&v[2]); }
                if v.len() >= 4 { params.dir0       = parse_f64(&v[3]); }
                if v.len() >= 5 { params.taudir     = parse_f64(&v[4]); }
                if v.len() >= 6 { params.forcescale = parse_i32(&v[5]); }
                if v.len() >= 7 { params.vdir       = parse_i32(&v[6]); }
                if v.len() >= 8 { params.ntwdir     = parse_i32(&v[7]); }
            }
        },
        "PERTURBATION" => {
            // gromosXX format (one data line):
            //   NTG NRDGL RLAM DLAMT ALPHLJ ALPHC NLAM NSCALE
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if v.len() >= 1 { params.ntg    = parse_i32(&v[0]); }
                if v.len() >= 2 { params.nrdgl  = parse_i32(&v[1]); }
                if v.len() >= 3 { params.rlam   = parse_f64(&v[2]); }
                if v.len() >= 4 { params.dlamt  = parse_f64(&v[3]); }
                if v.len() >= 5 { params.alphlj = parse_f64(&v[4]); }
                if v.len() >= 6 { params.alphc  = parse_f64(&v[5]); }
                if v.len() >= 7 { params.nlam   = parse_i32(&v[6]); }
                if v.len() >= 8 { params.nscale = parse_i32(&v[7]); }
            }
        },
        "ENERGYMIN" => {
            // gromosXX format:
            //   Line 0: NTEM NCYC DELE DX0 DXM NMIN FLIM
            if let Some(line) = data_lines.first() {
                let v = parse_values(line);
                if v.len() >= 1 {
                    params.ntem = parse_i32(&v[0]);
                }
                // NCYC at v[1] - conjugate gradient cycles, skip
                if v.len() >= 3 {
                    params.dele = parse_f64(&v[2]);
                }
                if v.len() >= 4 {
                    params.dx0 = parse_f64(&v[3]);
                }
                if v.len() >= 5 {
                    params.dxm = parse_f64(&v[4]);
                }
                if v.len() >= 6 {
                    params.nmin = parse_usize(&v[5]);
                }
                if v.len() >= 7 {
                    params.flim = parse_f64(&v[6]);
                }
            }
        },
        _ => {
            // Unknown block - stored in raw_blocks by caller
        },
    }

    Ok(())
}

impl ImdParameters {
    /// Compute the individual lambda and its derivative with respect to RLAM.
    ///
    /// Following gromosXX: individual_lambda = RLAM^NLAM,
    /// d(individual_lambda)/d(RLAM) = NLAM * RLAM^(NLAM-1).
    pub fn lambda_and_derivative(&self) -> (f64, f64) {
        let nlam = self.nlam as f64;
        let l = self.rlam.powf(nlam);
        let dl = if self.nlam <= 0 {
            0.0
        } else {
            nlam * self.rlam.powf(nlam - 1.0)
        };
        (l, dl)
    }
}

/// Split a line into whitespace-separated tokens
fn parse_values(line: &str) -> Vec<String> {
    line.split_whitespace().map(|s| s.to_string()).collect()
}

fn parse_f64(s: &str) -> f64 {
    s.parse::<f64>().unwrap_or(0.0)
}

fn parse_i32(s: &str) -> i32 {
    s.parse::<i32>().unwrap_or(0)
}

fn parse_usize(s: &str) -> usize {
    s.parse::<usize>().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_tmp(content: &str, suffix: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!("gromos_imd_test_{suffix}.tmp"));
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn test_default_parameters() {
        let params = ImdParameters::default();
        assert_eq!(params.nstlim, 1000);
        assert_eq!(params.dt, 0.002);
        assert_eq!(params.ntc, 1); // No constraints
    }

    #[test]
    fn test_parse_gromosxx_format() {
        // Actual gromosXX .in file format (positional)
        let content = "\
TITLE
Two Argon atoms - LJ pair reference
END
SYSTEM
#      NPM      NSM
         1        0
END
STEP
#   NSTLIM         T        DT
        10       0.0     0.002
END
BOUNDCOND
#      NTB    NDFMIN
         0         0
END
FORCE
#      NTF array
# bonds    angles    imp.     dihe     charge nonbonded
# H        H         H        H
     0        0         0        0     0  1
# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)
     2     1      2
END
PAIRLIST
#       ALGORITHM       NSNB    RCUTP   RCUTL   SIZE    TYPE
        standard        1       0.8     1.4     auto    chargegroup
END
NONBONDED
#   NLRELE
         0
#    APPAK      RCRF     EPSRF  NSLFEXCL
       0.0       1.4       1.0         1
END
CONSTRAINT
#       NTC
        1
#       NTCP
        shake
#       NTCP0(1)
        0.0001
#       NTCS
        shake
#       NTCS0(1)
        0.0001
END
WRITETRAJ
#     NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB
         1         0         0         1         1         0         0
END
PRINTOUT
#     NTPR      NTPP
         1         0
END
MULTIBATH
    weak-coupling
#   NBATHS
    1
#   TEMP0  TAU
    300.0  -1.0
#   DOFSET
    1
#   LAST   COM-BATH  IR-BATH
    2      1         1
END
INITIALISE
#  NTIVEL  NTISHK  NTINHT  NTINHB
        0       0       0       0
#  NTISHI  NTIRTC  NTICOM
        1       0       0
#  NTISTI
        0
#      IG   TEMPI
   210185     0.0
END
";
        let path = write_tmp(content, "gxx_format");
        let params = read_imd_file(&path).expect("Failed to parse gromosXX format");

        assert_eq!(params.npm, 1);
        assert_eq!(params.nsm, 0);
        assert_eq!(params.nstlim, 10);
        assert_eq!(params.t0, 0.0);
        assert_eq!(params.dt, 0.002);
        assert_eq!(params.ntb, 0);
        assert_eq!(params.ndfmin, 0);
        assert_eq!(params.ntf, [0, 0, 0, 0, 0, 1]);
        assert_eq!(params.negr, 2);
        assert_eq!(params.nre, vec![1, 2]);
        assert_eq!(params.nlrele, 0);
        assert_eq!(params.appak, 0.0);
        assert_eq!(params.rcrf, 1.4);
        assert_eq!(params.epsrf, 1.0);
        assert_eq!(params.nslfexcl, 1);
        assert_eq!(params.nsnb, 1);
        assert_eq!(params.rcutp, 0.8);
        assert_eq!(params.rcutl, 1.4);
        assert_eq!(params.ntc, 1);
        assert_eq!(params.shake_tol, 0.0001);
        assert_eq!(params.ntwx, 1);
        assert_eq!(params.ntwe, 1);
        assert!(params.ntwf);
        assert_eq!(params.ntpr, 1);
        assert_eq!(params.temp_bath[0].temp0, vec![300.0]);
        assert_eq!(params.temp_bath[0].tau, vec![-1.0]);
        assert_eq!(params.ntivel, 0);
        assert_eq!(params.ig, 210185);
        assert_eq!(params.tempi, 0.0);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn test_parse_water_box_format() {
        let content = "\
TITLE
216 SPC waters in rectangular box
END
SYSTEM
#      NPM      NSM
         1        0
END
STEP
#   NSTLIM         T        DT
       100       0.0     0.002
END
BOUNDCOND
#      NTB    NDFMIN
         1         0
END
FORCE
#      NTF array
# bonds    angles    imp.     dihe     charge nonbonded
# H        H         H        H
     1        1         1        1     1  1
# NEGR    NRE(1)
     1     648
END
NONBONDED
#   NLRELE
         1
#    APPAK      RCRF     EPSRF  NSLFEXCL
       0.0       0.9      62.0         1
END
PAIRLIST
#       ALGORITHM       NSNB    RCUTP   RCUTL   SIZE    TYPE
        standard        5       0.8     0.9     auto    chargegroup
END
";
        let path = write_tmp(content, "water_box");
        let params = read_imd_file(&path).expect("Failed to parse water box");

        assert_eq!(params.nstlim, 100);
        assert_eq!(params.ntb, 1);
        assert_eq!(params.ntf, [1, 1, 1, 1, 1, 1]);
        assert_eq!(params.negr, 1);
        assert_eq!(params.nre, vec![648]);
        assert_eq!(params.nlrele, 1);
        assert_eq!(params.rcrf, 0.9);
        assert_eq!(params.epsrf, 62.0);
        assert_eq!(params.nsnb, 5);
        assert_eq!(params.rcutl, 0.9);

        std::fs::remove_file(path).ok();
    }

    /// Verify MULTIBATH algorithm mapping matches gromosXX exactly:
    ///   "weak-coupling" → 0, "nose-hoover" → 1, "nose-hoover-chains N" → N (N ≥ 2)
    #[test]
    fn test_multibath_algorithm_mapping() {
        fn make_imd(algo_line: &str) -> String {
            format!("\
TITLE\ntest\nEND\nSYSTEM\n#      NPM      NSM\n         1        0\nEND\n\
MULTIBATH\n    {algo_line}\n#   NBATHS\n    1\n#   TEMP0  TAU\n    300.0  0.1\n#   DOFSET\n    1\n#   LAST   COM-BATH  IR-BATH\n    648    1         1\nEND\n")
        }

        // 0 = Berendsen / weak-coupling
        let p = read_imd_file(&write_tmp(&make_imd("weak-coupling"), "nhc_wc")).unwrap();
        assert_eq!(p.temp_bath[0].algorithm, 0, "weak-coupling should map to 0");

        // 1 = NHC single
        let p = read_imd_file(&write_tmp(&make_imd("nose-hoover"), "nhc_nh")).unwrap();
        assert_eq!(p.temp_bath[0].algorithm, 1, "nose-hoover should map to 1");

        // N = NHC chain of length N (N >= 2)
        let p =
            read_imd_file(&write_tmp(&make_imd("nose-hoover-chains  3"), "nhc_chain3")).unwrap();
        assert_eq!(
            p.temp_bath[0].algorithm, 3,
            "nose-hoover-chains 3 should map to 3"
        );
        assert_eq!(p.temp_bath[0].nhc_chain, 3, "nhc_chain should be 3");

        // Numeric fallback: "0" → 0, "1" → 1
        let p = read_imd_file(&write_tmp(&make_imd("0"), "nhc_num0")).unwrap();
        assert_eq!(p.temp_bath[0].algorithm, 0, "numeric 0 should map to 0");

        let p = read_imd_file(&write_tmp(&make_imd("1"), "nhc_num1")).unwrap();
        assert_eq!(p.temp_bath[0].algorithm, 1, "numeric 1 should map to 1");
    }
}
