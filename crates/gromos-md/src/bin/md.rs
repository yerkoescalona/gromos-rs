//! md - Molecular Dynamics simulation engine (gromosXX-compatible CLI)
//!
//! Usage: md @topo <topology> @conf <coordinates> @input <parameters> [@fin <final_conf>]
//!           [@trc <trajectory>] [@tre <energies>] [@trf <forces>] [@trv <velocities>]
//!           [@verb <level>]
//!
//! The main GROMOS-RS molecular dynamics simulation program.
//! Command-line interface matches gromosXX md++ conventions.
//! All simulation parameters are read from the @input (.imd/.in) file.

use gromos::{
    algorithm::{
        BerendsenBarostatParameters,
        ShakeParameters, NtcMode,
        AlgorithmSequence, SimulationState,
        Forcefield, LeapFrogVelocity, LeapFrogPosition,
        TemperatureCalculation, EnergyCalculation, ShakeAlgorithm, SettleAlgorithm, LincsAlgorithm,
        RemoveCOMMotion, BerendsenThermostat,
        PressureCalculation, VirialType,
        BerendsenBarostat, BerendsenBarostatParams,
        SteepestDescentAlgorithm,
    },
    configuration::{Box as SimBox, Configuration},
    interaction::{
        nonbonded::CRFParameters,
    },
    io::{
        coordinate::read_coordinates,
        energy::{EnergyFrame, EnergyWriter},
        imd::read_imd_file,
        posres::{read_posresspec, read_refpos, build_posres_entries},
        topology::{build_topology, read_topology_file},
        trajectory::TrajectoryWriter,
        force::ForceWriter,
        EdsBlock, EdsStatsWriter, EdsVrWriter, GamdBlock, GamdBoostWriter, GamdStatsWriter,
    },
    math::{Periodicity, Rectangular, Vacuum, Vec3},
    random::generate_velocities,
    pairlist::{PairlistContainer, StandardPairlistAlgorithm},
    validation::{
        validate_configuration, validate_coordinates, validate_energy, validate_topology,
    },
};
use std::env;
use std::process;
use std::time::Instant;

/// Simple drop-timer that logs elapsed time at debug level.
struct Timer {
    name: &'static str,
    start: Instant,
}
impl Timer {
    fn new(name: &'static str) -> Self {
        Self { name, start: Instant::now() }
    }
}
impl Drop for Timer {
    fn drop(&mut self) {
        log::debug!("{} took {:.3} ms", self.name, self.start.elapsed().as_secs_f64() * 1000.0);
    }
}

fn print_usage() {
    eprintln!("md - Molecular Dynamics simulation (gromosXX-compatible)");
    eprintln!();
    eprintln!("Usage: md @topo <topology> @conf <coordinates> @input <parameters>");
    eprintln!("          [@fin <final_conf>] [@trc <trajectory>] [@tre <energies>]");
    eprintln!("          [@trf <forces>] [@trv <velocities>] [@verb <level>]");
    eprintln!();
    eprintln!("Input files:");
    eprintln!("  @topo       Molecular topology file (.topo/.top)");
    eprintln!("  @pttopo     Perturbation topology file (.ptp)");
    eprintln!("  @conf       Initial coordinates and restart data (.cnf/.g96)");
    eprintln!("  @input      Input parameter file (.imd/.in) — all simulation settings");
    eprintln!("  @posresspec Position restraints specification");
    eprintln!("  @refpos     Position restraints reference positions");
    eprintln!("  @distrest   Distance restraints specification");
    eprintln!("  @angrest    Angle restraints specification");
    eprintln!("  @dihrest    Dihedral restraints specification");
    eprintln!("  @colvarres  Collective variable restraints");
    eprintln!("  @jval       J-value restraints specification");
    eprintln!("  @xray       X-ray restraints specification");
    eprintln!("  @sym        Symmetry restraints specification");
    eprintln!("  @rdc        RDC restraints specification");
    eprintln!("  @order      Order-parameter restraints specification");
    eprintln!("  @lud        Local elevation umbrella database");
    eprintln!("  @led        Local elevation coordinate definition");
    eprintln!("  @bsleus     Ball & Stick Local Elevation topology");
    eprintln!("  @friction   Atomic friction coefficients");
    eprintln!("  @qmmm       QM/MM specification file");
    eprintln!("  @gamd       GaMD restraints specification");
    eprintln!();
    eprintln!("Output files:");
    eprintln!("  @fin        Final configuration output (.cnf)");
    eprintln!("  @trc        Coordinate trajectory (.trc)");
    eprintln!("  @trv        Velocity trajectory (.trv)");
    eprintln!("  @trf        Force trajectory (.trf)");
    eprintln!("  @trs        Special trajectory (.trs)");
    eprintln!("  @tre        Energy trajectory (.tre)");
    eprintln!("  @trg        Free energy trajectory (.trg)");
    eprintln!("  @bae        Block averaged energy trajectory");
    eprintln!("  @bag        Block averaged free-energy trajectory");
    eprintln!();
    eprintln!("Control:");
    eprintln!("  @verb       Verbosity level (0=quiet, 1=info, 2=debug)");
    eprintln!("  @print      Print additional information (pairlist, force)");
    eprintln!("  @version    Print version information");
    eprintln!("  @develop    Run untested development code");
    eprintln!();
    eprintln!("All simulation parameters (timestep, cutoffs, thermostat, etc.)");
    eprintln!("are specified in the @input file, following gromosXX conventions.");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  md @topo system.topo @conf initial.cnf @input run.imd");
    eprintln!("  md @topo system.topo @conf initial.cnf @input run.imd @trc out.trc @tre out.tre");
    eprintln!("  md @topo system.topo @conf initial.cnf @input run.imd @verb 1");
}

/// gromosXX-compatible command-line arguments.
/// These are file paths and control flags — all simulation parameters
/// come from the @input file.
#[derive(Debug)]
struct MDArgs {
    // Required input files
    topo_file: String,
    conf_file: String,
    input_file: String,
    // Optional input files
    pttopo_file: Option<String>,
    posresspec_file: Option<String>,
    refpos_file: Option<String>,
    distrest_file: Option<String>,
    angrest_file: Option<String>,
    dihrest_file: Option<String>,
    gamd_file: Option<String>,
    // Output files
    fin_file: Option<String>,
    trc_file: Option<String>,
    trv_file: Option<String>,
    trf_file: Option<String>,
    trs_file: Option<String>,
    tre_file: Option<String>,
    trg_file: Option<String>,
    bae_file: Option<String>,
    bag_file: Option<String>,
    // Control
    verbose: usize,
    print_flags: Vec<String>,
    develop: bool,
}

impl Default for MDArgs {
    fn default() -> Self {
        Self {
            topo_file: String::new(),
            conf_file: String::new(),
            input_file: String::new(),
            pttopo_file: None,
            posresspec_file: None,
            refpos_file: None,
            distrest_file: None,
            angrest_file: None,
            dihrest_file: None,
            gamd_file: None,
            fin_file: None,
            trc_file: None,
            trv_file: None,
            trf_file: None,
            trs_file: None,
            tre_file: None,
            trg_file: None,
            bae_file: None,
            bag_file: None,
            verbose: 0,
            print_flags: Vec::new(),
            develop: false,
        }
    }
}

fn parse_args(args: Vec<String>) -> Result<MDArgs, String> {
    let mut md_args = MDArgs::default();

    let mut i = 1;
    while i < args.len() {
        let arg = args[i].as_str();

        // Helper: get the next argument value
        macro_rules! next_val {
            ($name:expr) => {{
                i += 1;
                if i >= args.len() {
                    return Err(format!("Missing value for {}", $name));
                }
                args[i].clone()
            }};
        }

        match arg {
            // Required input files
            "@topo" => md_args.topo_file = next_val!("@topo"),
            "@conf" => md_args.conf_file = next_val!("@conf"),
            "@input" | "@imd" => md_args.input_file = next_val!("@input"),
            // Optional input files
            "@pttopo" => md_args.pttopo_file = Some(next_val!("@pttopo")),
            "@posresspec" => md_args.posresspec_file = Some(next_val!("@posresspec")),
            "@refpos" => md_args.refpos_file = Some(next_val!("@refpos")),
            "@distrest" => md_args.distrest_file = Some(next_val!("@distrest")),
            "@angrest" => md_args.angrest_file = Some(next_val!("@angrest")),
            "@dihrest" => md_args.dihrest_file = Some(next_val!("@dihrest")),
            "@gamd" => md_args.gamd_file = Some(next_val!("@gamd")),
            // Output files
            "@fin" => md_args.fin_file = Some(next_val!("@fin")),
            "@trc" | "@traj" => md_args.trc_file = Some(next_val!(arg)),
            "@trv" => md_args.trv_file = Some(next_val!("@trv")),
            "@trf" => md_args.trf_file = Some(next_val!("@trf")),
            "@trs" => md_args.trs_file = Some(next_val!("@trs")),
            "@tre" | "@ene" => md_args.tre_file = Some(next_val!(arg)),
            "@trg" => md_args.trg_file = Some(next_val!("@trg")),
            "@bae" => md_args.bae_file = Some(next_val!("@bae")),
            "@bag" => md_args.bag_file = Some(next_val!("@bag")),
            // Control
            "@verb" | "@verbose" | "@v" => {
                md_args.verbose = next_val!("@verb")
                    .parse()
                    .map_err(|_| format!("Invalid value for @verb: {}", args[i]))?;
            },
            "@print" => {
                md_args.print_flags.push(next_val!("@print"));
            },
            "@version" => {
                println!("GROMOS-RS md {}", env!("CARGO_PKG_VERSION"));
                process::exit(0);
            },
            "@develop" => {
                md_args.develop = true;
            },
            // Catch unrecognized
            _ if arg.starts_with('@') => {
                return Err(format!("Unknown argument: {}", arg));
            },
            _ => {
                return Err(format!("Unexpected argument: {} (did you forget @?)", arg));
            },
        }
        i += 1;
    }

    if md_args.topo_file.is_empty() {
        return Err("Missing required argument @topo".to_string());
    }
    if md_args.conf_file.is_empty() {
        return Err("Missing required argument @conf".to_string());
    }
    if md_args.input_file.is_empty() {
        return Err("Missing required argument @input".to_string());
    }

    Ok(md_args)
}

fn main() {
    let t_total = Instant::now();
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args.contains(&"--help".to_string()) || args.contains(&"-h".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    // Parse command-line arguments (file paths only, gromosXX style)
    let md_args = match parse_args(args) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        },
    };

    // Set up logging
    let filter = match md_args.verbose {
        0 => "info",
        _ => "debug",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(filter))
        .init();

    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║                   GROMOS-RS MD Engine                        ║");
    println!("║           Rust Implementation of GROMOS MD                   ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    log::info!("GROMOS-RS MD simulation starting");
    log::debug!("Verbose level: {}", md_args.verbose);

    // === Read simulation parameters from @input file ===
    println!("Loading input parameters: {}", md_args.input_file);
    let imd = match read_imd_file(&md_args.input_file) {
        Ok(p) => p,
        Err(e) => {
            log::error!("Failed to read input file: {}", e);
            eprintln!("Error reading input file: {}", e);
            process::exit(1);
        },
    };
    log::debug!("Input parameters: steps={}, dt={}, ntb={}", imd.nstlim, imd.dt, imd.ntb);

    // Extract simulation parameters from IMD file
    let n_steps = imd.nstlim;
    let dt = imd.dt;
    let cutoff = imd.rcutl;
    let epsilon = 1.0; // CRF interior dielectric (always 1 in GROMOS)
    let rf_epsilon = imd.epsrf;
    let rf_kappa = imd.appak;
    let pairlist_update = imd.nsnb;
    let ntf = imd.ntf;
    let ntf_bond = ntf[0] != 0;
    let ntf_angle = ntf[1] != 0;
    let ntf_improper = ntf[2] != 0;
    let ntf_dihedral = ntf[3] != 0;
    // Constraint algorithm selection: gromosXX dispatches the solute algorithm
    // on NTCP (only relevant when NTC>1) and the solvent algorithm on NTCS
    // (only relevant when solvent molecules exist), independently.
    let solute_constrained = imd.ntc > 1;
    let solute_lincs = solute_constrained && imd.ntcp == 2;
    let solute_shake = solute_constrained && !solute_lincs;

    let solvent_constrained = imd.ntcs > 0 && imd.nsm > 0;
    let solvent_settle = solvent_constrained && imd.ntcs == 3;
    let solvent_lincs = solvent_constrained && imd.ntcs == 2;
    let solvent_shake = solvent_constrained && !solvent_settle && !solvent_lincs;

    let shake_enabled = solute_shake || solvent_shake;
    let settle_enabled = solvent_settle;
    let lincs_enabled = solute_lincs || solvent_lincs;
    let shake_tolerance = imd.shake_tol;
    let nstxout = imd.ntwx;
    let nstener = imd.ntwe;
    let nstlog = imd.ntpr;
    let temperature = if !imd.temp_bath.is_empty() && !imd.temp_bath[0].temp0.is_empty() {
        imd.temp_bath[0].temp0[0]
    } else {
        300.0
    };
    let thermostat_tau = if !imd.temp_bath.is_empty() && !imd.temp_bath[0].tau.is_empty() {
        imd.temp_bath[0].tau[0]
    } else {
        -1.0
    };
    let thermostat_on = thermostat_tau > 0.0;

    // Derive output file paths (from @args or defaults)
    let trc_file = md_args.trc_file.clone().unwrap_or_else(|| "md.trc".to_string());
    let tre_file = md_args.tre_file.clone().unwrap_or_else(|| "md.tre".to_string());

    // === gromosXX initialization order: parameters → topology → coordinates ===

    // === 1. Load topology (read_topology) ===
    println!("Loading topology: {}", md_args.topo_file);
    log::debug!("Reading topology file: {}", md_args.topo_file);
    let _timer = Timer::new("Topology loading");

    let topo_data = match read_topology_file(&md_args.topo_file) {
        Ok(data) => data,
        Err(e) => {
            log::error!("Failed to read topology: {}", e);
            eprintln!("Error reading topology: {}", e);
            process::exit(1);
        },
    };

    // gromosXX convention: read_topology() then topo.solvate(0, nsm)
    let mut topo = build_topology(topo_data);

    // === 2. Load coordinates (read_configuration) ===
    // Load coordinates BEFORE solvation so we can determine actual NSM
    println!("Loading coordinates: {}", md_args.conf_file);
    log::debug!("Reading coordinate file: {}", md_args.conf_file);
    let _timer = Timer::new("Coordinate loading");

    let coord_data = match read_coordinates(&md_args.conf_file) {
        Ok(data) => data,
        Err(e) => {
            log::error!("Failed to read coordinates: {}", e);
            eprintln!("Error reading coordinates: {}", e);
            process::exit(1);
        },
    };

    let positions = coord_data.positions;
    let velocities = coord_data.velocities;
    let box_dims = coord_data.box_dims;

    // Determine actual NSM from coordinate file if solvent template exists
    let actual_nsm = if imd.nsm > 0 && !topo.solvent_atom_template.is_empty() {
        let atoms_per_solvent = topo.solvent_atom_template.len();
        let n_solute = topo.solute.num_atoms();
        let remaining = positions.len().saturating_sub(n_solute);
        if remaining % atoms_per_solvent != 0 {
            eprintln!("Error: ({} coords - {} solute) is not divisible by {} atoms/solvent",
                positions.len(), n_solute, atoms_per_solvent);
            process::exit(1);
        }
        let nsm_from_coords = remaining / atoms_per_solvent;
        if nsm_from_coords != imd.nsm {
            println!("  Adjusting NSM: {} (imd) -> {} (from coordinates)", imd.nsm, nsm_from_coords);
        }
        nsm_from_coords
    } else {
        imd.nsm
    };
    topo.solvate(actual_nsm);

    println!("  Solute atoms: {}", topo.solute.num_atoms());
    if !topo.solvents.is_empty() {
        let sv = &topo.solvents[0];
        println!("  Solvent: {} molecules × {} atoms = {} atoms",
            sv.num_molecules, sv.atoms_per_molecule(), sv.total_atoms());
    }
    println!("  Total atoms: {}", topo.num_atoms());
    println!("  Bonds: {}", topo.solute.bonds.len());
    println!("  Angles: {}", topo.solute.angles.len());
    println!("  Dihedrals: {}", topo.solute.proper_dihedrals.len());
    println!("  Impropers: {}", topo.solute.improper_dihedrals.len());
    println!("  Chargegroups: {}", topo.chargegroups.len());
    println!();

    // Validate topology
    log::debug!("Validating topology");
    let topo_validation = validate_topology(&topo);
    if topo_validation.has_errors() {
        topo_validation.print();
        topo_validation.print_summary();
        if topo_validation.has_fatal() {
            log::error!("Fatal errors in topology - cannot continue");
            process::exit(1);
        }
        log::warn!("Topology has errors, but continuing");
    } else if !topo_validation.warnings.is_empty() {
        topo_validation.print();
        log::debug!("{} warnings in topology", topo_validation.warnings.len());
    } else {
        log::debug!("Topology validation passed");
    }

    println!("  Positions loaded: {}", positions.len());
    if !velocities.is_empty() {
        println!("  Velocities loaded: {}", velocities.len());
    }
    println!(
        "  Box: ({:.4}, {:.4}, {:.4}) nm",
        box_dims.x, box_dims.y, box_dims.z
    );
    println!("  Box type: {}", match coord_data.box_type {
        0 => "vacuum",
        1 => "rectangular",
        2 => "triclinic",
        3 => "truncated octahedron",
        _ => "unknown",
    });
    println!();

    // check_configuration: atom count must match
    if positions.len() != topo.num_atoms() {
        log::error!(
            "Atom count mismatch: topology={}, coordinates={}",
            topo.num_atoms(),
            positions.len()
        );
        eprintln!(
            "Error: Number of atoms in topology ({}) != coordinates ({})",
            topo.num_atoms(),
            positions.len()
        );
        process::exit(1);
    }

    // Parse input file for GAMD/EDS blocks
    let gamd_block = {
        log::debug!("Parsing input file for GAMD block: {}", md_args.input_file);
        match GamdBlock::parse_file(&md_args.input_file) {
            Ok(block) => block,
            Err(e) => {
                log::warn!("Failed to parse GAMD block: {}", e);
                None
            },
        }
    };

    let eds_block = {
        log::debug!("Parsing input file for EDS block: {}", md_args.input_file);
        match EdsBlock::parse_file(&md_args.input_file) {
            Ok(block) => block,
            Err(e) => {
                log::warn!("Failed to parse EDS block: {}", e);
                None
            },
        }
    };

    if let Some(ref block) = gamd_block {
        println!();
        println!("GAMD Parameters detected:");
        println!("  Search mode:  {:?}", block.search_mode);
        println!("  Boost form:   {:?}", block.boost_form);
        println!("  Threshold:    {:?}", block.threshold_type);
        println!("  Sigma0 dih:   {:.2}", block.sigma0_dih);
        println!("  Sigma0 tot:   {:.2}", block.sigma0_tot);
        if let (Some(k), Some(e)) = (block.k_tot, block.e_tot) {
            println!("  K_tot:        {:.6}", k);
            println!("  E_tot:        {:.2}", e);
        }
        println!();
    }

    if let Some(ref block) = eds_block {
        println!();
        println!("EDS Parameters detected:");
        println!("  Num states:   {}", block.num_states);
        println!("  Form:         {:?}", block.form);
        println!("  S values:     {:?}", block.s_values);
        println!("  E offsets:    {:?}", block.e_offsets);
        println!("  Temperature:  {:.1} K", block.temperature);
        if block.search_enabled {
            println!(
                "  AEDS enabled: E_max={:.2}, E_min={:.2}",
                block.e_max, block.e_min
            );
        }
        println!();
    }

    // Check for conflicting modes
    if gamd_block.is_some() && eds_block.is_some() {
        log::error!("Cannot enable both GAMD and EDS simultaneously");
        eprintln!("Error: Both GAMD and EDS blocks found in input file");
        eprintln!("       These methods cannot be used together in the same simulation");
        process::exit(1);
    }

    // Create GAMD parameters if enabled
    let mut gamd_params = gamd_block.as_ref().map(|block| {
        log::info!("Creating GAMD parameters from input block");
        block.to_parameters()
    });

    // Validate coordinates
    log::debug!("Validating coordinates");
    let coord_validation = validate_coordinates(&positions, Some(box_dims));
    if coord_validation.has_errors() {
        coord_validation.print();
        coord_validation.print_summary();
        if coord_validation.has_fatal() {
            log::error!("Fatal errors in coordinates - cannot continue");
            process::exit(1);
        }
        log::warn!("Coordinates have errors, but continuing");
    } else if !coord_validation.warnings.is_empty() {
        coord_validation.print();
        log::debug!(
            "{} warnings in coordinates",
            coord_validation.warnings.len()
        );
    } else {
        log::debug!("Coordinate validation passed");
    }

    // Create configuration
    log::debug!("Creating configuration");
    let mut conf = Configuration::new(topo.num_atoms(), 1, 1);
    conf.current_mut().pos = positions.clone();
    conf.current_mut().vel = if imd.ntivel == 1 {
        log::info!(
            "Generating initial velocities (NTIVEL=1): T={:.2} K, seed={}",
            imd.tempi,
            imd.ig
        );
        generate_velocities(imd.tempi, imd.ig as u32, &topo.mass)
    } else if velocities.len() == topo.num_atoms() {
        velocities.clone()
    } else {
        vec![Vec3::ZERO; topo.num_atoms()]
    };
    conf.current_mut().box_config = SimBox::rectangular(box_dims.x, box_dims.y, box_dims.z);
    conf.copy_current_to_old();

    // Create EDS parameters if enabled (now that we have num_atoms)
    let mut eds_params = if let Some(ref block) = eds_block {
        log::info!("Creating EDS parameters from input block");
        if block.search_enabled {
            match block.to_aeds_parameters(topo.num_atoms()) {
                Ok(aeds) => Some(aeds),
                Err(e) => {
                    log::error!("Failed to create AEDS parameters: {}", e);
                    eprintln!("Error: Failed to create AEDS parameters: {}", e);
                    process::exit(1);
                },
            }
        } else {
            match block.to_parameters(topo.num_atoms()) {
                Ok(eds) => {
                    // Wrap in AEDS for uniform handling
                    Some(gromos::eds::AEDSParameters::new(eds, 0.0, 0.0, false))
                },
                Err(e) => {
                    log::error!("Failed to create EDS parameters: {}", e);
                    eprintln!("Error: Failed to create EDS parameters: {}", e);
                    process::exit(1);
                },
            }
        }
    } else {
        None
    };

    // Validate configuration
    log::debug!("Validating configuration (topology + coordinates)");
    let conf_validation = validate_configuration(&topo, &conf);
    if conf_validation.has_errors() {
        conf_validation.print();
        conf_validation.print_summary();
        if conf_validation.has_fatal() {
            log::error!("Fatal errors in configuration - cannot continue");
            process::exit(1);
        }
        log::warn!("Configuration has errors, but continuing");
    } else if !conf_validation.warnings.is_empty() {
        conf_validation.print();
        log::debug!(
            "{} warnings in configuration",
            conf_validation.warnings.len()
        );
    } else {
        log::debug!("Configuration validation passed");
    }

    // Setup nonbonded interactions
    println!("Setting up nonbonded interactions:");
    println!("  Cutoff:      {:.3} nm", cutoff);
    println!("  Epsilon:     {:.2}", epsilon);
    println!("  RF epsilon:  {:.2}", rf_epsilon);
    println!("  RF kappa:    {:.3} nm^-1", rf_kappa);
    println!();

    log::debug!("Calculating CRF parameters");
    let crf_params = CRFParameters::new(
        cutoff,
        epsilon,
        rf_epsilon,
        rf_kappa,
    );
    log::debug!(
        "CRF parameters: crf_2cut3i={:.6}, crf_cut3i={:.6}",
        crf_params.crf_2cut3i,
        crf_params.crf_cut3i
    );

    log::debug!("Converting LJ parameter matrix");
    let lj_params = Forcefield::convert_lj_parameters(&topo);
    log::debug!(
        "LJ parameter matrix: {}x{} atom types",
        lj_params.len(),
        if lj_params.is_empty() {
            0
        } else {
            lj_params[0].len()
        }
    );

    log::debug!("Initializing pairlist");
    let mut pairlist = PairlistContainer::new(
        imd.rcutp, // short range cutoff (RCUTP)
        cutoff,    // long range cutoff (RCUTL)
        0.0,       // skin (no extra distance)
    );
    pairlist.update_frequency = pairlist_update;
    log::debug!(
        "Pairlist update frequency: {} steps",
        pairlist.update_frequency
    );

    let use_chargegroups = !topo.chargegroups.is_empty();
    let pairlist_algorithm = StandardPairlistAlgorithm::new(use_chargegroups);
    let periodicity = if box_dims.x == 0.0 && box_dims.y == 0.0 && box_dims.z == 0.0 {
        Periodicity::Vacuum(Vacuum)
    } else {
        Periodicity::Rectangular(Rectangular::new(box_dims))
    };

    // Initial pairlist generation
    log::debug!("Generating initial pairlist");
    pairlist_algorithm.update(&topo, &conf, &mut pairlist, &periodicity);
    println!("  Initial pairlist: {} pairs", pairlist.total_pairs());
    println!();

    // Load position restraints if requested
    let posres = if imd.ntpor > 0 {
        if let Some(ref por_file) = md_args.posresspec_file {
            use gromos::interaction::restraints::{PositionRestraint, PositionRestraints};

            // Read restrained atom indices from POSRESSPEC
            let restrained_atoms = read_posresspec(por_file).unwrap_or_else(|e| {
                eprintln!("Error reading position restraint spec: {}", e);
                process::exit(1);
            });

            // Get reference positions:
            // NTPORB=1: from @refpos file
            // NTPORB=0: from startup configuration (positions)
            let ref_positions = if imd.ntporb >= 1 {
                if let Some(ref rpr_file) = md_args.refpos_file {
                    read_refpos(rpr_file).unwrap_or_else(|e| {
                        eprintln!("Error reading reference positions: {}", e);
                        process::exit(1);
                    })
                } else {
                    eprintln!("NTPORB={} but no @refpos file specified", imd.ntporb);
                    process::exit(1);
                }
            } else {
                // Use startup positions as reference
                positions.clone()
            };

            // Build entries combining atom indices with reference positions
            let entries = build_posres_entries(&restrained_atoms, &ref_positions);

            let mut pr = PositionRestraints::new();
            for entry in &entries {
                pr.add_restraint(PositionRestraint::new(
                    entry.atom,
                    entry.reference_pos,
                    imd.cpor,
                ));
            }
            println!("  Position restraints: {} atoms, CPOR={:.1} kJ/(mol·nm²)",
                entries.len(), imd.cpor);
            Some(pr)
        } else {
            eprintln!("NTPOR={} but no @posresspec file specified", imd.ntpor);
            process::exit(1);
        }
    } else {
        None
    };

    // === Build Algorithm Sequence (gromosXX pattern) ===
    let is_minimization = imd.ntem > 0;

    if is_minimization {
        println!("Setting up algorithm sequence: Steepest Descent Energy Minimization");
        println!("  NTEM:  {} (steepest descent)", imd.ntem);
        println!("  DELE:  {:.6} kJ/mol", imd.dele);
        println!("  DX0:   {:.6} nm", imd.dx0);
        println!("  DXM:   {:.6} nm", imd.dxm);
        println!("  NMIN:  {}", imd.nmin);
        println!("  FLIM:  {:.6}", imd.flim);
        println!();
    } else {
        println!("Setting up algorithm sequence: Leap-Frog");
    }
    let mut md_sequence = AlgorithmSequence::new();

    if is_minimization {
        // Energy minimization sequence:
        //   1. Forcefield (compute forces)
        //   2. SteepestDescent (exchange + position update along force direction)
        //   3. EnergyCalculation (finalize energies)

        // Forcefield (bonded + nonbonded forces)
        let mut forcefield = Forcefield::new(
            lj_params,
            crf_params,
            periodicity,
            pairlist,
            pairlist_algorithm,
        );
        forcefield.ntf_bond = ntf_bond;
        forcefield.ntf_angle = ntf_angle;
        forcefield.ntf_improper = ntf_improper;
        forcefield.ntf_dihedral = ntf_dihedral;
        forcefield.parallel_nonbonded = topo.num_atoms() > 100;
        if !topo.solvent_atom_template.is_empty() {
            forcefield.atoms_per_solvent = topo.solvent_atom_template.len();
        }
        if let Some(pr) = posres.clone() {
            forcefield.position_restraints = pr;
        }
        md_sequence.push(Box::new(forcefield));

        // Steepest Descent minimizer
        let sd = SteepestDescentAlgorithm::new()
            .with_tolerance(imd.dele)
            .with_step_sizes(imd.dx0, imd.dxm)
            .with_min_steps(imd.nmin)
            .with_force_limit(imd.flim);
        md_sequence.push(Box::new(sd));

        // SHAKE constraints (gromosXX applies constraints even during minimization)
        if shake_enabled {
            let ntc_mode = if solute_shake {
                match imd.ntc {
                    3 => NtcMode::AllBonds,
                    2 => NtcMode::HydrogenBonds,
                    _ => NtcMode::SolventOnly,
                }
            } else {
                NtcMode::SolventOnly
            };
            let mut shake_alg = ShakeAlgorithm::new(ShakeParameters {
                tolerance: shake_tolerance,
                max_iterations: 1000,
                ntc: ntc_mode,
            });
            shake_alg.include_solvent = solvent_shake;
            if imd.ntishk >= 1 {
                shake_alg.shake_initial_positions = true;
            }
            if imd.ntishk >= 2 {
                shake_alg.shake_initial_velocities = true;
            }
            md_sequence.push(Box::new(shake_alg));
        }

        // SETTLE constraints (analytical rigid-water solver, NTCS=settle)
        if settle_enabled {
            md_sequence.push(Box::new(SettleAlgorithm::new()));
        }

        // LINCS constraints (linear constraint solver, NTCP/NTCS=lincs)
        if lincs_enabled {
            let ntc_mode = match imd.ntc {
                3 => NtcMode::AllBonds,
                2 => NtcMode::HydrogenBonds,
                _ => NtcMode::SolventOnly,
            };
            md_sequence.push(Box::new(LincsAlgorithm::new(
                ntc_mode,
                imd.lincs_order_solute,
                imd.lincs_order_solvent,
                solute_lincs,
                solvent_lincs,
            )));
        }

        // Energy calculation (finalize totals)
        // Note: TemperatureCalculation is NOT included for EM — gromosXX EM
        // does not compute kinetic energy (E_total = E_pot).
        md_sequence.push(Box::new(EnergyCalculation::new()));
    } else {
        // Standard MD sequence

        // 1. COM motion removal (gromosXX: first in sequence, before forcefield)
        if imd.nticom >= 1 || imd.nscm != 0 {
            md_sequence.push(Box::new(RemoveCOMMotion::new(imd.nticom, imd.nscm)));
        }

        // 2. Forcefield (bonded + nonbonded forces)
        let mut forcefield = Forcefield::new(
            lj_params,
            crf_params,
            periodicity,
            pairlist,
            pairlist_algorithm,
        );
        forcefield.ntf_bond = ntf_bond;
        forcefield.ntf_angle = ntf_angle;
        forcefield.ntf_improper = ntf_improper;
        forcefield.ntf_dihedral = ntf_dihedral;
        forcefield.parallel_nonbonded = topo.num_atoms() > 100;
        if !topo.solvent_atom_template.is_empty() {
            forcefield.atoms_per_solvent = topo.solvent_atom_template.len();
        }
        if imd.couple_pressure {
            forcefield.virial_type = match imd.pressure_parameters.as_ref().map(|p| p.virial).unwrap_or(0) {
                2 => VirialType::Molecular,
                1 => VirialType::Atomic,
                _ => VirialType::None,
            };
        }
        if let Some(pr) = posres {
            forcefield.position_restraints = pr;
        }
        md_sequence.push(Box::new(forcefield));

    // 3. Leap-Frog velocity step (exchange_state + v update)
        md_sequence.push(Box::new(LeapFrogVelocity::new()));

        // 3b. Berendsen thermostat (between velocity and position update, gromosXX convention)
        if thermostat_on {
            // Compute degrees of freedom: 3N - N_constraints - NDFMIN
            let n_atoms = topo.num_atoms();
            let n_solute = topo.num_solute_atoms();
            let atoms_per_solvent = if !topo.solvent_atom_template.is_empty() {
                topo.solvent_atom_template.len()
            } else {
                1
            };
            let n_solvent_molecules = if atoms_per_solvent > 0 && n_atoms > n_solute {
                (n_atoms - n_solute) / atoms_per_solvent
            } else {
                0
            };
            let solvent_constraint_dof = if shake_enabled || settle_enabled || lincs_enabled {
                n_solvent_molecules * topo.solvent_constraint_template.len()
            } else {
                0
            };
            // TODO: solute constraint DOF for NTC=2,3
            let solute_constraint_dof = 0usize;
            let total_dof = (3 * n_atoms - solvent_constraint_dof - solute_constraint_dof) as f64
                - imd.ndfmin as f64;
            println!("  Thermostat DOF: {:.0} (3*{} - {} solvent_constr - {} NDFMIN)",
                total_dof, n_atoms, solvent_constraint_dof, imd.ndfmin);
            md_sequence.push(Box::new(BerendsenThermostat::new_single_bath(
                temperature, thermostat_tau, total_dof, n_atoms,
            )));
        }

        // 4. Leap-Frog position step (r update)
        md_sequence.push(Box::new(LeapFrogPosition::new()));

        // 5. SHAKE constraints (if enabled)
        if shake_enabled {
            let ntc_mode = if solute_shake {
                match imd.ntc {
                    3 => NtcMode::AllBonds,
                    2 => NtcMode::HydrogenBonds,
                    _ => NtcMode::SolventOnly,
                }
            } else {
                NtcMode::SolventOnly
            };
            let mut shake_alg = ShakeAlgorithm::new(ShakeParameters {
                tolerance: shake_tolerance,
                max_iterations: 1000,
                ntc: ntc_mode,
            });
            shake_alg.include_solvent = solvent_shake;
            // gromosXX: NTISHK controls initial position/velocity shaking
            // NTISHK=1: shake positions, NTISHK=2: shake positions + velocities
            if imd.ntishk >= 1 {
                shake_alg.shake_initial_positions = true;
            }
            if imd.ntishk >= 2 {
                shake_alg.shake_initial_velocities = true;
            }
            md_sequence.push(Box::new(shake_alg));
        }

        // 5b. SETTLE constraints (analytical rigid-water solver, NTCS=settle)
        if settle_enabled {
            md_sequence.push(Box::new(SettleAlgorithm::new()));
        }

        // 5c. LINCS constraints (linear constraint solver, NTCP/NTCS=lincs)
        if lincs_enabled {
            let ntc_mode = match imd.ntc {
                3 => NtcMode::AllBonds,
                2 => NtcMode::HydrogenBonds,
                _ => NtcMode::SolventOnly,
            };
            md_sequence.push(Box::new(LincsAlgorithm::new(
                ntc_mode,
                imd.lincs_order_solute,
                imd.lincs_order_solvent,
                solute_lincs,
                solvent_lincs,
            )));
        }

        // 6. Temperature/kinetic energy calculation
        md_sequence.push(Box::new(TemperatureCalculation::new()));

        // 7. Pressure calculation and barostat (if pressure coupling is on)
        if imd.couple_pressure {
            let virial_type = match imd.pressure_parameters.as_ref().map(|p| p.virial).unwrap_or(0) {
                2 => VirialType::Molecular,
                1 => VirialType::Atomic,
                _ => VirialType::None,
            };
            md_sequence.push(Box::new(PressureCalculation::new(virial_type)));

            let pp = imd.pressure_parameters.as_ref();
            md_sequence.push(Box::new(BerendsenBarostat::new(BerendsenBarostatParams {
                pressure0: pp.map(|p| p.pressure0[0][0]).unwrap_or(1.0),
                compressibility: pp.map(|p| p.compressibility[0][0]).unwrap_or(4.575e-4),
                tau: pp.map(|p| p.tau_p).unwrap_or(0.5),
            })));
        }

        // 8. Energy finalization
        md_sequence.push(Box::new(EnergyCalculation::new()));
    } // end of MD vs minimization branch

    println!("  Sequence: {}", md_sequence.algorithm_names().join(" → "));
    println!();

    // Initialize the sequence
    let mut sim_state = SimulationState::new(dt, n_steps);
    md_sequence.init(&topo, &mut conf, &sim_state).unwrap_or_else(|e| {
        eprintln!("Error initializing algorithm sequence: {}", e);
        process::exit(1);
    });

    // Thermostat is now wired into the algorithm sequence (above)
    if thermostat_on {
        println!("Setting up thermostat: Berendsen");
        println!("  Target temp:   {:.1} K", temperature);
        println!("  Coupling time: {:.3} ps", thermostat_tau);
        println!();
    }

    // Setup barostat (from PRESSURESCALE block)
    let barostat_params = if imd.couple_pressure {
        let pp = imd.pressure_parameters.as_ref();
        Some(BerendsenBarostatParameters {
            target_pressure: pp.map(|p| p.pressure0[0][0]).unwrap_or(1.0),
            coupling_time: pp.map(|p| p.tau_p).unwrap_or(0.5),
            compressibility: pp.map(|p| p.compressibility[0][0]).unwrap_or(4.575e-4),
            isotropic: true,
        })
    } else {
        None
    };

    if let Some(ref params) = barostat_params {
        println!("Setting up barostat: Berendsen");
        println!("  Target pres:   {:.1} bar", params.target_pressure);
        println!("  Coupling time: {:.3} ps", params.coupling_time);
        println!();
    }

    // Setup SHAKE constraints
    let shake_params = if shake_enabled {
        let ntc_mode = match imd.ntc {
            3 => NtcMode::AllBonds,
            2 => NtcMode::HydrogenBonds,
            _ => NtcMode::SolventOnly,
        };
        Some(ShakeParameters {
            tolerance: shake_tolerance,
            max_iterations: 1000,
            ntc: ntc_mode,
        })
    } else {
        None
    };

    if let Some(ref params) = shake_params {
        println!("Setting up constraints: SHAKE");
        println!("  NTC mode:      {:?}", params.ntc);
        println!("  Tolerance:     {:.6}", params.tolerance);
        println!("  Max iter:      {}", params.max_iterations);
        println!();
    }
    if settle_enabled {
        println!("Setting up constraints: SETTLE (analytical rigid water)");
        println!();
    }
    if lincs_enabled {
        println!("Setting up constraints: LINCS (solute={}, solvent={}, order={}/{})",
            solute_lincs, solvent_lincs, imd.lincs_order_solute, imd.lincs_order_solvent);
        println!();
    }

    // Setup trajectory writer
    let mut traj_writer = match TrajectoryWriter::new(
        &trc_file,
        "GROMOS-RS MD trajectory",
        false, // velocities
        false, // forces
    ) {
        Ok(w) => w,
        Err(e) => {
            eprintln!("Error creating trajectory file: {}", e);
            process::exit(1);
        },
    };

    // Setup energy writer
    let mut ene_writer = match EnergyWriter::new(&tre_file, "GROMOS-RS MD energies") {
        Ok(w) => w,
        Err(e) => {
            eprintln!("Error creating energy file: {}", e);
            process::exit(1);
        },
    };

    // Setup force trajectory writer (if @trf given)
    let mut force_writer = if let Some(ref trf_file) = md_args.trf_file {
        match ForceWriter::new(trf_file, "GROMOS-RS MD forces", true) {
            Ok(w) => {
                println!("  Force output:  {}", trf_file);
                Some(w)
            },
            Err(e) => {
                eprintln!("Error creating force file: {}", e);
                process::exit(1);
            },
        }
    } else {
        None
    };

    // Setup GaMD writers if enabled
    let mut gamd_stats_writer = if gamd_params.is_some() {
        let stats_file = "gamd_stats.dat";
        match GamdStatsWriter::new(stats_file, "GROMOS-RS GaMD Statistics") {
            Ok(mut w) => {
                w.set_write_interval(10); // Write every 10 steps
                println!("  GaMD stats:   {}", stats_file);
                Some(w)
            },
            Err(e) => {
                eprintln!("Warning: Failed to create GaMD stats file: {}", e);
                None
            },
        }
    } else {
        None
    };

    let mut gamd_boost_writer = if gamd_params.is_some() {
        let boost_file = "gamd_boost.dat";
        match GamdBoostWriter::new(boost_file, "GROMOS-RS GaMD Boost Potential") {
            Ok(w) => {
                println!("  GaMD boost:   {}", boost_file);
                Some(w)
            },
            Err(e) => {
                eprintln!("Warning: Failed to create GaMD boost file: {}", e);
                None
            },
        }
    } else {
        None
    };

    // Setup EDS writers if enabled
    let mut eds_stats_writer = if eds_params.is_some() {
        let stats_file = "eds_stats.dat";
        match EdsStatsWriter::new(stats_file, "GROMOS-RS EDS Statistics") {
            Ok(mut w) => {
                w.set_write_interval(10); // Write every 10 steps
                println!("  EDS stats:    {}", stats_file);
                Some(w)
            },
            Err(e) => {
                eprintln!("Warning: Failed to create EDS stats file: {}", e);
                None
            },
        }
    } else {
        None
    };

    let mut eds_vr_writer = if eds_params.is_some() {
        let vr_file = "eds_vr.dat";
        match EdsVrWriter::new(vr_file, "GROMOS-RS EDS Reference Energy") {
            Ok(w) => {
                println!("  EDS V_R:      {}", vr_file);
                Some(w)
            },
            Err(e) => {
                eprintln!("Warning: Failed to create EDS V_R file: {}", e);
                None
            },
        }
    } else {
        None
    };

    // Parameters summary
    if is_minimization {
        println!("Energy Minimization Parameters:");
        println!("  Max steps:     {}", n_steps);
        println!("  Tolerance:     {} kJ/mol", imd.dele);
        println!("  Step size:     {} nm (max {})", imd.dx0, imd.dxm);
        println!("  Traj output:   {}", trc_file);
        println!("  Energy output: {}", tre_file);
        println!();

        println!("╔══════════════════════════════════════════════════════════════╗");
        println!("║              Starting Energy Minimization                    ║");
        println!("╚══════════════════════════════════════════════════════════════╝");
        println!();

        log::info!(
            "Starting energy minimization: max {} steps, dele={} kJ/mol",
            n_steps,
            imd.dele
        );
    } else {
        println!("MD Parameters:");
        println!("  Steps:         {}", n_steps);
        println!("  Time step:     {} ps", dt);
        println!(
            "  Total time:    {} ps",
            n_steps as f64 * dt
        );
        println!("  Temperature:   {} K", temperature);
        println!("  Traj output:   {}", trc_file);
        println!("  Energy output: {}", tre_file);
        println!();

        println!("╔══════════════════════════════════════════════════════════════╗");
        println!("║                   Starting MD Simulation                     ║");
        println!("╚══════════════════════════════════════════════════════════════╝");
        println!();

        log::info!(
            "Starting MD simulation: {} steps, dt={} ps",
            n_steps,
            dt
        );
    }

    let start_time = Instant::now();
    let init_elapsed = t_total.elapsed();
    log::info!("Initialization wall time: {:.3} s", init_elapsed.as_secs_f64());
    let mut energy_history: Vec<(f64, f64, f64)> = Vec::new();
    let mut minimization_converged = false;
    let mut prev_min_energy = f64::MAX;
    let mut actual_steps: usize = 0;

    // Main simulation loop
    for step in 0..=n_steps {
        let time = step as f64 * dt;

        log::debug!("Step {}: time = {:.6} ps", step, time);

        // Run the algorithm sequence for this step
        md_sequence.run_step(&topo, &mut conf, &sim_state).unwrap_or_else(|e| {
            eprintln!("Error at step {}: {}", step, e);
            process::exit(1);
        });

        // Debug: dump forces at step 0
        if step == 0 && md_args.verbose >= 2 {
            log::debug!("=== Forces at step 0 (old state, after exchange) ===");
            for i in 0..topo.num_atoms() {
                let f = conf.old().force[i];
                log::debug!("  Atom {:2}: ({:18.9}, {:18.9}, {:18.9})", i+1, f.x, f.y, f.z);
            }
        }

        // Apply GAMD boost if enabled
        if let Some(ref mut gamd) = gamd_params {
            let dihedral_energy = 0.0; // TODO: Separate dihedral energy from bonded
            let total_potential = conf.current().energies.potential_total;

            // Update GAMD statistics
            gamd.update_statistics(dihedral_energy, total_potential);

            // Update parameters if in GaMD search mode
            if gamd.search_mode == gromos::gamd::SearchMode::GamdSearch {
                gamd.update_all_parameters();
            }

            // Apply boost potential
            // For now, use simplified version with only total potential boost
            // TODO: Separate dihedral forces for dual boost
            let dihedral_forces = vec![Vec3::ZERO; topo.num_atoms()];
            let other_forces: Vec<Vec3> = conf.current().force.clone();

            let boost = gamd.apply_boost(
                &mut conf,
                dihedral_energy,
                &dihedral_forces,
                total_potential,
                &other_forces,
            );

            // Add boost to potential energy
            conf.current_mut().energies.potential_total += boost;

            log::debug!(
                "GAMD boost applied: ΔV={:.4}, V_new={:.4}",
                boost,
                conf.current().energies.potential_total
            );

            // Write GaMD statistics
            if let Some(ref mut writer) = gamd_stats_writer {
                if let Err(e) = writer.write_frame(step, dihedral_energy, total_potential, gamd) {
                    log::warn!("Failed to write GaMD stats at step {}: {}", step, e);
                }
            }

            // Write GaMD boost potential
            if let Some(ref mut writer) = gamd_boost_writer {
                // For now, assume boost is only from total potential
                let boost_dih = 0.0;
                let boost_pot = boost;
                if let Err(e) = writer.write_frame(step, time, boost, boost_dih, boost_pot) {
                    log::warn!("Failed to write GaMD boost at step {}: {}", step, e);
                }
            }
        }

        // Apply EDS if enabled
        if let Some(ref mut aeds) = eds_params {
            let current_potential = conf.current().energies.potential_total;

            // For EDS, we need to:
            // 1. Update each state's energy (for simplicity, using current potential as state 0)
            // 2. Calculate reference energy V_R
            // 3. Apply EDS forces

            // TODO: This is a simplified implementation
            // Full EDS requires multiple state energies calculated simultaneously
            // For now, we'll just set up the framework

            let eds = &mut aeds.eds;

            // Update state 0 with current potential
            if eds.num_states > 0 {
                eds.states[0].energy = current_potential;
            }

            // Calculate reference energy based on form
            match eds.form {
                gromos::eds::EDSForm::SingleS => eds.calculate_reference_energy_single_s(),
                gromos::eds::EDSForm::MultiS => eds.calculate_reference_energy_multi_s(),
                gromos::eds::EDSForm::PairS => {
                    log::warn!("PairS EDS form not yet fully implemented");
                    eds.calculate_reference_energy_single_s();
                },
            }

            // Apply EDS forces (modifies configuration)
            eds.apply_forces(&mut conf);

            // Replace potential energy with reference energy
            conf.current_mut().energies.potential_total = eds.reference_energy;

            log::debug!(
                "EDS applied: V_R={:.4}, original V={:.4}",
                eds.reference_energy,
                current_potential
            );

            // Update AEDS parameters if search enabled
            if aeds.search_enabled {
                // TODO: Implement AEDS parameter updates
                log::debug!("AEDS search mode active");
            }

            // Write EDS statistics
            if let Some(ref mut writer) = eds_stats_writer {
                if let Err(e) = writer.write_frame(step, eds) {
                    log::warn!("Failed to write EDS stats at step {}: {}", step, e);
                }
            }

            // Write EDS reference energy
            if let Some(ref mut writer) = eds_vr_writer {
                if let Err(e) = writer.write_frame(step, time, eds.reference_energy) {
                    log::warn!("Failed to write EDS V_R at step {}: {}", step, e);
                }
            }
        }

        // Validate energy
        // gromosXX convention: after the algorithm sequence, energies are in old()
        // (Forcefield wrote to current(), exchange_state moved it to old(),
        //  Temperature_Calculation and Energy_Calculation also write to old())
        let state = conf.old();
        let n_dof = topo.inverse_mass.len() * 3;
        let temp = state.temperature(n_dof);
        let ene_validation = validate_energy(
            state.energies.kinetic_total,
            state.energies.potential_total,
            state.energies.total(),
            temp,
        );
        if ene_validation.has_errors() && md_args.verbose > 0 {
            ene_validation.print();
            log::warn!("Energy validation failed at step {}", step);
        }

        // Store energy for drift check
        energy_history.push((
            state.energies.kinetic_total,
            state.energies.potential_total,
            state.energies.total(),
        ));

        actual_steps = step;

        // Track per-step energy change for minimization display
        let min_de = if is_minimization {
            let current_pot = state.energies.potential_total + state.energies.special_total;
            let de = if prev_min_energy < f64::MAX / 2.0 {
                current_pot - prev_min_energy
            } else {
                0.0
            };
            prev_min_energy = current_pot;
            de
        } else {
            0.0
        };

        // Minimization convergence check (mirrors gromosXX: compare potential_total + special_total)
        if is_minimization && step > imd.nmin {
            if min_de.abs() < imd.dele {
                minimization_converged = true;
                println!();
                println!("*** Energy minimization CONVERGED at step {} ***", step);
                println!("    dE = {:.6e} kJ/mol (tolerance: {:.6e})", min_de.abs(), imd.dele);
                println!("    E_pot = {:.10e} kJ/mol", state.energies.potential_total);
                println!();
                // Write final frame before breaking
                if let Err(e) = traj_writer.write_frame(step, time, &conf) {
                    eprintln!("Error writing final trajectory frame: {}", e);
                }
                break;
            }
        }

        // Log progress
        if step % nstlog == 0 {
            if is_minimization {
                println!("Step {:6}  E_pot: {:18.10e}  dE: {:12.4e}",
                    step, state.energies.potential_total, min_de);
            } else {
                println!("Step {:6}  Time: {:8.3} ps  E_pot: {:18.10e}  E_kin: {:18.10e}  E_tot: {:18.10e}  T: {:6.1} K",
                    step, time, state.energies.potential_total, state.energies.kinetic_total,
                    state.energies.total(), temp);
            }
            log::debug!(
                "  Bond: {:.10e}  LJ: {:.10e}  CRF: {:.10e}",
                state.energies.bond_total,
                state.energies.lj_total,
                state.energies.crf_total
            );
        }

        // Write trajectory
        if step % nstxout == 0 {
            if let Err(e) = traj_writer.write_frame(step, time, &conf) {
                eprintln!("Error writing trajectory: {}", e);
            }
        }

        // Write energies
        if step % nstener == 0 {
            let volume = conf.old().box_config.volume();
            let pressure = conf.old().pressure();
            let energies = conf.old().energies.clone();

            let ene_frame = EnergyFrame {
                time,
                kinetic: energies.kinetic_total,
                potential: energies.potential_total,
                total: energies.total(),
                temperature: temp,
                volume,
                pressure,
                bond: energies.bond_total,
                angle: 0.0,
                improper: 0.0,
                dihedral: 0.0,
                lj: energies.lj_total,
                coul_real: energies.crf_total, // CRF energy
                coul_recip: 0.0,
                coul_self: 0.0,
                shake: 0.0,
                restraint: 0.0,
                extra: Vec::new(),
            };

            if let Err(e) = ene_writer.write_frame(&ene_frame) {
                eprintln!("Error writing energy: {}", e);
            }
        }

        // Write forces (same frequency as trajectory)
        if let Some(ref mut fw) = force_writer {
            if step % nstxout == 0 {
                let forces = &conf.old().force;
                let constraint_forces = if shake_enabled || settle_enabled || lincs_enabled {
                    Some(conf.old().constraint_force.as_slice())
                } else {
                    None
                };
                if let Err(e) = fw.write_frame(step, time, forces, constraint_forces) {
                    eprintln!("Error writing forces: {}", e);
                }
            }
        }

        // Advance simulation state for next step
        sim_state.advance();
    }

    if is_minimization {
        if minimization_converged {
            log::info!("Energy minimization converged at step {}", actual_steps);
        } else {
            log::warn!("Energy minimization did NOT converge within {} steps", n_steps);
            println!();
            println!("*** WARNING: Energy minimization did NOT converge ***");
            println!("    Max steps reached: {}", n_steps);
            println!("    Consider increasing NSTLIM or adjusting DELE");
            println!();
        }
    } else {
        log::info!("MD loop completed - {} steps", n_steps);
    }

    // Check energy drift
    log::debug!("Checking energy drift over trajectory");
    use gromos::validation::check_energy_drift;
    let drift_report = check_energy_drift(&energy_history);
    if drift_report.has_errors() || !drift_report.warnings.is_empty() {
        drift_report.print();
        drift_report.print_summary();
    } else {
        log::debug!("Energy drift check passed");
    }

    // Finalize output files
    log::debug!("Finalizing output files");
    if let Err(e) = traj_writer.flush() {
        log::error!("Failed to flush trajectory: {}", e);
        eprintln!("Error flushing trajectory: {}", e);
    } else {
        log::debug!("Trajectory file finalized: {}", trc_file);
    }

    if let Err(e) = ene_writer.finalize() {
        log::error!("Failed to finalize energy file: {}", e);
        eprintln!("Error finalizing energy file: {}", e);
    } else {
        log::debug!("Energy file finalized: {}", tre_file);
    }

    if let Some(ref mut fw) = force_writer {
        if let Err(e) = fw.flush() {
            log::error!("Failed to flush force file: {}", e);
            eprintln!("Error flushing force file: {}", e);
        } else if let Some(ref trf_file) = md_args.trf_file {
            log::debug!("Force file finalized: {}", trf_file);
        }
    }

    // Finalize GaMD writers if enabled
    if let Some(ref mut writer) = gamd_stats_writer {
        if let Err(e) = writer.finalize() {
            log::error!("Failed to finalize GaMD stats file: {}", e);
            eprintln!("Error finalizing GaMD stats file: {}", e);
        } else {
            log::debug!("GaMD stats file finalized");
        }
    }

    if let Some(ref mut writer) = gamd_boost_writer {
        if let Err(e) = writer.finalize() {
            log::error!("Failed to finalize GaMD boost file: {}", e);
            eprintln!("Error finalizing GaMD boost file: {}", e);
        } else {
            log::debug!("GaMD boost file finalized");
        }
    }

    // Finalize EDS writers if enabled
    if let Some(ref mut writer) = eds_stats_writer {
        if let Err(e) = writer.finalize() {
            log::error!("Failed to finalize EDS stats file: {}", e);
            eprintln!("Error finalizing EDS stats file: {}", e);
        } else {
            log::debug!("EDS stats file finalized");
        }
    }

    if let Some(ref mut writer) = eds_vr_writer {
        if let Err(e) = writer.finalize() {
            log::error!("Failed to finalize EDS V_R file: {}", e);
            eprintln!("Error finalizing EDS V_R file: {}", e);
        } else {
            log::debug!("EDS V_R file finalized");
        }
    }

    // Write final configuration if @fin was specified
    if let Some(ref fin_path) = md_args.fin_file {
        let positions = &conf.current().pos;
        let velocities = &conf.current().vel;
        let box_vec = conf.current().box_config.dimensions();
        let vels = if velocities.iter().any(|v| *v != Vec3::ZERO) {
            Some(velocities.as_slice())
        } else {
            None
        };
        let box_opt = if box_vec.x > 0.0 { Some(box_vec) } else { None };
        let title = format!("Final configuration after {} steps", actual_steps);
        if let Err(e) = gromos::io::g96::write_g96(fin_path, &title, positions, vels, box_opt, Some(&topo)) {
            eprintln!("Error writing final configuration: {}", e);
        } else {
            log::info!("Final configuration written to: {}", fin_path);
        }
    }

    let elapsed = start_time.elapsed();
    let total_elapsed = t_total.elapsed();
    log::info!("Simulation wall time: {:.2} s", elapsed.as_secs_f64());
    log::info!("Total wall time: {:.3} s (init: {:.3} s, sim: {:.3} s)",
        total_elapsed.as_secs_f64(), init_elapsed.as_secs_f64(), elapsed.as_secs_f64());

    println!();
    println!("╔══════════════════════════════════════════════════════════════╗");
    if is_minimization {
        println!("║              Energy Minimization Complete                    ║");
    } else {
        println!("║                   Simulation Complete                        ║");
    }
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();
    println!("Statistics:");
    println!("  Total steps:     {}", actual_steps);
    if is_minimization {
        println!("  Converged:       {}", if minimization_converged { "YES" } else { "NO" });
    } else {
        println!(
            "  Simulation time: {:.3} ps",
            n_steps as f64 * dt
        );
    }
    println!("  Wall time:       {:.2} s", elapsed.as_secs_f64());
    println!(
        "  Performance:     {:.1} ns/day",
        (n_steps as f64 * dt * 1e-3) / elapsed.as_secs_f64() * 86400.0
    );
    println!();
    println!("Output files:");
    if let Some(ref fin_path) = md_args.fin_file {
        println!("  Final conf: {}", fin_path);
    }
    println!("  Trajectory: {}", trc_file);
    println!("  Energies:   {}", tre_file);
    println!();
    println!("Done!");
}
