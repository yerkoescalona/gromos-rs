//! GROMOS-RS CLI: Unified command-line interface
//!
//! This binary provides a single entry point for all GROMOS tools
//! using clap's multicall feature for BusyBox-style operation.
//!
//! ## Usage
//!
//! ```bash
//! # Run as subcommand
//! gromos md -c input.imd -t system.top
//!
//! # Or via symlink (multicall mode)
//! md -c input.imd -t system.top
//!
//! # Install symlinks for all tools
//! gromos --install ~/.local/bin
//! ```

use clap::{Command, Parser, Subcommand};
use std::path::Path;

mod commands;

#[derive(Parser)]
#[command(
    name = "gromos",
    author = "GROMOS Developers",
    version,
    about = "GROMOS-RS molecular dynamics toolkit",
    long_about = "Unified CLI for GROMOS molecular dynamics simulation and analysis tools.\n\n\
                  This single binary provides access to 100+ tools that were previously \
                  separate executables. Use --install to create symlinks for backward compatibility."
)]
struct Cli {
    /// Install symlinks for all tools to the specified directory
    #[arg(long, value_name = "PATH")]
    install: Option<std::path::PathBuf>,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    // ============================================================
    // Simulation commands
    // ============================================================
    /// Run molecular dynamics simulation
    #[cfg(feature = "cmd-md")]
    Md(commands::simulation::MdArgs),

    /// Energy minimization
    #[cfg(feature = "cmd-minimize")]
    Minimize(commands::simulation::MinimizeArgs),

    // ============================================================
    // Analysis commands
    // ============================================================
    /// Calculate radial distribution function
    #[cfg(feature = "cmd-rdf")]
    Rdf(commands::analysis::RdfArgs),

    /// Calculate RMSD
    #[cfg(feature = "cmd-rmsd")]
    Rmsd(commands::analysis::RmsdArgs),

    /// Hydrogen bond analysis
    #[cfg(feature = "cmd-hbond")]
    Hbond(commands::analysis::HbondArgs),

    /// Radius of gyration
    #[cfg(feature = "cmd-gyrate")]
    Gyrate(commands::analysis::GyrateArgs),
}

fn main() {
    // Check if called via symlink (multicall mode)
    let args: Vec<String> = std::env::args().collect();
    let program_name = Path::new(&args[0])
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("gromos");

    // If not called as "gromos", treat program name as command
    if program_name != "gromos" {
        dispatch_multicall(program_name, &args[1..]);
        return;
    }

    // Normal CLI parsing
    let cli = Cli::parse();

    if let Some(install_path) = cli.install {
        install_symlinks(&install_path).expect("Failed to install symlinks");
        return;
    }

    match cli.command {
        Some(cmd) => dispatch_command(cmd),
        None => {
            // Print help if no command given
            Cli::parse_from(["gromos", "--help"]);
        }
    }
}

fn dispatch_command(cmd: Commands) {
    match cmd {
        #[cfg(feature = "cmd-md")]
        Commands::Md(args) => commands::simulation::run_md(args),

        #[cfg(feature = "cmd-minimize")]
        Commands::Minimize(args) => commands::simulation::run_minimize(args),

        #[cfg(feature = "cmd-rdf")]
        Commands::Rdf(args) => commands::analysis::run_rdf(args),

        #[cfg(feature = "cmd-rmsd")]
        Commands::Rmsd(args) => commands::analysis::run_rmsd(args),

        #[cfg(feature = "cmd-hbond")]
        Commands::Hbond(args) => commands::analysis::run_hbond(args),

        #[cfg(feature = "cmd-gyrate")]
        Commands::Gyrate(args) => commands::analysis::run_gyrate(args),
    }
}

fn dispatch_multicall(cmd: &str, args: &[String]) {
    // Map legacy/short names to commands
    match cmd {
        "md" | "gromos-md" => {
            #[cfg(feature = "cmd-md")]
            {
                let full_args = std::iter::once("gromos".to_string())
                    .chain(std::iter::once("md".to_string()))
                    .chain(args.iter().cloned());
                let cli = Cli::parse_from(full_args);
                if let Some(Commands::Md(args)) = cli.command {
                    commands::simulation::run_md(args);
                }
            }
        }
        "rdf" | "g_rdf" => {
            #[cfg(feature = "cmd-rdf")]
            {
                let full_args = std::iter::once("gromos".to_string())
                    .chain(std::iter::once("rdf".to_string()))
                    .chain(args.iter().cloned());
                let cli = Cli::parse_from(full_args);
                if let Some(Commands::Rdf(args)) = cli.command {
                    commands::analysis::run_rdf(args);
                }
            }
        }
        // Add more tool mappings here
        _ => {
            eprintln!("Unknown command: {}", cmd);
            std::process::exit(1);
        }
    }
}

/// Install symlinks for all tools
fn install_symlinks(install_path: &Path) -> std::io::Result<()> {
    let current_exe = std::env::current_exe()?;
    
    // Create directory if it doesn't exist
    std::fs::create_dir_all(install_path)?;

    // List of all tools to create symlinks for
    let tools = [
        // Simulation
        "md", "minimize",
        // Analysis
        "rdf", "rmsd", "hbond", "gyrate",
        // Add all 104 tools here...
    ];

    for tool in tools {
        let link_path = install_path.join(tool);
        if link_path.exists() {
            std::fs::remove_file(&link_path)?;
        }
        #[cfg(unix)]
        std::os::unix::fs::symlink(&current_exe, &link_path)?;
        println!("Created symlink: {}", link_path.display());
    }

    // Legacy GROMOS name aliases
    let legacy = [
        ("g_rdf", "rdf"),
        ("g_rmsd", "rmsd"),
        ("g_hbond", "hbond"),
    ];

    for (old_name, _new_name) in legacy {
        let link_path = install_path.join(old_name);
        if link_path.exists() {
            std::fs::remove_file(&link_path)?;
        }
        #[cfg(unix)]
        std::os::unix::fs::symlink(&current_exe, &link_path)?;
        println!("Created legacy symlink: {}", link_path.display());
    }

    println!("\nSymlinks installed to: {}", install_path.display());
    println!("Add to PATH: export PATH=\"{}:$PATH\"", install_path.display());

    Ok(())
}
