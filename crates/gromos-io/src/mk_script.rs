//! IMD file text-level manipulation for mk_script.
//!
//! mk_script reads a base .imd file and overrides specific parameter values
//! per job. Parameters are identified by name (matching the GROMOS convention)
//! and located within their IMD block.
//!
//! The override mechanism works at the text level: each IMD block has a known
//! layout mapping parameter names to (data_line_index, field_index) positions.
//! This avoids needing a full IMD serializer while supporting the flexible
//! .jobs column scheme.

use crate::jobs::JobSpec;
use std::collections::HashMap;
use std::fmt::Write as FmtWrite;

/// Mapping of IMD parameter names to their location in the IMD text.
///
/// Each entry: (block_name, data_line_index (1-based), field_index (0-based))
/// This covers the most commonly overridden parameters. Parameters that use
/// array notation like TEMP0[1] are handled specially.
fn param_location(name: &str) -> Option<(&'static str, usize, usize)> {
    // (block, data_line, field_index)
    Some(match name {
        // INITIALISE block: 1 data line with 10 fields
        "NTIVEL" => ("INITIALISE", 1, 0),
        "NTISHK" => ("INITIALISE", 1, 1),
        "NTINHT" => ("INITIALISE", 1, 2),
        "NTINHB" => ("INITIALISE", 1, 3),
        "NTISHI" => ("INITIALISE", 1, 4),
        "NTIRTC" => ("INITIALISE", 1, 5),
        "NTICOM" => ("INITIALISE", 1, 6),
        "NTISTI" => ("INITIALISE", 1, 7),
        "IG" => ("INITIALISE", 1, 8),
        "TEMPI" => ("INITIALISE", 1, 9),

        // STEP block: 1 data line
        "NSTLIM" => ("STEP", 1, 0),
        "T" => ("STEP", 1, 1),
        "DT" => ("STEP", 1, 2),

        // BOUNDCOND block: 1 data line
        "NTB" => ("BOUNDCOND", 1, 0),
        "NDFMIN" => ("BOUNDCOND", 1, 1),

        // POSITIONRES block: 1 data line
        "NTPOR" => ("POSITIONRES", 1, 0),
        "NTPORB" => ("POSITIONRES", 1, 1),
        "NTPORS" => ("POSITIONRES", 1, 2),
        "CPOR" => ("POSITIONRES", 1, 3),

        // ENERGYMIN block: 1 data line
        "NTEM" => ("ENERGYMIN", 1, 0),
        "NCYC" => ("ENERGYMIN", 1, 1),
        "DELE" => ("ENERGYMIN", 1, 2),
        "DX0" => ("ENERGYMIN", 1, 3),
        "DXM" => ("ENERGYMIN", 1, 4),
        "NMIN" => ("ENERGYMIN", 1, 5),
        "FLIM" => ("ENERGYMIN", 1, 6),

        // COMTRANSROT block: 1 data line
        "NSCM" => ("COMTRANSROT", 1, 0),

        // CONSTRAINT block
        "NTC" => ("CONSTRAINT", 1, 0),

        // NONBONDED block: 1 data line
        "NLRELE" => ("NONBONDED", 1, 0),

        // WRITETRAJ block: 1 data line
        "NTWX" => ("WRITETRAJ", 1, 0),
        "NTWSE" => ("WRITETRAJ", 1, 1),
        "NTWV" => ("WRITETRAJ", 1, 2),
        "NTWF" => ("WRITETRAJ", 1, 3),
        "NTWE" => ("WRITETRAJ", 1, 4),
        "NTWG" => ("WRITETRAJ", 1, 5),
        "NTWB" => ("WRITETRAJ", 1, 6),

        // PRINTOUT block: 1 data line
        "NTPR" => ("PRINTOUT", 1, 0),

        // PRESSURESCALE/BAROSTAT
        "NTP" => ("PRESSURESCALE", 1, 0),
        "COUPLE" => ("PRESSURESCALE", 1, 0), // alias

        _ => return None,
    })
}

/// Apply job-specific overrides to a base IMD file (as text lines).
///
/// Parameters are matched by name to their block/field location.
/// Array parameters like TEMP0[N] and PRSBTH[N] are handled specially.
///
/// Returns the modified IMD content as a string.
pub fn apply_job_overrides(base_imd: &str, job: &JobSpec) -> String {
    // Build a per-block, per-data-line override index
    // Key: (block_name, data_line_idx) → Vec<(field_idx, new_value)>
    let mut overrides: HashMap<(String, usize), Vec<(usize, String)>> = HashMap::new();

    // Collect TEMP0[N] overrides for special MULTIBATH handling
    let mut temp0_overrides: HashMap<usize, String> = HashMap::new();

    for (name, value) in &job.params {
        // Handle TEMP0[N] array parameters
        if let Some(rest) = name.strip_prefix("TEMP0[") {
            if let Some(idx_str) = rest.strip_suffix(']') {
                if let Ok(idx) = idx_str.parse::<usize>() {
                    temp0_overrides.insert(idx, value.clone());
                    continue;
                }
            }
        }

        // Handle standard parameters via the location table
        if let Some((block, line_idx, field_idx)) = param_location(name) {
            overrides
                .entry((block.to_string(), line_idx))
                .or_default()
                .push((field_idx, value.clone()));
        }
    }

    let mut output = String::with_capacity(base_imd.len());
    let mut current_block = String::new();
    let mut data_line_idx: usize = 0;

    for line in base_imd.lines() {
        let trimmed = line.trim();

        if trimmed == "END" {
            current_block.clear();
            data_line_idx = 0;
            writeln!(output, "{}", line).unwrap();
            continue;
        }

        // Detect block start
        if !trimmed.is_empty()
            && !trimmed.starts_with('#')
            && trimmed
                .chars()
                .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit() || c == '_')
        {
            current_block = trimmed.to_string();
            data_line_idx = 0;
            writeln!(output, "{}", line).unwrap();
            continue;
        }

        if trimmed.starts_with('#') || trimmed.is_empty() {
            writeln!(output, "{}", line).unwrap();
            continue;
        }

        data_line_idx += 1;

        // Special handling for MULTIBATH TEMP0/TAU line
        if current_block == "MULTIBATH" && data_line_idx == 3 && !temp0_overrides.is_empty() {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            let n_baths = parts.len() / 2;
            let mut new_line = String::new();
            for i in 0..n_baths {
                let temp = if let Some(v) = temp0_overrides.get(&(i + 1)) {
                    v.clone()
                } else {
                    parts[i * 2].to_string()
                };
                let tau = parts[i * 2 + 1];
                write!(new_line, "{:>10}{:>8}", temp, tau).unwrap();
            }
            writeln!(output, "{}", new_line).unwrap();
            continue;
        }

        // Check if this data line has overrides
        let key = (current_block.clone(), data_line_idx);
        if let Some(field_overrides) = overrides.get(&key) {
            let mut parts: Vec<String> = trimmed
                .split_whitespace()
                .map(|s| s.to_string())
                .collect();
            for (field_idx, value) in field_overrides {
                if *field_idx < parts.len() {
                    parts[*field_idx] = value.clone();
                }
            }
            // Reconstruct with reasonable spacing
            let new_line: String = parts.iter().map(|p| format!("{:>12}", p)).collect();
            writeln!(output, "{}", new_line).unwrap();
            continue;
        }

        writeln!(output, "{}", line).unwrap();
    }

    output
}

/// Configuration for mk_script generation.
#[derive(Debug, Clone)]
pub struct MkScriptConfig {
    /// System name (used in filename templates)
    pub system: String,
    /// Path to MD binary
    pub bin: String,
    /// Working directory
    pub dir: String,
    /// Input files: key → path
    pub files: HashMap<String, String>,
}

/// Generate a shell run script for a single job.
pub fn generate_run_script(
    config: &MkScriptConfig,
    template: &crate::script_template::ScriptTemplate,
    job: &JobSpec,
    prev_job: Option<&JobSpec>,
) -> String {
    let num = job.job_id;
    let sys = &config.system;

    let input_name = template.expanded_filename("input", sys, num);
    let output_name = template.expanded_filename("output", sys, num);
    let fin_name = template.expanded_filename("coord", sys, num);
    let trx_name = template.expanded_filename("outtrx", sys, num);
    let tre_name = template.expanded_filename("outtre", sys, num);

    // Initial coordinates: first job uses @files coord, rest use previous job's output coords
    let coord = if let Some(prev) = prev_job {
        template.expanded_filename("coord", sys, prev.job_id)
    } else {
        config
            .files
            .get("coord")
            .cloned()
            .unwrap_or_else(|| "input.cnf".to_string())
    };

    let topo = config
        .files
        .get("topo")
        .cloned()
        .unwrap_or_else(|| "system.top".to_string());

    let mut script = String::new();
    writeln!(script, "#!/bin/sh").unwrap();
    writeln!(script, "# mk_script generated run script").unwrap();
    writeln!(script, "# Job {} of system {}", num, sys).unwrap();
    writeln!(script).unwrap();

    // MD command
    writeln!(script, "{} \\", config.bin).unwrap();
    writeln!(script, "  @topo {} \\", topo).unwrap();
    writeln!(script, "  @conf {} \\", coord).unwrap();
    writeln!(script, "  @fin {} \\", fin_name).unwrap();
    writeln!(script, "  @input {} \\", input_name).unwrap();
    writeln!(script, "  @trc {} \\", trx_name).unwrap();
    writeln!(script, "  @tre {} \\", tre_name).unwrap();

    // Position restraint files (check if NTPOR is set and > 0)
    let has_posres = job
        .get_i32("NTPOR")
        .map(|v| v > 0)
        .unwrap_or(false);
    if has_posres {
        if let Some(por) = config.files.get("posresspec") {
            writeln!(script, "  @posresspec {} \\", por).unwrap();
        }
        if let Some(rpr) = config.files.get("refpos") {
            writeln!(script, "  @refpos {} \\", rpr).unwrap();
        }
    }

    writeln!(script, "  > {}", output_name).unwrap();

    // Chain to next job
    if let Some(lastcmd) = template.misc.get("lastcommand") {
        let next_num = num + 1;
        let next_script = crate::script_template::ScriptTemplate::expand(lastcmd, sys, next_num);
        writeln!(script).unwrap();
        writeln!(script, "# chain to next job").unwrap();
        writeln!(script, "{}", next_script).unwrap();
    }

    script
}
