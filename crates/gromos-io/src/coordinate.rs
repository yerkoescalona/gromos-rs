//! Parser for GROMOS coordinate files (.conf/.cnf)
//!
//! Format example:
//! ```text
//! TITLE
//!   System description
//! END
//! POSITION
//!     1 RES    ATOM      1    x        y        z
//! END
//! VELOCITY (optional)
//!     1 RES    ATOM      1    vx       vy       vz
//! END
//! BOX
//!     lx   ly   lz
//! END
//! ```

use gromos_core::configuration::{Box as SimBox, Configuration};
use crate::IoError;
use gromos_core::math::Vec3;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Read GROMOS coordinate file
pub fn read_coordinate_file<P: AsRef<Path>>(
    path: P,
    num_temp_groups: usize,
    num_energy_groups: usize,
) -> Result<Configuration, IoError> {
    let file = File::open(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    let reader = BufReader::new(file);

    let mut positions = Vec::new();
    let mut velocities = Vec::new();
    let mut box_dims = Vec3::ZERO;
    let mut in_position = false;
    let mut in_velocity = false;
    let mut in_box = false;

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        // Skip comments and empty lines
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        // Check for block markers
        if trimmed == "POSITION" {
            in_position = true;
            in_velocity = false;
            in_box = false;
            continue;
        } else if trimmed == "VELOCITY" {
            in_position = false;
            in_velocity = true;
            in_box = false;
            continue;
        } else if trimmed == "BOX" {
            in_position = false;
            in_velocity = false;
            in_box = true;
            continue;
        } else if trimmed == "END" {
            in_position = false;
            in_velocity = false;
            in_box = false;
            continue;
        }

        // Parse data
        if in_position {
            // CRITICAL: Skip first 24 characters (residue/atom metadata)
            // Format: "    1 RES    ATOM      1"  (24 chars) then x y z
            if line.len() < 24 {
                continue;
            }

            let coords = &line[24..];
            let parts: Vec<&str> = coords.split_whitespace().collect();

            if parts.len() >= 3 {
                let x: f64 = parts[0].parse().map_err(|_| {
                    IoError::ParseError(format!("Invalid x coordinate: {}", parts[0]))
                })?;
                let y: f64 = parts[1].parse().map_err(|_| {
                    IoError::ParseError(format!("Invalid y coordinate: {}", parts[1]))
                })?;
                let z: f64 = parts[2].parse().map_err(|_| {
                    IoError::ParseError(format!("Invalid z coordinate: {}", parts[2]))
                })?;

                positions.push(Vec3::new(x, y, z));
            }
        } else if in_velocity {
            // Same format as POSITION
            if line.len() < 24 {
                continue;
            }

            let coords = &line[24..];
            let parts: Vec<&str> = coords.split_whitespace().collect();

            if parts.len() >= 3 {
                let vx: f64 = parts[0]
                    .parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid vx: {}", parts[0])))?;
                let vy: f64 = parts[1]
                    .parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid vy: {}", parts[1])))?;
                let vz: f64 = parts[2]
                    .parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid vz: {}", parts[2])))?;

                velocities.push(Vec3::new(vx, vy, vz));
            }
        } else if in_box {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();

            if parts.len() >= 3 {
                let lx: f64 = parts[0]
                    .parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid box x: {}", parts[0])))?;
                let ly: f64 = parts[1]
                    .parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid box y: {}", parts[1])))?;
                let lz: f64 = parts[2]
                    .parse()
                    .map_err(|_| IoError::ParseError(format!("Invalid box z: {}", parts[2])))?;

                box_dims = Vec3::new(lx, ly, lz);
            }
        }
    }

    // Validate
    if positions.is_empty() {
        return Err(IoError::FormatError("No positions found".to_string()));
    }

    let n_atoms = positions.len();

    // Create configuration
    let mut conf = Configuration::new(n_atoms, num_temp_groups, num_energy_groups);

    // Set positions
    conf.current_mut().pos = positions;

    // Set velocities (if present)
    if !velocities.is_empty() {
        if velocities.len() != n_atoms {
            return Err(IoError::FormatError(format!(
                "Velocity count ({}) doesn't match atom count ({})",
                velocities.len(),
                n_atoms
            )));
        }
        conf.current_mut().vel = velocities;
    }

    // Set box
    if box_dims != Vec3::ZERO {
        conf.current_mut().box_config = SimBox::rectangular(box_dims.x, box_dims.y, box_dims.z);
    }

    Ok(conf)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Minimal 4-atom GROMOS .conf content (cg16-style).
    /// Position lines: first 24 chars are residue/atom metadata, then x y z.
    const CG16_CONF: &str = "\
TITLE
  4-atom test system (cg16)
END
POSITION
    1 RES    C         1   4.197156491  0.505921049  2.679733124
    1 RES    H         2   4.300000000  0.600000000  2.700000000
    1 RES    H         3   4.100000000  0.400000000  2.650000000
    1 RES    H         4   4.200000000  0.500000000  2.800000000
END
BOX
  10.0  10.0  10.0
END
";

    fn write_tmp(content: &str, suffix: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!("gromos_test_{suffix}.tmp"));
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn test_parse_cg16() {
        let path = write_tmp(CG16_CONF, "cg16_conf");
        let conf = read_coordinate_file(&path, 1, 1).expect("Failed to parse inline conf");

        assert_eq!(conf.current().pos.len(), 4);

        let pos0 = conf.current().pos[0];
        assert!((pos0.x - 4.197156491_f32).abs() < 1e-5);
        assert!((pos0.y - 0.505921049_f32).abs() < 1e-5);
        assert!((pos0.z - 2.679733124_f32).abs() < 1e-5);

        let box_dims = conf.current().box_config.dimensions();
        assert!((box_dims.x - 10.0_f32).abs() < 1e-5);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn test_parse_missing_file_returns_error() {
        let result = read_coordinate_file("/nonexistent/path/file.conf", 1, 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_conf_with_velocities() {
        let content = "\
TITLE
  test
END
POSITION
    1 RES    C         1   1.000000000  2.000000000  3.000000000
END
VELOCITY
    1 RES    C         1   0.010000000  0.020000000  0.030000000
END
BOX
  5.0  5.0  5.0
END
";
        let path = write_tmp(content, "conf_vel");
        let conf = read_coordinate_file(&path, 1, 1).expect("Failed to parse conf with velocities");

        assert_eq!(conf.current().pos.len(), 1);
        assert_eq!(conf.current().vel.len(), 1);
        assert!((conf.current().vel[0].x - 0.01_f32).abs() < 1e-5);

        std::fs::remove_file(path).ok();
    }
}
