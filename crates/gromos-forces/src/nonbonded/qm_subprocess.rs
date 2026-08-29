//! Shared subprocess mechanics for file-in/file-out QM engine providers ([`super::xtb`],
//! [`super::mopac`]). Crate-private — not a new public trait. Each engine's actual physics (input
//! format, output parsing, units) stays engine-specific in its own module; only the mechanical
//! parts (work-dir setup, stale-output removal, spawn, exit-status check) are common enough to
//! share. See `architecture.md` rule 3 ("a small, stable trait taxonomy — resist proliferation")
//! for why this is a helper module and not a `QmEngine` trait mirroring gromosXX's `QM_Worker`
//! base class.

use crate::provider::ProviderError;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Create `work_dir` if missing (idempotent — the directory is reused, overwritten, on every
/// call, no temp-directory dependency).
pub(crate) fn ensure_work_dir(work_dir: impl Into<PathBuf>) -> Result<PathBuf, ProviderError> {
    let work_dir = work_dir.into();
    std::fs::create_dir_all(&work_dir)
        .map_err(|e| ProviderError::ComputationFailed(format!("failed to create work_dir: {e}")))?;
    Ok(work_dir)
}

/// Remove a possibly-stale output file before running, so a crashed engine can't leave a prior
/// result silently readable as this call's answer. Best-effort: absence is not an error.
pub(crate) fn remove_stale(path: &Path) {
    let _ = std::fs::remove_file(path);
}

/// Run `binary` with `args` in `work_dir`; a spawn failure or nonzero exit status becomes a
/// `ProviderError` carrying stderr.
pub(crate) fn run_subprocess(
    binary: &str,
    work_dir: &Path,
    args: &[&str],
) -> Result<(), ProviderError> {
    let output = Command::new(binary)
        .current_dir(work_dir)
        .args(args)
        .output()
        .map_err(|e| ProviderError::ComputationFailed(format!("failed to run {binary}: {e}")))?;
    if !output.status.success() {
        return Err(ProviderError::ComputationFailed(format!(
            "{binary} exited with {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr)
        )));
    }
    Ok(())
}
