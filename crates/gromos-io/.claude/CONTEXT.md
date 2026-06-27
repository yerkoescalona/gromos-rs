# gromos-io — stage contract

## Job
L0 I/O service. All file format parsing and writing. The only place in the workspace where
GROMOS file formats are read or written.

## Inputs (consumes from)
`gromos-core` types only.

## Outputs (public API)
Parsers + writers for every GROMOS file format used in the workspace.

## Status
- Topology parser ✓: all blocks (SOLUTEATOM, SOLVENTCONSTR, LJPARAMETERS, PRESSUREGROUPS, ...)
- Coordinate reader/writer ✓: POSITION, POSITIONRED, VELOCITYRED, GENBOX
- Trajectory ✓: `TrajectoryReader` (3-column and legacy 7-column POSITIONRED; GENBOX optional for vacuum),
  `TrajectoryWriter::write_frame` (via `write_trc_frame`), `write_trc_frame(step, time, positions, box)`
- Energy/force writers ✓: ENERTRJ, FREEFORCERED/CONSFORCERED
- PDB: `write_pdb_positions(path, title, positions, box, topology)` ✓
- G96: `write_g96(path, title, positions, velocities, box, topology)` ✓
- Free energy: `FreeEnergyWriter` (FREEENERGY03), `read_free_energy_trajectory` ✓
- MTB, IFP, IMD, PTP parsers ✓; position restraints (.por/.rpr) ✓

## Key files
```
src/topology.rs    — topology file parser (~1200 LOC; P4 split candidate)
src/trajectory.rs  — TrajectoryReader + TrajectoryWriter (write_trc_frame)
src/g96.rs         — write_g96, G96Writer
src/pdb.rs         — write_pdb_positions
src/free_energy.rs — FreeEnergyWriter, read_free_energy_trajectory
src/imd.rs         — IMD parameter file parser
```

## Decisions logged
- **POSITIONRED format:** standard is 3-column (x y z); reader accepts both 3-col and legacy 7-col.
  `write_trc_frame` always writes 3-column (standard GROMOS format).
- **GENBOX optional:** vacuum trajectories omit GENBOX; reader returns `Vec3::ZERO` when missing.
- **PERTURBATION block is multi-line:** splits NTG/NRDGL/RLAM/DLAMT on line 1 and ALPHLJ/ALPHC/NLAM/NSCALE
  on line 2. Always combine `data_lines` before tokenising. Bug: `data_lines.first()` silently dropped line 2.
- **Energy output:** full f64 scientific notation for exact comparison.

## Crate-specific rules
- **ALL parsing/writing lives here.** No file format code in binaries or other crates.
- **All analysis binaries** use `write_g96`, `write_pdb_positions`, `TrajectoryWriter::write_trc_frame`
  — never hand-roll `BufWriter`/`writeln!` for GROMOS formats in binary entry points.
