# gromos-io — stage contract

## Job
L0 I/O service. All file format I/O. The only place in the workspace where GROMOS file formats are parsed or written.

## Inputs (consumes from)
`gromos-core` types only.

## Outputs (public API)
Parsed structs for every GROMOS file format used in the workspace.

## Status
- Topology parser ✓: all blocks (SOLUTEATOM, SOLVENTCONSTR, LJPARAMETERS, PRESSUREGROUPS, ...)
- Coordinate reader ✓: POSITION, POSITIONRED, VELOCITY, VELOCITYRED, GENBOX
- Energy/trajectory/force writers ✓: ENERTRJ, POSITIONRED, FREEFORCERED/CONSFORCERED
- MTB parser ✓: MTBUILDBLSOLUTE, MTBUILDBLEND, MTBUILDBLSOLVENT (with continuation lines)
- IFP parser ✓: all type codes + SINGLEATOMLJPAIR + MIXEDATOMLJPAIR
- Position restraints parser ✓ (.por/.rpr)

## Key files
```
src/topology.rs    — topology file parser (~1200 LOC; P4 split candidate)
src/coordinate.rs  — coordinate/GENBOX parser
src/imd.rs         — IMD parameter file parser
```

## Decisions logged
- **Energy output:** full f64 scientific notation for exact comparison.
- **BONDANGLETYPE / HARMBONDANGLETYPE / DIHEDRALTYPE:** intentionally unsupported. These are GROMOS96-era back-compat shims for data already carried by BONDANGLEBENDTYPE/TORSDIHEDRALTYPE in unified form. No .top in the repo corpus uses them; silently ignored by design. Locked via `test_parse_dihedral_and_improper_type_conversions`.
- **Unit-conversion audit:** all conversions in topology.rs verified line-by-line against `in_topology.cc` (lines 854, 928, 1055). No bugs found.
- **Doc style:** Rust → KaTeX + `[^label]` footnotes.
- **PERTURBATION block is multi-line.** The gromosXX `.imd` PERTURBATION block splits its 8 parameters across two non-comment lines (NTG NRDGL RLAM DLAMT on line 1; ALPHLJ ALPHC NLAM NSCALE on line 2). The parser originally used `data_lines.first()` and silently dropped ALPHLJ/ALPHC/NLAM/NSCALE (they defaulted to 0). Fixed by combining all data lines before tokenising. **Pattern to watch:** any multi-value block whose comment says "one data line" but whose test inputs show two — always combine `data_lines` before splitting values.

## Crate-specific rules
- **ALL parsing/writing lives here.** No file format code in binaries or other crates.
