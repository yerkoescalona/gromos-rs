# gromos-tools — stage contract

## Job
L3 system builder. System construction binaries: the gromos++-style tools for preparing simulation inputs.

## Inputs (consumes from)
`gromos-core` + `gromos-io` (+ `gromos-forces` for single-point energy — not yet).

## Outputs
Standalone binaries: make_top, com_top, check_top, pdb2g96, sim_box, ion, mk_script,
make_pt_top, prep_posres, build_box.

## Status
- All system-building tools built ✓
- make_top tested: GB3 (56 res, 457 atoms) + Na+ with 54A7 ✓
- **No reference tests yet** — target of P2 roadmap (hand-craft minimal systems, mirror `gromosXX_references/` harness)

## Key files
```
src/bin/topology/   — make_top, com_top, check_top, etc.
src/bin/box/        — sim_box, build_box, ion, etc.
```

## Crate-specific rules
- **Refactor make_top / com_top into library functions** (not binary-only `main`s) so pyo3-gromos can call them — P3 target.
- **File parsing always delegates to gromos-io.** Never duplicate format code in tool binaries.
- **AtomSpecifier:** build the gromos++ grammar as a facade over gromos-core primitives (P2.1), not as a standalone reimplementation.
