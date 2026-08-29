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
- `sim_box` follows gromos++'s algorithm exactly (box from the longest solute atom–atom distance,
  template centred and replicated `int(box/template)+1` times, hydrogen-blind clash filter, GENBOX
  output); `tests/sim_box_reference.rs` requires gromos++'s box, count and positions (2026-08-30).
- Other tools: no reference tests yet — target of P2 roadmap (hand-craft minimal systems, mirror
  `gromosXX_references/` harness). Adding one: put gromos++'s output under `tests/data/<tool>/` with a
  README that records the exact command, and compare in `tests/<tool>_reference.rs`.

## Key files
```
src/bin/topology/   — make_top, com_top, check_top, etc.
src/bin/box/        — sim_box, build_box, ion, etc.
tests/              — reference tests against gromos++ output (`tests/data/<tool>/README.md`)
```

## Crate-specific rules
- **Refactor make_top / com_top into library functions** (not binary-only `main`s) so pyo3-gromos can call them — P3 target.
- **File parsing always delegates to gromos-io.** Never duplicate format code in tool binaries.
- **AtomSpecifier:** build the gromos++ grammar as a facade over gromos-core primitives (P2.1), not as a standalone reimplementation.
