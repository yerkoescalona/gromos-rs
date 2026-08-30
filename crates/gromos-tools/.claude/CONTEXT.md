# gromos-tools — stage contract

## Job
L3 system builder. System construction binaries: the gromos++-style tools for preparing simulation inputs.

## Inputs (consumes from)
`gromos-core` + `gromos-io` (+ `gromos-forces` for single-point energy — not yet).

## Outputs
Standalone binaries: make_top, com_top, check_top, pdb2g96, sim_box, ion, mk_script,
make_pt_top, prep_posres, build_box.

## Status (audit 2026-08-30, every tool run on real inputs against gromos++ where a counterpart exists)

**Verified against gromos++ output** (reference tests under `tests/`, fixtures + commands in
`tests/data/<tool>/README.md`; gromos++ = `.local/gromosPlsPls`):

| tool | what is compared | test |
|------|------------------|------|
| `make_top` | topology reads back identical: atoms, exclusions, 1-4, every bonded term and parameter, solvent — peptide `NH3+ ALA GLY COO-` and a two-MTB methanol; gromos++'s own reader accepts ours | `topology_references.rs` |
| `com_top` | same, aladip + CH4 | `topology_references.rs` |
| `red_top` | same, atoms 1–6 of aladip | `topology_references.rs` |
| `make_pt_top` | `PERTATOMPARAM` reads back identical | `topology_references.rs` |
| `sim_box` | box, solvent count, every position (352 SPC around methanol) | `sim_box_reference.rs` |
| `copy_box` | every position and label, `@dir x` | `box_tool_references.rs` |
| `inbox` | every position | `box_tool_references.rs` |
| `ion` | the water gromos++ replaces (`@potential`, its energy rule: all atoms of other molecules, shifted Coulomb) | `box_tool_references.rs` |
| `build_box` | positions and box (methanol, `@nsm 3 @dens 790`) — checked by hand, identical | — |
| `check_top` | runs; gromos++ additionally cross-checks dihedral types against the building blocks (it flags one in the tutorial's methanol topology that we accept) | — |
| `atominfo` | gromos++'s columns, 1-based IAC | — |

**Run, but not gromos++'s program** — same name, simplified or different semantics; results were
not cross-checked. Treat as utilities, not ports: `gch` (checks constraint distances; gromos++
regenerates hydrogen coordinates), `check_box` (box statistics; gromos++ reports close atom pairs),
`bin_box` (bins a trajectory; gromos++ builds a binary mixture), `unify_box` (rescales to a box;
gromos++ converts between box types), `pert_top` (compares two topologies; gromos++ modifies
selected atoms), `con_top` (concatenates; gromos++ converts force-field parameters), `link_top`,
`addvirt_top`, `prep_eds`, `prep_noe`, `amber2gromos`, `ran_box`, `ran_solvation` (random
placement with other argument models), `mk_script` (needs `@template`).

**No gromos++ counterpart:** `traj2pdb`, `top2build`, `prep_posres`, `prep_xray` (gromos++'s
needs clipper), `mpi_scaling`.

Conventions: `make_top`, `com_top`, `check_top`, `pdb2g96`, `ion`, `build_box`, `prep_posres` take
`--flag` (clap) *and* `@flag` (via `gromos_args`); the others take `@flag` only.

Adding a reference test: put gromos++'s output under `tests/data/<tool>/` with a README that records
the exact command, and compare in a `tests/*_references.rs` (`read_topology_file` for topologies —
formatting-independent — or the `parse` helper of `box_tool_references.rs` for coordinates).

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
