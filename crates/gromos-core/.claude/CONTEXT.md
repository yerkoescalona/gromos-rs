# gromos-core — stage contract

## Job
L0 data core. Owns the one `(Topology, Configuration)` that everything else shares. No physics, no I/O.

## Inputs (consumes from)
None — foundation crate. All other crates depend on this one.

## Outputs (public API)
`Topology`, `Configuration`, `Vec3`, `BoundaryCondition` (Vacuum / Rectangular / Triclinic),
`Algorithm` trait + `AlgorithmSequence`, `AtomSelection`,
`StandardPairlistAlgorithm`, `CellListPairlistAlgorithm`.

## Status
- Boundary: Vacuum ✓, Rectangular ✓, Triclinic ✓ (NTB=-1 truncated-octahedron, wired)
- Topology: solvent expansion, chargegroups, exclusions ✓
- AtomSelection: **underbuilt** — `parse_solvent` ignores spec; `parse_molecule` rejects non-first molecules; name/residue search solute-only → P2.1
- Pairlist: `CellListPairlistAlgorithm` (O(N) for rectangular boxes) implemented and validated by set-equality vs `StandardPairlistAlgorithm`; **not yet wired into md.rs** — system-size heuristic is the remaining step (P1.6 deferred)

## Key files
```
src/algorithm.rs      — Algorithm trait, AlgorithmSequence
src/math.rs           — Vec3, BoundaryCondition, truncoct_triclinic_*
src/topology.rs       — Topology struct
src/selection.rs      — AtomSelection (P2.1: extend to AtomSpecifier facade)
src/pairlist.rs       — StandardPairlistAlgorithm, CellListPairlistAlgorithm
```

## Crate-specific rules
- **No physics here.** Energy/force computation belongs in gromos-forces (L1), never in this crate.
- **Solvent-ness is an attribute, not a partition.** Do AtomSelection work in a way that doesn't calcify around the old solute/solvent index split before the Dim 10 instancing refactor lands.
- **BoundaryCondition:** vacuum if box_dims = (0,0,0), rectangular otherwise; GENBOX `box_type` field is parsed but only the dimensions matter — except NTB=-1 which routes to Triclinic via `truncoct_triclinic_box`.
