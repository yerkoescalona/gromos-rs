# gromos-core — stage contract

## Job
L0 data core. Owns `(Topology, Configuration)` shared by everything. No physics, no I/O.

## Inputs (consumes from)
None — foundation crate.

## Outputs (public API)
`Topology`, `Configuration`, `Vec3`/`Mat3`, `Periodicity` (Vacuum/Rectangular/Triclinic),
`Algorithm` trait + `AlgorithmSequence`, `AtomSelection`, `Stat`,
`PairlistAlgorithm` (enum: Standard/CellList), `StandardPairlistAlgorithm`, `CellListPairlistAlgorithm`,
`gather_chain`, `gather_bond`, `gather_molecules` (PBC gathering).

## Status
- Boundary: Vacuum ✓, Rectangular ✓, Triclinic ✓ (NTB=-1 truncated-octahedron)
- Topology: Dim 10 instancing model complete — `moltypes[0]` = solute, `moltypes[1..]` = solvent types,
  `instances[k].role` = Solute/Solvent; flat arrays derived from instances; `Solvent` struct removed;
  `Solute` is a ZST shell (atoms moved to moltypes). Direct field access: `topo.moltypes[0].bonds`.
- AtomSelection: full gromos-rs grammar ✓ — `a:name`, `1:res(N:CA)`, `not()`, `minus()`, `;`-union,
  `m:` (solute by role), `s:` (solvent by role), syntax-error hints
- PBC gathering: `gather.rs` — `gather_chain`, `gather_bond`, `gather_molecules` ✓
- Statistics: `stat.rs` — `Stat` with `ave()`, `rmsd()`, `ee()` (block averaging) ✓
- Pairlist: `PairlistAlgorithm` enum wired everywhere (9a-0 ✓); `CellListPairlistAlgorithm` bit-identical to Standard (9a-1 ✓, margin=0); `water_1000_spc_gridcell` reference validates against gromosXX Grid_Cell_Pairlist

## Key files
```
src/topology.rs   — Topology, MoleculeType, MoleculeInstance, Role
src/selection.rs  — AtomSelection, full grammar + error hints
src/gather.rs     — PBC molecule gathering
src/stat.rs       — Stat, block-averaging ee()
src/pairlist.rs   — StandardPairlistAlgorithm, CellListPairlistAlgorithm
src/math.rs       — Vec3, Periodicity, truncoct_triclinic_*
```

## Crate-specific rules
- **No physics.** Energy/force belongs in gromos-forces.
- **Dim 10 convention:** `moltypes[0]` is always the solute (GROMOS `moltype[0]` convention).
  Direct field access: `topo.moltypes[0].bonds.push(...)` — no accessor methods.
- **Role-based dispatch:** `is_solvent_atom(i)`, `role_of_atom(i)` — never check `i >= num_solute_atoms()`.
- **gather.rs:** `nearest_image(ri, rj)` returns `ri − rj` (displacement); gathered position = `rj + nearest_image(ri, rj)`.
