# gromos-forces — stage contract

## Job
L1 pure compute. Every force and energy calculation for classical terms.
Single source of physics — both MD engine and analysis facade call this layer.

## Inputs (consumes from)
`gromos-core`: Topology, Configuration, BoundaryCondition, pairlist algorithms.

## Outputs (public API)
Bonded + nonbonded force/energy; pairlist; restraint forces; `single_point_energy()`.

## Status
- LJ + CRF nonbonded ✓; 1-4 interactions ✓; twin-range (RCUTP/RCUTL) ✓
- All bonded types ✓: quartic/harmonic bonds, cos-harmonic/harmonic angles, dihedrals, impropers, cross-dihedrals
- Perturbed nonbonded ✓: soft-core LJ+CRF, pairlist correction, self/excluded/1-4/PERTATOMPAIR corrections
- Position restraints ✓; distance restraints ✓
- `energy.rs` ✓ — `single_point_energy(topo, positions, box_dims, params)` — used by `ener` analysis binary
- ⚠️ Perturbed RF self-term needs second-sourcing from GROMOS book (`perturbed_nonbonded_term.cc:596,749,1444`)
- Dihedral restraints: **not yet** — next after CellList wiring (P1.6)

## Key files
```
src/bonded/          — bonded force calculations (split from bonded.rs)
src/nonbonded/       — LJ+CRF, rf_excluded, pairlist loops (split from nonbonded.rs)
src/energy.rs        — single_point_energy() entry point for analysis
src/restraints.rs    — position + distance restraints
```

## Crate-specific rules
- **Single source of physics.** The moment `ener.rs` or any analysis program re-implements LJ, you've recreated gromos++.
- **Bonded force vectors:** `v = pos(i) - pos(j)` (GROMOS convention).
- **Do NOT port the Martina solute/solvent misclassification bug** (`extended_grid_pairlist_algorithm.cc:1309`).
- **Perturbed pairlist approach:** gromos-rs uses one combined pairlist + correction (subtract state-A, add perturbed) — mathematically equivalent to separate perturbed pairlists when `alpha_lj` is correct.
- **`perturbed_lj_crf_interaction` LJ is cutoff-independent.** `crf.cutoff_sq` only affects CRF soft-core.
- Grep `interaction/` for `bug|fixme|wrong|hack` before porting any new subsystem.
