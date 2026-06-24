# gromos-forces — stage contract

## Job
L1 pure compute. Every force and energy calculation for classical terms. The single source of
physics — both the MD engine and the analysis facade call this layer; nothing re-implements it.

## Inputs (consumes from)
`gromos-core`: Topology, Configuration, BoundaryCondition, pairlist algorithms.

## Outputs (public API)
Bonded + nonbonded force/energy computation; pairlist construction; restraint forces.

## Status
- LJ + CRF nonbonded ✓; 1-4 interactions (cs6/cs12 + scaled CRF) ✓
- All bonded types ✓: quartic/harmonic bonds, cos-harmonic/harmonic angles, dihedrals, impropers, cross-dihedrals
- NTF flag control, RF excluded interactions (forces + energy + self-terms) ✓
- **Perturbed nonbonded** ✓: `perturbed_pairlist_correction` + `perturbed_self_energy_correction` + `perturbed_excluded_correction` + `perturbed_one_four_correction` + `perturbed_atom_pair_correction`; all wired in `Forcefield`; `ch4_water_fep` reference test passes.
- Pairlist: chargegroup-based, atom-based, twin-range (RCUTP/RCUTL with force caching) ✓
  - ⚠️ md.rs still uses `StandardPairlistAlgorithm` (O(N²)); `CellListPairlistAlgorithm` lives in gromos-core but is not wired into the md binary (P1.6)
- Position restraints ✓
- Distance / dihedral / angle restraints: **not yet** (P1.6)
- Distance-field / local elevation: **not yet** (P1.6)

## Key files
```
src/bonded.rs          — all bonded force calculations (~1300 LOC, P4 split candidate)
src/nonbonded.rs       — LJ+CRF, rf_excluded, pairlist loops (~1500 LOC, P4 split candidate)
src/electrostatics.rs  — CRF/PME parameters
```

## Crate-specific rules
- **This is the single source of physics.** The moment `ener.rs` or any analysis program re-implements LJ, you have recreated gromos++.
- **Bonded force vectors:** `v = pos(i) - pos(j)` (gromosXX convention).
- **Do NOT port the Martina solute/solvent misclassification bug** (`extended_grid_pairlist_algorithm.cc:1309`); classify pairs by both atoms' roles.
- **Second-source the perturbed RF self-term** before porting FEP. gromosXX authors flagged it as untrusted: `perturbed_nonbonded_term.cc:596,749` and `:1444`. Derive from the GROMOS book independently; diff against C++; investigate disagreements.
- **gromosXX uses separate perturbed pairlists** (`m_perturbed_pairlist.solute_short` / `.solute_long`) distinct from the regular pairlist. gromos-rs uses a single combined pairlist and applies a correction (subtract state-A, add perturbed) — this is mathematically equivalent when the effective `alpha_lj` is correct. The correction approach is in `perturbed_pairlist_correction`.
- **`perturbed_lj_crf_interaction` LJ is cutoff-independent.** The `crf.cutoff_sq` only affects CRF soft-core terms; LJ soft-core uses only `alpha_lj`, `b_ljs_lambda2`, and `c126 = c12/c6`. Confusion here is a common wrong turn during debugging.
- Grep `interaction/` for `bug|fixme|wrong|hack` before porting any new subsystem.
