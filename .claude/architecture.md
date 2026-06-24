# Architecture guidelines — how to shape code

These are the durable shapes every refactor should honor. They are *guidelines, not tasks* — the
"why" lives in `FUTURE.md` (Dimensions 1, 2, 9, 10, 12, 14). When a change would violate one of
these, that's a signal to stop and reconsider, not a rule to route around quietly.

## The one principle everything else serves

**gromos-rs is one machine, not two.** GROMOS's fatal split is that its engine (gromosXX) and its
analysis suite (gromos++) never shared a core, so the physics is implemented twice. Our entire
architecture exists to make the engine and the analysis facade *the same machine* over a shared core.
Every other guideline is in service of this.

## The layered model

Think in layers; keep them from leaking into each other.

- **L0 — Data core.** One `(Topology, State)`, owned in `gromos-core`. Storage is data-oriented
  (Structure-of-Arrays, not Array-of-Structs) so it's cache-friendly and SIMD/GPU-ready. Repeated
  molecules are *instanced* (stored once, referenced N times), not expanded. "Solvent-ness", QM-zone
  membership, constraint/temperature/energy group are **attributes**, never structural partitions.
- **L1 — Pure compute.** Forces/energies (`PotentialProvider`) and analysis quantities (`Observable`)
  are *pure functions* over `(Topology, State, SpatialIndex)`. **This layer is shared by the engine
  and the analysis facade.** It is the single source of physics.
- **L2 — Evolve.** Mutating steppers (integrators, constraints, thermostats, barostats) that own and
  advance state. MD-specific.
- **L3 — Orchestration.** A run is a composition. An MD step composes L0+L1+L2 and loops over *time*;
  an analysis composes L0+L1 and loops over *frames*. Default orchestration is an explicit, ordered
  sequence (GROMOS-faithful, debuggable). Don't foreclose a lazy compute-graph orchestration later.
- **L4 — Facade & bindings.** The gromos++-style analysis API and py-gromos are thin over L0–L3 and
  contain **zero physics**.

The payoff: `ener` calls the same provider an MD step calls; `rmsd` calls shared geometry; gathering
is one L0 service. The ~79k lines gromos++ spent re-implementing physics simply never get written.

## The "lock now" invariants

Settle these on paper before the refactor; they cost almost nothing now and are very expensive to
retrofit (see FUTURE.md Dim 12.5 / 14.5):

1. **One owned `(Topology, State)` core — never grow a second `System`/`gcore`.** This is the single
   thing gromos++ got wrong.
2. **L1 purity.** Every energy/force/observable is computed in exactly one place, called by both
   engine and analysis. If you're about to write a physics term inside an analysis program, stop —
   it belongs in L1.
3. **A small, stable trait taxonomy — resist proliferation.** Aim for about five contracts:
   `PotentialProvider`, `Observable`, `Stepper`, `SpatialIndex`, `Reader`/`Writer`. Everything else is
   data or a composition of these.
4. **One state-history model.** Leapfrog needs previous state, analysis needs frame windows,
   checkpointing needs snapshots — pick one representation so all three agree.

## The provider pattern (the keystone)

Express force/energy computation as **additive providers over shared state**, classical terms
included. The shape that matters (exact signature will evolve):

- **data-in / scattered-forces-out** — takes positions/charges + a region, returns energy and
  per-atom force/virial contributions onto an *arbitrary atom subset* (embedding forces land on MM
  atoms too).
- **stateful** (`&mut self`) — SCF loops, wavefunction reuse, running stats need it; don't assume
  pure functions.
- **fallible** (`Result`) and **extensible return** — external engines fail; ML providers return
  uncertainty. Design for both from the start.

If classical LJ+CRF is *already* a provider, then QM/MM, ML potentials, and new force terms are just
more providers — additive, touching no existing physics. That is the whole point.

## Supporting shapes

- **Spatial index is a *service*, not an MD pairlist.** Give it a query API ("neighbors of these
  atoms within r, as pairs / as a graph") so one structure serves the MD pairlist, an ML model's
  radial graph, and the QM-zone gather. (FUTURE Dim 9e.)
- **Typed units boundary.** GROMOS works in kJ/mol·nm; QM in Hartree/Bohr; ML in eV/Å. Make unit
  conversions explicit/typed so mismatches are caught, not silent force bugs.
- **Representation vs role are orthogonal.** Instancing answers "stored once or many?"; role/attributes
  answer "treated how?". Keep them independent so one water can be promoted into the QM region (or to
  solute) without a data-model migration.

## Performance philosophy (don't chase C++ micro-opts)

C++ gromos has no secret weapon: templates are Rust monomorphization, OpenMP is `rayon`, intrinsics
are `wide`/SIMD — gromos-rs already matches these. **Don't try to out-tune the inner loop; you'll
tie.** Win by refusing the architecture C++ was forced into: gromosXX forks the same physics per
backend (serial/omp/mpi/cuda) and per variant (perturbed). Write the physics **once** and let a
*scheduler* choose the backend. Prefer SoA, one kernel many schedulers, and a spatial-index service
over four pairlist implementations.

## The trap to watch for

The failure mode is **partial sharing drifting into two cores** — reusing some primitives while
re-implementing others "just for this program." That is exactly how gromos++ began. When adding an
analysis or tool, the default is to consume L0/L1; re-implementing is the exception that needs a
reason.
