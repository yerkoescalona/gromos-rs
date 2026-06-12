# gromos-rs — FUTURE.md (where Rust can overtake gromosXX + gromos++)

Companion to `PLAN.md`. PLAN.md is about **parity** (reproduce the C++ engine, bit-for-bit,
feature by feature). FUTURE.md is about **divergence** — the dimensions where a Rust rewrite can
do things the two C++ codebases (gromosXX/md++ ≈ 190k LOC, gromos++ ≈ 79k LOC) structurally
cannot, or can only do with enormous duplication and risk.

Read PLAN.md first. Nothing here matters until the reference suite is green and stays green; every
idea below must be gated behind "the bit-for-bit tests still pass."

---

## 0. The premise, answered directly

> *"I am suspicious that C++ gromos exploits the language for optimizations I cannot reach in Rust."*

It does not. I read the hot paths on both sides. Here is the actual scorecard:

| C++ gromos optimization | What it is | gromos-rs equivalent | Status |
|---|---|---|---|
| Template specialization of the nonbonded term (`Nonbonded_Term`, boundary/interaction-function as template params) | Compile-time dispatch, no vtable in the inner loop | Generics `fn …<BC: BoundaryCondition>` → **monomorphization** (identical codegen) | **Already have it** (`nonbonded.rs`) |
| Template `bool` flags (virial on/off, perturbation on/off) | One source, many specialized loops | `const VIRIAL: bool` const generics | **Already have it** (`process_pairs<BC, const VIRIAL: bool>`) |
| `#pragma omp parallel` over the pairlist (54 pragmas) | Thread-parallel outer loop | `rayon` `par_chunks` | **Already have it** (`*_parallel` innerloops) |
| Hand-written SSE/AVX in places, `-O3 -march=native` | SIMD | `wide` (`f64x4`) + autovec under `target-cpu=native` | **Partially** (`lj_crf_interaction_simd_x4`) |
| Separate MPI master/slave nonbonded sets (67 files touch MPI) | Distributed pairlist | not yet — but see §1 | open |
| Hand-written CUDA kernels (`cukernel/`) | GPU offload | not yet — but see §4 | open |

**The uncomfortable truth:** C++ doesn't have a secret weapon. Templates are monomorphization,
OpenMP is rayon, intrinsics are `wide`/`std::simd`. Where C++ *wins today* is only that those
paths are already written and tuned. Where C++ **structurally loses** is everything that the
language makes expensive to evolve: it pays for its performance with **three-to-five parallel
copies of the same physics** (serial / OMP / MPI / CUDA / perturbed) that must be kept in sync by
hand. That duplication is the seam Rust pries open.

**So stop trying to beat C++ at its own micro-optimizations.** You'll match them and that's it.
Win on the axes below, where the rewrite is not a rewrite of the *code* but of the *architecture
the code was forced into by C++'s constraints.*

---

## The reframe

gromosXX's performance architecture is **"one algorithm, hand-forked N times per backend."**
Look at `interaction/nonbonded/interaction/`: `nonbonded_set`, `perturbed_nonbonded_set`,
`cuda_nonbonded_set`, `mpi_nonbonded_master`, `mpi_nonbonded_slave`, `omp_nonbonded_interaction`,
`solvent_innerloop`, `spc_table` — these are largely the *same loop* re-expressed for each
execution context. Every physics change is an N-way edit. That is the tax C++ levies, and it caps
how fast the engine can actually improve.

Rust's lever is not "a faster inner loop." It is **"one expression of the physics that targets
many backends without forking."** Get that right and you don't beat gromos by 10% on a benchmark —
you beat it on the derivative: features land 5× faster, GPU comes nearly free, and the analysis
suite (gromos++) stops being a second 79k-LOC codebase and becomes a thin façade. That is the
overperformance that compounds.

The dimensions, highest-leverage first:

---

## Dimension 1 — Data-oriented core (SoA), the one micro-architectural win that *is* worth taking

Both codebases store positions Array-of-Structs: gromosXX `math::VArray` = `vector<GenericVec<double>>`
(`d_v[3]` inline), gromos-rs `Vec<Vec3>`. AoS wastes cache and blocks clean vectorization: loading
`x` for 4 atoms touches 4 cache lines and interleaves the `y,z` you didn't want.

**The lever:** Rust makes a Structure-of-Arrays core *ergonomic and safe* in a way C++ won't,
because in Rust you can keep the `Vec3` ergonomics at the API boundary (an index type + accessor)
while the storage is `struct Soa { x: Vec<f64>, y: Vec<f64>, z: Vec<f64> }`. C++ gromos cannot
adopt SoA without rewriting `VArray` and every one of the thousands of `pos(i)` call sites by hand,
across all N forked loops — so it never will.

**Plan**
- Introduce `gromos-core::soa::PositionStore` (and vel/force) behind the existing `State` API.
  Keep `pos[i] -> Vec3` working via accessor; migrate the nonbonded/pairlist hot path to consume
  raw `&[f64]` lanes.
- The pairlist already chunks by charge group — feed SoA lanes into `wide::f64x4` directly instead
  of gathering from `Vec<Vec3>`.
- Gate behind a criterion benchmark: target ≥1.5× on `water_216` single-step nonbonded vs the
  current AoS path before committing.

**Why Rust wins:** the migration is *local* (storage + accessors) because the borrow checker and
the single non-forked loop mean there's one place to change, not five. **Payoff:** real,
measurable single-thread speedup *and* it unlocks Dimensions 3/4 (SIMD/GPU want SoA anyway).
**Risk:** medium; touches the bit-for-bit hot path. Mitigation: SoA must be exercised by the full
reference suite at every step.

---

## Dimension 2 — One physics kernel, many backends (kill the N-way fork)

This is the big one — the structural defeat of C++ gromos.

**The lever:** express each interaction (LJ+CRF, bonded, SHAKE) **once** as a pure function over
SoA lanes + parameters, with backend (serial / rayon / SIMD / GPU) chosen by the *caller*, not
baked into a forked file. gromos-rs already half-does this with `const VIRIAL: bool` and `BC`
generics; the goal is to finish it so that `perturbed`, `mpi`, and `cuda` are *not* separate
implementations of the loop but separate *drivers* over the same kernel.

**Plan**
- Define a `Kernel` boundary: `fn lj_crf(pair_block, params, &mut accumulators)` with no knowledge
  of how pairs were scheduled.
- Backends become schedulers: `serial`, `rayon`, (later) `gpu` all call the identical kernel body.
- Perturbation/FEP becomes a *parameter* (`λ`-interpolated params) fed to the same kernel, not
  `perturbed_nonbonded_innerloop.cc`. (PLAN.md P1.5 FEP lands here for free.)
- Delete the conceptual need for `omp_*` / `mpi_*` / `cuda_*` triplication before it's ever written.
- **Shape the kernel as an additive `PotentialProvider` (data-in / scattered-forces-out, stateful,
  fallible) from the start** — even classical LJ+CRF implements it. This is the keystone that makes
  Dimension 12 (QM/MM + ML) *additive* rather than a retrofit: a QM/ML engine is then just another
  provider. See Dim 12.5 for the exact invariants — they are decisions for *this* refactor, not later.

**Why Rust wins:** traits + monomorphization give zero-cost backend selection; C++ would need the
same but its existing code has already committed to the forked layout and won't be untangled.
**Payoff:** the compounding one — every future physics feature is written once. **Risk:** high
design cost up front; this is an architecture decision to make *early*, before EDS/GaMD/REMD
(PLAN.md P1.5) calcify into their own forks the way they did in C++.

---

## Dimension 3 — SIMD as the default, not the exception

gromos-rs has `lj_crf_interaction_simd_x4` but most innerloops are scalar. gromosXX's vectorization
is `-O3` autovec plus a few hand spots — fragile, compiler-version-dependent.

**The lever:** with the SoA core (Dim 1) and the single kernel (Dim 2), `wide::f64x4` (→ `std::simd`
when stable) becomes the *natural* expression of the inner loop: process 4 pairs/lanes at once,
masked tails. Portable across x86 AVX2 and ARM NEON from one source — gromosXX's hand SIMD is not.

**Plan**
- Once SoA lands, make the SIMD path the primary `lj_crf` kernel; scalar becomes the remainder loop.
- Vectorize the reciprocal/`rsqrt` heavy part (the LJ+CRF core is `1/r`, `1/r^6` — SIMD-friendly).
- Benchmark on both AVX2 and (if available) NEON to prove the portability claim.

**Payoff:** 2–4× on the nonbonded core, portably. **Risk:** low–medium; numerics must stay within
the 1e-6 force / 1e-8 energy tolerances — SIMD reassociation can perturb the last bits, so verify
against references and possibly pin reduction order.

---

## Dimension 4 — GPU from one codebase (`cubecl`/`wgpu`), not a second language

gromosXX's GPU path is hand-written CUDA in `cukernel/` — NVIDIA-only, a separate language, a
separate copy of the physics, and a maintenance sink most builds don't even compile.

**The lever:** Rust GPU compute (`cubecl`, `wgpu`, or `rust-gpu`) lets the **same kernel body**
(Dim 2) target GPU. `cubecl` in particular compiles Rust to CUDA/ROCm/WGPU — you write the LJ+CRF
loop once and it runs on CPU SIMD *and* GPU. That is precisely the thing C++ gromos pays a whole
`.cu` tree to fake, and only for NVIDIA.

**Plan (later — after Dims 1–3 prove the single-kernel model)**
- Prototype the nonbonded kernel under `cubecl`; validate against `water_216` reference.
- GPU pairlist (cell list) as a second kernel.
- Keep it feature-gated (`--features gpu`); CPU path remains the correctness ground truth.

**Payoff:** vendor-portable GPU MD from one source — a capability gromosXX cannot match without a
second CUDA rewrite per vendor. **Risk:** high; bit-for-bit parity with a CPU reference is hard on
GPU (FP reductions). Treat GPU as "fast approximate, validated against CPU," not bit-identical.
**Sequence late** — this is the reward for getting Dim 2 right, not a starting point.

---

## Dimension 5 — Collapse gromosXX + gromos++ into zero-copy, no second codebase

gromos++ is 79k LOC that **re-implements** geometry, PBC, topology reading, energy evaluation —
because in the C++ world the analysis tools couldn't link the MD engine's internals cleanly, so
they were copied. PLAN.md already names this ("seamless merge, no duplication"); FUTURE.md states
the *performance* consequence: the duplication isn't just a maintenance cost, it's a correctness
and speed cost (two LJ implementations that can silently disagree; analysis that re-reads and
re-gathers trajectories naively).

**The lever:** Rust's module/crate system + the single physics kernel (Dim 2) make the gromos++
façade a *thin* API over the engine. `ener` calls the *same* `lj_crf` kernel the MD step uses —
not a hardcoded copy (today `gromos-analysis/ener.rs` hardcodes σ/ε; PLAN.md P2.4). Single-point
energy on a trajectory frame is then automatically SoA + SIMD + (eventually) GPU.

**Plan** (this is PLAN.md P2, viewed through the performance lens)
- `gromos-analysis` depends on `gromos-forces`/`gromos-integrators`; delete every reimplemented
  primitive (energy, geometry, fit, stats).
- Trajectory analysis runs the production kernel → analysis inherits every Dim 1–4 speedup for free.

**Payoff:** one physics implementation total (vs three in C++ land), and analysis tools that are as
fast as the engine. **Risk:** low, mostly mechanical; it's already the stated architecture.

---

## Dimension 6 — Streaming / out-of-core trajectory analysis

gromos++ tools largely load frames and process serially. Modern trajectories are larger than RAM.

**The lever:** Rust's iterator + `rayon` model makes a **streaming, parallel, back-pressured**
trajectory pipeline natural: `frames().par_bridge().map(analyze).reduce(...)` with memory-mapped
or chunked I/O, never materializing the whole trajectory. Combined with the zero-copy façade
(Dim 5), each frame's energy/geometry uses the SIMD kernel.

**Plan**
- A `TrajectoryStream` iterator over `gromos-io` readers (mmap where possible).
- Analysis programs (rmsd, rdf, ene_ana) become `fold`/`reduce` over the stream.
- Block-averaging error estimates (PLAN.md P2.2) fit naturally as online reducers.

**Payoff:** analyze trajectories that don't fit in RAM, in parallel, at engine speed. C++ gromos++
would need a substantial I/O and threading rework to match. **Risk:** low.

---

## Dimension 7 — Determinism & reproducibility as a guaranteed feature

C++ MD reproducibility is famously fragile: OMP reduction order, `-ffast-math`, compiler version,
and `-march=native` all perturb the last bits. gromosXX manages it by discipline, not by guarantee.

**The lever:** Rust gives you *control* over this. You can offer a `--deterministic` mode with
fixed reduction order (no `fast-math` reassociation by default in Rust — FP is strict unless you
opt in), pinned thread-count reduction trees, and a portable RNG (you already did this:
`GslMt19937` matches GSL bit-for-bit, PLAN.md 1.1). Make "same input → same bytes, any machine,
any thread count" a *tested, advertised* property.

**Plan**
- A deterministic parallel reduction (fixed-order tree) for energies/virial, selectable at runtime.
- CI matrix: run the reference suite at `RAYON_NUM_THREADS=1,2,4,8` and assert identical output.
- Document the guarantee. This is a genuine scientific selling point gromosXX cannot honestly make.

**Payoff:** a correctness/reproducibility guarantee no mainstream C++ MD code offers. **Risk:** low;
mostly a matter of choosing reduction strategy and adding CI. Note tension with Dim 3 SIMD — the
deterministic mode may forgo some SIMD reassociation; offer both modes.

---

## Dimension 8 — Distributed without the MPI fork (optional, far out)

gromosXX's MPI domain decomposition is 67 files of master/slave bookkeeping — a third copy of the
nonbonded path. If distributed scaling is ever needed, the single-kernel model (Dim 2) means the
domain decomposition is a *scheduler* over the same kernel, not a reimplementation. Modern Rust
options (`mpi` crate, or a custom async transport) keep the physics untouched.

**Status:** deliberately last. Single-node SoA+SIMD+GPU covers the vast majority of GROMOS-scale
systems; MPI is only worth it once you're GPU-saturated. Listed for completeness so it's designed
*for* (one kernel) and not *against* (a fork) from day one.

---

## Dimension 9 — Pairlist & charge groups: from O(N²) to cache-coherent O(N)

This is the dimension I under-weighted in the first draft, and on reflection it's the **single
biggest scaling gap in gromos-rs today** — more urgent than SIMD or GPU, because it's an
*asymptotic* gap, not a constant-factor one.

**The honest current state.** gromos-rs has exactly one pairlist algorithm:
`StandardPairlistAlgorithm` (`gromos-core/src/pairlist.rs`), which is **O(N_cg²)** — it loops over
every pair of charge-group centers-of-geometry and tests the cutoff. That's perfectly fine for the
reference suite (216 waters = ~648 CGs), and it's bit-for-bit correct. But it does not scale: at
50k–100k atoms the pairlist construction itself dominates the step. gromosXX already moved past
this — it ships **four** pairlist implementations: `standard_pairlist_algorithm`,
`standard_pairlist_algorithm_atomic`, `grid_cell_pairlist`, and `extended_grid_pairlist_algorithm`
(the grid/cell ones are O(N) via spatial decomposition). **So here gromos-rs is behind on the
asymptotics, not ahead.** Catching up is non-negotiable; the question is how to catch up in a way
that then *overtakes*.

### 9a — The catch-up: a charge-group-aware cell list (O(N))

A cell (linked-cell) list bins particles into a grid of cutoff-sized cells; each particle only
checks its own + 26 neighbor cells → O(N). The GROMOS subtlety that the naïve textbook version gets
wrong: **GROMOS's cutoff unit is the charge group, not the atom.** The pairlist tests distance
between CG centers-of-geometry (`chargegroup.center_of_geometry`, already in `topology.rs`), and the
reaction-field self-term correctness depends on whole neutral charge groups being inside or outside
the cutoff together. So the cell list must **bin charge groups by their COG**, not atoms — then
expand accepted CG pairs into atom pairs. Get this wrong and RF energies drift.

**Plan**
- `CellListPairlistAlgorithm` behind the existing `update<BC>()` interface — drop-in alongside
  `StandardPairlistAlgorithm`, chosen by system size or `.imd` flag, so the O(N²) path stays as the
  always-correct reference oracle for tests.
- Bin CG COGs into a grid sized to `max(cutoff_long + skin)`; build half-neighbor lists (i<j) to
  avoid double-counting, matching the standard algorithm's pair set exactly.
- Validate by asserting the cell-list pair set is *identical* to the standard algorithm's on every
  reference system (a set-equality test, not just an energy test) before trusting it.

This is parity work. The next two are where Rust overtakes.

### 9b — The overtake: spatial reordering of the atom arrays (the cache lever C++ won't pull)

A cell list naturally produces a *spatial bucketing* of atoms. The big win is to then **renumber the
atoms in memory by cell** (Morton / Z-order or simple cell-major order) so that atoms close in space
are close in the SoA arrays (Dim 1). Then the inner loop's neighbor accesses are cache-coherent
instead of scattered — often a larger speedup than SIMD, because the nonbonded loop is
memory-latency-bound, not compute-bound.

gromosXX **cannot easily do this**: renumbering atoms would invalidate the thousands of places that
index `pos(i)` by absolute atom number across its forked loops, its I/O, its analysis. So it leaves
atoms in input order and eats the cache misses. gromos-rs, with the single non-forked loop (Dim 2)
and an indexed SoA core (Dim 1), can keep an internal permutation transparently: the *engine* works
in spatially-sorted order, and the permutation is applied only at I/O boundaries (read/write
trajectories in original order). **This is a structural overperform** — the same property (one loop,
owned data layout) that powers Dimensions 1–2 pays off again here.

**Plan**
- Maintain `perm: Vec<u32>` (sorted → original) and its inverse; sort on each pairlist rebuild.
- Engine loops consume sorted SoA; `gromos-io` writers un-permute at the boundary.
- Benchmark `water_216` and a synthetic 50k-atom box: expect the gap vs the O(N²) path to widen
  dramatically with size, and the sorted-vs-unsorted gap to grow with system size too.

### 9c — Parallel + verlet-skin rebuild criterion

- **Parallel cell build** with rayon: cells are independent buckets → lock-free parallel binning and
  per-cell neighbor enumeration. One implementation, the scheduler (Dim 2) decides serial vs
  parallel — not gromosXX's separate `omp`/`mpi`/`cuda` pairlist forks.
- **Displacement-triggered rebuild.** Today gromos-rs rebuilds on a fixed step interval (NSNB,
  matching gromosXX). A more robust and often faster criterion: track max atom displacement since the
  last build and rebuild only when `2·max_disp > skin`. This is *safer* (never misses a pair that
  crossed the skin) and *faster* (skips rebuilds when nothing moved). Offer it as an opt-in mode,
  keep NSNB as the bit-for-bit default for the reference suite.

### 9d — Charge groups as a first-class, queryable structure

Charge groups sit at the intersection of three concerns: the pairlist cutoff unit (9a), the RF
self-term correctness, and **topology composition** (next section — combining two topologies must
renumber CG boundaries and exclusion offsets consistently). They're currently a `Vec` of atom-index
lists on `Topology`. Promoting them to a queryable primitive (CG-of-atom, atoms-of-CG, COG cache)
serves the pairlist, the RF term, the AtomSpecifier work (PLAN.md P2.1), and the builder below — one
structure, four consumers. This is the no-duplication principle applied to the CG concept itself.

### 9e — Design the cell list as a spatial-index *service* (a Dim 12 prerequisite)

Don't build a single-purpose MD pairlist. The MD pairlist, an ML potential's radial graph (different
cutoff, edge lists, sometimes angles/triplets), and the QM-zone gather (MM charges within cutoff of
the QM region, rebuilt each step) are **three different spatial queries** — and QM/MM/ML (Dim 12)
needs all three from *one* spatial structure. Give the cell list a **query-based public API**
("neighbors of these atoms within radius r, as pairs / as a graph") rather than baking in the MD
pairlist's assumptions. This costs nothing extra now and is the difference between Dim 12 reusing the
spatial index and Dim 12 forking its own (the gromosXX disease). Lock this in when you write 9a.

**Payoff:** O(N) scaling (catch-up) *plus* cache-coherent inner loops (overtake) *plus* a robust
rebuild criterion — and one cell-list abstraction instead of gromosXX's four pairlist files.
**Risk:** medium. 9a must be set-equality-validated against the O(N²) oracle; 9b's permutation must
be invisible at I/O (un-permute on write, or trajectories come out scrambled). Sequence 9a early
(it's the real scaling blocker), 9b after the SoA core lands.

---

## Compositional topology construction in py-gromos

The second thing you asked about. This is less "beat C++ on FLOPs" and more "beat the C++ *workflow*
so decisively that it changes who can use GROMOS" — and it's the natural home for the Python API
(PLAN.md P3).

### The C++ workflow, and why it's the thing to beat

Building a system in gromos today is a **batch pipeline of stringly-typed CLI programs over files**:
`make_top` (MTB building blocks + IFP force field → `.top`) → `com_top` (combine topologies) →
`check_top` (validate, prints to stdout) → `sim_box`/`build_box` (solvate) → `ion` (add ions) →
hand-edited `.imd`. Each stage serializes to a file, errors surface late as text, nothing composes
in memory, and you cannot build-then-simulate without a file round-trip. gromos-rs already
re-implements all these as tools (`gromos-tools`) — but porting the *pipeline* would inherit its
ergonomics. The opportunity is to expose the **objects** the pipeline manipulates and let Python
compose them.

### The core insight: a topology is an algebra

Molecules compose. MTB building blocks are reusable parameterized units; a force field is the
resolution context that gives them meaning; a system is built by combining units and patching the
seams (end groups, links). That's a small algebra, and the API should *be* that algebra rather than
hiding it behind file I/O. Three primitives:

- **`ForceField`** — the resolution context (one IFP + its MTB building blocks). Every block and
  atom type resolves against exactly one force field, so a type mismatch is caught at *composition
  time*, not when md later fails to parse.
- **`BuildingBlock`** — an immutable template (an MTB entry): atoms, partial charges, bonds,
  charge-group boundaries, link atoms for head/tail.
- **`Topology`** — a value type that composes: combining two topologies (`+`) is `com_top` as an
  operator; the Rust core renumbers atoms, **charge groups (9d), exclusions and 1-4 pairs** so the
  result is always consistent. *This is why the CG/exclusion bookkeeping must be a first-class Rust
  primitive — composition is where it gets stressed.*

### The design rule that makes it safe

**Python expresses intent (verbs); the Rust core maintains every invariant.** Python never
hand-edits an exclusion list or a charge-group boundary — it calls builder verbs (`add`, `link`,
`patch`, `+`, `solvate`) and the Rust layer (the same `build_exclusions` / `build_14_pairs` /
`build_lj_matrix` the engine already uses, `topology.rs`) regenerates the derived data. So an
in-progress topology is *always* valid by construction, and there is exactly one implementation of
those invariants — shared with the MD engine. (Same no-duplication theme as the physics kernel.)

### What the API should look like

A fluent builder for the common polymer case, operator composition for mixing, methods for
solvation/ions, eager structured validation, and — the real payoff — a direct hand-off to the
simulation with no file round-trip:

```python
import gromos as gx

ff = gx.ForceField.load("54a7")                  # IFP + MTB blocks; the resolution context

# make_top as a first-class call: sequence -> linked, patched, exclusions/1-4 built
pep = gx.Topology.from_sequence(
    ff, ["ALA", "GLY", "SER"], n_term="NH3+", c_term="COO-"
)

# ...or compose residue-by-residue when you need control
pep = (gx.Builder(ff)
        .start("NH3+").add("ALA").add("GLY").add("SER").end("COO-")
        .build())

# com_top as `+`; Rust renumbers atoms, charge groups, exclusions, 1-4
system = pep + gx.Topology.from_block(ff, "NA+") * 3      # add 3 sodium ions

# sim_box / build_box behind a method
system = system.solvate("spc", box="cubic", min_dist=0.9)

# check_top as a structured result, not stdout text
report = system.check()
if not report.ok:
    raise ValueError(report.errors)          # typed errors, caught early

# the payoff: build -> minimize -> simulate in one process, zero files
sim = gx.Simulation(system, gx.imd.nvt(temp=300, steps=10_000))
traj = sim.minimize().run()
```

### The Polars-lazy option (worth prototyping)

PLAN.md P3 already says to study Polars' pyo3 patterns. The deepest version of this idea is a
**`LazyTopology`**: builder verbs record an operation plan (an AST) instead of materializing, and
`.collect()`/`.build()` validates and realizes it once. Benefits mirror Polars: the whole build is
validated as a unit (catch "ALA links to a block with no head atom" before allocating anything),
the plan is introspectable and serializable (a reproducible, diffable system recipe — far better
than a 50k-line `.top`), and expensive steps (exclusion generation over a big polymer) run once at
collect. A system becomes a short, version-controllable *recipe*, not a giant generated artifact.
This is a genuine overtake: gromos++ has no equivalent of "the system as a composable, lazy,
introspectable program."

### Why this overperforms gromos++ (not just "is nicer")

1. **One address space, no round-trips** — build → solvate → minimize → simulate → analyze without
   serializing `.top`/`.cnf`/`.trc` between every step. The C++ pipeline *must* round-trip.
2. **Eager, typed validation** — errors at composition time with structured data, not stdout text
   from `check_top` after the fact.
3. **Composability** — `+`, `*`, `.solvate()` are real operations on real objects; the C++ tools
   can't compose, they regenerate.
4. **Same invariants as the engine** — the topology you build is the exact Rust struct md consumes,
   so "it built" and "it runs" can't diverge. In C++ land make_top and md are different programs
   that can disagree.
5. **Reproducible recipes** — a 20-line Python script (or a lazy plan) replaces a generated
   megabyte `.top`; reviewable, diffable, parameterizable.

**Plan (extends PLAN.md P3)**
- Expose `ForceField`, `BuildingBlock`, `Topology` (with `+`/`*`/`solvate`/`check`/`from_sequence`)
  in `pyo3-gromos`, backed by the existing `gromos-tools` make_top/com_top logic refactored into
  reusable `gromos-core`/`gromos-tools` library functions (not binary-only `main`s).
- Promote charge groups + exclusions to first-class primitives (9d) so composition is consistent.
- `Builder` fluent API first (concrete, testable); `LazyTopology` as a second iteration once the
  eager path is solid.
- Round-trip parity test: a topology built compositionally must serialize to a `.top` byte-identical
  to the one `make_top` produces for the same system — the reference-test discipline, applied to the
  builder.

**Risk:** low on the physics (it reuses validated engine code); medium on API design — get the
core algebra (`ForceField` context, `Topology` as a composable value, Rust owns invariants) right
*first*, because the fluent/lazy sugar is cheap to change but the object model is not.

---

## Dimension 10 — Dissolve the solute/solvent split: separate *representation* from *role*

You flagged this one yourself, and it's sharp: the solute/solvent division is *simultaneously* one
of GROMOS's best ideas and one of its most annoying constraints. The reason it's both is that it
**fuses two orthogonal concerns into one rigid boundary** — and once you see the two concerns
separately, the fix is obvious and it's exactly the kind of thing the unification should deliver.

### How gromosXX and gromos++ evaluate it (they differ — and gromos++ is the worse one)

- **gromosXX** keeps a single contiguous atom array: solute occupies `[0, num_solute_atoms)`,
  solvent is appended after, and `is_solvent(i)` is literally `i >= num_solute_atoms()`
  (`topology.h`). The win: the solvent *topology* (charges, LJ types, bonds/constraints) is stored
  **once** per solvent type and logically instanced `num_solvent_molecules` times — you don't repeat
  the water topology 216×. It even exploits the regularity for speed: a dedicated `solvent_innerloop`
  + hard-coded `spc_table.h` fast path. So gromosXX's split is an *index threshold over one array* —
  rigid, but at least uniform in memory.
- **gromos++** is harder and is the source of most of the pain you feel. `System` stores solute and
  solvent in **separate containers of different classes**: `mol(i)` returns a `Molecule`, `sol(i)`
  returns a `Solvent`, with separate `numMolecules()` / `numSolvents()` (`gcore/System.h`). Every
  analysis program therefore has *two code paths*, and `AtomSpecifier` bifurcates into `m:` (solute
  molecule) vs `s:` (solvent) namespaces. There is no clean way to say "this one water is part of my
  analysis selection like a solute" without special-casing.
- **gromos-rs today** is a hybrid (and currently the worst of both): `Topology::solvate()` keeps a
  `Solvent` *template* (good — stored once) **but also physically expands** the flat
  `mass`/`charge`/`iac` arrays per solvent atom (`topology.rs:580`), *and* carries the
  solute/solvent boundary at `n_solute`. So it pays the duplication it was trying to avoid *and*
  inherits the structural split. This is fixable precisely because it's young.

### The core insight: two concerns are being fused

The label "solvent" is silently asserting **two independent things at once**:

1. **Representation** — "store this molecule's topology once and *instance* it N times" (a
   storage/deduplication concern; the genuine, valuable idea).
2. **Role / semantics** — "treat these atoms as rigid bulk: SETTLE/SHAKE-constrained, in the cutoff
   solvent fast-path, in the `s:` analysis namespace, usually excluded from fitting/restraints" (a
   labeling concern).

Because the two are fused, you *cannot* have a molecule that is stored compactly **but** treated as
solute (your "promote one catalytic water" case), nor an ion treated as solvent for one analysis.
Every painful situation you described is a place where you wanted concern 1 and concern 2 to point
different directions, and the data model wouldn't let them.

### The fix: instancing for representation, attributes for role

Decouple them. This is the data-model change the unification should make once, in `gromos-core`, so
that *both* the engine and the analysis façade inherit it:

- **Representation = a `MoleculeType` registry + instances.** A molecule type is stored once; the
  system holds *instances* (`{ moltype_id, atom_offset }`). Per-atom **parameters** (charge, IAC,
  mass, bonds, constraints) are looked up through the moltype — stored once, never expanded —
  exactly gromosXX's solvent trick, but **generalized to any repeated molecule**, not only water.
  You get the memory win for 1000 lipids or 100 copies of a protein, not just solvent. (This is also
  GROMACS's `moltype`/`molblock` model — proven at scale.) Positions/velocities/forces stay explicit
  per atom in the contiguous SoA arrays (Dim 1) — they must, every atom moves independently — but
  the *parameter* indirection kills the duplication gromos-rs currently has.
- **Role = per-atom / per-molecule attributes**, queried uniformly. "Solvent-ness", constraint
  group, temperature/energy/pressure group (already present), restraint membership, "analyzed" — all
  become *labels*, not a partition of the arrays. `is_solvent(i)` becomes `role(i) == Solvent`, a
  *filter*, not an index threshold. The `s:`/`m:` namespace split collapses into one AtomSpecifier
  grammar over attributes (this is the PLAN.md P2.1 selection work, done right from the start).

### Why this is a real overtake, and what it buys each consumer

- **Analysis (the façade, Dim 5):** programs iterate **one** atom array and filter by attribute — no
  more dual `mol()`/`sol()` code paths, no `m:`/`s:` bifurcation. The single thing that makes
  gromos++ analysis tedious disappears structurally.
- **"One water as solute, another as solvent":** trivial. It's two instances with different `role`
  labels, or a cheap **de-instancing** operation — `system.promote(water_42)` materializes that one
  instance as an explicit, editable molecule (so you can restrain it, perturb its LJ, give it special
  treatment) while the other 215 waters stay instanced and compact. In gromosXX this requires
  rebuilding the topology; here it's a local operation. **This is the ergonomic win that answers your
  exact frustration.**
- **Pairlist (Dim 9):** the moltype carries an optional "regular/rigid + fast-kernel" hint, so the
  solvent fast-path (gromosXX's `solvent_innerloop`/SPC table) is *preserved* — you still know
  "instance N is an SPC water, use the vectorized 3-site kernel" — but now it's keyed on *moltype*,
  not on the rigid solute/solvent boundary. Speed without the rigidity.
- **Composition (the topology builder):** instancing *is* how `water.instanced(216)` and `na * 3`
  should be represented — this dimension is the data-model foundation under the `*` operator in the
  py-gromos builder. The two sections are the same idea from two ends.
- **Memory/SoA:** parameters stored once per moltype, positions contiguous → smaller, more
  cache-friendly than today's expanded arrays.

### The cost to be honest about

Instancing adds one indirection (`atom → instance → moltype → parameter`) on every parameter
lookup. Mitigations: parameters that are truly per-atom-hot (charge, IAC) can still be materialized
into a flat array *as a cache* derived from the moltype (best of both — dedup as source of truth,
flat array as a rebuildable accelerator), while bonds/constraints/exclusions stay in the moltype
where the dedup matters most. De-instancing must correctly renumber charge groups and exclusions —
which is exactly why Dimension 9d (charge groups as a first-class primitive) is a prerequisite.

**Plan**
- In `gromos-core`: introduce `MoleculeType` + `MoleculeInstance`; make `Topology` a list of
  instances over a moltype registry. Keep `is_solvent`/`num_solute_atoms` working as *derived*
  attribute queries during migration so the reference suite never breaks.
- Replace `solvate()`'s array expansion with instancing; the flat `charge`/`iac` arrays become a
  rebuildable cache, not the source of truth.
- Add `role` as a per-instance attribute; route SHAKE/SETTLE, the solvent cutoff fast-path, and the
  AtomSpecifier `s:`/`m:` grammar through it.
- Implement `promote(instance)` (de-instancing) with CG/exclusion renumbering (needs 9d).
- Validate: every reference system must produce byte-identical output before/after the model change
  — this is a pure refactor of representation, the physics is untouched.

**Payoff:** the solute/solvent rigidity that annoys you stops being structural; the storage win
generalizes beyond water; analysis loses its dual code paths; and the builder's `*` operator gets
its foundation. **Risk:** medium-high — it's surgery on the central data structure — but it's a
*representation* change with an exact bit-for-bit oracle (the current output), so it's verifiable.
Do it **before** the analysis façade (Dim 5) and the topology builder calcify around the old split.

---

## Dimension 11 — Don't port the bugs: differential auditing of gromosXX

Your fear is well-founded and specific: a *faithful* port faithfully reproduces whatever is in the
C++ — including its bugs — and the 190k-LOC complexity hides them. I went through the physics paths
looking for this, and the result is nuanced enough to need its own discipline. **Three categories,
and you must treat them differently:**

- **(A) Genuine bugs you must NOT port** — defects in gromosXX. If you reproduce them, you've
  inherited a bug; if you "fix" them silently, your bit-for-bit reference test breaks and you can't
  tell a fix from a regression.
- **(B) Quirks you MUST port** — intentional-or-not GROMOS conventions that define "GROMOS
  behavior." Reproduce them exactly *and document them as deliberate*, or your output silently
  disagrees with every published GROMOS result.
- **(C) Latent divergences** — places gromos-rs already differs from gromosXX but **no reference
  test exercises that path yet**, so nothing has caught it. These are the dangerous ones: they look
  fine until you wire the feature, then everything downstream is subtly wrong.

The C++ developers flag many of their own issues in comments (`grep -rniE 'bug|fixme|wrong|hack'
interaction/ algorithm/ math/`). That grep is the single highest-yield audit tool — start there.

### Concrete findings from this pass

| # | Location | Finding | Class | Action for gromos-rs |
|---|----------|---------|-------|----------------------|
| 1 | `math.rs:90` vs `boundary_implementation.cc:285-318` | **gromos-rs already diverges on triclinic nearest-image.** gromos-rs uses fractional-coordinate `frac - frac.round()`; gromosXX uses iterative `while`-loop reduction in z→y→x order over the lower-triangular box. These are **not equivalent for strongly triclinic cells** — the fractional `round()` does not always yield the Cartesian nearest image. | **C (latent)** | Triclinic is "defined but not wired" (PLAN P1.4), so no test catches this *yet*. When you wire it, a triclinic reference test **will fail**. Port the `while`-loop z→y→x version *before* trusting triclinic; don't assume the textbook `round()` matches GROMOS. |
| 2 | `extended_grid_pairlist_algorithm.cc:1309` | Dev-flagged (`//BUG …--martina`): a pair is bucketed **solvent-solvent** when only `a2 >= num_solute_atoms`, without checking `a1`. A solute–solvent pair can be misclassified. | **A (real bug)** | You don't have a grid pairlist yet (Dim 9a). When you build it, classify pairs by **both** atoms' roles. **Root cause is the rigid solute/solvent index split — exactly what Dim 10 removes.** Fixing the data model dissolves the whole bug class. |
| 3 | `perturbed_nonbonded_term.cc:596,749` | Dev comments: *"Chris: CHECK! I'm not sure if the self-term correction is not wrong…"* on the perturbed RF self-term. | **A? (unverified)** | FEP/TI is not ported yet (PLAN P1.5). When you do, **derive the perturbed RF self-term independently** from the GROMOS book, don't transcribe the C++ — the original authors weren't sure. |
| 4 | `perturbed_nonbonded_term.cc:1444` | Dev comment: *"there is a bug here from the previous version."* | **A? (unverified)** | Same as #3 — flag the perturbed term as untrusted; second-source it when porting FEP. |
| 5 | `berendsen_thermostat.cc:105`, `nosehoover_thermostat.cc:129,185` | *"small flexible constraints hack!"* — special-casing in the thermostat DOF/coupling path. | **B (quirk)** | Nosé-Hoover is unported (PLAN P1.3). Reproduce the hack *if* you support flexible constraints; document it as deliberate so it doesn't look like a porting mistake later. |
| 6 | `latticesum.h:351,372,380` | Multiple `// wrong!!!` lines — but these are **commented-out** wrong formulas next to the corrected ones (Ewald/lattice-sum γ̂). | **benign** | Informational. When porting PME/lattice-sum, use the *active* formula; the `wrong!!!` lines document a trap the authors already fell into — don't re-derive into the same hole. |

The headline is **#1**: it's not a hypothetical, it's in your tree *right now*, and the only reason
it hasn't bitten is that triclinic is unwired. That is the exact shape of the problem you're
worried about — except here you didn't port a C++ bug, you "improved" past a GROMOS convention and
will lose bit-for-bit parity. Both directions of divergence are caught by the same discipline below.

### The discipline (this is the real deliverable)

The reference suite (PLAN.md, 28/28) is your bug oracle — but **only for wired paths.** Every
*unwired* path (triclinic, FEP, grid pairlist, Nosé-Hoover, EDS/GaMD, lattice-sum) is an
un-audited surface where category-C divergences and category-A inherited bugs both hide. So:

1. **No unwired path gets trusted without a reference test first.** Before wiring triclinic, write
   a triclinic gromosXX reference; it will immediately surface finding #1. Make this a rule: the
   reference test comes *before* the wiring, not after.
2. **Second-source the physics, don't transcribe it.** For anything the devs flagged as uncertain
   (#3, #4) or anything subtle (RF self-terms, virial, Ewald), derive it from the GROMOS book /
   primary paper independently, implement to *that*, then diff against the C++. Where they disagree,
   you've found either a C++ bug (A) or your own — investigate, don't assume the C++ is right.
3. **Grep the C++ for self-flagged issues per subsystem before porting it** — `bug|fixme|wrong|hack`
   over the files you're about to translate. Two minutes, catches the ones the authors already knew.
4. **When you must reproduce a quirk (B), encode it as a named, documented decision** (like PLAN.md's
   "BONDANGLETYPE intentionally unsupported" entries) so a future reader can't mistake a deliberate
   GROMOS-ism for a porting bug — and so you can later offer a "corrected mode" behind a flag.
5. **A `--gromos-compat` vs `--corrected` split, eventually.** Where gromosXX is genuinely wrong
   (e.g. #2 if it survives into a release), you may want bit-for-bit compat mode *and* a corrected
   mode. Rust makes this a clean feature/flag; design for it so "faithful" and "correct" can coexist
   instead of being forced to choose.

**Payoff:** you stop inheriting 30 years of accreted C++ defects *and* you stop silently diverging
from GROMOS where you meant to match it — with a repeatable method, not luck. **Risk:** low; it's
process, not code. The only cost is discipline: tests before wiring, derivation before transcription.

---

## Dimension 12 — QM/MM and ML potentials: kill the interop tax (the flagship dimension)

This is the one with the most upside — it's where the field is *moving* and where C++ gromos is
*structurally* weakest. It is also the dimension that most justifies starting the architecture
rethink **now**, before the classical engine calcifies, because the cost of retrofitting a clean
QM/MM/ML seam later is enormous and the cost of designing for it today is nearly zero. Treat this
section as the reason Dimensions 1, 2, 9 and 10 are worth doing carefully: they are the substrate
that makes this one possible.

### 12.1 — How gromosXX actually does it (grounded in the source)

gromosXX already has the *hooks*: `interaction/qmmm/` has ten workers — external-engine drivers
(`orca`, `gaussian`, `dftb`, `mndo`, `mopac`, `turbomole`, `xtb`), a `ghost` worker, and a
**neural-network worker** (`nn_worker.cc`). Reading the code, two architectural facts define the
weakness:

1. **The `QM_Worker` base class is shaped around file I/O.** Its virtual interface is
   `write_input` / `write_qm_atom` / `write_mm_atom` / `write_charge` / `open_input` /
   `run_calculation` / `open_output` / `process_output` (`qm_worker.h`). The *abstraction itself*
   assumes "format a text input file → fork/exec an external program → parse its text output." Every
   one of the ten workers is a hand-written instance of that, fragile to each package's output-format
   drift. A native, in-process engine **does not fit this shape** — you'd be forcing arrays through a
   file-shaped hole.
2. **The ML path embeds CPython.** `nn_worker.cc` uses `pybind11/embed.h` and drives SchNetPack
   *through an embedded Python interpreter* (`py::str model_path`, `py::list val_models_paths`). So
   every force evaluation crosses the C++↔Python boundary — GIL, per-step tensor marshalling — and
   the deployment couples a whole Python/SchNetPack environment to the MD binary. For ML-driven MD
   **that boundary is frequently the dominant cost**, and the deployment coupling is its own tax.

The genuine complexity lives in `QM_Zone` (`qm_zone.h`) and is worth respecting, not underestimating:
link/capping atoms (covalent bonds crossing the QM/MM boundary need capping H), MM point charges
within a cutoff of the QM region **rebuilt every step**, charge-on-spring (polarisable) embedding via
an SCF iteration, and distance-scaled charges. `QMMM_Interaction` is itself an `Interaction` subclass
— so even in gromosXX, QM/MM is conceptually "a force contributor." That validates the target design;
the problem is purely that the *worker abstraction* and the *ML runtime coupling* are wrong.

### 12.2 — The Rust overtake (why this is categorically better, not just faster)

ML interatomic potentials are the fastest-growing area in molecular simulation, and the incumbents
(ASE calculators, OpenMM-ML, the GROMACS NN bridges, gromosXX's `nn_worker`) **all pay the
Python-boundary tax** — it's the universal pattern, not a gromos failing. A native-Rust MD engine
that runs ML potentials *in-process on its own data layout* is a capability none of those stacks
have structurally:

- **Native in-process inference, zero language boundary.** `ort` (ONNX Runtime), `tch` (libtorch
  C++ API, no Python), or pure-Rust `candle`/`burn` run the model in the same process, reading the
  same SoA arrays (Dim 1), writing the same force accumulators (Dim 2). No interpreter, no GIL, no
  per-step serialize/deserialize. This is the layer that makes Python-coupled ML-MD slow, deleted.
- **QM/MM as an additive force provider, not a fork** — over the single-kernel data model (Dim 2),
  so embedding reads the same charges/positions everything else does.
- **One trait, many engines** instead of ten forked workers + an embedded interpreter.

### 12.3 — The architectural problems to overcome (be honest about these)

This is not free; the hard parts are real and most are *interface-design* problems you want to solve
on paper before writing engine code:

- **(P1) The abstraction must be data-in/data-out, not file-in/file-out.** Invert gromosXX: the
  trait takes arrays (positions, charges, region) and returns energy + gradients. File-based external
  engines are *one implementation* that internally formats/parses; the file dance must never leak
  into the trait shape. Getting this inversion right is the whole ballgame.
- **(P2) The model-export boundary.** Most trained potentials live in PyTorch/Python. The on-ramp is
  *train in Python → export → load in Rust*. But not every architecture exports cleanly: equivariant
  message-passing models (NequIP/Allegro/MACE, built on `e3nn`) use custom ops and dynamic neighbor
  graphs that **may not lower to ONNX**. TorchScript via `tch` is more permissive but heavier. Decide
  the on-ramp early and test export on the *specific* target architectures — this gates which models
  are first-class.
- **(P3) Multiple, differently-parameterized neighbor structures.** The MD pairlist, an ML model's
  radial graph (different/larger cutoff, edge lists, sometimes triplets/angles), and the QM-zone
  gather (MM charges within cutoff of QM, rebuilt each step) are *three different* spatial queries.
  The cell list (Dim 9) must therefore be a **general spatial-index service** that can emit all
  three, not a single-purpose MD pairlist. **This is a Dim 9 design constraint to lock in now.**
- **(P4) Units and conventions, typed.** QM works in Hartree/Bohr, ML models in eV/Å or kcal/mol,
  GROMOS in kJ/mol·nm. gromosXX manages this with `print_unit_factors` and discipline. Use Rust
  newtypes (or an explicit boundary conversion struct) so unit mismatches are *compile/return-time*
  errors, not silent 23×/27.2× force bugs — a classic, brutal ML-MD failure mode.
- **(P5) Forces flow onto MM atoms too.** Electrostatic embedding means the QM density polarizes and
  the *MM* point charges feel a force. So a provider returns gradients on **both** the QM region and
  the embedding MM atoms — i.e. contributions on an *arbitrary atom subset*, not "the QM region
  only." The accumulator interface must accept partial, scattered contributions (also what Dim 2's
  classical providers want).
- **(P6) Stateful, iterative coupling.** Polarisable embedding is an SCF loop (mutual QM↔MM
  polarization); QM workers reuse the previous wavefunction as an initial guess. The trait must
  allow `&mut self` (stateful workers) and iteration — it cannot assume pure functions.
- **(P7) Uncertainty / committee / active learning.** `nn_worker` already loads *validation models*
  (ensemble). Modern ML-MD wants ensemble disagreement → uncertainty → trigger a DFT recompute /
  active learning ("learn on the fly"). Let a provider optionally return **uncertainty** alongside
  energy/forces, and let the scheduler swap/escalate providers at runtime (ML → DFT on high
  uncertainty). Design the return type for this from the start.
- **(P8) Determinism boundary.** QM/ML are not bit-deterministic, so the Dim 11 reference discipline
  must classify these paths as "validated against reference energies/forces within tolerance," a
  *separate test tier* from the bit-for-bit classical suite.
- **(P9) Process lifecycle & failure** for external engines: crashes, timeouts, NaNs. Rust `Result`
  + supervised subprocess + fallback; async (Dim 6) to overlap the (slow) QM call with other work.

### 12.4 — The model landscape: SchNet today, but design model-agnostic

gromosXX bet on **SchNet** (schnetpack). The field has moved fast, and the interface must not be
SchNet-shaped:
- **SchNet** — continuous-filter convolutions; simple, relatively export-friendly. The incumbent.
- **NequIP / Allegro** — E(3)-equivariant message passing; high accuracy, but `e3nn` ops are an
  export-friction risk (P2).
- **MACE** — higher-order equivariant messages; current SOTA accuracy/efficiency. Crucially,
  **MACE-OFF / MACE-MP foundation models** target organic/biomolecular chemistry off-the-shelf —
  directly relevant to a biomolecular engine; a foundation potential you load and run, no training.
- **ANI** (Behler-Parrinello-style) — simpler descriptors, very ONNX-friendly; a pragmatic first target.
- **TorchMD-NET** (equivariant transformer), and emerging **universal/foundation potentials**.

What varies across them and therefore must be *parameters of the interface*, not assumptions:
cutoff radius; graph structure (pairs vs triplets/angles); equivariance (→ export path); whether
forces come analytically or via autograd; units; element/embedding conventions. **Design the trait so
swapping SchNet → MACE-OFF is a config change, not an interface change.**

### 12.5 — What to decide NOW (the "start immediately" deliverable)

The user's instinct is right and it's the single most valuable takeaway in this document: you can
make QM/MM/ML *seamless* later **for free** if you adopt a handful of invariants now, while the
classical engine is still malleable. These are not QM/ML code — they are shape constraints on the
classical refactors you're already going to do:

1. **Force computation = additive providers over shared state (Dim 2), contributing to arbitrary
   atom subsets.** Express even classical LJ+CRF / bonded as `Provider`s that read SoA state (Dim 1)
   and scatter into shared accumulators. If the *classical* engine is built this way, QM/MM/ML is
   literally just another provider — no retrofit. **This is the keystone decision; everything else
   follows.**
2. **A provider trait that is data-in/data-out and allows state + uncertainty.** Sketch it now even
   if only LJ implements it:
   ```rust
   trait PotentialProvider {
       fn contribute(
           &mut self,                      // &mut: stateful workers (SCF, wavefunction reuse) — P6
           region: &AtomSet,               // arbitrary subset, incl. embedding MM atoms — P5
           state: &State,                  // shared SoA positions/charges — Dim 1
           neigh: &dyn SpatialIndex,       // query the shared cell list — P3
           out: &mut Accumulators,         // scatter energy + per-atom forces + virial — Dim 2
       ) -> Result<ProviderExtra>;         // ProviderExtra: optional uncertainty, etc. — P7
   }
   ```
   The exact signature will change; the *shape* (subset in, scattered forces out, stateful, fallible,
   extensible return) is what must be right early.
3. **The cell list (Dim 9) is a spatial-index *service*, not an MD pairlist.** Its public API must be
   query-based ("give me neighbors of these atoms within r") so it can serve the MD pairlist, an ML
   graph, and the QM-zone gather from one spatial structure (P3).
4. **Region/role membership is a first-class data-model attribute (Dim 10).** "In the QM zone" is an
   attribute, and **promoting a solvent water into the QM region is exactly the Dim 10 de-instancing
   operation** — so Dim 10 is a *prerequisite* for clean QM/MM, not an independent track. Coupling
   and decoupling parts of the system is the instancing/role model doing its job.
5. **A typed units boundary** (P4) introduced with the engine, not bolted on at the QM seam.
6. **A second test tier** (non-bit-for-bit, tolerance-based, validated against reference
   energies/forces) stood up alongside the classical bit-for-bit suite (P8, Dim 11).

Adopt 1–6 and the eventual QM/MM/ML work is *additive*: new `Provider` impls + adapters, touching no
classical code. Skip them and you'll be doing gromosXX's retrofit — a forked subsystem bolted onto a
data model that never expected it.

**Sequence:** the QM/MM/ML *engines* are consumers of Dim 1+2+9+10 and land after them. But the
**invariants in 12.5 are decisions for the very next refactor**, not later — they cost almost
nothing now and are very expensive to retrofit. This is the concrete sense in which "this affects
the plan now": the classical force-evaluation refactor should be written as providers from the
start. **Risk:** medium-high on the engines (export friction P2, determinism P8, native-inference
parity); **low** on the invariants — they make the classical code cleaner regardless of whether QM/ML
ever ships.

---

## Dimension 13 — Coarse-grained as an ecosystem bridge (Martini)

gromosXX has coarse-grained support (`is_coarse_grained`, the `c16_cg` regression). The forward bet
isn't *just* supporting CG — it's making gromos-rs CG **compatible with Martini**, the dominant CG
force field, whose ecosystem lives in GROMACS. That turns CG support from a feature into a *bridge*:
import Martini topologies/parameters and run them in gromos-rs.

**The lever:** "overtake by interoperability." gromosXX is a closed GROMOS-force-field world. The
unified gromos-rs data model — instancing (Dim 10: CG beads are *massively* repeated, so the
moltype-registry dedup shines) + a pluggable force-field loader (the `ForceField` resolution context
from the compositional-topology section) — makes gromos-rs a candidate **polyglot engine**: GROMOS
IFP *and* Martini, behind one loader, on one engine. A `martini` import path sits alongside the
GROMOS IFP/MTB path.

**The honest caveat (this is the "I hope so" part):** Martini does **not** share GROMOS's nonbonded
conventions — different shifted-LJ treatment, combination rules, its own water/RF-vs-cutoff scheme,
specific bead-type interaction matrix. The bridge is contingent on whether Martini's interaction
treatment maps onto (or can be added to) the gromos-rs nonbonded kernel. **This needs a focused
investigation before it's a commitment** — it may require a Martini-specific nonbonded mode, not
just a parameter import.

**Plan (prototype, gated on the investigation)**
- Investigate Martini's nonbonded conventions vs the GROMOS LJ+CRF kernel; decide whether it's a
  parameter-import or needs a `martini` interaction mode.
- If viable: a Martini topology/parameter loader producing the same `Topology` the engine consumes;
  validate energies against a GROMACS-Martini reference on a small bead system.

**Payoff:** access to the entire Martini ecosystem from a Rust engine — a reach no GROMOS build has.
**Risk:** high/uncertain — conventions may not align; this is a research bet, not a port. (Note:
**polarisable/charge-on-spring force fields are explicitly out of scope** — low relevance for the
target use cases.)

---

## Dimension 14 — The unifying architecture (the meta-dimension behind all the others)

You liked `PotentialProvider` because it's the right *shape*. This dimension is the realization that
**`PotentialProvider` is one instance of a single, more general pattern** — and that getting the
whole-system architecture right is the decision that saves the most future reimplementation, because
it's the one that determines whether gromos-rs becomes one codebase or repeats GROMOS's two.

### 14.1 — gromosXX and gromos++ as wholes (the honest holistic read)

The two C++ codebases are **two architectures that do not share a core**, and that single fact is
the origin of the 190k+79k-LOC duplication:

- **gromosXX (the engine) is a pipeline of mutating algorithms over a shared state triad.** Three
  long-lived objects — `Topology` (fixed parameters), `Configuration` (mutable state, with
  `current()`/`old()` for leapfrog history), `Simulation` (run parameters) — are threaded through an
  ordered `Algorithm_Sequence`. `create_md_sequence.cc` literally `push_back`s the steps in order:
  COM-removal → Forcefield → EDS/GaMD → barostat → leapfrog-position/velocity → thermostat →
  constraints → energy/pressure calc. Each `Algorithm::apply(topo, conf, sim)` mutates the shared
  state in place. The `Forcefield` is itself an Algorithm wrapping a sub-pipeline of `Interaction`s.
  Clean, comprehensible — but the physics is *trapped inside mutating algorithm objects*, not exposed
  as reusable pure functions, and the backend forks (omp/mpi/cuda) live inside each.
- **gromos++ (the toolbox) is ~107 standalone frame-loop `main()`s over a *separate* core
  (`gcore::System`).** The canonical program (e.g. `rmsd.cc`): parse `Arguments` (@-flags), build a
  `System`, set boundary + gather method, open an `InG96` trajectory, `while (!ic.eof()) { ic >> sys;
  gather(); …analyze…; cout << result; }`, with `AtomSpecifier` for selection. Because it could not
  cleanly link gromosXX's internals, it **re-implements geometry, PBC, energy, statistics** in its
  own `gcore`/`gmath`/`bound`. That re-implementation *is* the 79k LOC.

So: **a pipeline architecture and a frame-loop architecture, bolted to two different data cores.**
The deep waste isn't any one duplicated function — it's that the two halves were never the same
machine.

### 14.2 — Where gromos-rs + py-gromos are heading right now

gromos-rs currently **mirrors gromosXX**: an `AlgorithmSequence` of `Box<dyn Algorithm>` over
`Topology`/`Configuration`/`State` (`gromos-core/src/algorithm.rs`). py-gromos/pyo3 is heading toward
a **compositional Python API**: `Simulation` plus per-step pyclasses (`Forcefield`,
`LeapFrogIntegrator`, `LeapFrogVelocity`, …) you assemble into a sequence
(`pyo3-gromos/src/{simulation,algorithm_sequence}.rs`). PLAN.md P2 already names the goal — make
gromos-analysis a *facade* over engine primitives, not a second core. **The trajectory is correct.**
The risk is doing it by *partial* sharing (some primitives reused, some re-implemented) and drifting
back into two cores by accident — exactly how gromos++ started. The architecture below is how you
make the sharing *total* and structural rather than incidental.

### 14.3 — What architecture could be "perfect"? Four candidates, honest trade-offs

The question deserves real options, not one decree:

- **(A) Pipeline-of-Algorithms** — the status quo (gromosXX/gromos-rs). `Vec<Box<dyn Algorithm>>`,
  each mutates shared state. *Pro:* simple, matches GROMOS, ordering is explicit. *Con:* physics is
  trapped inside mutating objects (not reusable by analysis); `dyn` dispatch at the coarse level;
  hard to parallelize/GPU/differentiate because "mutate in place" is the only contract.
- **(B) Provider/Operator split + explicit state** — the `PotentialProvider` idea generalized: cleanly
  separate **pure compute** (forces, energies, observables: data-in → data-out) from **mutating
  steppers** (integrators, constraints, thermostats that own and evolve state). *Pro:* the pure layer
  is **shared by engine and analysis** — this is the duplication-killer; it's parallel/GPU/diff
  friendly and testable. *Con:* you must draw the pure/mutating line deliberately (SHAKE mutates;
  is it a stepper? yes).
- **(C) ECS / data-oriented (Entity-Component-System flavor)** — atoms = entities; position, velocity,
  charge, role, qm-membership, charge-group = *components* (columns); providers/steppers = *systems*.
  *Pro:* maps **perfectly** onto SoA (Dim 1), instancing + role attributes (Dim 10), and adding a
  property (polarizability, QM flag) becomes "add a column," not struct surgery; archetype storage is
  cache-coherent. *Con:* full ECS frameworks are heavy and can over-abstract MD's tight numeric loops
  — take the *data-oriented column model*, not a game-engine framework wholesale.
- **(D) Dataflow / compute graph** (JAX-MD / TorchMD-NET style) — express the simulation as a DAG of
  ops the runtime schedules, fuses, lowers to CPU/GPU, and can **autodiff**. *Pro:* GPU (Dim 4),
  *differentiable simulation* (in-loop ML force-field training, sensitivity/inverse design), kernel
  fusion; the `LazyTopology` idea generalized to the whole run. *Con:* highest complexity; bit-for-bit
  parity is harder; premature to build now.

### 14.4 — The recommended synthesis: a layered architecture, where **analysis and MD are the same machine**

Don't pick one; layer them. Take B's compute/evolve split, store it on a C-style data-oriented core,
keep A's explicit sequence as the *default* orchestration, and leave D as an optional future
orchestration you don't foreclose:

```
L4  Facade & bindings   gromos++-style analysis API + py-gromos — thin over L0–L3, zero physics
L3  Orchestration       a run = a composition. MD run = loop over time; analysis = loop over frames.
                        Default: explicit Sequence (eager, GROMOS-faithful).  Optional later: compute graph (D).
L2  Evolve (mutating)   Steppers: integrators, constraints, thermostats, barostats — compose, own state history
L1  Compute (pure)      PotentialProvider (forces/energy) + Observable (analysis quantity) — pure over (Topo, State, SpatialIndex)
L0  Data core           SoA + instancing + role attributes (Dim 1/10) + SpatialIndex service (Dim 9e); ONE (Topology, State)
```

**The keystone insight — the thing that collapses gromos++ into a facade:**

> An **MD step** is an orchestration over L0 + L1 + L2 (loop over *time*, evolving state).
> An **analysis** is an orchestration over L0 + L1 only (loop over *trajectory frames*, no steppers).
> **They share L0 and L1 entirely.** `ener` calls the *same* `PotentialProvider` the MD step calls;
> `rmsd` calls the *same* geometry/fit primitives; gathering is one L0 service. The 79k LOC of
> gromos++ that re-implemented physics simply *does not get written* — it becomes L4 sugar over L1.

`PotentialProvider` was L1 for forces. `Observable` is the same shape for analysis quantities
(RMSD, RDF, energy decomposition). Both are pure functions over the shared core — which is why both
the engine and the toolbox can call them, and why there is exactly one implementation of each.

### 14.5 — The decisions that save the most reimplementation (ranked)

1. **One owned data core (L0), forever.** `Topology` + `State` defined once in `gromos-core`;
   readers/writers, engine, analysis, and Python *all* consume it. **Never grow a second
   `System`/`gcore`.** This is the single decision gromos++ got wrong; it is the whole ballgame.
2. **The pure compute layer (L1) is the only source of physics.** Every energy/force/observable is a
   provider/observable over L0; analysis and MD both call them. The moment analysis re-implements an
   energy term (today: `ener.rs` hardcodes LJ σ/ε), you've re-created gromos++. PLAN.md P2.2's
   "single-point energy entry point" is literally this rule.
3. **A small, *stable* trait taxonomy — resist proliferation.** Aim for ~five contracts:
   `PotentialProvider`, `Observable`, `Stepper`, `SpatialIndex`, `Reader`/`Writer`. Everything else
   is data or a composition of these. Decide it early; it's the contract all crates implement.
4. **Decide the state-history model once** (the `current()`/`old()` question). Leapfrog needs
   previous state; analysis needs frame windows; checkpointing needs snapshots. Pick one
   representation (e.g. a small ring buffer of `State`) so steppers, analysis, and IO all agree —
   a subtle decision that, gotten wrong, forces reimplementation across L2/L3/L4.
5. **Trajectory = a stream of `State`s (Dim 6).** One mechanism serves analysis input, checkpointing,
   and live MD output. Don't build separate "read a frame" and "write a frame" worlds.
6. **Python binds L0 + trait *objects*, never reimplements.** Every L1/L2 trait gets a thin pyclass;
   Python composition = building an L3 orchestration. The current pyo3 `Forcefield`/`LeapFrogIntegrator`
   classes are exactly this — formalize it as the rule so py-gromos can never drift into its own physics.
7. **Keep L1/L2 pure/composable enough that a compute-graph orchestration (D) can be *added* later**
   without rewriting physics. Build eager now (GROMOS-faithful, debuggable); don't foreclose lazy/GPU/
   differentiable. This is cheap optionality if you respect the L1 purity boundary, very expensive to
   retrofit if you don't.

### 14.6 — Lock now vs keep optional

- **Lock now** (they shape every near-term refactor): L0 single core (#1), L1 purity (#2), the
  five-trait taxonomy (#3), state-history model (#4). These are the "decide on the next refactor"
  items — same urgency as the Dim 12 prerequisites, of which they are the generalization.
- **Keep optional / don't build yet:** the dataflow graph (D), differentiable simulation, full ECS
  framework adoption. Design so they're *additive*, build them only when a concrete need (GPU at
  scale, ML-FF training in-loop) justifies the complexity.

**Bottom line:** the "perfect" architecture is **B (provider/operator) on a C-flavored data core (L0),
orchestrated by A (explicit sequence) today and optionally D (graph) later — with L0 and L1 shared by
the engine and the analysis facade so they are literally the same machine.** That sharing is the
decision that turns GROMOS's two codebases into one, and it is the highest-leverage call in this
entire document. `PotentialProvider` is where it starts; L0+L1 is where it pays off.

---

## Sequencing — how to actually get the win (and prove it)

**Dimension 14 (the unifying architecture) is the frame around this whole list, not a step in it.**
Its "lock now" decisions — one data core (L0), pure compute layer (L1) shared by engine + analysis,
the five-trait taxonomy, the state-history model — are the contract that the steps below implement.
Settle them first on paper; they cost nothing and prevent the gromos++ two-cores drift.

The dimensions are ordered by leverage, but the *safe build order* is:

1. **Benchmark harness first.** You can't claim overperformance without a baseline. Extend the
   existing `nonbonded_bench`/`math_bench`: add end-to-end MD-step, pairlist, and SHAKE benches;
   `cargo bench -- --save-baseline c++parity`; ideally run the gromosXX binary on the same systems
   and record wall-clock side-by-side. **This is PLAN.md P4's benchmarking item, pulled to the
   front** — without it everything below is vibes.
2. **Dimension 9a (cell-list pairlist)** — the real *asymptotic* gap; pull this early, it gates any
   system bigger than the reference suite. Validate by set-equality against the O(N²) oracle.
3. **Dimension 1 (SoA)** — the one real micro-win, and the prerequisite for 3/4/9b.
4. **Dimension 2 (single kernel)** — the architectural decision; make it *before* EDS/GaMD/REMD
   (PLAN.md P1.5) get written, or they'll fork like they did in C++.
5. **Dimensions 3 (SIMD default) + 9b (spatial reorder)** — fall out of 1+2; both are cache/compute
   wins on the now-unified loop.
6. **Dimension 10 (instancing + role attributes)** — the solute/solvent refactor; land it *before*
   the façade and the builder calcify around the old split. Pure representation change, bit-for-bit
   oracle. Depends on 9d (CG primitives).
7. **Dimension 5 (zero-copy façade)** — already the plan; now it inherits 1–3/9/10 for free, and
   loses its dual solute/solvent code paths.
8. **Compositional topology (py-gromos)** — extends PLAN.md P3; can proceed in parallel with the
   engine work since it reuses validated `gromos-tools` logic. Depends on 9d (CG/exclusion
   primitives) and 10 (instancing underpins the `*` operator).
9. **Dimensions 6, 7** — low-risk, high-credibility, can go in parallel.
10. **Dimension 12 (QM/MM + ML potentials)** — highest *strategic* value; a consumer of 1+2, so it
    lands after them, **but prototype the `PotentialProvider` trait early** so QM/MM/ML is designed as
    a force provider, not bolted on as a fork.
11. **Dimensions 4 (GPU), 8 (MPI)** — the payoff for having done 2 correctly; late, feature-gated.
12. **Dimension 13 (Martini/CG bridge)** — a research bet gated on a nonbonded-conventions
    investigation; not on the critical path, revisit once the engine + instancing model are solid.

**Cross-cutting throughout: Dimension 11 (differential audit).** Not a phase — a rule applied to
every step: grep the C++ for self-flagged bugs before porting a subsystem; write the reference test
*before* wiring an unwired path; second-source uncertain physics instead of transcribing it. Fix
the triclinic divergence (#1) the moment you wire triclinic.

**Non-negotiable invariant:** every step keeps the 28/28 reference suite green. Performance work
that breaks bit-for-bit parity is a regression, not an optimization — except the explicitly
"fast/approximate" GPU path, which is validated against the CPU reference rather than bit-matched.

---

## The honest bottom line

You will **not** beat gromosXX by out-optimizing its inner loop — its templates already compile to
roughly what your generics do, and you've matched its OMP with rayon and its intrinsics with `wide`.

You **will** beat it by refusing the architecture C++ forced on it: the N-way forked physics. Build
**one** SoA-backed, SIMD-first kernel that every backend (serial, parallel, GPU, distributed) and
every consumer (MD engine, analysis façade) calls — and the C++ codebase's 190k+79k lines of
hand-synchronized duplication becomes your competitive advantage, because each of those forks is a
place gromos can't move and you can. The overperformance isn't a faster benchmark on day one; it's
that on day 365 you've shipped GPU, portable SIMD, streaming analysis, and reproducibility
guarantees from a single physics implementation, while the C++ trees are still keeping five copies
in sync.
