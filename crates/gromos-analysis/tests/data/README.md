# gromos-analysis reference outputs (gromos++, `.local/gromosPlsPls`, 2026-08-30)

`shared/` is `crates/gromos-md/tests/gromosXX_references/shared/`, `aladip.trc` the trajectory
`gromosXX_references/aladip_solvated/expected/trajectory.trc`. gromos++'s gather/Boltzmann warnings
are stripped from the stored outputs.

| dir | gromos++ command |
|-----|------------------|
| tser | `tser @topo shared/aladip.topo @pbc r @prop "d%1:1,2 a%1:1,2,3 t%1:1,2,3,4 tp%1:2,3,4,5%60" @traj aladip.trc @dist 5` |
| mdf | `mdf @topo shared/aladip.topo @pbc r @centre 1:1 1:6 @with 1:3-12 @traj aladip.trc` (writes `MIN_1:1.dat`, `MIN_1:6.dat`) |
| dg_ener | `dg_ener @temp 300 @stateA eA.dat @stateB eB.dat` and `… @col 3 3` (synthetic series, `eA`/`eB` from a seeded generator) |
| dfmult | `dfmult @temp 300 @stateR eR.dat @endstates eX1.dat eX2.dat` |
| matrix_overlap | `matrix_overlap @m1 m1.dat @m2 m2.dat @dim 3` |

`tests/gromospp_references.rs` runs ours with the same arguments and compares every number
(1e-6 relative; `-nan` matched literally) and every label.
| bilayer_dist | `bilayer_dist @topo shared/aladip.topo @pbc r @atoms 1:a @selection s:OW @grid 10 @traj aladip.trc` and `… @density` |
| bilayer_oparam | `bilayer_oparam @topo shared/aladip.topo @pbc r @atoms 1:2-11 @refvec 0 0 1 @traj aladip.trc` |
| jval | `jval @topo shared/aladip.topo @pbc r @jval aladip.jval [@timeseries or @rmsd] @traj aladip.trc` |
| edyn | `edyn @topo shared/aladip.topo @pbc r @atoms 1:a @eigenvalues 1 2 @traj aladip.trc` — `EIVAL`, `EIFLUC`, `COVAR`, `COVATOM` compared; gromos++'s projections (`EIVEC`, `EVPRJ`, `ESSDYN`, `PRJMAX/MIN`) use the covariance columns instead of the eigenvectors (`edyn.cc` discards the return value of `diagonaliseSymmetric`) and are not compared — ours use the eigenvectors, checked by an invariant test |
| gca | `gca @topo shared/aladip.topo @pbc r @prop "t%1:<dihedral.txt>%30%-60%60" @traj shared/aladip.conf`, `… @prop "d%1:1,2%0.2 a%1:1,2,3%100" @mobile first`, `… @prop "a%1:1,2,3%100"` — POSITION blocks compared (gromos++ also copies the input VELOCITY block) |
