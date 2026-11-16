# Scriven bubble growth

A vapour bubble of 50 µm grows in liquid 1.25 K above saturation, in a cube of
300 µm. The validation of Section 3.2 and 3.3, and all studies of the
four-fold anisotropy of Section 4 and Appendices B and C.

Mesh levels are named by the node count per edge of the `.dom` block: `075`
is 74³ = 405 224 cells. `st` and `aniso` are structured meshes, `poly`,
`run-` and `hcpoly` polyhedral ones.

| folder | what | paper | source state |
|---|---|---|---|
| `setup/` | control file, mesh definition, initial temperature profile and hooks of the 40³ runs of the revision (`solver-settings/`) | | |
| `submitted/` | the case as published with the submitted version; see the note below | | |
| `mesh-convergence/polyhedral/` | polyhedral 75³ to 150³ at 2 µs, with the job scripts | Section 3.2 | a state just before `9bdd6e62e8`, with the projected interface distance not yet committed; its record equals that of `9bdd6e62e8` over 600 steps, verified |
| `mesh-convergence/structured/` | structured 75³ to 150³ at 2 µs (`aniso-*-off`), radius records only: no control file survives | Section 3.2 | branch `extended-gradient-stencil`, at or after `a41ccae804`, stated |
| `time-step/` | structured 75³ at 2, 1 µs (`st75-dt2000`, `st75-dt1000`), 100³ and 125³ at 1 µs; polyhedral 75³ at 1 and 10 µs (`run-75-dt1`, `run-75-dt10`) | Section 3.2 | structured: `9bdd6e62e8`, verified; polyhedral: as `mesh-convergence/polyhedral` |
| `budgets/` | budget runs of the submitted version (`scriven`, `scriven75`, `scriven125`) and the polyhedral 75³ budget run `run-75-budget` | Section 3.3 | `bdf9354578`, verified; `run-75-budget`: `9bdd6e62e8`, verified |
| `solver-settings/` | 40³, 0.4 ms, one setting changed per run: `base`, `urfT1` (no under-relaxation of the energy), `upwindT2` (upwind energy convection), `outer30`; `SC_*` vary the iteration cap of the linear solvers | Section 3.3; response on the linear solvers | `6adb80f599`, verified |
| `anisotropy/uniform-mdot/` | the computed mass flux replaced by a uniform one (`UNIFORM_MASS_TRANSFER yes`), with and without transport | Section 4.3, Appendix C.3 | `9bdd6e62e8`, verified |
| `anisotropy/transport-removed/` | `st75-inflate-unif`: the bubble grown by the source alone, with the volume-fraction transport removed through the phase properties (see its `control`); `sigmatch-075`: the same with σ/ρ_l matched; `st75-eqrho*`, `st75-noadv-unif`: the equal-density attempts, which stop the growth and are not used | Section 4.3, Appendix C.5 | uncommitted states of 3 August 2026 between `9bdd6e62e8` and `6adb80f599`, verified |
| `anisotropy/surface-tension/` | surface tension ×0.25, 1, 2, 4 at 0.1 µs | Appendix C.6 | `9bdd6e62e8`, verified |
| `anisotropy/blend-weight/compressive-ladder/` | blend weight frozen at the compressive limit (`VOF_BLEND_GAMMA 1.0`), 75³ to 150³, and the 150³ baseline | Section 4.7, Appendix C.10 | the state of `6adb80f599`, built ten hours before that commit, verified |
| `anisotropy/blend-weight/limits-075/` | 75³ with the blend weight at 1 (`st75-hyperc`), at 0 (`st75-uq`), and the reference (`st75-gamma-off`) | Section 4.5, Appendix C.7 | an uncommitted state of 4 August 2026 just before `6adb80f599`, verified |
| `anisotropy/blend-weight/polyhedral/` | the compressive limit on polyhedral meshes; the boundedness corrector does not converge there | not used | as `limits-075` |
| `anisotropy/blend-weight/superseded-2us/` | the first ladder at 2 µs, replaced by `compressive-ladder` | not used | as `limits-075` |
| `anisotropy/constant-beta/` | one constant acceptor weight on every face (`ISOTROPIC_VOF_BETA`); it breaks boundedness, see the `# VOF_CORR` lines | not used | uncommitted states of 3 August 2026 |
| `anisotropy/extended-stencil/` | the coupled extended gradient stencil, 75³ to 150³ (`-extonly`, `-extw1`), radius records only | Appendix B.3 | branch `extended-gradient-stencil`, at or after `a41ccae804`, stated |

`bench-data.dat` has one row per time step: time, interfacial area, vapour
volume, volume-equivalent radius; later builds add a column. The paper reports
the radius of the front, computed from the front files, which are not
included.

**`submitted/` is not the set-up of the convergence runs.** Its `control` is
that of the 10 µs runs, its structured `bench-data-struct-*.dat` are records
of 10 µs runs (`bench-data-struct-100.dat` opens with a row of NaN), its
polyhedral records were made before the interface distance was projected, and
its `User_Mod/Vof_Initialization_Ellipsoid.f90` initialises no bubble on a
polyhedral mesh (it loops over `Grid % cells_n_nodes(c)`, which is negative
for a polyhedral cell; the hooks in `setup/` and `budgets/` take the absolute
value). Its `RUNS.md` gives figure and table numbers of the submitted version.

**`solver-settings/upwindT` did not run with upwind convection**: its control
file has `ADVECTION_SCHEME_FOR_ENERGY` twice and the solver reads the first,
`minmod`; its record equals that of `base`. `upwindT2` repeats it with the
keyword replaced: the growth over the first 0.4 ms is 5.9 % lower than in
`base` (2.9 × 10⁻⁴ lower in `urfT1`). `outer30` stopped after four rows.
