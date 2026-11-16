# T-Flows: solver and records of the runs of the paper

This branch holds, as they were used, the solver and the saved records behind
the results of

> J. Kren, B. Ničeno, Y. Sato, *Sharp-interface VOF method for phase-change
> simulations on unstructured meshes*, Computers & Fluids (in revision, 2026).

For a documented version of the same solver with cases that run and are
checked from the tree, see the branch `release`.

## Layout

One folder per case, one subfolder per study. Each case folder has a
`README.md` that says which part of the paper each study belongs to and which
build of the solver made it.

| folder | case | paper |
|---|---|---|
| [`Stefan/`](Stefan) | Stefan problem | Sections 3.1, 3.3 |
| [`Sucking/`](Sucking) | sucking problem | Sections 3.1, 3.3 |
| [`Scriven/`](Scriven) | bubble growth in superheated liquid; the anisotropy studies | Sections 3.2, 3.3, 4; Appendices B, C |
| [`StaticDroplet/`](StaticDroplet) | static droplet, spurious currents | Appendix D |
| [`Kinematic/`](Kinematic) | kinematic expansion test and standard advection tests | Section 4.4, Appendix A.2 |
| [`scripts/`](scripts) | analysis and plotting scripts | |
| [`T-Flows/`](T-Flows) | the solver source and the sign patch | |

The same names recur inside the case folders:

| subfolder | what it holds |
|---|---|
| `setup/` | the control file, mesh definition and hooks (`User_Mod/`) from which the runs of the revision were made |
| `submitted/` | the case as published with the submitted version |
| `budgets/` | the mass and energy budget runs of the submitted version |
| `revision/` | the runs made for the revision |
| any other | one study; each run is one subfolder with its `control` and its records |

A run folder holds what the run read and wrote, except meshes, field output
and records larger than 1 MB (listed in `OMITTED.tsv`).

## The solver

`T-Flows/Sources` is the source of commit `e36541a802` of the development branch
`unstructured-boiling-release`: commit `2c095ac070` (2026-11-11) of
`T-Flows-jan-review`, with which the runs of the revision were made, without
the parts of the development version that the paper does not use. It is
reduced as in the branch `release` to `Generate`, `Convert`, `Divide`,
`Process/Cpu`, `Shared` and `Libraries`. The correction of the liquid source
term of the volume-fraction equation is `T-Flows/patches/vof_source_sign.patch`
and is not applied.

The runs of the submitted version were made between November 2025 and August
2026 with at least seven earlier states of the source, the runs of the
revision with two. The case READMEs name the state of every study. Every
commit they name is an ancestor of `2c095ac070`. T-Flows does not print the
commit it was built from, so a run was tied to a source state afterwards, from
build and job scripts, the md5 of the `Process` kept in its folder, and the
set of control keywords the binary reads: *verified* means that a script, a
log or a binary shows it, *stated* that only the authors' notes say so.

Between the earlier states and `2c095ac070` lie the corrections of the VOF
module made up to October 2026. Those that touch code these cases run:

| commit | change | cases it can affect |
|---|---|---|
| `3b620baae7` | floor on the distance from a cell centre to the interface, 1 % of the cell size | all phase-change cases |
| `47ff80a316`, `bd051e76db` | one heat capacity per cut face and no advection across it in the energy equation | all phase-change cases |
| `9bdd6e62e8`, `8d36302bb7` | interface distance as a projection on the interface normal | polyhedral phase-change runs made before 31 July 2026 |
| `2a86af041a` | fluid properties in boundary cells (walls were free slip before) | static droplet |
| `7fc14fcc00`, `ada705e1e6` | volume-fraction gradient next to symmetry planes | Stefan, Sucking |
| `4e76d2fe20` | boundary faces kept off the front | cases whose interface meets a boundary |
| `9e4b1f2482` | pivot of the incomplete-Cholesky preconditioner for non-symmetric matrices | all |
| `dbf23475f4` | compiler option for double precision | all |

What was measured between builds: the source of this branch gives the same
records as `2c095ac070`, byte for byte, for the Stefan run `R_urf1` and the
kinematic run `st075-dt2-step`; the review solver and the solver of
`6adb80f599` are compared term by term for the Stefan problem in
`Stefan/revision`; the static droplet of the revision reproduces the submitted
traces to 7 × 10⁻¹⁰ (hexahedral) and 4 × 10⁻⁷ (polyhedral).

## Not included

- **Meshes and field output.** The polyhedral meshes are 0.7 to 5.8 GB each,
  the field output of the runs several hundred gigabytes. The structured
  meshes are made by `Generate` from the `.dom` files; the polyhedral ones
  cannot be regenerated from this branch.
- **Records larger than 1 MB**: full solver logs and long traces, listed in
  `OMITTED.tsv`.
- **The annular demonstration** (Section 5 of the paper).
- **Runs and scripts that are not of the paper**: a few folders of the
  archive belong to solver options that the paper does not use.
- **The reference data of Bureš & Sato (2021)** used in the growth-rate
  figure, and **the manuscript sources**.
- **Copies.** Records that were copied out of a run folder for the analysis
  are kept once, in their run folder.

## Notes for readers

- The files are kept as they were used, with one substitution: the home
  directory of the workstation, its host name and the user name on the
  cluster are replaced by `/home/user`, `workstation` and `user`. Nothing
  else is changed. Absolute paths in scripts and records therefore name
  folders of the machines on which the runs were made, and comments refer to
  notes that are not included.
- The analysis scripts read the field output through paths or environment
  variables that point outside this branch; those marked **S** in
  `scripts/FIGURES-submitted.md` need nothing but NumPy and SciPy. Figure and
  table numbers in that file are those of the submitted version.
- The hooks of a case are compiled into `Process`. Build one binary per case
  (`Cases/build_case.sh` of the branch `release` does this safely).

## Licence

GPL-3.0 (`LICENSE`), the licence of T-Flows since 30 August 2026; see
`NOTICE.md`.

Paul Scherrer Institute PSI, Laboratory for Simulation and Modelling.
