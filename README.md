# T-Flows: sharp-interface phase change on unstructured meshes

**The version of the flow solver T-Flows behind the paper**

> J. Kren, B. Ničeno, Y. Sato, *Sharp-interface VOF method for phase-change
> simulations on unstructured meshes*, Computers & Fluids (in revision, 2026).

[T-Flows](https://github.com/DelNov/T-Flows) is a finite-volume solver for
incompressible flows on unstructured meshes of arbitrary polyhedral cells,
written in modern Fortran. This branch holds the part of it that the paper
uses: the pre-processors, the CPU solver with the phase-change model, and the
cases with which the model can be run and checked.

What the phase-change model does:

- transports the volume fraction with an **algebraic VOF scheme** (CICSAM or
  STACS, with a boundedness corrector and sub-stepping) on hexahedral and
  polyhedral cells;
- reconstructs the interface as the surface of volume fraction 1/2 (the
  *front*);
- computes the temperature gradient on each side of the interface with a
  **least-squares stencil in which the arm that crosses the interface ends at
  the interface, at the saturation temperature**;
- obtains the **mass transfer rate from the jump of the conductive heat flux**
  at the interface (Stefan condition);
- imposes the saturation temperature in the energy equation through the flux
  over the cut distance (a Dirichlet condition at the interface);
- couples this to one-fluid momentum with continuum surface force, smoothed
  curvature and Rhie–Chow interpolation.

## Quick start

Requirements: gfortran and an MPI library (the cases of the paper were built
with gfortran 7.5 and 14.2 and Open MPI), GNU make, Python 3 with NumPy for
the checks. Linux, 64 bit.

```bash
Cases/build_tools.sh                 # Generate, Divide, Convert      (1 min)
Cases/build_case.sh Cases/Stefan     # Process with the Stefan hooks  (2 min)

cd Cases/Stefan
python3 run_stefan.py first 100 5.0e-3 1800 SIMPLE_UNDERRELAXATION_FOR_ENERGY=1.0
python3 check_stefan.py runs/first   # PASS: equal to the reference record
```

The run takes three to four minutes on one core. Building by hand, and the
traps of building for several cases in one tree, are described in
[doc/building.md](doc/building.md).

## Cases

| case | what it is | status in this tree |
|---|---|---|
| [Cases/Stefan](Cases/Stefan) | one-dimensional Stefan problem with the mass and energy budgets of Section 3.3 | run from this tree, equal to the record of the paper's revision bit for bit |
| [Cases/Kinematic](Cases/Kinematic) | radial expansion of a sphere with prescribed fluxes (Section 4.4) | run from this tree on a mesh made here, equal to the record of the paper's revision bit for bit |

The Sucking, Scriven, static-droplet and annular cases are not staged in this
branch yet, and the polyhedral meshes of the paper (0.7 to 5.8 GB each) are
not distributed. The set-ups and the saved records of the runs of the paper
are in the branch `paper-runs`.

## Documentation

| | |
|---|---|
| [doc/phase-change.md](doc/phase-change.md) | the model, where each part of it is in the source, and the control keywords |
| [doc/building.md](doc/building.md) | building the four programs, and building for a case |
| [Cases/README.md](Cases/README.md) | how the cases are organised |
| [KNOWN_ISSUES.md](KNOWN_ISSUES.md) | what is known to be wrong or untested in this version |
| [NOTICE.md](NOTICE.md) | authors, licence, third-party code |

The T-Flows manual, with the general use of `Generate`, `Convert`, `Divide`
and `Process`, is the `readme.md` of the
[T-Flows repository](https://github.com/DelNov/T-Flows).

## What this branch is, exactly

- **Source.** `Sources/` is the tree of commit `0ce2122640` of the branch
  `unstructured-boiling-release` of the authors' development repository,
  reduced to `Generate`, `Convert`, `Divide`, `Process/Cpu`, `Shared` and
  `Libraries`. That commit is commit `2c095ac070` (2026-11-11) of the
  branch `T-Flows-jan-review`, with which the runs of the paper's revision
  were made, without the parts of the development version that the paper
  does not use (listed below). What the paper uses is unchanged: both cases
  of this branch reproduce the records of those runs bit for bit, and so do
  the Sucking, Scriven, static-droplet and restart tests made before and
  after each removal. Two more things differ from the commit:
  symbolic links that pointed outside these folders are removed, and the one
  of them that the build depends on
  (`Sources/Process/Cpu/User_Mod/Vof_Initialization_Ellipsoid.f90`) is
  replaced by the file it pointed to. `MANIFEST.json` lists every file with
  its sha256 and everything that was left out.
- **Not included:** the GPU port, the T-Flows test suite and manual, the
  utilities; and from the solver, the features of the development version
  that no case of the paper uses: the geometric VOF (with the surfactant
  transport, the curvature fit at walls and the interface viscosity that
  only it could switch on), nucleation with the microlayer and near-wall
  models, the coupling of several domains (conjugate heat transfer), mesh
  adaptation, variable surface tension with the CSS and Marangoni terms,
  and the Lee model of mass transfer. A control file that asks for any of
  them stops the run with a message that names the keyword.
- **Kept, although the paper does not use them:** the turbulence models,
  particle tracking, porosity and STL input, the PETSc interface, passive
  scalars and the other parts of T-Flows itself, as they are in the
  development branch. They are not documented or checked here. The
  synthetic inlet eddies and the surface exporter are kept for the annular
  case of the paper.
- **Which results it reproduces.** The runs of the paper were made with
  several builds of the solver; this source is the latest of them, with one
  correction made after them, the sign of the phase-change source of the
  volume-fraction equation (`KNOWN_ISSUES.md`, item 1). The branch
  `paper-runs` states which build produced which result.

## Layout

```
Sources/            Generate, Convert, Divide, Process/Cpu, Shared, Libraries
Binaries/           where the programs are put by their makefiles
Cases/              cases, build scripts, default user hooks
doc/                documentation of the phase-change model and of the build
patches/            the sign correction of KNOWN_ISSUES.md item 1 (patch -R undoes it)
MANIFEST.json       every file with its sha256; what was left out
```

## Licence

GPL-3.0 (`LICENSE`), the licence of T-Flows since 30 August 2026. Third-party
code and attributions: [NOTICE.md](NOTICE.md).

Paul Scherrer Institute PSI, Laboratory for Simulation and Modelling.
