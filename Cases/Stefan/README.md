# Stefan problem

A planar interface between vapour and saturated liquid moves away from a
heated wall. The vapour conducts heat from the wall (10 K above saturation)
to the interface; the liquid is at the saturation temperature and is pushed
out through the far boundary. Water and steam at 1 bar, 2 mm domain, interface
at 0.5 mm at the start. Section 3.1 and 3.3 of the paper.

## Run

```bash
../build_tools.sh            # once
../build_case.sh .           # Process with the hooks of this case
python3 run_stefan.py first 100 5.0e-3 1800 SIMPLE_UNDERRELAXATION_FOR_ENERGY=1.0
python3 check_stefan.py runs/first
```

`run_stefan.py <name> <N> <dt> <steps> [KEY=VALUE ...]` makes `runs/<name>`,
writes the mesh with `N` cells along the domain, replaces or appends the
given keywords in `control`, and runs on one rank. `--cubic` makes the cells
cubes; `--stage` only prepares the folder.

The check prints the position of the front against the Stefan solution and
compares `caf-terms.dat` with `expected/caf-terms.dat`. That file is the
record of the run `R_urf1` of the review series (the command above, made
with the full development tree). Run from this tree with gfortran 7.5, the
result is equal to it in every digit.

## Files

| file | |
|---|---|
| `control` | the control file of the submitted paper: under-relaxation of the energy equation 0.95, 6 to 12 outer iterations |
| `stefan.dom`, `box_parameters.ini` | mesh definition for `Generate`; the box that holds the vapour at the start |
| `User_Mod/` | initial condition, and the terms of the mass and energy budgets |
| `expected/caf-terms.dat` | reference record |

## What the hooks write

`caf-terms.dat`, one row per time step. Rates are written, not integrals.

| column | quantity |
|---|---|
| 1, 2 | time step, time [s] |
| 3 | mass in the domain [kg] |
| 4 | vapour volume [m³] |
| 5, 6 | mass and volume flux through the outflow [kg/s], [m³/s] |
| 7 | sum of the mass transfer rate, last estimate of the step [kg/s] |
| 8 | sum of the mass transfer rate that the volume-fraction step used [kg/s] |
| 9 | heat entering through the wall [W] |
| 10, 11 | heat removed by the interface coupling of the energy equation, vapour and liquid side [W] |
| 12, 13 | sensible enthalpy above saturation, vapour and liquid cells [J] |
| 14 | mean position of the points where the front cuts the cell connections [m] |
| 15 | number of faces the front intersects |
| 16 | largest sum of face fluxes over the cells [m³/s] |
| 17, 18 | largest deviation of the temperature from saturation, smallest temperature [K] |
| 19, 20 | conductive heat flux at the front elements, vapour and liquid side [W] |

`caf-vof.dat` holds the terms of the volume-fraction step. With a file
`caf_mode.ini` (`run_stefan.py --mode 1 <flux>`) the hooks impose a uniform
temperature and a given mass flux, the check of the convective term of the
energy equation.

## The settings of the paper's revision

The submitted run is `python3 run_stefan.py base 100 5.0e-3 1800`. Its energy
balance is off by 3.5 %, which comes from two settings: the under-relaxation
of the energy equation with six outer iterations, and the higher-order
convective correction on the faces cut by the front. The series of the
revision removes both with
`SIMPLE_UNDERRELAXATION_FOR_ENERGY=1.0 ADVECTION_SCHEME_FOR_ENERGY=upwind`.
That series was run with an earlier build of the solver; the branch
`paper-runs` has its job lists and records.
