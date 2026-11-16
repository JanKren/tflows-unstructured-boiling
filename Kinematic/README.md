# Kinematic tests

The volume-fraction step alone, with face fluxes prescribed by the hooks:
momentum, pressure and energy are not solved, and there is no mass transfer
and no surface tension. Section 4.4 and Appendix A.2 of the paper.

| folder | what | paper | source state |
|---|---|---|---|
| `setup/` | the hooks (`User_Mod/`) and the control templates `control.expansion` and `kinematic.in.expansion` | | |
| `expansion/` | radial expansion of a sphere with the exact flux of a point source, on the Scriven cube: `st*` structured, `po*` polyhedral, `diag*` with the field output kept for diagnosis; with the staging scripts (`setup_expansion.sh`, `run_case.sh`, `queue*.sh`), `scripts/` and `logs/` | Section 4.4 | `2c095ac070`, verified |
| `advection/` | translation of a disc and of a sphere, and the three-dimensional deformation field, on cubes and polyhedral meshes at several Courant numbers, with their staging scripts | Appendix A.2 | `2c095ac070`, verified |

Run names in `expansion/`: `st075` is the structured 75³ mesh, `dt2` a time
step of 2 µs; `-step` keeps the flux only where the mean volume fraction of a
face's cells exceeds 1/2 (the velocity jump), `-sol` keeps it everywhere
(solenoidal), `-mask` has no flux inside the exact sphere, `-mix` scales the
flux by the mean volume fraction, and `-samp` starts from a sub-sampled volume
fraction. `kinematic.in` of each run states its settings.

`bench-data.dat` has one row per step: time, front area, volume of the sphere
phase, radius from the volume, radius of the front, minimum and maximum of
the volume fraction, the error norm E₁ against the initial field, and the
volume of the cells with a volume fraction between 0.01 and 0.99.
`kin-flux.dat` records the prescribed flux and its divergence, `kin-final.dat`
the error norms at the end.

The hooks in `setup/User_Mod` are dated after the binaries that made these
runs. With these hooks and the source in `T-Flows/`, the run
`expansion/st075-dt2-step` is reproduced bit for bit (branch `release`,
`Cases/Kinematic`); the advection runs were not repeated.
