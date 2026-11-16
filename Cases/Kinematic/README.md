# Kinematic tests

The volume-fraction step alone, with prescribed face fluxes. Momentum,
pressure and energy are not solved (zero outer iterations), and there is no
mass transfer and no surface tension. The hook
`User_Mod/Beginning_Of_Compute_Vof.f90` sets the face fluxes before every
step; each field is solenoidal on the mesh by construction.

| mode | field | paper |
|---|---|---|
| 1 | radial expansion: the exact flux of a point source through each face, with the source strength following the analytical radius of the Scriven bubble | Section 4.4 |
| 2 | uniform translation, optionally reversed | Appendix A.2 |
| 3 | three-dimensional deformation field, reversing in time | Appendix A.2 |

## Radial expansion

```bash
../build_tools.sh                  # once
../build_case.sh .                 # Process with the hooks of this case
./make_mesh.sh 75 8                # the Scriven cube, 74^3 cells, 8 subdomains
./setup_expansion.sh st075-dt2-step meshes/st075 2.0e-6 500 25 3 0
./run_case.sh st075-dt2-step 8     # 4 min on 8 cores
python3 check_expansion.py st075-dt2-step
```

The sixth argument of `setup_expansion.sh` selects the velocity inside the
interface band:

| value | velocity | result on the structured mesh |
|---|---|---|
| 0 | the liquid velocity continued through the interface; solenoidal | the sphere stays round |
| 3 | the flux is kept only on faces whose two cells have a mean volume fraction above 1/2, so the velocity vanishes behind the computed interface | the four-fold deformation of the coupled runs |
| 1, 2 | no flux inside the exact sphere; flux scaled by the mean volume fraction | variants |

`expected/st075-dt2-step` is the record of the run of that name behind
Section 4.4 (structured 75³, time step 2 µs, velocity jump).

## Files

| file | |
|---|---|
| `bubble.dom` | the 300 µm cube of the Scriven case for `Generate` |
| `control.expansion`, `kinematic.in.expansion` | templates that `setup_expansion.sh` fills in |
| `control.adv` | control template of the translation and deformation tests |
| `User_Mod/` | the prescribed fluxes, the initial sphere, the records |
| `make_mesh.sh`, `setup_expansion.sh`, `run_case.sh` | staging and running |
| `check_expansion.py` | radii against the exact radius; comparison with the record |

`kinematic.in`, read by the hooks:

```
1                        mode: 1 expansion, 2 translation, 3 deformation
0.0 0.0 0.0  50.0e-6     sphere centre and radius
4.0647 1.680e-7 0.0 0.0  mode 1: growth constant and liquid diffusivity
3                        velocity inside the interface band (table above)
0                        initial volume fraction: 0 production, 1 sub-sampled
```

## Records

`bench-data.dat`, one row per step: time; area of the front; volume of the
sphere phase; radius from the volume; radius of the front (area-weighted
distance of the front elements from the centre); minimum of the volume
fraction; maximum minus one; the error norm E₁ against the initial field;
volume of the cells with a volume fraction between 0.01 and 0.99.

`kin-flux.dat`, one row per step: step; time; source strength; exact radius;
exact volume; largest divergence of the prescribed fluxes as a fraction of a
cell volume per step; number of cells that hold the source; their flux sum.

`kin-final.dat`: the error norms and the volume change at the end.

The shape measure of the paper (front radius on the body diagonals minus that
on the axes) is computed from the front files `bubble-front-ts*.pvtu` by
`analyze_kinematic_expansion.py` of the branch `paper-runs`.

## Translation and deformation tests

`control.adv` and the same hooks, with meshes of unit cubes and a plane mesh
for the disc. The runs of Appendix A.2, with their mesh scripts and records,
are in the branch `paper-runs` (`cases/review/kinematic/adv`). They are not
staged by a script here.
