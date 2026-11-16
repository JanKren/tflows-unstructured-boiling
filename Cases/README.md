# Cases

Each folder holds one case as it was run for the paper: the `control` file,
the hooks in `User_Mod/`, the mesh definition, and scripts that stage and run
it from this tree.

```bash
Cases/build_tools.sh              # once: Generate, Divide, Convert
Cases/build_case.sh Cases/<case>  # Process with the hooks of that case
```

`build_case.sh` leaves the binary in the case folder as `Process`. Build again
when you change case: the hooks are compiled into the solver
(`doc/building.md`).

| folder | case | mesh | cost |
|---|---|---|---|
| `Stefan` | Stefan problem; mass and energy budgets | `Generate`, 100 × 5 × 5 cells | 4 min, 1 core |
| `Kinematic` | radial expansion with prescribed fluxes; translation and deformation tests | `Generate`, 74³ cells | 4 min, 8 cores |

The Sucking, Scriven, static-droplet and annular cases are not staged here
yet. Their set-ups and records are in the branch `paper-runs`
(`cases/submitted`, `cases/review`).

`default_hooks/` is a copy of the empty hooks of the solver, which
`build_case.sh` restores before every build.

Records that the checks compare with are in each case's `expected/`. They are
the records of the runs behind the paper's revision, made with the full
development tree; "equal" means every printed digit.
