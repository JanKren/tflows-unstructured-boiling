# Static droplet

A droplet at rest, in equilibrium between surface tension and pressure. The
velocities it develops (spurious currents) and the error of the pressure jump
measure the balance of the surface-tension and pressure forces. Appendix D of
the paper.

| folder | what | paper | source state |
|---|---|---|---|
| `submitted/` | the runs of the table of the submitted version, R/h = 12.8 and 25.6 on hexahedral and polyhedral meshes (`hex64`, `hex128`, `poly64`, `poly128`), with their hooks (`User_Mod/`); `sweeps/` holds the control files of the smoothing and Rhie–Chow variants | Appendix D (submitted) | not identified: the binary of 12 August 2026 reads the control keywords of `5061d3d8bf` |
| `revision/` | the runs of the revision: `stage/` (hooks, control files, build and job scripts), `traces/` and `traces_pilot/` (the spurious-current records), `smoke/` (short local checks), `analysis/` (comparison with the submitted traces) | Appendix D; response on the static drop | `2c095ac070` (`stage/COMMIT`), verified |

The long runs of the revision at R/h = 12.8 over several viscous times are not
included: they were on the cluster when this branch was made.
