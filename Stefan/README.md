# Stefan problem

A planar interface moves away from a heated wall while the vapour conducts
heat to it. Water and steam at 1 bar. Sections 3.1 and 3.3 of the paper.

| folder | what | paper | source state |
|---|---|---|---|
| `setup/` | control file, mesh definition (`stefan.dom`, `box_parameters.ini`), hooks with the mass and energy budget terms, and `run_stefan.py`, which stages and runs one case | | |
| `submitted/` | the case as published with the submitted version: set-up, interface position and the exact solution | Section 3.1 (submitted) | builds of January to March 2026, not identified |
| `budgets/` | the budget runs of the submitted version: `stefan` (5 ms), `stefan_dt2` (2.5 ms), `stefan_tol`, `stefan_mean`, `stefan_vol`; `records/` holds the extracted budget record; `build_case.sh` built all budget cases of Stefan, Sucking and Scriven | Section 3.3 (submitted) | `bdf9354578`, verified |
| `revision/` | the series of the revision on the energy budget, the mass error and the linear solvers; the job lists are `series*_jobs.txt`, each line one run | Section 3.3, Table `tab:stefan_energy` | see below |

Runs in `revision/`:

| runs | solver | |
|---|---|---|
| `R_*`, `base_review`, `A_review*` | `2c095ac070` | verified |
| `base_paper`, `P_*`, `Pfix_*`, `A_paper*`, `paper_*` | `6adb80f599` | verified |
| `C_*`, `R_sign`, `R_urf1_sign` | `2c095ac070` with `T-Flows/patches/vof_source_sign.patch` | stated |
| `P_sign` | `6adb80f599` with the same patch | stated |

`caf-terms.dat` of each run has one row per time step; its columns are listed
at the top of `setup/User_Mod/End_Of_Time_Step.f90`. `A_*` runs impose a
uniform temperature (file `caf_mode.ini`) to test the convective term.
