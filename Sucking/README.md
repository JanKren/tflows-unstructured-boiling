# Sucking problem

A planar interface moves into superheated liquid, which is pushed out of the
domain. Water and steam at 1 bar. Sections 3.1 and 3.3 of the paper.

| folder | what | paper | source state |
|---|---|---|---|
| `setup/` | control file, mesh definition, initial temperature profile and hooks with the budget terms, from which the runs of the revision were made | | |
| `submitted/` | the case as published with the submitted version: set-up, interface position and velocity, exact solution | Section 3.1 (submitted) | builds of January to March 2026, not identified |
| `budgets/` | the budget run of the submitted version (`sucking`); it was built with `../Stefan/budgets/build_case.sh` | Section 3.3 (submitted) | `bdf9354578`, verified |
| `revision/` | `S_base` (the production settings) and `S_tight` (linear solvers converged to 1e-10) | Section 3.3; response on the linear solvers | `2c095ac070`, verified |
