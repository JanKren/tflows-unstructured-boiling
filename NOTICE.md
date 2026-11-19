# Notices

T-Flows is developed by the T-Flows Development Team. This branch is a reduced
copy of commit `0ce2122640` of the branch `unstructured-boiling-release`, which is
commit `2c095ac070` of the branch `T-Flows-jan-review` without the parts
that the paper does not use.

It is distributed under the GNU General Public License, version 3 (`LICENSE`).
The licence text is the file `license` of the T-Flows main branch at commit
`f289644d1d` (30 August 2026). The commit this branch is taken from still
carried the earlier proprietary licence text in its own `license` file; that
file is not included here.

## Authors of the paper

Jan Kren, Bojan Ničeno, Yohei Sato. Paul Scherrer Institute PSI.

## Contributions

- T-Flows was written by **Bojan Ničeno** and the contributors listed in the
  history of the [T-Flows repository](https://github.com/DelNov/T-Flows).
- The phase-change model of the paper, the corrections to the VOF module made
  in the review of October 2026, the cases and their hooks were developed by
  **Jan Kren**, with Bojan Ničeno and Yohei Sato.

## Third-party software

- **METIS 5.1.0** (G. Karypis, University of Minnesota), Apache-2.0. The static
  libraries in `Sources/Libraries/Metis_5.1.0_Linux_64` are linked into
  `Divide` and `Process`.
- **Isoap** (J. López and co-workers), the iso-surface extraction on
  polyhedral cells behind `Sources/Shared/Isoap_Mod`:
  <https://data.mendeley.com/datasets/4rcf98s74c>.
- **MPI.** An MPI library is needed to build `Process` for parallel runs and
  is not included.
- **PETSc** is optional in T-Flows and is not used by the cases of this branch.
