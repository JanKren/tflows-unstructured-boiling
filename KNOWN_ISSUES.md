# Known issues of this version

Each item states what is wrong or untested, how large the effect is where it
was measured, and what the paper's runs did about it.

## 1. Sign of the liquid source term of the volume-fraction equation

`Sources/Process/Cpu/Vof_Mod/Utilities/Mass_Transfer_Vof_Source.f90` adds the
source term of the volume-fraction equation, the mass transfer rate divided by
the liquid density, with the wrong sign. The mass that leaves the domain then differs from the
mass that evaporated by 2 ρ_v/ρ_l of the expelled mass: 1.25 × 10⁻³ for water
and steam at 1 bar, on every record of the Stefan problem. With the sign
corrected the error is −2 × 10⁻⁶. The interface speed changes by 0.12 %.

The correction is `patches/vof_source_sign.patch`. **It is not applied in this
tree**, because the runs of the paper were made without it, except the Stefan
runs that measure it (branch `paper-runs`). Apply it with
`patch -p1 < patches/vof_source_sign.patch` from the root of this branch.

## 2. Boundedness of the volume fraction

After the corrector sweeps of the CICSAM and STACS face values
(`MAX_CORRECTION_CYCLES_BETA_VOF`, default 2), values of the volume fraction
below 0 and above 1 are set to 0 and 1 (`Vof_Mod/Core/Compute_Vof.f90`). The
solver prints the number of corrector sweeps of every step on the lines
`# VOF_CORR:`. A run in which the corrector is exhausted on most steps
completes without a warning; check these lines.

With the blend weight frozen at the compressive limit (`VOF_BLEND_GAMMA 1.0`)
the corrector does not converge on polyhedral meshes: 66 to 99 % of the steps
exhaust it, and raising the cap from 2 to 20 changes nothing.

## 3. Settings of the diagnostic switches

`VOF_BLEND_GAMMA`, `ISOTROPIC_VOF_BETA`, `UNIFORM_MASS_TRANSFER`,
`EXTENDED_GRADIENT_STENCIL` and `GRADIENT_STENCIL_WEIGHT_EXPONENT` exist for
the analysis of Section 4 and Appendices B and C of the paper. Their defaults
leave the solver unchanged. They are not production options:
`ISOTROPIC_VOF_BETA` with a value the scheme does not stay below breaks
boundedness on nearly every step, and the extended stencil is untested on
polyhedral meshes.

## 4. Parts of the solver outside the paper

The turbulence models, particle tracking, porosity, the PETSc interface,
passive scalars and the GPU-related modules in `Sources/Shared` are part of
T-Flows and are in the source as they are in the development branch. The
paper does not use them and nothing in this branch checks them.

Features of the development version that no case of the paper uses are not in
this branch (README.md). A control file that asks for one of them stops the
run with a message (`Sources/Shared/Control_Mod/Basic_Functions/
Stop_On_Removed_Features.f90`):

| keyword | stops the run when |
|---|---|
| `GeometricVOF`, `NUCLEATION`, `ADAPT_GRID`, `VARIABLE_SURFACE_TENSION` | set to anything but `no` |
| `MASS_TRANSFER_MODEL` | set to anything but `none` or `temperature_gradients` |
| `SURFACE_TENSION_MODEL` | set to anything but `csf` |
| `NUMBER_OF_DOMAINS` | larger than 1 |

T-Flows ignores keywords it does not know, so a keyword of a removed feature
that is not in this list is ignored without a word.

## 5. Cases compiled into the solver

The hooks of a case are compiled into `Process`. A `Process` built for one
case runs another case without its initialisation and without an error.
Build with `Cases/build_case.sh`, which restores the default hooks first, and
keep one binary per case.

## 6. Mesh size and block size

A block of `N` nodes in a `.dom` file of `Generate` gives `N − 1` cells in
that direction. The "structured 75³" mesh of the paper has 74³ = 405 224
cells.
