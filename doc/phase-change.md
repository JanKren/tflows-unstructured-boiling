# The phase-change model in the source

Section numbers are those of the revised manuscript (draft of 16 November
2026). Paths are relative to `Sources/Process/Cpu` unless they start with
`Shared/`.

## One time step

`Main_Pro.f90` does, in this order:

1. **Volume fraction** (`Vof_Mod/Core/Main_Vof.f90`, `Compute_Vof.f90`). The
   assembly (`Discretize.f90`) first computes the mass transfer rate
   (`Vof_Mod/Utilities/Mass_Transfer_Estimate.f90`): the temperature
   gradients of both phases with the saturation temperature at the front,
   and from them the rate of every cell that holds a front element. The
   volume fraction is then advanced with the face fluxes of the last
   pressure solution and the source of the phase change, and the front is
   rebuilt.
2. **Outer iterations**, each one: momentum, pressure with the volume source
   of the phase change, velocity correction, and energy
   (`Process_Mod/Compute_Energy.f90`). The energy step computes the mass
   transfer rate again and imposes the saturation temperature on the faces
   that the front intersects.

## Where each part is

| part of the method | paper | source |
|---|---|---|
| Volume-fraction equation, assembly and solution | §2.1 | `Vof_Mod/Core/Compute_Vof.f90`, `Vof_Mod/Core/Discretize.f90` |
| CICSAM and STACS face values: the weight of the acceptor cell, the blend of the compressive and the bounded value on the angle between the interface normal and the face | §2.1, App. A.1 | `Vof_Mod/Core/Predict_Beta.f90` |
| Boundedness corrector of the face values; clipping to [0, 1] after the sweeps | §2.1, App. A.1 | `Vof_Mod/Core/Correct_Beta.f90`, `Vof_Mod/Core/Compute_Vof.f90` |
| Courant number of the volume-fraction step and sub-stepping above `MAX_COURANT_VOF` | App. A.2 | `Vof_Mod/Core/Max_Courant_Number.f90`, `Vof_Mod/Core/Compute_Vof.f90` |
| Source of the volume-fraction equation from the phase change | §2.1 | `Vof_Mod/Utilities/Mass_Transfer_Vof_Source.f90` (see `KNOWN_ISSUES.md`, item 1) |
| The front: the surface of volume fraction 1/2, its elements, and where it cuts the connection of two cell centres | §2.1 | `Front_Mod/` (`Place_Front_At_Value.f90`, `Find_Front_Elements.f90`, `Mark_Cells_And_Faces.f90`), `Shared/Isoap_Mod` |
| Mixture properties | §2.1 | `Vof_Mod/Utilities/Update_Physical_Properties.f90` |
| Volume source of the phase change in the pressure equation | §2.2 | `Vof_Mod/Utilities/Mass_Transfer_Pressure_Source.f90` |
| Least-squares gradient with the arm that crosses the interface cut at the front and set to the saturation temperature | §2.3, App. A.3 | `Vof_Mod/Utilities/Calculate_Grad_Matrix_With_Front.f90`, `Grad_Component_No_Refresh_With_Front.f90`, `Grad_Variable_With_Front.f90` |
| Extrapolation of the gradient along the interface normal into the cells that hold the front (six calls per step: three components, two phases) | §2.3 | `Vof_Mod/Utilities/Extrapolate_Normal_To_Front.f90` |
| Mass transfer rate from the jump of the conductive heat flux over the area of the front element (Stefan condition) | §2.3 | `Vof_Mod/Utilities/Mass_Transfer_Estimate.f90` |
| Energy equation in the advective form; on a face that the front intersects, the conductive flux to the neighbour is replaced by the flux to the saturation temperature over the distance to the front, on the matrix diagonal | §2.4 | `Process_Mod/Compute_Energy.f90` (coefficients `Vof % a12`, `Vof % a21` from `Mass_Transfer_Estimate.f90`) |
| Curvature, its smoothing, continuum surface force | §2.5.2 | `Vof_Mod/Curvature/Curvature_Csf.f90`, `Smooth_Curvature.f90`, `Smooth_Scalar.f90`, `Vof_Mod/Utilities/Surface_Tension_Force.f90` |
| Rhie–Chow interpolation with the corrections of Choi and Gu | §2.5.3 | `Process_Mod/Rhie_And_Chow.f90` |
| Outer iterations | §2.6 | `Process_Mod/` (`Compute_Momentum.f90`, `Compute_Pressure.f90`, `Compute_Energy.f90`), `Main_Pro.f90` |
| Extended gradient stencil with inverse-distance weights | App. B.3 | `Vof_Mod/Utilities/Find_Extended_Stencil.f90` |

The interface-modified gradient that the solver saves in the result files
(`Results_Mod/Save_Vtu_Fields.f90`, field `Temperature Gradients`) is the
gradient before the extrapolation.

## Control keywords of the model

Defaults are those of the source. Two-phase keywords take one value per
phase; phase 0 is the vapour and phase 1 the liquid in the cases here.

| keyword | default | meaning |
|---|---|---|
| `INTERFACE_TRACKING` | `no` | switches the VOF model on |
| `TRACK_FRONT` | `no` | builds the front; needed by the phase-change model |
| `MASS_TRANSFER_MODEL` | `none` | `TEMPERATURE_GRADIENTS` selects the model of the paper |
| `LATENT_HEAT` | 1.0 | J/kg |
| `SATURATION_TEMPERATURE` | 100 | in the unit of the temperature field |
| `PHASE_DENSITIES`, `PHASE_VISCOSITIES`, `PHASE_CONDUCTIVITIES`, `PHASE_CAPACITIES` | | one value per phase |
| `SURFACE_TENSION` | 0.0 | N/m |
| `ADVECTION_SCHEME_FOR_VOF` | `upwind` | `cicsam` or `stacs` |
| `MAX_COURANT_VOF` | 0.25 | above it the volume-fraction step is sub-stepped |
| `MAX_SUBSTEP_CYCLES_VOF` | 100 | cap of the sub-steps |
| `MAX_CORRECTION_CYCLES_BETA_VOF` | 2 | sweeps of the boundedness corrector |
| `MAX_ITERATIONS_FOR_VOF_SOLVER`, `TOLERANCE_FOR_VOF_SOLVER` | 6, 1e-6 | linear solver of the volume fraction |
| `MAX_SMOOTHING_CYCLES_CURVATURE_VOF` | 2 | smoothing cycles of the volume fraction before the curvature (the spelling `MAX_SMOOTH_CYCLES_CURVATURE_VOF` is accepted too) |
| `GRADIENT_METHOD_FOR_VOF` | `least_squares` | `gauss_theorem` in the cases of the paper |
| `GU_CORRECTION`, `CHOI_CORRECTION` | | corrections of the Rhie–Chow interpolation |
| `BLEND_SYSTEM_MATRICES` | `yes` | keeps the upwind part of the convection in the matrix; with it the energy equation is solved in the advective form |
| `MIN_SIMPLE_ITERATIONS`, `MAX_SIMPLE_ITERATIONS` | 3, 12 | outer iterations per time step |

Switches of the analysis in Section 4 and Appendices B and C. Their defaults
leave the solver unchanged:

| keyword | default | effect |
|---|---|---|
| `UNIFORM_MASS_TRANSFER` | `no` | gives every front element the same mass flux per unit area, with the total of the domain unchanged |
| `VOF_BLEND_GAMMA` | −1 (off) | replaces the blend weight of the face value by a constant: 1 is the compressive limit (Hyper-C), 0 the bounded limit (ULTIMATE-QUICKEST) |
| `ISOTROPIC_VOF_BETA` | −1 (off) | replaces the weight of the acceptor cell by one constant on every interface face |
| `EXTENDED_GRADIENT_STENCIL` | `no` | adds the neighbours over edges and corners to the gradient stencil |
| `GRADIENT_STENCIL_WEIGHT_EXPONENT` | 0.0 | exponent of the inverse-distance weight of the gradient stencil |

## What the solver prints

- `# VOF_CORR:` at every volume-fraction step: the sweeps of the boundedness
  corrector, and the number of cells left below 0 and above 1 after them.
- The hooks of the cases write their own records (`bench-data.dat`,
  `caf-terms.dat`, `kin-flux.dat`); each case's `README.md` gives the columns.

## The hooks of a case

A case brings its initial condition and its output as Fortran files in
`User_Mod/`, compiled into `Process` (see `doc/building.md`). The hooks that
the cases here use:

| hook | called | used for |
|---|---|---|
| `Initialize_Variables.f90` | once | interface and temperature at the start |
| `Beginning_Of_Compute_Vof.f90` | before each volume-fraction step | prescribed face fluxes of the kinematic tests |
| `End_Of_Compute_Vof.f90`, `End_Of_Compute_Energy.f90`, `End_Of_Compute_Pressure.f90`, `Beginning_Of_Iteration.f90` | in each step | terms of the mass and energy budgets |
| `End_Of_Time_Step.f90` | after each step | the records |
| `Types.f90` | | variables the hooks share |
