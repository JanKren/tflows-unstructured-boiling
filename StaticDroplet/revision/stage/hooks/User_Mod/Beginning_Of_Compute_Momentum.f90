!==============================================================================!
  subroutine User_Mod_Beginning_Of_Compute_Momentum(Flow, Turb, Vof, Sol)
!------------------------------------------------------------------------------!
!   Exact-curvature test of the balanced-force property (Francois et al.).     !
!                                                                              !
!   EXACT_CURVATURE_TEST 0   off (default): the solver's curvature is used     !
!   EXACT_CURVATURE_TEST 1   curvature = 2/R0 in every cell, so the surface    !
!                            force is sigma kappa grad(fun), a pure gradient   !
!   EXACT_CURVATURE_TEST 2   curvature = 2/R0 on the solver's own support      !
!                            (where it computed a non-zero value), elsewhere 0 !
!                                                                              !
!   The curvature is computed once per time step in Main_Vof and consumed in   !
!   Surface_Tension_Force (from Compute_Momentum) and in Rhie_And_Chow (from   !
!   Compute_Pressure).  This hook runs before both, in every outer iteration.  !
!   The sign is taken from the solver's own field.                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
!------------------------------[Local parameters]------------------------------!
  real, parameter :: R0 = 0.4
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer, save            :: mode  = 0
  logical, save            :: first = .true.
  integer                  :: c
  real                     :: s, kap
!==============================================================================!

  if(first) then
    call Control % Read_Int_Item('EXACT_CURVATURE_TEST', 0, mode, .false.)
    first = .false.
    if(First_Proc()) print '(a,i2)', ' # EXACT_CURVATURE_TEST mode: ', mode
  end if
  if(mode .eq. 0) return

  Grid => Flow % pnt_grid

  s = 0.0
  do c = Cells_In_Domain()
    s = s + Vof % curv(c)
  end do
  call Global % Sum_Real(s)
  kap = sign(2.0 / R0, s)

  if(mode .eq. 1) then
    Vof % curv(:) = kap
  else
    do c = lbound(Vof % curv, 1), ubound(Vof % curv, 1)
      if(abs(Vof % curv(c)) > FEMTO) Vof % curv(c) = kap
    end do
  end if

  end subroutine
