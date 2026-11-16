!==============================================================================!
  subroutine Surface_Tension_Force(Vof, i)
!------------------------------------------------------------------------------!
!   Computes surface tension forces for momentum                               !
!                                                                              !
!   Recall how momentum equations look in Compute_Momentum procedure:          !
!                                                                              !
!     /             /              /               /             /             !
!    |     du      |              |               |             |              !
!    | rho -- dV + | rho u u dS = | mu DIV u dS - | GRAD p dV + | f dV         !
!    |     dt      |              |               |             |              !
!   /             /              /               /             /               !
!                                                                              !
!   Here, f will be the surface tension force.  If we use material             !
!   derivative for the sake of shortness, the equation will read:              !
!                                                                              !
!     /             /               /             /                            !
!    |     Du      |               |             |                             !
!    | rho -- dV = | mu DIV u dS - | GRAD p dV + | f dV                        !
!    |     Dt      |               |             |                             !
!   /             /               /             /                              !
!                                                                              !
!   and with definition of surface tension force:                              !
!                                                                              !
!     /             /               /                        /                 !
!    |     Du      |               |                         |                 !
!    | rho -- dV = | mu DIV u dS - | GRAD p dV + sigma kappa | GRAD c dV       !
!    |     Dt      |               |                         |                 !
!   /             /               /                         /                  !
!                                                                              !
!   In order to balance pressure with surface tension force, last two terms    !
!   should always be treated in exactly the same manner!  One could even       !
!   write the equation like this:                                              !
!                                                                              !
!     /             /               /                                          !
!    |     Du      |               |                                           !
!    | rho -- dV = | mu DIV u dS - | GRAD ( p + sigma kappa c ) dV             !
!    |     Dt      |               |                                           !
!   /             /               /       |                   |                !
!                                         |<---- together --->|                !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  integer, intent(in)     :: i
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: col
  real, contiguous,  pointer :: surf_fx(:), surf_fy(:), surf_fz(:)
  integer                    :: c
!==============================================================================!

  ! Get out of here if surface tension is neglected
  if(Vof % surface_tension < TINY) return

  ! Take aliases
  Flow    => Vof % pnt_flow
  Grid    => Vof % pnt_grid
  surf_fx => Vof % surf_fx
  surf_fy => Vof % surf_fy
  surf_fz => Vof % surf_fz
  ! The force acts on the sharp colour field.  With the smoothed field
  ! spurious currents grew on fine hexahedral meshes until the run failed
  ! (static Scriven bubble, h = 3.125 um, real water viscosity: NaN at
  ! 48 us); on the sharp field the Courant number stays flat.
  ! Rhie_And_Chow.f90 must use the same field.  Curvature still comes from
  ! the smoothed field.
  col => Vof % fun

  ! Units here:
  ! N/m * 1/m * 1/m = N / m^3
  select case(i)
    case(1)
      do c = Cells_In_Domain_And_Buffers()
        surf_fx(c) = Vof % surf_tens(c)     &
                   * Vof % curv(c)          &
                   * col % x(c)
      end do
    case(2)
      do c = Cells_In_Domain_And_Buffers()
        surf_fy(c) = Vof % surf_tens(c)     &
                   * Vof % curv(c)          &
                   * col % y(c)
      end do
    case(3)
      do c = Cells_In_Domain_And_Buffers()
        surf_fz(c) = Vof % surf_tens(c)     &
                   * Vof % curv(c)          &
                   * col % z(c)
      end do

  end select

  end subroutine
