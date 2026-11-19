!==============================================================================!
  subroutine Mass_Transfer_Vof_Source(Vof, b)
!------------------------------------------------------------------------------!
!   Phase-change source of the volume-fraction equation                        !
!                                                                              !
!   Called from Vof_Mod/Core/Discretize                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  real, intent(inout)     :: b(Vof % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  integer                   :: c, e, vapour, liquid
!==============================================================================!

  ! Take aliases
  Grid => Vof % pnt_grid
  Flow => Vof % pnt_flow

  ! If not a problem with mass transfer, get out of here
  if(Flow % mass_transfer_model .eq. NO_MASS_TRANSFER) return

  ! Distinguish between liquid and vapor
  call Vof % Get_Vapour_And_Liquid_Phase(vapour, liquid)

  !-------------------!
  !   Volume source   !
  !-------------------!

  ! The colour function is the volume fraction of phase 1, the liquid in
  ! every case with phase change (the denser phase listed second, so that
  ! liquid - vapour = 1).  Evaporation (m_dot > 0) consumes liquid: the
  ! equation of the liquid fraction carries -m_dot / rho_l, beside the
  ! expansion m_dot (1/rho_v - 1/rho_l) that Mass_Transfer_Pressure_Source
  ! puts into the divergence of the face fluxes.  With this sign the mass
  ! in the domain plus the mass through the outflow stays constant.  With
  ! the opposite one, which this routine had until November 2026, each
  ! kilogram evaporated left two kilograms of liquid behind, a mass error
  ! of 2 rho_v / rho_l of the mass expelled and a relative 2 rho_v / rho_l
  ! in the interface velocity.
  !
  ! Were the vapour phase 1, the colour would be the vapour fraction and
  ! its source +m_dot / rho_v, which the expression below does not give.
  do c = Cells_In_Domain_And_Buffers()
    e = Vof % Front % elem_in_cell(c)  ! front element

    if(e .ne. 0) then
      b(c) = b(c)                              &
           - Vof % m_dot(c) * (liquid-vapour)  &
           * (1.0/Vof % phase_dens(liquid))
    end if
  end do

  end subroutine
