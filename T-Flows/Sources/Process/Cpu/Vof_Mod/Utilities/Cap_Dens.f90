!==============================================================================!
  real function Cap_Dens(Vof, Flow, c)
!------------------------------------------------------------------------------!
!   Density times heat capacity of cell c, as the energy equation sees it.     !
!                                                                              !
!   Without a mass transfer model the properties follow the colour function    !
!   (Flow % density and Flow % capacity, blended in the cut cells).  With one, !
!   the energy equation takes the sharp values of the phase that the cell      !
!   belongs to: phase 1 where the colour function is above one half, phase 0   !
!   where it is below, and the blended values at exactly one half.             !
!                                                                              !
!   Everything that transports or stores enthalpy has to use this one rule:    !
!   advection, the upwind part of the matrix and inertia.  Taking the blended  !
!   values in one place and the sharp ones in another moves enthalpy with one  !
!   set of properties and stores it with another, which at 1 MPa water and     !
!   steam (rho*cp of 3.9e6 and 1.4e4) is a factor of 280.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),  target :: Vof
  type(Field_Type), target :: Flow
  integer,      intent(in) :: c
!==============================================================================!

  Cap_Dens = Flow % capacity(c) * Flow % density(c)

  if(Flow % mass_transfer_model .ne. NO_MASS_TRANSFER) then
    if(Vof % fun % n(c) > 0.5) then
      Cap_Dens = Vof % phase_capa(1) * Vof % phase_dens(1)
    else if(Vof % fun % n(c) < 0.5) then
      Cap_Dens = Vof % phase_capa(0) * Vof % phase_dens(0)
    end if
  end if

  end function
