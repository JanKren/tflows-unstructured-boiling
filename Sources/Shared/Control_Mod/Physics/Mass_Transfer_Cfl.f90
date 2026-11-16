!==============================================================================!
  subroutine Mass_Transfer_Cfl(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the maximum allowed Stefan-source dilatation Courant number used by
!>  the mass-transfer crash-guard clamp (Mass_Transfer_Estimate); falls back to
!>  0.25 when the keyword is absent.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control
  real,   intent(out) :: val
  logical, optional   :: verbose
!==============================================================================!

  call Control % Read_Real_Item('MASS_TRANSFER_CFL', 0.25, val, verbose)

  end subroutine
