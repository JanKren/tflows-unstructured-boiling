!==============================================================================!
  subroutine Extract_All_Surfaces(Exporter, Grid, Flow, Vof)
!------------------------------------------------------------------------------!
!   Extract all active cutting plane surfaces                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surface_Exporter_Type), intent(inout) :: Exporter
  type(Grid_Type),             intent(in)    :: Grid
  type(Field_Type),            intent(in), optional :: Flow
  type(Vof_Type),              intent(in), optional :: Vof
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  do i = 1, Exporter % n_cuts
    if (Exporter % cuts(i) % active) then
      call Exporter % Extract_Surface(i, Grid, Flow, Vof)
    end if
  end do

  end subroutine