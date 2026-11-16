!==============================================================================!
  subroutine Export_All_Surfaces(Exporter)
!------------------------------------------------------------------------------!
!   Export all extracted surfaces to VTU files                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surface_Exporter_Type), intent(inout) :: Exporter
!-----------------------------------[Locals]-----------------------------------!
  integer :: i
!==============================================================================!

  do i = 1, Exporter % n_cuts
    if (Exporter % cuts(i) % active) then
      call Exporter % Export_Surface_Vtu(i)
    end if
  end do

  end subroutine