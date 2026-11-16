!==============================================================================!
  subroutine Set_Save_Interval(Exporter, interval)
!------------------------------------------------------------------------------!
!   Set the save interval for surface extraction (in time steps)             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surface_Exporter_Type), intent(inout) :: Exporter
  integer,                     intent(in)    :: interval
!==============================================================================!

  Exporter % save_interval = max(1, interval)
  print *, 'Surface export save interval set to: ', Exporter % save_interval, ' time steps'

  end subroutine