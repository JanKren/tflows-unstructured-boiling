!==============================================================================!
  function Should_Save_Now(Exporter, current_time_step) result(should_save)
!------------------------------------------------------------------------------!
!   Check if surfaces should be saved at the current time step               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surface_Exporter_Type), intent(inout) :: Exporter
  integer,                     intent(in)    :: current_time_step
  logical                                     :: should_save
!==============================================================================!

  should_save = .false.

  ! Check if enough time steps have passed since last save
  if (current_time_step - Exporter % last_saved_step >= Exporter % save_interval) then
    should_save = .true.
    Exporter % last_saved_step = current_time_step
  end if

  end function