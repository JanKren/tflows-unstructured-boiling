!==============================================================================!
  subroutine End_Parallel_Abort(Global)
!------------------------------------------------------------------------------!
!>  Sequential counterpart of the parallel MPI_Abort terminator.  Exits the
!>  program immediately with a non-zero status so the scheduler registers the
!>  failure (the plain `stop` in Message % Error exits with status 0).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in) :: Global  !! global comunication object
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
!==============================================================================!

  error stop 1

  end subroutine
