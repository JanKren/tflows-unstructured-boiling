!==============================================================================!
  subroutine Abort_Parallel(Global)
!------------------------------------------------------------------------------!
!>  Nothing to tear down when there is only one processor.
!------------------------------------------------------------------------------!
!   The caller follows this with a stop, which on one rank is the whole job.    !
!   See the parallel version for why the two exits have to be told apart.       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type), intent(in) :: Global  !! global communication object
!==============================================================================!

  Unused(Global)

  flush(output_unit)

  end subroutine
