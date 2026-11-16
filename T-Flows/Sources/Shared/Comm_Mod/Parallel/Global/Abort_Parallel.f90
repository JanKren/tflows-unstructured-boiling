!==============================================================================!
  subroutine Abort_Parallel(Global)
!------------------------------------------------------------------------------!
!>  Tears the whole job down from a single processor, for use when something
!>  has gone wrong and there is no reason to believe the others will agree.
!------------------------------------------------------------------------------!
!   Why this exists next to End_Parallel                                       !
!                                                                              !
!   End_Parallel is the ORDERLY exit: a barrier, then Mpi_Finalize.  Both are   !
!   collectives, so it only works when every rank arrives - which is the right  !
!   thing at the end of a run, where they all do.                              !
!                                                                              !
!   It is exactly the wrong thing on an error.  A rank that stops on a fault    !
!   the others have not met blocks in the barrier waiting for them, and they    !
!   are still computing and will never come; the job hangs instead of failing,  !
!   and hangs with the diagnosis already printed and nobody reading it.  It     !
!   then has to be killed by a wall-clock limit, which on a cluster means the   !
!   error is reported as a timeout.                                            !
!                                                                              !
!   This was harmless until the grid began to adapt.  The errors T-Flows had    !
!   before - a missing keyword, an unreadable file - are met by every rank at   !
!   the same instant, so the collectives are satisfied and the exit is clean.   !
!   Refinement is the first thing here where a rank can reach a fault ALONE,    !
!   because each one judges its own cells; and it duly did, twice in one day.   !
!                                                                              !
!   Mpi_Abort is the right primitive for that because it is NOT collective:     !
!   one rank calling it takes the job down.                                     !
!                                                                              !
!   The flush matters.  Mpi_Abort does not unwind anything, so a message still  !
!   sitting in a buffer is lost, and the whole point of aborting here is to be  !
!   able to read what went wrong.                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Global  !! global communication object
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!==============================================================================!

  Unused(Global)

  flush(output_unit)

  call Mpi_Abort(MPI_COMM_WORLD, 1, error)

  end subroutine
