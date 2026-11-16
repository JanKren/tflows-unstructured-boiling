!==============================================================================!
  subroutine End_Parallel_Abort(Global)
!------------------------------------------------------------------------------!
!>  Forcibly terminates ALL ranks from a single caller via MPI_Abort, with a
!>  non-zero exit code.  Used on fatal errors that may be raised by only some
!>  ranks (e.g. a locally corrupt front), where the collective End_Parallel
!>  (Mpi_Barrier + Mpi_Finalize) would deadlock the ranks that never error.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Comm_Type) :: Global  !! global communication object
!-----------------------------------[Locals]-----------------------------------!
  integer :: error
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Global)
!==============================================================================!

  call Mpi_Abort(MPI_COMM_WORLD, 1, error)

  end subroutine
