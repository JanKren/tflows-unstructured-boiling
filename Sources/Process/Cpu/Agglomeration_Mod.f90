#include "../../Shared/Assert.h90"
#include "../../Shared/Unused.h90"

!==============================================================================!
  module Agglomeration_Mod
!------------------------------------------------------------------------------!
!>  An aggregation multigrid V-cycle, used as a preconditioner for the native
!>  Conjugate Gradient solver.  It exists to remove the mesh dependence of the
!>  iteration count, which is what the Incomplete Cholesky preconditioner
!>  cannot do and what dominates the cost of the pressure solve.
!------------------------------------------------------------------------------!
!   Why this, and why not the refinement hierarchy                             !
!                                                                              !
!   Measured on the rising bubble case, the average residual contraction of     !
!   the ICC preconditioned CG is                                                !
!                                                                              !
!      8192 cells (h=1/16)   rho = 0.718     41.7 iterations per 1e6           !
!     16000 cells (h=1/20)   rho = 0.757     49.5                              !
!     65536 cells (h=1/32)   rho = 0.836     77.2                              !
!                                                                              !
!   which is the textbook O(1/h) of any single level preconditioner.  The      !
!   adapted grid, 15332 cells carrying two refinement levels, measured         !
!   rho = 0.752, against 0.757 for a UNIFORM grid of 16000 cells.  The two are  !
!   the same to within the scatter, so the cells of two sizes meeting cost      !
!   nothing in conditioning: an adapted grid is exactly as hard as a uniform    !
!   grid with as many cells, no harder.                                        !
!                                                                              !
!   That is the reason this module ignores Grid % cell_level and the           !
!   parent/child links entirely.  There is no hierarchy-specific difficulty     !
!   for a hierarchy-aware coarsening to exploit; the difficulty is the plain    !
!   mesh dependence that multigrid of any flavour addresses.  Aggregating on    !
!   the matrix graph also keeps the preconditioner usable on grids that were    !
!   never refined, which is most of them.                                      !
!------------------------------------------------------------------------------!
!   What it achieves, and what it does not                                     !
!                                                                              !
!   Same case, same measure, this cycle in place of Incomplete Cholesky:        !
!                                                                              !
!                          iterations per 1e6 reduction                        !
!                            icc      amg     gain                              !
!      8192 cells           41.7     16.8     2.5                              !
!     16000 cells           49.5     18.4     2.7                              !
!     65536 cells           76.3     24.7     3.1                              !
!     15332, adapted        48.4     19.1     2.5                              !
!                                                                              !
!   The gain grows with the grid, which is the expected shape: the cycle is     !
!   removing a mesh dependence, so it wins by more the more mesh there is.      !
!                                                                              !
!   It does NOT remove that dependence outright - 16.8 to 24.7 over an eight    !
!   fold increase in cells is flatter than 41.7 to 76.3, but it is not flat.    !
!   The prolongation here is piecewise constant, and unsmoothed aggregation is  !
!   known to leave a slow growth behind.  Smoothing the prolongator, or a K-    !
!   cycle in place of the V-cycle, is where the rest of it is, and neither is   !
!   done here.                                                                 !
!------------------------------------------------------------------------------!
!   Parallel behaviour                                                         !
!                                                                              !
!   Each rank aggregates only the rows it owns and drops couplings to cells    !
!   held by other ranks, so the V-cycle is a rank-local operator.  This is      !
!   deliberately the SAME approximation the existing Incomplete Cholesky        !
!   preconditioner already makes (its substitutions skip buffer entries), so    !
!   the parallel semantics of the solver do not change: the preconditioner is   !
!   block Jacobi over ranks either way, the answer is unaffected, and only the  !
!   iteration count depends on the decomposition.                              !
!                                                                              !
!   The cost of that choice is real and worth stating.  On the adapted grid,   !
!   iterations per 1e6 reduction:                                              !
!                                                                              !
!               1 rank    2 ranks   4 ranks                                    !
!     icc        48.4      51.8      58.0                                      !
!     amg        19.1      27.6      36.0                                      !
!                                                                              !
!   Both degrade, but the cycle degrades faster, for a reason particular to    !
!   multigrid: its whole business is moving error across the domain in one     !
!   sweep, and a rank boundary the coarse levels cannot see is a wall the      !
!   error has to diffuse through instead.  Incomplete Cholesky was never doing !
!   that, so it has less to lose.  The cycle is still ahead at four ranks, by  !
!   1.6 rather than 2.5, but the trend says the advantage keeps narrowing.     !
!   Recovering it needs coarse levels that span ranks, which is a parallel     !
!   coarsening problem and is not attempted here.                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Matrix_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !> Coarsen no further than this, and never build more levels than this.
  integer, parameter :: AGG_MIN_UNKNOWNS =  40
  integer, parameter :: AGG_MAX_LEVELS   =  24

  !> Sweeps of symmetric Gauss-Seidel used to finish off the coarsest level.
  integer, parameter :: AGG_COARSE_SWEEPS = 12

  !--------------------------------------------!
  !   One level of the aggregation hierarchy   !
  !--------------------------------------------!
  type Agg_Level_Type

    integer :: n = 0    !! unknowns on this level

    ! The level's operator, in the same CSR layout Matrix_Mod uses, except
    ! that the diagonal is located through dia(:) rather than assumed first.
    integer, allocatable :: row(:)
    integer, allocatable :: col(:)
    integer, allocatable :: dia(:)
    real,    allocatable :: val(:)
    real,    allocatable :: d_inv(:)  !! reciprocal of the diagonal

    ! Which aggregate of the NEXT (coarser) level each unknown belongs to.
    ! This is the whole of the prolongation: it is piecewise constant, so
    ! restriction is a sum over the aggregate and prolongation a broadcast.
    integer, allocatable :: agg(:)

    ! Work space, so a cycle allocates nothing
    real,    allocatable :: x(:), b(:), r(:)

  end type

  !----------------------------------!
  !   The hierarchy as a whole       !
  !----------------------------------!
  type Agglomeration_Type

    integer :: n_levels = 0
    type(Agg_Level_Type), allocatable :: lev(:)

    contains
      procedure :: Setup      !! builds the hierarchy from a system matrix
      procedure :: V_Cycle    !! applies one V-cycle, the preconditioner itself
      procedure, private :: Aggregate     !! pairwise matching on one level
      procedure, private :: Coarse_Level  !! Galerkin product for one level
      procedure, private :: Smooth        !! Gauss-Seidel sweep

  end type

  !---------------------------------------------------------------------!
  !   The hierarchy is rebuilt from scratch at the start of every        !
  !   solve, because the pressure matrix changes with the density        !
  !   field every time step.  Since setup and use are therefore never    !
  !   separated by another solve, one instance is enough and it can be   !
  !   held here rather than threaded through the solver's interfaces.    !
  !---------------------------------------------------------------------!
  type(Agglomeration_Type), target :: Agglomeration

  contains

#   include "Agglomeration_Mod/Setup.f90"
#   include "Agglomeration_Mod/Aggregate.f90"
#   include "Agglomeration_Mod/Coarse_Level.f90"
#   include "Agglomeration_Mod/Smooth.f90"
#   include "Agglomeration_Mod/V_Cycle.f90"

  end module
