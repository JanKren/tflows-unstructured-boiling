!==============================================================================!
  subroutine Setup(Agg, ni, A)
!------------------------------------------------------------------------------!
!>  Builds the whole aggregation hierarchy from the system matrix A, using only
!>  its first 'ni' rows and the couplings among them.
!------------------------------------------------------------------------------!
!   Only the rows this rank owns are taken, and within them only the columns     !
!   that are also owned.  Entries reaching into the buffer are dropped, which    !
!   makes the hierarchy, and therefore the whole cycle, rank local.  See the      !
!   note on parallel behaviour in Agglomeration_Mod.f90 for why that is the       !
!   right approximation to make here rather than a compromise.                    !
!                                                                                !
!   The hierarchy is rebuilt on every call.  The pressure matrix changes with     !
!   the density field at every time step, and a stale hierarchy would be a       !
!   preconditioner for a matrix that no longer exists - still convergent, since  !
!   any symmetric positive definite operator is, but progressively less useful.   !
!   Rebuilding costs one pass over the matrix per level, geometrically           !
!   decreasing, which is a couple of matrix-vector products in total against     !
!   the tens of iterations the cycle saves.                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Agglomeration_Type), intent(inout) :: Agg  !! parent class
  integer,                   intent(in)    :: ni   !! unknowns owned here
  type(Matrix_Type),         intent(in)    :: A    !! system matrix
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, k, l, nnz, n_coarse
!==============================================================================!

  !--------------------------------------------------------!
  !                                                        !
  !   Level one is the owned part of the system matrix     !
  !                                                        !
  !--------------------------------------------------------!
  if(allocated(Agg % lev)) deallocate(Agg % lev)
  allocate(Agg % lev(AGG_MAX_LEVELS))

  ! Count the entries that survive the restriction to owned columns
  nnz = 0
  do i = 1, ni
    do k = A % row(i), A % row(i+1) - 1
      if(A % col(k) .le. ni) nnz = nnz + 1
    end do
  end do

  Agg % lev(1) % n = ni
  allocate(Agg % lev(1) % row  (ni + 1))
  allocate(Agg % lev(1) % col  (nnz))
  allocate(Agg % lev(1) % val  (nnz))
  allocate(Agg % lev(1) % dia  (ni))
  allocate(Agg % lev(1) % d_inv(ni))
  allocate(Agg % lev(1) % agg  (ni))
  allocate(Agg % lev(1) % x    (ni))
  allocate(Agg % lev(1) % b    (ni))
  allocate(Agg % lev(1) % r    (ni))

  nnz = 0
  do i = 1, ni
    Agg % lev(1) % row(i) = nnz + 1
    do k = A % row(i), A % row(i+1) - 1
      j = A % col(k)
      if(j .le. ni) then
        nnz = nnz + 1
        Agg % lev(1) % col(nnz) = j
        Agg % lev(1) % val(nnz) = A % val(k)
        if(j .eq. i) Agg % lev(1) % dia(i) = nnz
      end if
    end do
  end do
  Agg % lev(1) % row(ni + 1) = nnz + 1

  do i = 1, ni
    Assert(Agg % lev(1) % dia(i) > 0)
    if(abs(Agg % lev(1) % val(Agg % lev(1) % dia(i))) > TINY) then
      Agg % lev(1) % d_inv(i) = 1.0 / Agg % lev(1) % val(Agg % lev(1) % dia(i))
    else
      Agg % lev(1) % d_inv(i) = 0.0
    end if
  end do

  !----------------------------------------------------------!
  !                                                          !
  !   Coarsen until the levels stop getting usefully smaller  !
  !                                                          !
  !----------------------------------------------------------!
  Agg % n_levels = 1

  do l = 1, AGG_MAX_LEVELS - 1

    if(Agg % lev(l) % n .le. AGG_MIN_UNKNOWNS) exit

    call Agg % Aggregate(l, n_coarse)

    !----------------------------------------------------------------!
    !   Refuse a level that buys nothing.  Pairwise matching can      !
    !   fail to pair anything if the strong couplings all point at    !
    !   cells already taken, and a level the same size as its parent  !
    !   would make the cycle recurse without ever reaching the        !
    !   bottom.                                                       !
    !----------------------------------------------------------------!
    if(n_coarse .ge. Agg % lev(l) % n) exit

    call Agg % Coarse_Level(l, n_coarse)

    Agg % n_levels = l + 1

  end do

  end subroutine
