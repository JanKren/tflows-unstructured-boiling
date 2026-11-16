!==============================================================================!
  subroutine Coarse_Level(Agg, l, n_coarse)
!------------------------------------------------------------------------------!
!>  Forms the operator of level l+1 as the Galerkin product of the operator of
!>  level l with the piecewise constant prolongation held in Agg % lev(l) % agg.
!------------------------------------------------------------------------------!
!   With a piecewise constant prolongation the triple product P^T A P collapses  !
!   into a sum: the coarse entry (I,J) is the sum of every fine entry (i,j) with  !
!   i in aggregate I and j in aggregate J.  No interpolation weights appear, so   !
!   this is a single pass over the fine matrix.                                   !
!                                                                                !
!   Because the same map is used on both sides, a symmetric fine operator gives   !
!   a symmetric coarse one, which the Conjugate Gradient method above depends on. !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Agglomeration_Type), intent(inout) :: Agg       !! parent class
  integer,                   intent(in)    :: l         !! fine level
  integer,                   intent(in)    :: n_coarse  !! unknowns on l+1
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i, j, k, n, ic, jc, nnz, m
  integer, allocatable :: first(:), next(:)  !! aggregate -> its members
  integer, allocatable :: mark(:)            !! column -> slot in current row
  integer, allocatable :: touched(:)         !! columns used by current row
  integer              :: n_touched
!==============================================================================!

  n = Agg % lev(l) % n

  !----------------------------------------------------!
  !   Chain the members of each aggregate together,     !
  !   so the rows of one coarse row can be walked       !
  !----------------------------------------------------!
  allocate(first(n_coarse));  first(:) = 0
  allocate(next (n));         next (:) = 0
  do i = n, 1, -1                    ! backwards, so chains come out ascending
    ic        = Agg % lev(l) % agg(i)
    next(i)   = first(ic)
    first(ic) = i
  end do

  allocate(mark   (n_coarse));  mark(:) = 0
  allocate(touched(n_coarse))

  !----------------------------------------------!
  !   First pass: how many non-zeros will there   !
  !   be?  The row structure is not known until   !
  !   the aggregates have been walked once.       !
  !----------------------------------------------!
  nnz = 0
  do ic = 1, n_coarse
    n_touched = 0
    i = first(ic)
    do while(i .ne. 0)
      do k = Agg % lev(l) % row(i), Agg % lev(l) % row(i+1) - 1
        jc = Agg % lev(l) % agg(Agg % lev(l) % col(k))
        if(mark(jc) .eq. 0) then
          n_touched          = n_touched + 1
          touched(n_touched) = jc
          mark(jc)           = 1
        end if
      end do
      i = next(i)
    end do
    nnz = nnz + n_touched
    do m = 1, n_touched
      mark(touched(m)) = 0
    end do
  end do

  !-----------------------------------!
  !   Room for the coarse operator    !
  !-----------------------------------!
  Agg % lev(l+1) % n = n_coarse
  allocate(Agg % lev(l+1) % row  (n_coarse + 1))
  allocate(Agg % lev(l+1) % col  (nnz))
  allocate(Agg % lev(l+1) % val  (nnz))
  allocate(Agg % lev(l+1) % dia  (n_coarse))
  allocate(Agg % lev(l+1) % d_inv(n_coarse))
  allocate(Agg % lev(l+1) % agg  (n_coarse))
  allocate(Agg % lev(l+1) % x    (n_coarse))
  allocate(Agg % lev(l+1) % b    (n_coarse))
  allocate(Agg % lev(l+1) % r    (n_coarse))
  Agg % lev(l+1) % dia(:) = 0

  !----------------------------------------!
  !   Second pass: accumulate the entries  !
  !----------------------------------------!
  nnz = 0
  do ic = 1, n_coarse

    Agg % lev(l+1) % row(ic) = nnz + 1
    n_touched = 0

    i = first(ic)
    do while(i .ne. 0)
      do k = Agg % lev(l) % row(i), Agg % lev(l) % row(i+1) - 1
        jc = Agg % lev(l) % agg(Agg % lev(l) % col(k))
        if(mark(jc) .eq. 0) then
          nnz                        = nnz + 1
          n_touched                  = n_touched + 1
          touched(n_touched)         = jc
          mark(jc)                   = nnz
          Agg % lev(l+1) % col(nnz)  = jc
          Agg % lev(l+1) % val(nnz)  = 0.0
        end if
        Agg % lev(l+1) % val(mark(jc)) = Agg % lev(l+1) % val(mark(jc))  &
                                       + Agg % lev(l) % val(k)
      end do
      i = next(i)
    end do

    ! Find this row's diagonal and clear the markers for the next row
    do m = 1, n_touched
      jc = touched(m)
      if(jc .eq. ic) Agg % lev(l+1) % dia(ic) = mark(jc)
      mark(jc) = 0
    end do

  end do
  Agg % lev(l+1) % row(n_coarse + 1) = nnz + 1

  !-------------------------------------------------------------!
  !   The reciprocal of the diagonal, used by the smoother.      !
  !   A vanishing diagonal cannot arise from a sum over an        !
  !   aggregate of an M-matrix, but guard it rather than divide   !
  !   by it, because a stray zero here would poison the cycle     !
  !   silently rather than stop it.                               !
  !-------------------------------------------------------------!
  do ic = 1, n_coarse
    Assert(Agg % lev(l+1) % dia(ic) > 0)
    if(abs(Agg % lev(l+1) % val(Agg % lev(l+1) % dia(ic))) > TINY) then
      Agg % lev(l+1) % d_inv(ic) = 1.0 / Agg % lev(l+1) % val(  &
                                         Agg % lev(l+1) % dia(ic))
    else
      Agg % lev(l+1) % d_inv(ic) = 0.0
    end if
  end do

  deallocate(first, next, mark, touched)

  end subroutine
