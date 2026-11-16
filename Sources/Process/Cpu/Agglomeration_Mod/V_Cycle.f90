!==============================================================================!
  recursive subroutine V_Cycle(Agg, l)
!------------------------------------------------------------------------------!
!>  Applies one multigrid V-cycle on level 'l', taking the right hand side from
!>  Agg % lev(l) % b and leaving the approximate solution in Agg % lev(l) % x.
!------------------------------------------------------------------------------!
!   The cycle is arranged to be symmetric, so that the operator it represents    !
!   is a valid preconditioner for Conjugate Gradient:                            !
!                                                                                !
!     - the initial guess is zero, which makes the whole cycle linear in b;      !
!     - the pre-smoothing runs forwards and the post-smoothing backwards,        !
!       the two being transposes of one another;                                  !
!     - restriction is the transpose of prolongation, both being the same        !
!       piecewise constant aggregate map;                                         !
!     - the coarse operator is the Galerkin product, so it is symmetric if the   !
!       fine one is;                                                              !
!     - the coarsest level is finished with the same number of forward and       !
!       backward sweeps.                                                          !
!                                                                                !
!   Drop any one of those and CG will still appear to run, but its residual will !
!   stop falling monotonically and it will converge to the wrong answer or not   !
!   at all.                                                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Agglomeration_Type), intent(inout) :: Agg  !! parent class
  integer,                   intent(in)    :: l    !! level to work on
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, k, n, ic, sweep
!==============================================================================!

  n = Agg % lev(l) % n

  !-------------------------------------------------------!
  !                                                       !
  !   The coarsest level: just smooth it hard and return  !
  !                                                       !
  !-------------------------------------------------------!
  if(l .eq. Agg % n_levels) then
    Agg % lev(l) % x(1:n) = 0.0
    do sweep = 1, AGG_COARSE_SWEEPS
      call Agg % Smooth(l, .true.)
      call Agg % Smooth(l, .false.)
    end do
    return
  end if

  !------------------------!
  !                        !
  !   Pre-smooth, forward  !
  !                        !
  !------------------------!
  Agg % lev(l) % x(1:n) = 0.0
  call Agg % Smooth(l, .true.)

  !--------------------------!
  !   Residual, r = b - A x  !
  !--------------------------!
  do i = 1, n
    Agg % lev(l) % r(i) = Agg % lev(l) % b(i)
    do k = Agg % lev(l) % row(i), Agg % lev(l) % row(i+1) - 1
      j = Agg % lev(l) % col(k)
      Agg % lev(l) % r(i) = Agg % lev(l) % r(i)  &
                          - Agg % lev(l) % val(k) * Agg % lev(l) % x(j)
    end do
  end do

  !-------------------------------------------------!
  !   Restrict: the coarse right hand side is the    !
  !   sum of the residual over each aggregate        !
  !-------------------------------------------------!
  Agg % lev(l+1) % b(1:Agg % lev(l+1) % n) = 0.0
  do i = 1, n
    ic = Agg % lev(l) % agg(i)
    Agg % lev(l+1) % b(ic) = Agg % lev(l+1) % b(ic) + Agg % lev(l) % r(i)
  end do

  !----------------------------!
  !   Solve the coarse problem !
  !----------------------------!
  call Agg % V_Cycle(l+1)

  !-------------------------------------------------!
  !   Prolongate: every member of an aggregate       !
  !   picks up that aggregate's correction           !
  !-------------------------------------------------!
  do i = 1, n
    ic = Agg % lev(l) % agg(i)
    Agg % lev(l) % x(i) = Agg % lev(l) % x(i) + Agg % lev(l+1) % x(ic)
  end do

  !----------------------------!
  !                            !
  !   Post-smooth, backward    !
  !                            !
  !----------------------------!
  call Agg % Smooth(l, .false.)

  end subroutine
