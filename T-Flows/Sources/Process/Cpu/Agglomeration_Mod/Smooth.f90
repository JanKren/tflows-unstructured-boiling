!==============================================================================!
  subroutine Smooth(Agg, l, forward)
!------------------------------------------------------------------------------!
!>  One Gauss-Seidel sweep over level 'l', improving that level's x against
!>  that level's b, running forwards or backwards through the unknowns.
!------------------------------------------------------------------------------!
!   The direction is an argument because the cycle has to sweep forwards on the  !
!   way down and backwards on the way up.  A single direction would make the     !
!   V-cycle non-symmetric, and a non-symmetric preconditioner is not merely      !
!   slower for Conjugate Gradient, it breaks the short recurrence the method is   !
!   built on, so the iteration stops converging to the right answer.  Pairing     !
!   the two directions makes the cycle its own transpose.                        !
!                                                                                !
!   The vectors are reached through Agg rather than passed in, because passing    !
!   Agg % lev(l) % x into a dummy argument of a procedure that also receives Agg  !
!   would let the same storage be seen through two names, which the compiler is   !
!   entitled to assume cannot happen.                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Agglomeration_Type), intent(inout) :: Agg      !! parent class
  integer,                   intent(in)    :: l        !! level to smooth
  logical,                   intent(in)    :: forward  !! sweep direction
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, k, n
  real    :: sum
!==============================================================================!

  n = Agg % lev(l) % n

  if(forward) then

    do i = 1, n
      sum = Agg % lev(l) % b(i)
      do k = Agg % lev(l) % row(i), Agg % lev(l) % row(i+1) - 1
        j = Agg % lev(l) % col(k)
        if(j .ne. i) sum = sum - Agg % lev(l) % val(k) * Agg % lev(l) % x(j)
      end do
      Agg % lev(l) % x(i) = sum * Agg % lev(l) % d_inv(i)
    end do

  else

    do i = n, 1, -1
      sum = Agg % lev(l) % b(i)
      do k = Agg % lev(l) % row(i), Agg % lev(l) % row(i+1) - 1
        j = Agg % lev(l) % col(k)
        if(j .ne. i) sum = sum - Agg % lev(l) % val(k) * Agg % lev(l) % x(j)
      end do
      Agg % lev(l) % x(i) = sum * Agg % lev(l) % d_inv(i)
    end do

  end if

  end subroutine
