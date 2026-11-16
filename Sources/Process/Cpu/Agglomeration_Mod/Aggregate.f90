!==============================================================================!
  subroutine Aggregate(Agg, l, n_coarse)
!------------------------------------------------------------------------------!
!>  Groups the unknowns of level 'l' into aggregates by greedy pairwise
!>  matching along the strongest coupling, and returns how many aggregates
!>  were made.  The result is written into Agg % lev(l) % agg.
!------------------------------------------------------------------------------!
!   An unknown is paired with the unaggregated neighbour it is most strongly    !
!   coupled to, strength being the most negative off-diagonal entry.  Negative  !
!   is the sign an elliptic operator gives a coupling that actually carries      !
!   information; a positive off-diagonal is an artefact of mesh non-orthogonality !
!   and pairing along it would group cells the operator does not really connect. !
!                                                                               !
!   Whatever is left unmatched becomes an aggregate on its own.  That is not a   !
!   failure: a singleton simply passes through to the coarse level unchanged.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Agglomeration_Type), intent(inout) :: Agg       !! parent class
  integer,                   intent(in)    :: l         !! level to aggregate
  integer,                   intent(out)   :: n_coarse  !! aggregates made
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j, k, n, best
  real    :: best_val
!==============================================================================!

  n = Agg % lev(l) % n

  Agg % lev(l) % agg(1:n) = 0
  n_coarse = 0

  do i = 1, n

    if(Agg % lev(l) % agg(i) .ne. 0) cycle   ! already spoken for

    !-----------------------------------------------!
    !   Look for the strongest free partner of i    !
    !-----------------------------------------------!
    best     = 0
    best_val = 0.0
    do k = Agg % lev(l) % row(i), Agg % lev(l) % row(i+1) - 1
      j = Agg % lev(l) % col(k)
      if(j .eq. i)                     cycle  ! the diagonal is not a coupling
      if(Agg % lev(l) % agg(j) .ne. 0) cycle  ! partner already taken
      if(Agg % lev(l) % val(k) .lt. best_val) then
        best_val = Agg % lev(l) % val(k)
        best     = j
      end if
    end do

    !--------------------------------------------------!
    !   Open a new aggregate, with or without a mate   !
    !--------------------------------------------------!
    n_coarse = n_coarse + 1
    Agg % lev(l) % agg(i) = n_coarse
    if(best .gt. 0) Agg % lev(l) % agg(best) = n_coarse

  end do

  end subroutine
