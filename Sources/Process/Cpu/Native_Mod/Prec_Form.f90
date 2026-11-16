!==============================================================================!
  subroutine Prec_Form(Nat, ni, A, d, d_inv, prec)
!------------------------------------------------------------------------------!
!>  The Prec_Form subroutine is designed to form the preconditioning matrix 'd'
!>  from the provided matrix 'A'. This subroutine is a key component in
!>  preconditioned iterative solvers, enhancing convergence rates. It supports
!>  various preconditioning techniques, including Diagonal (Jacobi) and
!>  Incomplete Cholesky preconditioning. The choice of preconditioner is
!>  determined by the 'prec' parameter. For Diagonal preconditioning, only
!>  diagonal elements of 'A' are used. Incomplete Cholesky preconditioning
!>  involves a more complex calculation using both diagonal and lower
!>  triangular parts of 'A'. This subroutine also offers an option for no
!>  preconditioning, setting all elements of 'd' to 1.0. The 'Nat' object
!>  contains the grid and other related data, 'ni' specifies the size of the
!>  problem, and 'd' and 'd_inv' are output arrays representing the
!>  preconditioning matrix and its inverse, respectively.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type),        intent(in)  :: Nat       !! parent class
  integer,                   intent(in)  :: ni        !! number of unknowns
  type(Matrix_Type), target, intent(in)  :: A         !! system matrix
  real,                      intent(out) :: d(:)      !! preconditioned matrix
  real,                      intent(out) :: d_inv(:)  !! preconditioned inverse
                                                      !! matrix
  character(SL),             intent(in)  :: prec      !! preconditioner
!-----------------------------------[Locals]-----------------------------------!
  real                          :: sum, a_ki
  integer                       :: i, j, jt, k
  integer, contiguous,  pointer :: a_col(:), a_row(:), a_dia(:)
  real,    contiguous,  pointer :: a_val(:)
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Nat)
!==============================================================================!

  call Profiler % Start('Native_Prec_Form (all solvers)')

  ! Take some aliases
  a_col => A % col
  a_row => A % row
  a_dia => A % dia
  a_val => A % val

  !---------------------------------!
  !   1) Diagonal preconditioning   !
  !---------------------------------!
  if(prec .eq. 'jacobi') then
    !$omp parallel do private(i) shared (d, a_val, a_dia)
    do i = 1, ni
      d(i) = a_val(a_dia(i))
      d_inv(i) = 1.0 / d(i)
    end do
    !$omp end parallel do

  !--------------------------------------------!
  !   2) Incomplete Cholesky preconditioning   !
  !- - - - - - - - - - - - - - - - - - - - - - !
  !   Matrix storage requirements:             !
  !   1) First entry in each row must be the   !
  !      diagonal (a_col(a_row(i)) == i)       !
  !   2) Remaining column indices must be      !
  !      sorted in ascending order             !
  !                                            !
  !   This enables early loop termination when !
  !   processing lower/upper triangular parts  !
  !   during forward/backward substitution.    !
  !--------------------------------------------!
  else if(prec .eq. 'icc') then
    do i = 1, ni
      sum = a_val(a_dia(i))                ! take diagonal entry
      do j = a_row(i) + 1, a_row(i+1) - 1  ! row from the diagonal on
        k = a_col(j)
        if(k .gt. i) exit
        ! Pivot of M = (D+L) D^-1 (D+U), the form Prec_Solve applies:
        ! p_i = a_ii - sum_k a_ik a_ki / p_k, with d(k) = 1 / p_k here.
        ! a_ik**2 equals a_ik a_ki only for a symmetric A; upwinded
        ! advection makes the two differ by orders of magnitude, and
        ! a_ik**2 then drives pivots to zero or below.
        a_ki = a_val(j)
        do jt = a_row(k) + 1, a_row(k+1) - 1
          if(a_col(jt) .eq. i) then
            a_ki = a_val(jt)
            exit
          end if
        end do
        sum = sum - d(k) * a_val(j) * a_ki
      end do
      d_inv(i) = sum
      d(i)     = 1.0 / sum
    end do

  !------------------------------------------------!
  !   3) Aggregation multigrid                     !
  !- - - - - - - - - - - - - - - - - - - - - - - - +-----------------------!
  !   The whole preconditioner is the hierarchy, which is rebuilt here     !
  !   from the matrix as it stands now.  Nothing is left in d or d_inv,    !
  !   but they are still filled with the diagonal, because a caller that   !
  !   ever wanted to fall back on Jacobi would find them ready.            !
  !-----------------------------------------------------------------------!
  else if(prec .eq. 'amg') then
    do i = 1, ni
      d(i)     = a_val(a_dia(i))
      d_inv(i) = 1.0 / d(i)
    end do
    call Agglomeration % Setup(ni, A)

  !---------------------------!
  !   .) No preconditioning   !
  !---------------------------!
  else
    !$omp parallel do private(i) shared (d)
    do i = 1, ni
      d(i) = 1.0
    end do
    !$omp end parallel do
  end if

  call Profiler % Stop('Native_Prec_Form (all solvers)')

  end subroutine
