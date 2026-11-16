!==============================================================================!
  subroutine Preconditioner_For_Pressure(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads the preconditioner to use for the pressure equation alone, falling
!>  back on whatever PRECONDITIONER_FOR_SYSTEM_MATRIX asked for.
!------------------------------------------------------------------------------!
!   The pressure equation is the only elliptic one T-Flows solves, and it is     !
!   the only one whose iteration count grows as the mesh is refined.  Momentum   !
!   and VOF carry a strong time-derivative term on the diagonal, are solved to   !
!   convergence in a handful of iterations whatever the mesh, and have nothing   !
!   for a multigrid cycle to remove.  Building a hierarchy for them would be     !
!   setup cost against no saving at all.                                        !
!                                                                                !
!   Hence this second keyword.  Leaving it out changes nothing: the pressure     !
!   equation keeps the preconditioner every other equation is using.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
  character(SL)       :: val      !! preconditioner, on entry the default
  logical, optional   :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: def
!==============================================================================!

  ! Whatever the caller has already resolved becomes the default here
  def = val

  call Control % Read_Char_Item('PRECONDITIONER_FOR_PRESSURE',  &
                                 def,                           &
                                 val,                           &
                                 verbose)
  call String % To_Lower_Case(val)

  ! Check validity of the input
  if( val .ne. 'none'                     .and.  &
      val .ne. 'diagonal'                 .and.  &
      val .ne. 'jacobi'                   .and.  &
      val .ne. 'incomplete_cholesky'      .and.  &
      val .ne. 'icc'                      .and.  &
      val .ne. 'agglomeration_multigrid'  .and.  &
      val .ne. 'amg') then
    call Message % Error(72,                                          &
             'Unknown preconditioner for the pressure: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',                &
             file=__FILE__, line=__LINE__, one_proc=.true.)

  ! Use the same names as PETSc
  else
    if(val .eq. 'diagonal')                val = 'jacobi'
    if(val .eq. 'incomplete_cholesky')     val = 'icc'
    if(val .eq. 'agglomeration_multigrid') val = 'amg'
  end if

  end subroutine
