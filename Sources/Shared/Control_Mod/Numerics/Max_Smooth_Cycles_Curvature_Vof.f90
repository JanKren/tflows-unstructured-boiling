!==============================================================================!
  subroutine Max_Smooth_Cycles_Curvature_Vof(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  Reads maximum number of cycles for curvature estimation in VOF.
!>
!>  Accepts two spellings of the keyword:
!>    * MAX_SMOOTHING_CYCLES_CURVATURE_VOF  (canonical)
!>    * MAX_SMOOTH_CYCLES_CURVATURE_VOF     (deprecated alias, no "ING")
!>  Many historical test controls use the shorter spelling; without this alias
!>  it was silently ignored and the default (2) was used instead, hurting
!>  reproducibility.  If both are present the canonical spelling wins.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! max smoothing cycles for curvature in VOF
  logical, optional    :: verbose  !! controls output verbosity
!-----------------------------------[Locals]-----------------------------------!
  integer, parameter :: NOT_SET = -987654321  ! sentinel for "keyword absent"
  integer            :: val_canon, val_alias
!==============================================================================!

  ! Probe both spellings (no verbose -> no spurious "keyword not found" notes)
  call Control % Read_Int_Item('MAX_SMOOTHING_CYCLES_CURVATURE_VOF',  &
                                NOT_SET, val_canon)
  call Control % Read_Int_Item('MAX_SMOOTH_CYCLES_CURVATURE_VOF',     &
                                NOT_SET, val_alias)

  if(val_canon .ne. NOT_SET) then
    val = val_canon                       ! canonical spelling takes precedence
  else if(val_alias .ne. NOT_SET) then
    val = val_alias                       ! deprecated spelling accepted
    if(First_Proc()) then
      print '(a)', ' # NOTE! Keyword MAX_SMOOTH_CYCLES_CURVATURE_VOF is '  //  &
                   'deprecated.'
      print '(a)', ' #       Please use MAX_SMOOTHING_CYCLES_CURVATURE_VOF '//  &
                   'instead.'
    end if
  else
    val = 2                               ! default (as before)
  end if

  end subroutine
