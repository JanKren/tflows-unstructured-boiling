!==============================================================================!
  subroutine Gu_Correction(Control, corr, def, verbose)
!------------------------------------------------------------------------------!
!>  Reads from the control file if Gu's correction will be used.
!>
!>  The optional def argument sets what is assumed when the keyword is absent
!>  from the control file; it stays 'no' if def is not given.  Callers use it
!>  to turn the correction on by default for two-phase flow, where a
!>  non-gradient body force is always present.  An explicit GU_CORRECTION in
!>  the control file always wins over def.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  logical, intent(out) :: corr     !! output value, true or false
  logical, optional    :: def      !! value assumed if the keyword is absent
  logical, optional    :: verbose  !! verbosity of the output
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: val, def_val
!==============================================================================!

  def_val = 'no'
  if(present(def)) then
    if(def) def_val = 'yes'
  end if

  call Control % Read_Char_Item('GU_CORRECTION', def_val, val, verbose)
  call String % To_Upper_Case(val)

  if( val .eq. 'YES' ) then
    corr = .true.

  else if( val .eq. 'NO' ) then
    corr = .false.

  else
    call Message % Error(60,                                   &
             'Unknown state for Gu correction: '//trim(val)//  &
             '. \n This error is critical.  Exiting.',         &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
