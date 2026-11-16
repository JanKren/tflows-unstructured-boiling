!==============================================================================!
  subroutine Stop_On_Removed_Features(Control)
!------------------------------------------------------------------------------!
!>  This version of T-Flows has only the features that the cases of one paper
!>  need.  T-Flows does not object to keywords it does not read, so a control
!>  file asking for a removed feature would run without it and give no sign.
!>  This subroutine reads the control file once, before anything else is
!>  done, and stops the run with a message on any keyword that would switch
!>  on one of the removed features.  A keyword with a value that keeps its
!>  feature off (no, none, one domain, no scalars) passes.  Keywords are
!>  recognised in any mix of upper and lower case, and every line with the
!>  keyword is checked, not only the first one.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type) :: Control  !! parent class
!-----------------------------------[Locals]-----------------------------------!
  character(SL) :: key, val
  integer       :: i, n, i_val, ios
  logical       :: reached_end, feature_on
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: N_KEYS = 7
  ! Keywords which switch a removed feature on ...
  character(SL), parameter :: KEYS(N_KEYS) = [character(SL) ::       &
    'GEOMETRICVOF',                                                  &
    'NUCLEATION',                                                    &
    'MASS_TRANSFER_MODEL',                                           &
    'NUMBER_OF_DOMAINS',                                             &
    'ADAPT_GRID',                                                    &
    'VARIABLE_SURFACE_TENSION',                                      &
    'SURFACE_TENSION_MODEL']
  ! ... the values which keep it off, separated by commas, or, after "<=",
  ! the largest integer value which keeps it off ...
  character(SL), parameter :: OFF(N_KEYS) = [character(SL) ::        &
    'NO',                                                            &
    'NO',                                                            &
    'NONE,TEMPERATURE_GRADIENTS',                                    &
    '<=1',                                                           &
    'NO',                                                            &
    'NO',                                                            &
    'CSF']
  ! ... and the feature which is not in this version
  character(SL), parameter :: WHAT(N_KEYS) = [character(SL) ::       &
    'the geometric VOF',                                             &
    'the nucleation model (nor its microlayer, near-wall and '    // &
    'recession models)',                                             &
    'mass transfer models other than TEMPERATURE_GRADIENTS',         &
    'multiple domains (nor conjugate heat transfer)',                &
    'mesh adaptation',                                               &
    'a variable surface tension (nor Marangoni forces)',             &
    'the CSS model of surface tension']
!==============================================================================!

  rewind(Control % file_unit)

  do
    call File % Read_Line(Control % file_unit, reached_end)
    if(reached_end) exit

    key = Line % tokens(1)
    val = Line % tokens(2)
    call String % To_Upper_Case(key)
    call String % To_Upper_Case(val)

    do i = 1, N_KEYS
      if(key .ne. KEYS(i)) cycle

      ! Is the feature switched on by this value?
      if(OFF(i)(1:2) .eq. '<=') then
        read(OFF(i)(3:), *) n
        read(val, *, iostat=ios) i_val
        feature_on = (ios .ne. 0)
        if(.not. feature_on) feature_on = (i_val .gt. n)
      else
        feature_on = index(','//trim(OFF(i))//',', ','//trim(val)//',') .eq. 0
      end if

      if(feature_on) then
        call Message % Framed(72,                                             &
                 RED//'ERROR'//RESET//' in the control file',                 &
                 'The control file sets '//trim(Line % tokens(1))//' '     // &
                 trim(Line % tokens(2))//', but this version of T-Flows '  // &
                 'does not have '//trim(WHAT(i))//'.  Remove the '         // &
                 'keyword, or set it to a value which keeps the feature '  // &
                 'off, to run without it. \n Exiting!',                       &
                 one_proc=.true.)

        ! Every processor reads the same file and comes here as well, so
        ! they end the run together, after the first one flushed the message
        flush(6)
        call Global % End_Parallel
        stop 1
      end if
    end do
  end do

  end subroutine
