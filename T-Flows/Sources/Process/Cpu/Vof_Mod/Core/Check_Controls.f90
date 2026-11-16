!==============================================================================!
  subroutine Check_Controls(Vof, Flow, Turb)
!------------------------------------------------------------------------------!
!   Stops with a message on the combinations of controls that the VOF module   !
!   does not support, or would run with a result that is silently wrong.       !
!                                                                              !
!   Call it once per domain, after every control of the domain is read.        !
!   It decides from control values only, so all ranks reach the same verdict.  !
!                                                                              !
!   Each message says what is wrong and what to do about it.  A combination    !
!   that works but is known to be unreliable is reported with a NOTE and not   !
!   stopped, because production cases use it.                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),  target :: Vof
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  logical                  :: phase_change
!==============================================================================!

  if(.not. Flow % with_interface) return

  phase_change = (Flow % mass_transfer_model .ne. NO_MASS_TRANSFER)

  !---------------------------------------------!
  !                                             !
  !   Phase change needs what it is built on    !
  !                                             !
  !---------------------------------------------!
  if(phase_change .and. .not. Flow % heat_transfer) then
    call Message % Error(72,                                                  &
             'MASS_TRANSFER_MODEL is on, but HEAT_TRANSFER is not.  Phase ' // &
             'change is driven by the temperature field.  Set '             // &
             'HEAT_TRANSFER yes. \n Exiting!',                                 &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  if(phase_change .and. Flow % heat_transfer) then
    if(.not. Flow % t % blend_matrix) then
      call Message % Error(72,                                                &
               'BLEND_SYSTEM_MATRICES no is not supported with phase '      //&
               'change.  Without the upwind part in the matrix the energy ' //&
               'equation keeps a term that grows with the temperature '     //&
               'level itself, and with the saturation temperature of '      //&
               'water it diverges within a few steps.  Remove the '         //&
               'keyword (the default is yes). \n Exiting!',                   &
               file=__FILE__, line=__LINE__, one_proc=.true.)
    end if
  end if

  if(phase_change .and. Turb % model .ne. NO_TURBULENCE_MODEL .and.          &
                        Turb % model .ne. DNS) then
    call Message % Error(72,                                                  &
             'Phase change is not available with a turbulence model.  The ' //&
             'turbulent conductivity is dropped in the cells away from '    //&
             'the interface, so the heat flux that drives the phase '       //&
             'change would be wrong.  Use TURBULENCE_MODEL none (or dns) '  //&
             'or MASS_TRANSFER_MODEL none. \n Exiting!',                       &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  if(phase_change .and. .not. Vof % track_front) then
    call Message % Error(72,                                                  &
             'The algebraic VOF takes its phase-change source from the '    //&
             'front, which exists only with TRACK_FRONT yes.  Set '         // &
             'TRACK_FRONT yes. \n Exiting!',                                  &
             file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  end subroutine
