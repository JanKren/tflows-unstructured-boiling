!==============================================================================!
  subroutine Main_Results(Results, n_dom, Flow, Turb, Vof, Swarm, exit_now)
!------------------------------------------------------------------------------!
!>  The Main_Results subroutine serves as the central function for saving
!>  results in T-Flows. It orchestrates the process of saving numerical results
!>  for post-processing and backups, ensuring that data is outputted correctly
!>  across different domains and under various conditions.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Conditional result saving: Determines whether to save initial conditions !
!     and handles exit conditions for result saving.                           !
!   * Domain-specific processing: Iterates through domains to save results     !
!     and backups for each domain separately.                                  !
!   * Post-processing output: Saves results in .vtu format for visualization,  !
!     including detailed data like flow fields, turbulence modeling data,      !
!     VOF multiphase modeling, and particle swarms.
!   * User-customized output: Allows for integration of user-defined result    !
!     saving routines through calls to User_Mod_Save_Fields and                !
!     User_Mod_Save_Swarm.                                                     !
!   * Handles file-based triggers for saving (save_now) and exiting (exit_now) !
!     operations, allowing user to force saving or gracefult exit from the     !
!     Process at any time during the simulation.                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results    !! parent class
  integer             :: n_dom      !! number of domains
  type(Field_Type)    :: Flow(MD)   !! flow fields
  type(Turb_Type)     :: Turb(MD)   !! turbulence modelling variables
  type(Vof_Type)      :: Vof(MD)    !! multiphase modelling with VOF
  type(Swarm_Type)    :: Swarm(MD)  !! swarms of particles
  logical             :: exit_now   !! true for premature exit
!----------------------------------[Locals]------------------------------------!
  integer :: d         ! domain counter
  logical :: save_now   ! save everything, whatever the time step
  logical :: flag(3)    ! exit file, save file, wall clock: seen by this rank
  logical :: time_exit  ! wall clock limit reached
!==============================================================================!

  exit_now = .false.  ! set before any return below

  ! Get out of here if you don't want to save initial conditions all right
  if( .not. Results % initial .and.  &
      Time % First_Dt() .eq. Time % Curr_Dt() ) then
    return
  end if

  !-----------------------------------------------------------------------!
  !   One decision for all ranks.                                         !
  !                                                                       !
  !   Saving a backup or the results is collective, and so is leaving     !
  !   the time loop.  The trigger files, and the wall clock limit         !
  !   (every rank compares its own reading of the clock), are looked at   !
  !   by each rank on its own, so a file that appears, or a limit that    !
  !   is crossed, between the readings of two ranks sent some of them     !
  !   into the save and out of the loop and left the others in the next   !
  !   time step: the run hung there.  The ranks now agree on all three:   !
  !   if any of them sees one, all of them act on it.                     !
  !-----------------------------------------------------------------------!
  inquire(file='exit_now', exist=flag(1))
  inquire(file='save_now', exist=flag(2))
  flag(3) = Info % Time_To_Exit()
  call Global % Lor_Log_Array(3, flag)
  save_now  = flag(2)
  time_exit = flag(3)
  exit_now  = flag(1) .or. flag(3)

  ! Is it time to save the backup file?
  if(Time % Curr_Dt() .eq. Time % Last_Dt() .or.  &
     save_now                               .or.  &
     exit_now                               .or.  &
     Backup % Time_To_Save_Backup()         .or.  &
     time_exit) then
    do d = 1, n_dom
      call Control % Switch_To_Domain(d)
      call Backup % Save(Flow(d), Turb(d), Vof(d), Swarm(d), dom=d)
    end do
  end if

  ! Is it time to save results for post-processing?
  if(Time % Curr_Dt() .eq. Time % Last_Dt() .or.  &
     save_now                               .or.  &
     exit_now                               .or.  &
     Results % Time_To_Save_Results()       .or.  &
     time_exit) then

    do d = 1, n_dom
      call Control % Switch_To_Domain(d)
      call Results % Save_Vtu_Fields(Flow(d), Turb(d), Vof(d), Swarm(d),  &
                                     plot_inside=.true., domain=d)
      call Results % Save_Vtu_Fields(Flow(d), Turb(d), Vof(d), Swarm(d),  &
                                     plot_inside=.false., domain=d)

      if(Flow(d) % with_interface) then
        if(Vof(d) % track_front) then
          call Results % Save_Vtu_Front(Vof(d) % Front)
        end if
        if(Vof(d) % track_surface) then
          call Results % Save_Vtu_Surf(Vof(d) % surf)
        end if
      end if

      ! Write results in user-customized format
      call User_Mod_Save_Results(Flow(d), Turb(d), Vof(d), Swarm(d), domain=d)
      if(Flow(d) % with_particles) then
        call User_Mod_Save_Swarm(Flow(d), Turb(d), Vof(d), Swarm(d), domain=d)
      end if

    end do  ! through domains
  end if

  ! Is it time to save particles for post-processing?
  if(Time % Curr_Dt() .eq. Time % Last_Dt() .or.  &
     save_now                               .or.  &
     exit_now                               .or.  &
     Results % Time_To_Save_Swarm()         .or.  &
     time_exit) then

    do d = 1, n_dom
      call Control % Switch_To_Domain(d)
      if(Flow(d) % with_particles) then
        call Results % Save_Vtu_Swarm(Swarm(d), domain=d)

        call User_Mod_Save_Swarm(Flow(d), Turb(d), Vof(d), Swarm(d), domain=d)
      end if

    end do  ! through domains
  end if

  ! Remove the trigger files; the first process looks again, since a file
  ! can have been seen by another rank before it shows here
  if(save_now) call Remove_File('save_now')
  if(flag(1))  call Remove_File('exit_now')

  contains

  !============================================================================!
    subroutine Remove_File(name)
  !----------------------------------------------------------------------------!
  !   The first process deletes the named trigger file if it exists.           !
  !----------------------------------------------------------------------------!
    implicit none
  !--------------------------------[Arguments]---------------------------------!
    character(*), intent(in) :: name
  !---------------------------------[Locals]-----------------------------------!
    logical :: found
  !============================================================================!

    if(First_Proc()) then
      inquire(file=name, exist=found)
      if(found) then
        open (9, file=name, status='old')
        close(9, status='delete')
      end if
    end if

    end subroutine

  end subroutine
