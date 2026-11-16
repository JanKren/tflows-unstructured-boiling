!==============================================================================!
  subroutine Initialize(Monitor, Flow, restart, domain)
!------------------------------------------------------------------------------!
!>  This subroutine initializes the monitoring points for the Process sub-
!>  program in T-Flows. It reads the number of monitoring points and their
!>  coordinates from the control file, allocates memory for these points, and
!>  establishes file connections for each point. For each moonitoring point it
!>  finds the nearest computational cell and reports directly from that value,
!>  without any extrapolation. The monitoring files have extension 'monit-NNN'
!>  where NNN is point rank. The subroutine also handles the restart logic,
!>  ensuring that existing monitoring files are not over-written after restart.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Monitor_Type)      :: Monitor  !! parent class of Monitory_Type
  type(Field_Type), target :: Flow     !! flow field being monitored
  logical                  :: restart  !! true if current run is a restart
  integer,        optional :: domain   !! domain number
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: m, n, l, sc
  character(SL)            :: mon_file_name
  character(SL)            :: point_name
  real                     :: xyz(3), def(3)
!==============================================================================!

  ! Take alias for grid
  Grid => Flow % pnt_grid
  Monitor % pnt_grid => Grid

  ! Read number of Monitoring points from control file
  call Control % Read_Int_Item('NUMBER_OF_MONITORING_POINTS', 0, &
                                Monitor % n_points, .true.)

  if(Monitor % n_points .eq. 0) return

  ! Allocate memory for points
  allocate(Monitor % x(Monitor % n_points))
  allocate(Monitor % y(Monitor % n_points))
  allocate(Monitor % z(Monitor % n_points))

  ! Allocate memory accordingly
  allocate(Monitor % cell     (Monitor % n_points))
  allocate(Monitor % file_unit(Monitor % n_points))

  !----------------------------------------!
  !   Read Monitoring points coordinates   !
  !----------------------------------------!
  do n = 1, Monitor % n_points
    write(point_name, '(a,i3.3)') 'MONITORING_POINT_', n

    def = 0.  ! don't have a better idea what to set
    call Control % Read_Real_Vector(point_name, 3, def, xyz, .true.)

    Monitor % x(n) = xyz(1)
    Monitor % y(n) = xyz(2)
    Monitor % z(n) = xyz(3)
  end do

  !--------------------------------------------!
  !   Set the names for all Monitoring files   !
  !--------------------------------------------!
  call File % Set_Name(mon_file_name, extension='-monit.000', domain=domain)
  l = len_trim(mon_file_name)

  !-------------------------------!
  !   Find the Monitoring cells   !
  !-------------------------------!
  call Monitor % Find_Cells()

  !----------------------------------------------!
  !   Write first line in the Monitoring files   !
  !----------------------------------------------!
  do m = 1, Monitor % n_points

    if(Monitor % cell(m) > 0) then

      write(mon_file_name(l-2:l),'(I3.3)') m

      if(.not. restart) then
        call File % Open_For_Writing_Ascii(mon_file_name,  &
                                           Monitor % file_unit(m))
      else
        call File % Append_For_Writing_Ascii(mon_file_name,  &
                                             Monitor % file_unit(m))
      endif

      write(Monitor % file_unit(m), '(a20, 3f16.6)')   &
            '# Monitoring point:',                     &
            Grid % xc( Monitor % cell(m) ),            &
            Grid % yc( Monitor % cell(m) ),            &
            Grid % zc( Monitor % cell(m) )
      write(Monitor % file_unit(m), '(a11, 4a15)', advance='no')   &
            '#         ',                                          &
            ' U                 ',                                 &
            ' V                 ',                                 &
            ' W                 ',                                 &
            ' P                 '

      ! If heat transfer, write temperature too
      if(Flow % heat_transfer) then
        write(Monitor % file_unit(m), '(a15)', advance='no')   &
              ' T                 '
      end if

      do sc = 1, Flow % n_scalars
        write(Monitor % file_unit(m), '(a1,a4,a10)', advance='no')   &
              ' ', Flow % scalar(sc) % name, '          '
      end do

      ! End the line
      write(Monitor % file_unit(m),'(a)')  ' '


    end if

  end do

  end subroutine
