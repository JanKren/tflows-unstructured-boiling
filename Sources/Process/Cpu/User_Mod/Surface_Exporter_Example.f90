!==============================================================================!
  subroutine User_Mod_Surface_Exporter_Example(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   User-defined surface extraction and export functionality. Users can       !
!   define cutting planes here to extract specific surfaces from the domain   !
!   similar to OpenFOAM's surface extraction capability.                      !
!                                                                              !
!   USAGE:                                                                     !
!   1. Set surface_extraction_enabled to .true. to enable surface extraction  !
!   2. Uncomment and modify the cutting plane definitions below               !
!   3. Set the origin point (where the plane passes through)                 !
!   4. Set the normal vector (defines plane orientation)                      !
!   5. Set the save interval to control export frequency                      !
!                                                                              !
!   The extracted surfaces will be saved as .vtu files that can be           !
!   visualized in ParaView or VisIt.                                          !
!------------------------------------------------------------------------------!
  use Surface_Exporter_Mod
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
  real                     :: origin(3), normal(3)
  logical, save            :: first_call = .true.

  ! Enable/disable surface extraction
  logical, parameter       :: surface_extraction_enabled = .false.  ! Change to .true. to enable
!==============================================================================!

  ! Early exit if surface extraction is disabled
  if (.not. surface_extraction_enabled) return

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  u    => Flow % u
  v    => Flow % v
  w    => Flow % w
  t    => Flow % t

  !-------------------------------------------------------------!
  !   Define cutting planes (modify according to your needs)   !
  !-------------------------------------------------------------!
  if (first_call) then

    !-------------------------------------------------------------------!
    !   Set surface export frequency (independent of RESULTS_SAVE_INTERVAL)
    !-------------------------------------------------------------------!
    ! Set how often to save surfaces (in time steps)
    ! Examples:
    !   - Every time step: call Surface_Exporter % Set_Save_Interval(1)
    !   - Every 10 time steps: call Surface_Exporter % Set_Save_Interval(10)
    !   - Every 50 time steps: call Surface_Exporter % Set_Save_Interval(50)
    call Surface_Exporter % Set_Save_Interval(10)  ! Save every 10 time steps

    !-------------------------------------------------------------------!
    !   Example cutting plane definitions - MODIFY THESE AS NEEDED     !
    !-------------------------------------------------------------------!

    ! Example 1: Vertical plane at x = 0.0 (YZ plane) with enhanced extraction
    ! Uncomment the lines below to activate this cutting plane
    !
     origin(1) = 0.0   ! x-coordinate where plane passes through
     origin(2) = 0.0   ! y-coordinate (can be any value for YZ plane)
     origin(3) = 0.0   ! z-coordinate (can be any value for YZ plane)
     normal(1) = 1.0   ! normal pointing in x-direction
     normal(2) = 0.0
     normal(3) = 0.0
     ! Enhanced extraction with face centers for better density
     call Surface_Exporter % Add_Cut_Plane('x_plane_00', origin, normal, 5.0e-4, &
                                          extraction_method=3, search_radius=1.0e-5)

    ! Example 2: Horizontal plane at y = 0.25 (XZ plane)
    ! Uncomment the lines below to activate this cutting plane
    !
    ! origin(1) = 0.0   ! x-coordinate (can be any value for XZ plane)
    ! origin(2) = 0.25  ! y-coordinate where plane passes through
    ! origin(3) = 0.0   ! z-coordinate (can be any value for XZ plane)
    ! normal(1) = 0.0
    ! normal(2) = 1.0   ! normal pointing in y-direction
    ! normal(3) = 0.0
    ! call Surface_Exporter % Add_Cut_Plane('y_plane_025', origin, normal, 1.0e-6)

    ! Example 3: Plane at z = 0.1 (XY plane)
    ! Uncomment the lines below to activate this cutting plane
    !
    ! origin(1) = 0.0   ! x-coordinate (can be any value for XY plane)
    ! origin(2) = 0.0   ! y-coordinate (can be any value for XY plane)
    ! origin(3) = 0.1   ! z-coordinate where plane passes through
    ! normal(1) = 0.0
    ! normal(2) = 0.0
    ! normal(3) = 1.0   ! normal pointing in z-direction
    ! call Surface_Exporter % Add_Cut_Plane('z_plane_01', origin, normal, 1.0e-6)

    ! Example 4: Diagonal plane (45 degrees in XY plane)
    ! Uncomment the lines below to activate this cutting plane
    !
    ! origin(1) = 0.0   ! x-coordinate where plane passes through
    ! origin(2) = 0.0   ! y-coordinate where plane passes through
    ! origin(3) = 0.5   ! z-coordinate where plane passes through
    ! normal(1) = 1.0   ! normal components (will be automatically normalized)
    ! normal(2) = 1.0   ! This creates a 45-degree plane in XY
    ! normal(3) = 0.0
    ! call Surface_Exporter % Add_Cut_Plane('diagonal_xy', origin, normal, 1.0e-6)

    !-------------------------------------------------------------------!
    !   ADD YOUR OWN CUTTING PLANES HERE                               !
    !-------------------------------------------------------------------!
    !
    ! Template for adding new cutting planes:
    !
    ! origin(1) = ?     ! x-coordinate of point on plane
    ! origin(2) = ?     ! y-coordinate of point on plane
    ! origin(3) = ?     ! z-coordinate of point on plane
    ! normal(1) = ?     ! x-component of plane normal
    ! normal(2) = ?     ! y-component of plane normal
    ! normal(3) = ?     ! z-component of plane normal
    !
    ! Basic call:
    ! call Surface_Exporter % Add_Cut_Plane('my_plane_name', origin, normal, 1.0e-6)
    !
    ! Enhanced call with extraction options:
    ! call Surface_Exporter % Add_Cut_Plane('my_plane_name', origin, normal, tolerance, &
    !                                      extraction_method, search_radius)
    !
    ! Extraction methods:
    !   1 = Cell center intersections only (default, fastest)
    !   2 = Face centers within search radius (denser, slower)
    !   3 = Both methods combined (densest, slowest)
    !
    ! Search radius: Distance from plane to include points (only for methods 2 & 3)
    !   - Typical values: 0.001 to 0.1 depending on mesh size
    !   - Smaller = fewer points, larger = more points but potentially noisy

    first_call = .false.
  end if

  !---------------------------------------------------------!
  !   Extract and export surfaces (frequency controlled)   !
  !---------------------------------------------------------!

  ! Only process if cutting planes have been defined and it's time to save
  if (Surface_Exporter % n_cuts > 0 .and. &
      Surface_Exporter % Should_Save_Now(Time % Curr_Dt())) then

    ! Extract all surfaces based on current grid with field interpolation
    call Surface_Exporter % Extract_All_Surfaces(Grid, Flow, Vof)

    ! Export all extracted surfaces to VTU files
    call Surface_Exporter % Export_All_Surfaces()

    if (First_Proc()) then  ! Only print from first processor
      print *, 'Surface extraction completed for ', Surface_Exporter % n_cuts, &
               ' cutting planes at time step: ', Time % Curr_Dt()
    end if
  end if

  end subroutine
