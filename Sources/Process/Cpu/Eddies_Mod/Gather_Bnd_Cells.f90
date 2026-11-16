!==============================================================================!
  subroutine Gather_Bnd_Cells(Eddies)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to gather boundary cells corresponding to a
!>  specific boundary condition (where eddies are introduced) across all
!>  processors. This subroutine ensures that the data regarding eddy placement
!>  is consistent and continuous across processor boundaries, preventing
!>  discontinuities in the flow field.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Alias creation: Creates aliases for the grid (Grid) and the flow field   !
!     (Flow) for easier reference.                                             !
!   * Boundary cell counting: Counts the number of boundary cells in the       !
!     current processor that match the specified boundary condition.           !
!     This count is stored in eddies % n_bnd_cells.                            !
!   * Global boundary cell count: If running in a parallel environment, the    !
!     subroutine computes the total number of boundary cells across all        !
!     processors (eddies % n_bnd_cells_glo) and the start index of boundary    !
!     cells for each processor (s_bnd_cell_proc).                              !
!   * Boundary cell data gathering: Collects coordinates (x, y, z), wall       !
!     distance, and flow velocity (u, v, w) for each boundary cell that        !
!     matches the boundary condition. This data is stored in arrays within     !
!     the eddies structure.                                                    !
!   * Inlet plane position: Determines the average position of the inlet       !
!     plane based on the spread (difference between maximum and minimum        !
!     values) of the boundary cell coordinates. This helps in identifying      !
!     which plane (x, y, or z) the eddies are associated with.                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Eddies_Type), target :: Eddies  !! parent class; collection of eddies
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  integer, allocatable      :: n_bnd_cells_proc(:)
  integer, allocatable      :: s_bnd_cell_proc(:)
  integer                   :: c, n, cnt
  real                      :: x_dff, y_dff, z_dff, x_avg, y_avg, z_avg
  real                      :: x_min, x_max, y_min, y_max, z_min, z_max
  real                      :: radial_dist
!==============================================================================!

  ! Create an alias for Grid
  Grid => Eddies % pnt_grid
  Flow => Eddies % pnt_flow

  !------------------------------------------------------!
  !   Store boundary cells at given boundary condition   !
  !------------------------------------------------------!

  !--------------------------------------------------------------------!
  !   First pass: determine streamwise direction by finding bounds     !
  !--------------------------------------------------------------------!
  x_min = HUGE; x_max = -HUGE
  y_min = HUGE; y_max = -HUGE
  z_min = HUGE; z_max = -HUGE

  do c = -Grid % n_bnd_cells, -1
    if(Grid % Bnd_Cond_Name_At_Cell(c) .eq. Eddies % bc_name) then
      x_min = min(x_min, Grid % xc(c)); x_max = max(x_max, Grid % xc(c))
      y_min = min(y_min, Grid % yc(c)); y_max = max(y_max, Grid % yc(c))
      z_min = min(z_min, Grid % zc(c)); z_max = max(z_max, Grid % zc(c))
    end if
  end do

  ! Exchange min/max values across all processors
  call Global % Min_Real(x_min); call Global % Max_Real(x_max)
  call Global % Min_Real(y_min); call Global % Max_Real(y_max)
  call Global % Min_Real(z_min); call Global % Max_Real(z_max)

  ! Determine which direction has smallest spread (inlet plane)
  x_dff = x_max - x_min
  y_dff = y_max - y_min
  z_dff = z_max - z_min

  ! Count boundary cells in this processor (including radius filtering if specified)
  Eddies % n_bnd_cells = 0
  do c = -Grid % n_bnd_cells, -1
    if(Grid % Bnd_Cond_Name_At_Cell(c) .eq. Eddies % bc_name) then
      ! Apply radius filtering if specified
      if(Eddies % radius_limit < HUGE) then
        ! Calculate radial distance based on streamwise direction
        if(x_dff .eq. min(x_dff, y_dff, z_dff)) then
          ! x is streamwise, use y-z plane for radius
          radial_dist = sqrt(Grid % yc(c)**2 + Grid % zc(c)**2)
        else if(y_dff .eq. min(x_dff, y_dff, z_dff)) then
          ! y is streamwise, use x-z plane for radius
          radial_dist = sqrt(Grid % xc(c)**2 + Grid % zc(c)**2)
        else
          ! z is streamwise, use x-y plane for radius
          radial_dist = sqrt(Grid % xc(c)**2 + Grid % yc(c)**2)
        end if

        if(radial_dist <= Eddies % radius_limit) then
          Eddies % n_bnd_cells = Eddies % n_bnd_cells + 1
        end if
      else
        Eddies % n_bnd_cells = Eddies % n_bnd_cells + 1
      end if
    end if
  end do

  ! Estimate first boundary cell for each processor
  if(Parallel_Run()) then
    allocate(n_bnd_cells_proc(N_Procs()));  n_bnd_cells_proc(:) = 0;
    n_bnd_cells_proc(This_Proc()) = Eddies % n_bnd_cells
    call Global % Sum_Int_Array(N_Procs(), n_bnd_cells_proc)

    allocate(s_bnd_cell_proc(N_Procs())); s_bnd_cell_proc(:) = 0
    do n = 2, N_Procs()
      s_bnd_cell_proc(n) = s_bnd_cell_proc(n-1) + n_bnd_cells_proc(n-1)
    end do
  else
    allocate(s_bnd_cell_proc(0:0));  s_bnd_cell_proc(:) = 0;
  end if

  ! Gather coordinates from all processors
  Eddies % n_bnd_cells_glo = Eddies % n_bnd_cells
  call Global % Sum_Int(Eddies % n_bnd_cells_glo)
  if(First_Proc()) then
    print '(a,a,i6)', ' # Number of boundary cells at ',  &
                      trim(Eddies % bc_name),             &
                      Eddies % n_bnd_cells_glo
    if(Eddies % radius_limit < HUGE) then
      if(x_dff .eq. min(x_dff, y_dff, z_dff)) then
        print '(a,f8.4,a)', ' # Radius filtering applied: R <= ', Eddies % radius_limit, ' (x-streamwise)'
      else if(y_dff .eq. min(x_dff, y_dff, z_dff)) then
        print '(a,f8.4,a)', ' # Radius filtering applied: R <= ', Eddies % radius_limit, ' (y-streamwise)'
      else
        print '(a,f8.4,a)', ' # Radius filtering applied: R <= ', Eddies % radius_limit, ' (z-streamwise)'
      end if
    end if
  end if
  allocate(Eddies % bnd_xc(Eddies % n_bnd_cells_glo));  Eddies % bnd_xc(:) = 0.0
  allocate(Eddies % bnd_yc(Eddies % n_bnd_cells_glo));  Eddies % bnd_yc(:) = 0.0
  allocate(Eddies % bnd_zc(Eddies % n_bnd_cells_glo));  Eddies % bnd_zc(:) = 0.0
  allocate(Eddies % bnd_wd(Eddies % n_bnd_cells_glo));  Eddies % bnd_wd(:) = 0.0
  allocate(Eddies % bnd_u (Eddies % n_bnd_cells_glo));  Eddies % bnd_u (:) = 0.0
  allocate(Eddies % bnd_v (Eddies % n_bnd_cells_glo));  Eddies % bnd_v (:) = 0.0
  allocate(Eddies % bnd_w (Eddies % n_bnd_cells_glo));  Eddies % bnd_w (:) = 0.0

  cnt = 0
  do c = -Grid % n_bnd_cells, -1
    if(Grid % Bnd_Cond_Name_At_Cell(c) .eq. Eddies % bc_name) then
      ! Apply radius filtering if specified (same logic as counting phase)
      if(Eddies % radius_limit < HUGE) then
        ! Calculate radial distance based on streamwise direction
        if(x_dff .eq. min(x_dff, y_dff, z_dff)) then
          ! x is streamwise, use y-z plane for radius
          radial_dist = sqrt(Grid % yc(c)**2 + Grid % zc(c)**2)
        else if(y_dff .eq. min(x_dff, y_dff, z_dff)) then
          ! y is streamwise, use x-z plane for radius
          radial_dist = sqrt(Grid % xc(c)**2 + Grid % zc(c)**2)
        else
          ! z is streamwise, use x-y plane for radius
          radial_dist = sqrt(Grid % xc(c)**2 + Grid % yc(c)**2)
        end if

        if(radial_dist > Eddies % radius_limit) then
          cycle  ! Skip this cell if outside radius limit
        end if
      end if

      cnt = cnt + 1
      Eddies % bnd_xc(cnt + s_bnd_cell_proc(This_Proc())) = Grid % xc(c)
      Eddies % bnd_yc(cnt + s_bnd_cell_proc(This_Proc())) = Grid % yc(c)
      Eddies % bnd_zc(cnt + s_bnd_cell_proc(This_Proc())) = Grid % zc(c)
      ! Calculate wall distance for pipe inlet if not available
      if(Grid % wall_dist(c) <= 0.0) then
        ! For pipe inlet: wall distance = distance from pipe center
        ! Assuming pipe center is at y=0, z=0 (adjust if needed)
        Eddies % bnd_wd(cnt + s_bnd_cell_proc(This_Proc())) = &
          sqrt(Grid % yc(c)**2 + Grid % zc(c)**2)
      else
        Eddies % bnd_wd(cnt + s_bnd_cell_proc(This_Proc())) = Grid % wall_dist(c)
      end if
      Eddies % bnd_u (cnt + s_bnd_cell_proc(This_Proc())) = Flow % u % b(c)
      Eddies % bnd_v (cnt + s_bnd_cell_proc(This_Proc())) = Flow % v % b(c)
      Eddies % bnd_w (cnt + s_bnd_cell_proc(This_Proc())) = Flow % w % b(c)
    end if
  end do
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_xc)
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_yc)
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_zc)
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_wd)
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_u)
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_v)
  call Global % Sum_Real_Array(Eddies % n_bnd_cells_glo, Eddies % bnd_w)

  !-------------------------------------------------------!
  !           Find position of the inlet plane            !
  !   (Important for checking eddies leaving the plane)   !
  !-------------------------------------------------------!
  Eddies % x_plane = HUGE
  Eddies % y_plane = HUGE
  Eddies % z_plane = HUGE

  x_dff = maxval(Eddies % bnd_xc(:)) - minval(Eddies % bnd_xc(:))
  y_dff = maxval(Eddies % bnd_yc(:)) - minval(Eddies % bnd_yc(:))
  z_dff = maxval(Eddies % bnd_zc(:)) - minval(Eddies % bnd_zc(:))

  x_avg = 0.5*(maxval(Eddies % bnd_xc(:)) + minval(Eddies % bnd_xc(:)))
  y_avg = 0.5*(maxval(Eddies % bnd_yc(:)) + minval(Eddies % bnd_yc(:)))
  z_avg = 0.5*(maxval(Eddies % bnd_zc(:)) + minval(Eddies % bnd_zc(:)))

  if(x_dff .eq. min(x_dff, y_dff, z_dff)) Eddies % x_plane = x_avg
  if(y_dff .eq. min(x_dff, y_dff, z_dff)) Eddies % y_plane = y_avg
  if(z_dff .eq. min(x_dff, y_dff, z_dff)) Eddies % z_plane = z_avg

  end subroutine
