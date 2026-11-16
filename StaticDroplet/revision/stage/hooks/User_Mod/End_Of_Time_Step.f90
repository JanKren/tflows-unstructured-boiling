!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,    &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   Static droplet traces.                                                     !
!                                                                              !
!   spurious.dat      the six columns of the published runs, unchanged:        !
!                     time, volume of fun, u_rms (cell-count), u_max (global), !
!                     p_max - p_min, p_in - p_out                              !
!                                                                              !
!   spurious_ext.dat  where the largest velocity sits.  A polyhedral dual has  !
!                     small skewed cells in the box corners, and the global    !
!                     maximum of the published polyhedral run is there, in     !
!                     single-phase fluid 0.01 from a wall, not at the          !
!                     interface.  Columns:                                     !
!                      1 time                                                  !
!                      2 u_max, all cells                                      !
!                      3 u_max, interface cells (1e-3 < fun < 1 - 1e-3)        !
!                      4 u_max, band |r - R0| < 0.25 R0 about the box centre   !
!                      5 u_max, cells further than 0.2 from every wall         !
!                      6 u_max, sphere phase (fun > 0.5)                       !
!                      7 u_max, outer phase  (fun < 0.5)                       !
!                      8-10 position of the global maximum                     !
!                     11 kinetic energy, sum of rho u^2 vol / 2                !
!                     12 u_rms, volume weighted, all cells                     !
!                     13-15 centroid of fun                                    !
!                     16 min(fun)   17 max(fun)                                !
!                     18 p_in - p_out                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer, intent(in)      :: n_stat_t
  integer, intent(in)      :: n_stat_p
!------------------------------[Local parameters]------------------------------!
  real, parameter :: X0 = 1.0, Y0 = 1.0, Z0 = 1.0   ! box centre
  real, parameter :: R0 = 0.4, L_BOX = 2.0          ! sphere radius, box side
!--------------------------------[Locals]--------------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: fun
  integer                  :: c, fu, n_tot_cells
  real                     :: u_rms, a_in, a_out, p_in, p_out, a_vof
  real                     :: u_res, u_max, p_max, p_min
  real                     :: min_vfrac, max_vfrac
  real                     :: epsloc
  real                     :: um_if, um_band, um_core, um_in, um_out
  real                     :: xm, ym, zm, ke, vol_tot, u2v
  real                     :: cx, cy, cz, r, dw, f
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  fun  => Vof % fun

  epsloc = epsilon(epsloc)

  u_rms = 0.0
  a_in  = 0.0
  a_out = 0.0
  p_in  = 0.0
  p_out = 0.0
  a_vof = 0.0

  u_max   = 0.0
  um_if   = 0.0
  um_band = 0.0
  um_core = 0.0
  um_in   = 0.0
  um_out  = 0.0
  ke      = 0.0
  u2v     = 0.0
  vol_tot = 0.0
  cx      = 0.0
  cy      = 0.0
  cz      = 0.0
  p_max   = -HUGE
  p_min   = +HUGE
  min_vfrac = +HUGE
  max_vfrac = -HUGE

  ! Extremes of the volume fraction, which select the cells of each phase
  n_tot_cells = 0
  do c = Cells_In_Domain()
    min_vfrac = min(min_vfrac, fun % n(c))
    max_vfrac = max(max_vfrac, fun % n(c))
    n_tot_cells = n_tot_cells + 1
  end do
  call Global % Min_Real(min_vfrac)
  call Global % Max_Real(max_vfrac)
  call Global % Sum_Int(n_tot_cells)

  do c = Cells_In_Domain()
    u_res = sqrt(  Flow % u % n(c) ** 2  &
                 + Flow % v % n(c) ** 2  &
                 + Flow % w % n(c) ** 2)
    f = fun % n(c)

    u_rms = u_rms + u_res ** 2
    u_max = max(u_max, u_res)

    p_max = max(p_max, Flow % p % n(c))
    p_min = min(p_min, Flow % p % n(c))

    ! Published definition of the two pressure averages
    if(abs(max_vfrac - f) < epsloc) then
      a_in = a_in + Grid % vol(c)
      p_in = p_in + Flow % p % n(c) * Grid % vol(c)
    end if
    if(abs(f) < min_vfrac + epsloc) then
      a_out = a_out + Grid % vol(c)
      p_out = p_out + Flow % p % n(c) * Grid % vol(c)
    end if

    a_vof = a_vof + Grid % vol(c) * f

    ! Extended
    r  = sqrt(  (Grid % xc(c) - X0) ** 2  &
              + (Grid % yc(c) - Y0) ** 2  &
              + (Grid % zc(c) - Z0) ** 2)
    dw = min(Grid % xc(c), L_BOX - Grid % xc(c),  &
             Grid % yc(c), L_BOX - Grid % yc(c),  &
             Grid % zc(c), L_BOX - Grid % zc(c))

    if(f > 1.0e-3 .and. f < 1.0 - 1.0e-3) um_if   = max(um_if,   u_res)
    if(abs(r - R0) < 0.25 * R0)           um_band = max(um_band, u_res)
    if(dw > 0.2)                          um_core = max(um_core, u_res)
    if(f > 0.5) then
      um_in  = max(um_in,  u_res)
    else
      um_out = max(um_out, u_res)
    end if

    ke      = ke      + 0.5 * Flow % density(c) * u_res ** 2 * Grid % vol(c)
    u2v     = u2v     + u_res ** 2 * Grid % vol(c)
    vol_tot = vol_tot + Grid % vol(c)
    cx      = cx      + Grid % xc(c) * f * Grid % vol(c)
    cy      = cy      + Grid % yc(c) * f * Grid % vol(c)
    cz      = cz      + Grid % zc(c) * f * Grid % vol(c)
  end do

  call Global % Sum_Real(a_vof)
  call Global % Sum_Real(u_rms)
  u_rms = sqrt(u_rms / real(n_tot_cells))

  call Global % Max_Real(u_max)

  ! Where the global maximum is: the owner reports its position
  xm = -HUGE
  ym = -HUGE
  zm = -HUGE
  do c = Cells_In_Domain()
    u_res = sqrt(  Flow % u % n(c) ** 2  &
                 + Flow % v % n(c) ** 2  &
                 + Flow % w % n(c) ** 2)
    if(u_res >= u_max) then
      xm = Grid % xc(c)
      ym = Grid % yc(c)
      zm = Grid % zc(c)
      exit
    end if
  end do
  call Global % Max_Real(xm)
  call Global % Max_Real(ym)
  call Global % Max_Real(zm)

  call Global % Max_Real(p_max)
  call Global % Min_Real(p_min)
  call Global % Sum_Real(p_in)
  call Global % Sum_Real(p_out)
  call Global % Sum_Real(a_in)
  call Global % Sum_Real(a_out)
  p_in  = p_in  / a_in
  p_out = p_out / a_out

  call Global % Max_Real(um_if)
  call Global % Max_Real(um_band)
  call Global % Max_Real(um_core)
  call Global % Max_Real(um_in)
  call Global % Max_Real(um_out)
  call Global % Sum_Real(ke)
  call Global % Sum_Real(u2v)
  call Global % Sum_Real(vol_tot)
  call Global % Sum_Real(cx)
  call Global % Sum_Real(cy)
  call Global % Sum_Real(cz)

  ! Write the traces.  A restart continues them; a cold start replaces them.
  if(First_Proc()) then
    if(Time % Curr_Dt() .eq. 1) then
      call File % Delete('spurious.dat')
      call File % Delete('spurious_ext.dat')
    end if
    call File % Append_For_Writing_Ascii('spurious.dat', fu)
    write(fu,'(6(2X,E16.10E2))')  Time % Get_Time(), a_vof, u_rms, u_max, &
                                  p_max-p_min, p_in-p_out
    close(fu)
    call File % Append_For_Writing_Ascii('spurious_ext.dat', fu)
    write(fu,'(18(1X,ES15.7E3))')  Time % Get_Time(),                    &
                                   u_max, um_if, um_band, um_core,        &
                                   um_in, um_out, xm, ym, zm,             &
                                   ke, sqrt(u2v / vol_tot),               &
                                   cx / a_vof, cy / a_vof, cz / a_vof,    &
                                   min_vfrac, max_vfrac, p_in - p_out
    close(fu)
  end if

  end subroutine
