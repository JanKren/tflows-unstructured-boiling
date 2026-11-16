!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   Kinematic advection tests: per-step record of the advected sphere.         !
!                                                                              !
!   bench-data.dat columns:                                                    !
!    1 time            2 front area        3 volume of the sphere phase        !
!    4 R from volume   5 R_front (area-weighted centroid distance of front)    !
!    6 min(vof)        7 max(vof) - 1                                          !
!    8 E1 against the initial field                                            !
!    9 volume of cells with 0.01 < vof < 0.99                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer                  :: n_stat_t  ! 1st t.s. statistics turbulence
  integer                  :: n_stat_p  ! 1st t.s. statistics particles
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: fun
  type(Front_Type), pointer :: Front
  type(Elem_Type),  pointer :: Elem(:)
  integer                   :: c, e, fu
  real                      :: b_volume, surface, rad_volume, r_front
  real                      :: xg, yg, zg, f_min, f_max, e1, e1_ex, wid, t
  real, allocatable         :: frac(:)
!==============================================================================!

  ! Take aliases
  Grid  => Flow % pnt_grid
  fun   => Vof % fun
  Front => Vof % Front
  Elem  => Front % Elem

  b_volume = 0.0
  surface  = 0.0
  xg = 0.0;  yg = 0.0;  zg = 0.0
  f_min = +HUGE
  f_max = -HUGE
  e1  = 0.0
  wid = 0.0

  do c = Cells_In_Domain()
    b_volume = b_volume + Grid % vol(c) * (1.0 - fun % n(c))
    f_min = min(f_min, fun % n(c))
    f_max = max(f_max, fun % n(c))
    e1 = e1 + Grid % vol(c) * abs(fun % n(c) - kin_fun0(c))
    if(fun % n(c) > 0.01 .and. fun % n(c) < 0.99) wid = wid + Grid % vol(c)
    if(Vof % track_front) then
      e = Front % elem_in_cell(c)
      if(e > 0) then
        surface = surface + Elem(e) % area
        xg = xg + Elem(e) % area * Elem(e) % xe
        yg = yg + Elem(e) % area * Elem(e) % ye
        zg = zg + Elem(e) % area * Elem(e) % ze
      end if
    end if
  end do

  call Global % Sum_Real(b_volume)
  call Global % Sum_Real(surface)
  call Global % Sum_Real(xg)
  call Global % Sum_Real(yg)
  call Global % Sum_Real(zg)
  call Global % Sum_Real(e1)
  call Global % Sum_Real(wid)
  call Global % Min_Real(f_min)
  call Global % Max_Real(f_max)

  e1 = e1 / kin_vol0
  rad_volume = (3.0*b_volume/(4.0*PI))**(1.0/3.0)

  ! Area-weighted mean distance of the front elements from their centroid
  r_front = 0.0
  if(Vof % track_front .and. surface > 0.0) then
    xg = xg / surface;  yg = yg / surface;  zg = zg / surface
    do c = Cells_In_Domain()
      e = Front % elem_in_cell(c)
      if(e > 0) then
        r_front = r_front + Elem(e) % area                       &
                          * sqrt(  (Elem(e) % xe - xg)**2        &
                                 + (Elem(e) % ye - yg)**2        &
                                 + (Elem(e) % ze - zg)**2)
      end if
    end do
    call Global % Sum_Real(r_front)
    r_front = r_front / surface
  end if

  if(First_Proc()) then
    call File % Append_For_Writing_Ascii('bench-data.dat', fu)
    write(fu,'(9(1x,es18.10))') Time % Get_Time(), surface, b_volume,  &
                                rad_volume, r_front, f_min, f_max - 1.0,  &
                                e1, wid
    close(fu)
  end if

  !------------------------------------------------------------------!
  !   At the last step of a one-way translation, the error against   !
  !   the exact sphere at its final position                         !
  !------------------------------------------------------------------!
  if(Time % Curr_Dt() .eq. Time % Last_Dt()) then
    e1_ex = -1.0
    if(kin_mode .eq. 2 .and. kin_p(4) <= 0.0) then
      t = Time % Get_Time()
      allocate(frac(-Grid % n_bnd_cells:Grid % n_cells))
      call Kin_Sphere_Fraction(Grid, kin_c(1) + kin_p(1) * t,  &
                                     kin_c(2) + kin_p(2) * t,  &
                                     kin_c(3) + kin_p(3) * t,  &
                               kin_rad, kin_init, kin_cyl, frac)
      e1_ex = 0.0
      do c = Cells_In_Domain()
        e1_ex = e1_ex + Grid % vol(c) * abs((1.0 - fun % n(c)) - frac(c))
      end do
      call Global % Sum_Real(e1_ex)
      e1_ex = e1_ex / kin_vol0
      deallocate(frac)
    end if
    if(First_Proc()) then
      open(newunit=fu, file='kin-final.dat', action='write', status='replace')
      write(fu,'(a)') '# time  E1_vs_initial  E1_vs_exact_shifted  '  //  &
                      'dV_rel  min(vof)  max(vof)-1  dW_rel'
      write(fu,'(7(1x,es16.8))') Time % Get_Time(), e1, e1_ex,            &
                                 (b_volume - kin_vol0) / kin_vol0,        &
                                 f_min, f_max - 1.0, wid / kin_w0 - 1.0
      close(fu)
    end if
  end if

  end subroutine
