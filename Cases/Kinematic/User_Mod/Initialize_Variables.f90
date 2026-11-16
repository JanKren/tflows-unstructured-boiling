!==============================================================================!
  subroutine Kin_Sphere_Fraction(Grid, cx, cy, cz, rad, init, cyl, frac)
!------------------------------------------------------------------------------!
!   Volume fraction of a sphere in each cell (1 inside the sphere).            !
!                                                                              !
!   init = 0: the linear estimate (1 - d_min) / (d_max - d_min) on the node    !
!             distances, on every cell.  This is the production Scriven        !
!             initialisation, with abs() on the node count so that it also     !
!             works on polyhedral cells (negative node count).                 !
!   init = 1: axis-aligned hexahedra are sub-sampled on N_SUB**3 points;       !
!             other cells keep the linear estimate.                            !
!   cyl  = 1: a cylinder along z in place of the sphere (a disc in 2D).        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type), target :: Grid
  real,    intent(in)     :: cx, cy, cz, rad
  integer, intent(in)     :: init, cyl
  real                    :: frac(-Grid % n_bnd_cells:Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer, parameter :: N_SUB = 32
  integer            :: c, n, i_nod, i, j, k, n_in
  real               :: d, d_min, d_max, fz
  real               :: x_min, x_max, y_min, y_max, z_min, z_max, xs, ys, zs
!==============================================================================!

  frac(:) = 0.0

  ! Weight of the z distance: none for a cylinder along z
  fz = 1.0
  if(cyl .eq. 1) fz = 0.0

  do c = Cells_In_Domain_And_Buffers()

    d_min = +HUGE
    d_max = -HUGE
    do i_nod = 1, abs(Grid % cells_n_nodes(c))
      n = Grid % cells_n(i_nod, c)
      d = sqrt(  (Grid % xn(n) - cx)**2   &
               + (Grid % yn(n) - cy)**2   &
               + fz * (Grid % zn(n) - cz)**2) / rad
      d_min = min(d_min, d)
      d_max = max(d_max, d)
    end do
    Assert(d_max > -HUGE * 0.5)

    ! Cell is cut by the sphere
    if(d_min < 1.0 .and. d_max > 1.0) then

      if(init .eq. 1 .and. Grid % cells_n_nodes(c) .eq. 8) then
        x_min = +HUGE;  x_max = -HUGE
        y_min = +HUGE;  y_max = -HUGE
        z_min = +HUGE;  z_max = -HUGE
        do i_nod = 1, 8
          n = Grid % cells_n(i_nod, c)
          x_min = min(x_min, Grid % xn(n));  x_max = max(x_max, Grid % xn(n))
          y_min = min(y_min, Grid % yn(n));  y_max = max(y_max, Grid % yn(n))
          z_min = min(z_min, Grid % zn(n));  z_max = max(z_max, Grid % zn(n))
        end do
        n_in = 0
        do k = 1, N_SUB
          zs = z_min + (real(k) - 0.5) * (z_max - z_min) / real(N_SUB)
          do j = 1, N_SUB
            ys = y_min + (real(j) - 0.5) * (y_max - y_min) / real(N_SUB)
            do i = 1, N_SUB
              xs = x_min + (real(i) - 0.5) * (x_max - x_min) / real(N_SUB)
              if((xs-cx)**2 + (ys-cy)**2 + fz*(zs-cz)**2 < rad**2)  &
                n_in = n_in + 1
            end do
          end do
        end do
        frac(c) = real(n_in) / real(N_SUB**3)
      else
        frac(c) = (1.0 - d_min) / (d_max - d_min)
      end if

    ! Cell is entirely inside
    else if(d_max <= 1.0) then
      frac(c) = 1.0
    end if

  end do

  end subroutine

!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, Turb, Vof, Swarm, Sol)
!------------------------------------------------------------------------------!
!   Kinematic advection tests: a sphere of the dispersed phase (vof = 0)       !
!   in the continuous phase (vof = 1).  Velocities stay zero; the face fluxes  !
!   are imposed in User_Mod_Beginning_Of_Compute_Vof.                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: fun
  integer                  :: c, fu, inner, ios
  real, allocatable        :: frac(:)
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  fun  => Vof % fun

  ! Read the parameters of the test
  open(newunit=fu, file='kinematic.in', action='read', status='old')
  read(fu, *) kin_mode
  read(fu, *) kin_c(1), kin_c(2), kin_c(3), kin_rad
  read(fu, *) kin_p(1), kin_p(2), kin_p(3), kin_p(4)
  read(fu, *) inner
  read(fu, *) kin_init
  read(fu, *, iostat=ios) kin_cyl
  if(ios .ne. 0) kin_cyl = 0
  close(fu)

  allocate(frac    (-Grid % n_bnd_cells:Grid % n_cells));  frac(:)     = 0.0
  allocate(kin_fun0(-Grid % n_bnd_cells:Grid % n_cells));  kin_fun0(:) = 1.0

  call Kin_Sphere_Fraction(Grid, kin_c(1), kin_c(2), kin_c(3), kin_rad,  &
                           kin_init, kin_cyl, frac)

  do c = Cells_In_Domain_And_Buffers()
    fun % n(c) = 1.0 - frac(c)
  end do
  call Grid % Exchange_Cells_Real(fun % n)
  fun % o(:) = fun % n(:)

  ! Keep the initial state for the error norms
  kin_fun0(:) = fun % n(:)
  kin_vol0 = 0.0
  kin_w0   = 0.0
  do c = Cells_In_Domain()
    kin_vol0 = kin_vol0 + Grid % vol(c) * (1.0 - fun % n(c))
    if(fun % n(c) > 0.01 .and. fun % n(c) < 0.99) then
      kin_w0 = kin_w0 + Grid % vol(c)
    end if
  end do
  call Global % Sum_Real(kin_vol0)
  call Global % Sum_Real(kin_w0)

  if(First_Proc()) then
    open(newunit=fu, file='kin-init.dat', action='write', status='replace')
    write(fu, '(a)') '# mode init  volume_0  r_vol_0  w_0'
    write(fu, '(2i4,3es20.12)') kin_mode, kin_init, kin_vol0,  &
                                (3.0*kin_vol0/(4.0*PI))**(1.0/3.0), kin_w0
    close(fu)
  end if

  deallocate(frac)

  end subroutine
