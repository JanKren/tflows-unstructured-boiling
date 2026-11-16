!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   This function is computing benchmark for rising bubble.                    !
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
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: fun
  type(Front_Type), pointer:: Front
  type(Elem_Type), pointer :: Elem(:)
  integer                  :: c, e, fu
  real                     :: b_volume, surface, ne, rad_volume
  integer                  :: c1, c2, s
  real                     :: m_dom, out_rate, m_s
  real,               save :: out_cum  = 0.0    ! running egress [kg]
  real,               save :: out_prev = 0.0    ! previous step's egress rate
  real,               save :: m_zero   = 0.0    ! M_0, the invariant M_s
  logical,            save :: first    = .true.
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  fun  => Vof % fun
  Front => Vof % Front
  Elem => Front % Elem

  !-------------------!
  !   Bubble volume   !
  !-------------------!
  b_volume = 0.0
  surface = 0.0
  ne = 0
!  alpha_l = Vof%phase_cond(1)/(Vof%phase_capa(1)*Vof%phase_dens(1))

! do c = Cells_In_Domain()
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
    b_volume = b_volume + Grid % vol(c) * (1.0 - fun % n(c))
    e  = Vof % Front % elem_in_cell(c)
    if(e > 0) then
      ne = ne + 1
      surface = surface + Elem(e) % area
    end if
  end do

  call Global % Sum_Real(ne)
  call Global % Sum_Real(surface)
  call Global % Sum_Real(b_volume)
  rad_volume = (3.0*b_volume/(4*PI))**(1.0/3.0)

  ! Write to file

  !=========================================================================!
  !   Global mass budget, after Sato & Nicheno (JCP 2013, Fig. 8).           !
  !   M_d = sum rho V ; M_s = M_d + sum_t sum_outlet rho u.S dt              !
  !=========================================================================!
  m_dom = 0.0
  do c = 1, Grid % n_cells - Grid % Comm % n_buff_cells
    m_dom = m_dom + Flow % density(c) * Grid % vol(c)
  end do
  call Global % Sum_Real(m_dom)

  out_rate = 0.0
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 < 0) then
      if(Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW   .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. PRESSURE  .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. CONVECT) then
        out_rate = out_rate + Flow % density(c1) * Flow % v_flux % n(s)
      end if
    end if
  end do
  call Global % Sum_Real(out_rate)

  ! trapezoidal: Euler on this integral is first order and would dominate
  if(first) out_prev = out_rate
  out_cum  = out_cum + 0.5 * (out_rate + out_prev) * Flow % dt
  out_prev = out_rate
  m_s      = m_dom + out_cum

  ! M_0 is the invariant M_s, not M_d: this hook first runs at the END of
  ! step one, when one step's egress has already occurred
  if(first) then
    m_zero = m_s
    first  = .false.
  end if

  if(First_Proc()) then
    call File % Append_For_Writing_Ascii('mass-balance.dat', fu)
    write(fu,'(5(2X,E20.12E2))') Time % Get_Time(), m_dom, m_s,  &
                                 m_dom / m_zero, m_s / m_zero - 1.0
    close(fu)
  end if

  if(First_Proc()) then
    call File % Append_For_Writing_Ascii('bench-data.dat', fu)
    write(fu,'(4(2X,E16.8))') Time % Get_Time(), surface, b_volume, rad_volume
    close(fu)
  end if

  end subroutine
