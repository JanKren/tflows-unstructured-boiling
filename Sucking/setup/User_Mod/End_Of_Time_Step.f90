!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   CAF review instrumentation: every term of the global mass and energy       !
!   budgets, written once per time step to `caf-terms.dat`.  Rates are         !
!   written, not integrals, so that any quadrature can be applied afterwards.  !
!                                                                              !
!   Columns                                                                    !
!     1  time step                                                             !
!     2  time [s]                                                              !
!     3  M_d       mass in the domain, sum rho V                    [kg]       !
!     4  V_vap     vapour volume, sum (1 - alpha) V                 [m^3]      !
!     5  out_rate  mass flux through the outflow, rho(c1) * v_flux  [kg/s]     !
!     6  Q_out     volume flux through the outflow                  [m^3/s]    !
!     7  mdot_end  sum of m_dot, last estimate of the step          [kg/s]     !
!     8  mdot_vof  sum of m_dot the VOF source of this step used    [kg/s]     !
!     9  Q_wall    heat entering through the Dirichlet wall         [W]        !
!    10  Q_int_v   heat removed by the interface coupling, vapour   [W]        !
!    11  Q_int_l   heat removed by the interface coupling, liquid   [W]        !
!    12  H_v       sensible enthalpy above T_sat, cells alpha < 0.5 [J]        !
!    13  H_l       sensible enthalpy above T_sat, cells alpha > 0.5 [J]        !
!    14  x_front   mean x of the face intersections                 [m]        !
!    15  n_cut     number of faces the front intersects                        !
!    16  div_max   largest |sum of face fluxes| over the cells      [m^3/s]    !
!    17  dev_max   largest |T - T_sat| over the cells               [K]        !
!    18  T_min     smallest cell temperature                        [K]        !
!    19  q_0       sum of k_v grad T_0 . S over the front elements  [W]        !
!    20  q_1       sum of k_l grad T_1 . S over the front elements  [W]        !
!              (m_dot = (q_1 - q_0) / h_lv, last estimate of the step)         !
!                                                                              !
!   Q_int is the coupling of the energy equation, a_f (T_P - T_sat) summed     !
!   over the faces the front intersects, with the coefficients and the         !
!   temperature of the last outer iteration, which is what the last energy     !
!   solve removed from the cells.                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer,      intent(in) :: n_stat_t  ! start time step for Turb. stat.
  integer,      intent(in) :: n_stat_p  ! start time step for Swarm. stat.
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c, c1, c2, s, fu, n_cut
  real                     :: m_dom, v_vap, out_rate, q_out, mdot_end
  real                     :: q_wall, q_int_v, q_int_l, h_v, h_l, x_front
  real                     :: fc, q1, q2, div_max, dev_max, t_min, th
  real                     :: sum_q0, sum_q1
  real, allocatable        :: div(:)
!==============================================================================!

  Grid => Flow % pnt_grid

  !------------------------------------------------!
  !   Mass, vapour volume and sensible enthalpy    !
  !------------------------------------------------!
  m_dom    = 0.0
  v_vap    = 0.0
  mdot_end = 0.0
  sum_q0   = 0.0
  sum_q1   = 0.0
  h_v      = 0.0
  h_l      = 0.0
  dev_max  = 0.0
  t_min    = HUGE
  do c = Cells_In_Domain()
    m_dom    = m_dom + Flow % density(c) * Grid % vol(c)
    v_vap    = v_vap + (1.0 - Vof % fun % n(c)) * Grid % vol(c)
    mdot_end = mdot_end + Vof % m_dot(c)
    sum_q0   = sum_q0   + Vof % q_0(c)
    sum_q1   = sum_q1   + Vof % q_1(c)
    th = Flow % t % n(c) - Vof % t_sat
    if(Vof % fun % n(c) < 0.5) then
      h_v = h_v + Vof % phase_dens(0) * Vof % phase_capa(0) * th * Grid % vol(c)
    else
      h_l = h_l + Vof % phase_dens(1) * Vof % phase_capa(1) * th * Grid % vol(c)
    end if
    dev_max = max(dev_max, abs(th))
    t_min   = min(t_min, Flow % t % n(c))
  end do
  call Global % Sum_Real(m_dom)
  call Global % Sum_Real(v_vap)
  call Global % Sum_Real(mdot_end)
  call Global % Sum_Real(sum_q0)
  call Global % Sum_Real(sum_q1)
  call Global % Sum_Real(h_v)
  call Global % Sum_Real(h_l)
  call Global % Max_Real(dev_max)
  call Global % Min_Real(t_min)

  !---------------------------------------------------------!
  !   Boundary fluxes: outflow of mass, wall heat, and the  !
  !   divergence of the face fluxes in every cell           !
  !---------------------------------------------------------!
  allocate(div(-Grid % n_bnd_cells:Grid % n_cells));  div(:) = 0.0
  out_rate = 0.0
  q_out    = 0.0
  q_wall   = 0.0
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    div(c1) = div(c1) + Flow % v_flux % n(s)
    if(c2 > 0) div(c2) = div(c2) - Flow % v_flux % n(s)
    if(c2 < 0) then
      if(Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW   .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. PRESSURE  .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. CONVECT) then
        out_rate = out_rate + Flow % density(c1) * Flow % v_flux % n(s)
        q_out    = q_out    + Flow % v_flux % n(s)
      end if
      ! The heated wall is the one Dirichlet boundary of the temperature
      ! (an inflow at zero velocity in the paper's control).  The flux is
      ! the implicit two-point one the solver puts in the matrix.
      if(Grid % Bnd_Cond_Type(c2) .eq. INFLOW  .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. WALL) then
        fc = (  Grid % sx(s)**2 + Grid % sy(s)**2 + Grid % sz(s)**2 )   &
           / (  Grid % dx(s) * Grid % sx(s)                             &
              + Grid % dy(s) * Grid % sy(s)                             &
              + Grid % dz(s) * Grid % sz(s) )
        if(Vof % fun % n(c1) < 0.5) then
          q_wall = q_wall + Vof % phase_cond(0) * fc  &
                          * (Flow % t % n(c2) - Flow % t % n(c1))
        else
          q_wall = q_wall + Vof % phase_cond(1) * fc  &
                          * (Flow % t % n(c2) - Flow % t % n(c1))
        end if
      end if
    end if
  end do
  call Global % Sum_Real(out_rate)
  call Global % Sum_Real(q_out)
  call Global % Sum_Real(q_wall)
  div_max = 0.0
  do c = Cells_In_Domain()
    div_max = max(div_max, abs(div(c)))
  end do
  call Global % Max_Real(div_max)
  deallocate(div)

  !---------------------------------------------------------!
  !   Heat the interface coupling removed from the cells    !
  !---------------------------------------------------------!
  q_int_v = 0.0
  q_int_l = 0.0
  x_front = 0.0
  n_cut   = 0
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 > 0) then
      if(Vof % Front % intersects_face(s)) then
        q1 = Vof % a12(s) * (Flow % t % n(c1) - Vof % t_sat)
        q2 = Vof % a21(s) * (Flow % t % n(c2) - Vof % t_sat)
        if(Vof % fun % n(c1) < 0.5) then
          q_int_v = q_int_v + q1
        else
          q_int_l = q_int_l + q1
        end if
        if(Vof % fun % n(c2) < 0.5) then
          q_int_v = q_int_v + q2
        else
          q_int_l = q_int_l + q2
        end if
        x_front = x_front + Vof % Front % xs(s)
        n_cut   = n_cut + 1
      end if
    end if
  end do
  call Global % Sum_Real(q_int_v)
  call Global % Sum_Real(q_int_l)
  call Global % Sum_Real(x_front)
  call Global % Sum_Int (n_cut)
  if(n_cut > 0) x_front = x_front / real(n_cut)

  caf_dev_max = max(caf_dev_max, dev_max)

  if(First_Proc()) then
    call File % Append_For_Writing_Ascii('caf-terms.dat', fu)
    write(fu,'(i8,13(1x,es22.14e3),1x,i8,5(1x,es22.14e3))')             &
      Time % Curr_Dt(), Time % Get_Time(), m_dom, v_vap, out_rate,      &
      q_out, mdot_end, caf_mdot_vof, q_wall, q_int_v, q_int_l,          &
      h_v, h_l, x_front, n_cut, div_max, dev_max, t_min, sum_q0, sum_q1
    close(fu)
  end if

  end subroutine
