!==============================================================================!
  subroutine Mass_Transfer_Estimate(Vof)
!------------------------------------------------------------------------------!
!   Calculates pressure source due to phase change                             !
!                                                                              !
!   Called from Multiphase_Mod_Vof_Pressure_Correction                         !
!                                                                              !
!   Fills the interfacial mass transfer rate Vof % m_dot from the temperature  !
!   gradients extrapolated from each phase, and the interface conductances     !
!   Vof % a12 and Vof % a21.                                                   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Field_Type),    pointer :: Flow
  type(Front_Type),    pointer :: Front
  type(Var_Type),      pointer :: t
  real,    contiguous, pointer :: elem_sx(:), elem_sy(:), elem_sz(:)
  integer, contiguous, pointer :: elem_used(:)
  character(SL)                :: fname
  integer                      :: e, s, c, c1, c2, elem
  real                         :: phi_c1, phi_c2
  real                         :: dx1, dy1, dz1, dx2, dy2, dz2
  real                         :: d1, d2, st, st1, st2, face_s
  real                         :: sx1, sx2, sy1, sy2, sz1, sz2
  real                         :: nx_f, ny_f, nz_f, n_mag
  real                         :: r1_mag, r2_mag, cos_1, cos_2
  real                         :: q_0, q_1, sum_m_dot
  real                         :: sum_nab_t, sum_elem_ss
  integer                      :: n_floor_d, n_intf_face  ! floor diagnostic
  logical                      :: exists_file
  integer                      :: fu_mass_transfer
  real                         :: uni_flux, uni_m, uni_a
  integer                      :: n_clamped, vapour, liquid
  real                         :: mdot_cap, cfl, c_exp
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
  ! Floor on the cosine between the cell-centre connection and the ...
  ! ... interface normal.  It bounds the near-tangential limit, in which ...
  ! ... the normal distance vanishes and the coefficient would diverge.
  real,    parameter :: D_F_COS_MIN = 0.05
  ! Floor on the sub-cell distance to the interface, as a fraction of the ...
  ! ... local cell size.  Matches PSI-BOIL's Topology::close_to_cc = 1.0e-2.
  ! D_F_COS_MIN bounds the near-tangential limit; this bounds the case of ...
  ! ... the interface sweeping through a cell centre, where r itself vanishes.
  ! The two are independent and both are needed.
  real,    parameter :: CLOSE_TO_CC = 1.0e-2
  real,    parameter :: ONE_THIRD   = 1.0 / 3.0
!==============================================================================!

  ! Take aliases
  Grid  => Vof % pnt_grid
  Flow  => Vof % pnt_flow
  Front => Vof % Front
  t     => Flow % t

  ! If not a problem with mass transfer, get out of here
  if(Flow % mass_transfer_model == NO_MASS_TRANSFER) return

  if(Math % Approx_Real(Vof % latent_heat, 1.0) ) then
    if(First_Proc()) then
      print '(a,1pe12.4)', ' # Latent heat is: ', Vof % latent_heat
      print '(a)', 'You might forget to define LATENT_HEAT in control'
    end if
    call Message % Error(60,                                       &
                         'This error is critical.  Exiting.',      &
                         file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  call Work % Connect_Real_Cell(elem_sx, elem_sy, elem_sz)
  call Work % Connect_Int_Cell(elem_used)

  ! Initialize mass transfer term
  Vof % m_dot(:) = 0.0
  Vof % q_0(:) = 0.0
  Vof % q_1(:) = 0.0
  sum_elem_ss = 0.0

  !------------------------------------------------!
  !                                                !
  !   Compute gradients of temperature, imposing   !
  !    saturation temperature at the interface     !
  !                                                !
  !------------------------------------------------!
  call Vof % Grad_Variable_With_Front(t, Vof % t_sat)
  if(DEBUG) then
    call Grid % Save_Debug_Vtu("grad-t",                               &
                               scalar_cell = t % n,                    &
                               scalar_name = "t",                      &
                               vector_cell = (/t % x, t % y, t % z/),  &
                               vector_name = "grad-t")
  end if

  !----------------------------------------------!
  !                                              !
  !    Extrapolate gradients from each of the    !
  !   phases toward the cells at the interface   !
  !                                              !
  !----------------------------------------------!
  if(Flow % mass_transfer_model .eq. TEMPERATURE_GRADIENTS) then

    ! Initialize t_0 and t_1 ...
    Vof % t_0 % x(1:) = t % x(1:)
    Vof % t_0 % y(1:) = t % y(1:)
    Vof % t_0 % z(1:) = t % z(1:)
    Vof % t_1 % x(1:) = t % x(1:)
    Vof % t_1 % y(1:) = t % y(1:)
    Vof % t_1 % z(1:) = t % z(1:)

    ! ... and extrapolate to front
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_1 % x, towards=0)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_1 % y, towards=0)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_1 % z, towards=0)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_0 % x, towards=1)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_0 % y, towards=1)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_0 % z, towards=1)
  end if

  if(DEBUG) then
    fname = "grad-t-"
    write(fname(8:12), '(i5.5)') Time % Curr_Dt()
    call Grid % Save_Debug_Vtu(fname,                                  &
                               scalar_cell = t % n,                    &
                               scalar_name = "t",                      &
                               vector_cell = (/t % x, t % y, t % z/),  &
                               vector_name = "grad-t")
    fname = "grad-t0-"
    write(fname(8:12), '(i5.5)') Time % Curr_Dt()
    call Grid % Save_Debug_Vtu(fname,                                  &
                               scalar_cell = Vof % t_0 % n,            &
                               scalar_name = "t_0",                    &
                               vector_cell = (/Vof % t_0 % x,          &
                                               Vof % t_0 % y,          &
                                               Vof % t_0 % z/),        &
                               vector_name = "grad-t0")
    fname = "gra-t1-"
    write(fname(8:12), '(i5.5)') Time % Curr_Dt()
    call Grid % Save_Debug_Vtu(fname,                                  &
                               scalar_cell = Vof % t_1 % n,            &
                               scalar_name = "t_1",                    &
                               vector_cell = (/Vof % t_1 % x,          &
                                               Vof % t_1 % y,          &
                                               Vof % t_1 % z/),        &
                               vector_name = "grad-t1")
  end if

  !-------------------------------------------!
  !                                           !
  !   Computation of heat fluxes to or from   !
  !   interface towards the cells around it   !
  !                                           !
  !-------------------------------------------!

  !---------------------------------------------------!
  !   Distribute element areas to cell-based arrays   !
  !---------------------------------------------------!
  elem_sx(:) = 0.0
  elem_sy(:) = 0.0
  elem_sz(:) = 0.0
  do e = 1, Front % n_elems
    c = Front % Elem(e) % cell      ! in which cell it resides
    elem_sx(c) = Front % Elem(e) % sx
    elem_sy(c) = Front % Elem(e) % sy
    elem_sz(c) = Front % Elem(e) % sz
  end do
  call Grid % Exchange_Cells_Real(  &
         elem_sx(-Grid % n_bnd_cells:Grid % n_cells))
  call Grid % Exchange_Cells_Real(  &
         elem_sy(-Grid % n_bnd_cells:Grid % n_cells))
  call Grid % Exchange_Cells_Real(  &
         elem_sz(-Grid % n_bnd_cells:Grid % n_cells))

  if(DEBUG) then
    fname = "elem-s-"
    write(fname(8:12), '(i5.5)') Time % Curr_Dt()
    call Grid % Save_Debug_Vtu(fname,                      &
                               vector_cell = (/elem_sx,    &
                                               elem_sy,    &
                                               elem_sz/),  &
                               vector_name = "elem_s")
  end if

  !----------------------------------------------------!
  !   Count how many times will each element be used   !
  !----------------------------------------------------!
  elem_used(:) = 0
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 < 0) cycle  ! the front does not cross boundary faces

    phi_c1 = Vof % fun % n(c1)
    phi_c2 = Vof % fun % n(c2)

    ! If face crosses the 0.5 value, increase the ...
    ! ... number of times the element has been used
    if( (phi_c1 - 0.5) * (phi_c2 - 0.5) < 0 ) then
      elem_used(c1) = elem_used(c1) + 1
      elem_used(c2) = elem_used(c2) + 1
    end if
  end do
  call Grid % Exchange_Cells_Int(  &
         elem_used(-Grid % n_bnd_cells:Grid % n_cells))

  !------------------------------------------------------!
  !   Divide elements' surfaces by the number of times   !
  !     they are used to avoid excessive integration     !
  !------------------------------------------------------!
  !do c = Cells_In_Domain_And_Buffers()
  !  if(elem_used(c) > 0) then
  !    elem_sx(c) = elem_sx(c) / elem_used(c)
  !    elem_sy(c) = elem_sy(c) / elem_used(c)
  !    elem_sz(c) = elem_sz(c) / elem_used(c)
  !  end if
  !end do

  !----------------------------------------!
  !   Compute heat flux at the interface   !
  !----------------------------------------!
  Vof % a12(:) = 0.0
  Vof % a21(:) = 0.0
  n_floor_d   = 0
  n_intf_face = 0

  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 < 0) cycle  ! the front does not cross boundary faces

    phi_c1 = Vof % fun % n(c1)
    phi_c2 = Vof % fun % n(c2)

    Assert(phi_c1 .ne. 0.5)
    Assert(phi_c2 .ne. 0.5)

    if( (phi_c1 - 0.5) * (phi_c2 - 0.5) < 0 ) then

      ! Compute internal boundary coefficients on each side of the interface
      ! Units: W/(mK) * K/m * m^2 = W
      face_s = sqrt((Grid % sx(s))**2+(Grid % sy(s))**2+(Grid % sz(s))**2)
      if(face_s < 1.0e-30) cycle
      sx1 = (elem_sx(c1) * Grid % sx(s)) / face_s
      sy1 = (elem_sy(c1) * Grid % sy(s)) / face_s
      sz1 = (elem_sz(c1) * Grid % sz(s)) / face_s
      sx2 = (elem_sx(c2) * Grid % sx(s)) / face_s
      sy2 = (elem_sy(c2) * Grid % sy(s)) / face_s
      sz2 = (elem_sz(c2) * Grid % sz(s)) / face_s

      st1 = sqrt(sx1**2+sy1**2+sz1**2)
      st2 = sqrt(sx2**2+sy2**2+sz2**2)

      ! Simply sum the areas up.  Keep in mind that these
      ! have been corrected by the number of visits
      !sx = elem_sx(c1) + elem_sx(c2)
      !sy = elem_sy(c1) + elem_sy(c2)
      !sz = elem_sz(c1) + elem_sz(c2)
      !st = sqrt(sx**2 + sy**2 + sz**2)

      !st1 = sqrt(elem_sx(c1)**2+elem_sy(c1)**2+elem_sz(c1)**2)
      !st2 = sqrt(elem_sx(c2)**2+elem_sy(c2)**2+elem_sz(c2)**2)
      st = st1+st2

      dx1 = Grid % xc(c1) - Front % xs(s)
      dy1 = Grid % yc(c1) - Front % ys(s)
      dz1 = Grid % zc(c1) - Front % zs(s)

      dx2 = Grid % xc(c1) + Grid % dx(s) - Front % xs(s)
      dy2 = Grid % yc(c1) + Grid % dy(s) - Front % ys(s)
      dz2 = Grid % zc(c1) + Grid % dz(s) - Front % zs(s)

      ! Interface normal at the face.  The mean of two unit normals ...
      ! ... is not itself a unit vector, hence the renormalisation.
      nx_f = 0.5 * (Vof % nx(c1) + Vof % nx(c2))
      ny_f = 0.5 * (Vof % ny(c1) + Vof % ny(c2))
      nz_f = 0.5 * (Vof % nz(c1) + Vof % nz(c2))
      n_mag = norm2((/nx_f, ny_f, nz_f/))
      if(n_mag < TINY) cycle
      nx_f = nx_f / n_mag
      ny_f = ny_f / n_mag
      nz_f = nz_f / n_mag

      ! Normal distance from cell centre to interface, as the ...
      ! ... projection of the connection vector onto the normal:
      !
      !   d = |r| * max(|r_hat . n|, D_F_COS_MIN)
      !
      ! The projection is invariant under rotation of the frame.  The ...
      ! ... earlier form, norm2 of the componentwise product r*n, is ...
      ! ... not: for r = n = (1,1,0)/sqrt(2) it gives 0.707, and 1.000 ...
      ! ... for the same configuration turned 45 degrees about z.  The ...
      ! ... two agree exactly when r is axis-aligned, which is every ...
      ! ... face of a Cartesian grid, so only polyhedral results move.
      r1_mag = norm2((/dx1, dy1, dz1/))
      r2_mag = norm2((/dx2, dy2, dz2/))
      cos_1 = abs(dx1*nx_f + dy1*ny_f + dz1*nz_f) / max(r1_mag, TINY)
      cos_2 = abs(dx2*nx_f + dy2*ny_f + dz2*nz_f) / max(r2_mag, TINY)
      d1 = r1_mag * max(cos_1, D_F_COS_MIN)
      d2 = r2_mag * max(cos_2, D_F_COS_MIN)

      ! An assert that the two connection vectors oppose was dropped ...
      ! ... here: it failed on imperfections in periodic solutions.

      ! Guard against degenerate distances
      if(d1 < 1.0e-30 .or. d2 < 1.0e-30) cycle

      !--------------------------------------------------------------------!
      !   Floor the sub-cell distance at a fraction of the cell.           !
      !                                                                    !
      !   The conductances below are lambda*area/d.  As the interface      !
      !   sweeps across a cell centre d passes through arbitrarily small   !
      !   values, and nothing above stops it: the 1e-30 test catches an    !
      !   exact zero and nothing else.  On a 25 um cell a d of a           !
      !   nanometre inflates the conductance ~25000-fold, and that         !
      !   coefficient goes straight into the implicit energy matrix as a   !
      !   Robin term (Ac += a12, F += a12*t_sat).  An unbounded            !
      !   coefficient there does not stay well behaved in a segregated,    !
      !   under-relaxed solve - it throws the cell's temperature, cold     !
      !   or hot, and the excursions grow as the bubble does and more      !
      !   faces sweep more centres.                                        !
      !                                                                    !
      !   PSI-BOIL floors the same distance, at 1% of the cell size:       !
      !   Topology::distance_int_x returns                                 !
      !     std::max(dist, close_to_cc*clrold.dxc(i))                      !
      !   with close_to_cc = 1.0e-2 (Topology::topology.cpp).  Same        !
      !   fraction is used here; the resulting conductance is capped at    !
      !   100*lambda*area/dx, which is far above anything physical and     !
      !   far below anything that destabilises the matrix.                 !
      !--------------------------------------------------------------------!
      ! Diagnostic: the faces at which the floor acts.
      n_intf_face = n_intf_face + 1
      if(d1 < CLOSE_TO_CC * Grid % vol(c1)**ONE_THIRD .or.  &
         d2 < CLOSE_TO_CC * Grid % vol(c2)**ONE_THIRD)      &
        n_floor_d = n_floor_d + 1
      d1 = max(d1, CLOSE_TO_CC * Grid % vol(c1)**ONE_THIRD)
      d2 = max(d2, CLOSE_TO_CC * Grid % vol(c2)**ONE_THIRD)

      ! Cell c1 is in vapor and c2 is in liquid
      if(Vof % fun % n(c1) < 0.5) then
        Vof % a12(s) = Vof % a12(s) + Vof % phase_cond(0) * st / d1
        Vof % a21(s) = Vof % a21(s) + Vof % phase_cond(1) * st / d2
      ! Cell c1 is in liquid and c2 is in vapor
      else
        Vof % a12(s) = Vof % a12(s) + Vof % phase_cond(1) * st / d1
        Vof % a21(s) = Vof % a21(s) + Vof % phase_cond(0) * st / d2
      end if

    end if      ! face is at the front

  end do

  !-----------------------------------------------!
  !                                               !
  !   Mass transfer with extrapolated gradients   !
  !                                               !
  !-----------------------------------------------!
  sum_m_dot = 0.0

  if(Flow % mass_transfer_model .eq. TEMPERATURE_GRADIENTS) then
    do c = Cells_In_Domain_And_Buffers()
      if(elem_used(c) > 0) then
        e  = Front % elem_in_cell(c)
        if(e > 0) then
          ! Units: W/(mK) * K/m * m^2 = W
          q_0 = (  Vof % t_0 % x(c) * Front % Elem(e) % sx    &
                 + Vof % t_0 % y(c) * Front % Elem(e) % sy    &
                 + Vof % t_0 % z(c) * Front % Elem(e) % sz )  &
                 * Vof % phase_cond(0)
          ! Heat flux to the interface in cell c from phase 1
          ! Units: W/(mK) * K/m * m^2 = W
          q_1 = (  Vof % t_1 % x(c) * Front % Elem(e) % sx    &
                 + Vof % t_1 % y(c) * Front % Elem(e) % sy    &
                 + Vof % t_1 % z(c) * Front % Elem(e) % sz )  &
                 * Vof % phase_cond(1)
          st = sqrt(  (Front % Elem(e) % sx)**2    &
                    + (Front % Elem(e) % sy)**2    &
                    + (Front % Elem(e) % sz)**2)
          sum_elem_ss = sum_elem_ss + st
          ! Store q_0 and q_1 in Vof_Type
          Vof % q_0(c) = q_0
          Vof % q_1(c) = q_1

          Vof % m_dot(c) = (q_1-q_0) / Vof % latent_heat
          sum_m_dot = sum_m_dot + Vof % m_dot(c)
        end if
      end if
    end do

    ! Detailed diagnostics
    block
      real    :: sum_q0g, sum_q1g, sum_area_g
      real    :: sum_gradt0g, sum_gradt1g, sum_mdot_g
      integer :: n_mdot_g
      sum_q0g = 0.0
      sum_q1g = 0.0
      sum_area_g = 0.0
      sum_gradt0g = 0.0
      sum_gradt1g = 0.0
      sum_mdot_g = 0.0
      n_mdot_g = 0
      do c = Cells_In_Domain()
        if(elem_used(c) > 0) then
          e = Front % elem_in_cell(c)
          if(e > 0) then
            st = sqrt(Front % Elem(e) % sx**2  &
                    + Front % Elem(e) % sy**2  &
                    + Front % Elem(e) % sz**2)
            sum_area_g = sum_area_g + st
            sum_q0g = sum_q0g + Vof % q_0(c)
            sum_q1g = sum_q1g + Vof % q_1(c)
            sum_gradt0g = sum_gradt0g + sqrt(  Vof % t_0 % x(c)**2  &
                                             + Vof % t_0 % y(c)**2  &
                                             + Vof % t_0 % z(c)**2)
            sum_gradt1g = sum_gradt1g + sqrt(  Vof % t_1 % x(c)**2  &
                                             + Vof % t_1 % y(c)**2  &
                                             + Vof % t_1 % z(c)**2)
            sum_mdot_g = sum_mdot_g + Vof % m_dot(c)
            n_mdot_g = n_mdot_g + 1
          end if
        end if
      end do
      call Global % Sum_Real(sum_q0g)
      call Global % Sum_Real(sum_q1g)
      call Global % Sum_Real(sum_area_g)
      call Global % Sum_Real(sum_gradt0g)
      call Global % Sum_Real(sum_gradt1g)
      call Global % Sum_Real(sum_mdot_g)
      call Global % Sum_Int(n_mdot_g)
      if(First_Proc() .and. Iter % Current() == 1) then
        print '(a,e12.4,a,e12.4,a,e12.4,a,i0)',  &
          ' # AVOF_CMP sum_q0=', sum_q0g,        &
          ' sum_q1=', sum_q1g,                   &
          ' S_dom=', sum_area_g,                 &
          ' n_cells=', n_mdot_g
        print '(a,e12.4,a,e12.4,a,e12.4)',       &
          ' # AVOF_CMP sum_mdot=', sum_mdot_g,   &
          ' avg_gradt0=',                        &
          sum_gradt0g/max(n_mdot_g,1),           &
          ' avg_gradt1=',                        &
          sum_gradt1g/max(n_mdot_g,1)
      end if
    end block
  end if

  if(DEBUG) then
    fname = "m_dot0-"
    write(fname(8:12), '(i5.5)') Time % Curr_Dt()
    call Grid % Save_Debug_Vtu(fname,                      &
                               scalar_cell = Vof % m_dot0, &
                               scalar_name = "m_dot0")
    fname = "elem-s-"
    write(fname(8:12), '(i5.5)') Time % Curr_Dt()
    call Grid % Save_Debug_Vtu(fname,                      &
                               vector_cell = (/elem_sx,    &
                                               elem_sy,    &
                                               elem_sz/),  &
                               vector_name = "elem_s")
  end if

  ! Sum
  sum_nab_t = 0.0
  elem = 0
  do c = Cells_In_Domain()
    e  = Front % elem_in_cell(c)
    if(e > 0) then
      sum_nab_t = sum_nab_t + norm2((/t % x(c), t % y(c), t % z(c)/))
      elem = elem + 1
    end if
  end do
  call Global % Sum_Real(sum_m_dot)
  call Global % Sum_Real(sum_nab_t)
  call Global % Sum_Real(sum_elem_ss)

  !---------------------------------------------------------------------!
  !   Diagnostic: impose a spatially UNIFORM interfacial mass flux.     !
  !                                                                     !
  !   The computed m_dot is rescaled so that every interface element    !
  !   receives the same mass per unit area, while the domain total is   !
  !   left unchanged.  The bubble therefore grows at exactly the rate   !
  !   it would otherwise, but with no angular variation whatsoever in   !
  !   the source.  Any deformation that survives is attributable to     !
  !   the transport of the volume fraction alone.                       !
  !---------------------------------------------------------------------!
  if(Vof % uniform_m_dot) then

    ! The running sums above accumulate over domain AND buffer cells, so ...
    ! ... an interface cell on a partition boundary is counted twice once ...
    ! ... the global reduction is applied.  That is harmless for a ...
    ! ... diagnostic print but not for a normalisation: it would leave the ...
    ! ... domain total slightly altered, which is the one property this ...
    ! ... test relies on.  Re-sum over owned cells only, then reduce.
    uni_m = 0.0
    uni_a = 0.0
    do c = Cells_In_Domain()
      if(elem_used(c) > 0) then
        e = Front % elem_in_cell(c)
        if(e > 0) then
          st = norm2((/Front % Elem(e) % sx,   &
                       Front % Elem(e) % sy,   &
                       Front % Elem(e) % sz/))
          uni_m = uni_m + Vof % m_dot(c)
          uni_a = uni_a + st
        end if
      end if
    end do
    call Global % Sum_Real(uni_m)
    call Global % Sum_Real(uni_a)

    if(uni_a > TINY) then
      uni_flux = uni_m / uni_a                    ! [kg / (m^2 s)]

      ! Assignment covers buffers too, so the halo stays consistent
      do c = Cells_In_Domain_And_Buffers()
        if(elem_used(c) > 0) then
          e = Front % elem_in_cell(c)
          if(e > 0) then
            st = norm2((/Front % Elem(e) % sx,   &
                         Front % Elem(e) % sy,   &
                         Front % Elem(e) % sz/))
            Vof % m_dot(c) = uni_flux * st
          end if
        end if
      end do
    end if
  end if

  call Global % Sum_Int(n_floor_d)
  call Global % Sum_Int(n_intf_face)
  if(First_Proc()) then
    print '(a,i0,a,i0,a,f7.4)',                  &
      ' # ALG d_floor_hits=', n_floor_d,         &
      '  of_intf_faces=', n_intf_face,           &
      '  frac=', real(n_floor_d)/max(real(n_intf_face),1.0)
  end if

  if(DEBUG) then

    ! Check if the current process is the first one
    if(First_Proc()) then
      ! Inquire if 'statistics.dat' exists
      inquire(file='statistics.dat', exist=exists_file)

      ! If the file does not exist, create it and write headers
      if(.not. exists_file) then
        call File % Append_For_Writing_Ascii('statistics.dat',  &
                                             fu_mass_transfer)
        write(fu_mass_transfer, '(A, A, A, A)')    &
              "Time", "sum_m_dot", "sum_nab_t", "sum_elem_ss"
        close(fu_mass_transfer)  ! close using intrinsic close
      end if
    end if

    ! Append the mass transfer data to 'statistics.dat'
    if(First_Proc()) then
      call File % Append_For_Writing_Ascii('statistics.dat', fu_mass_transfer)
      write(fu_mass_transfer, '(4(2X, E16.8))')    &
            Time % Get_Time(), sum_m_dot, sum_nab_t, sum_elem_ss
      close(fu_mass_transfer)  ! close using intrinsic close
    end if
  end if

  !----------------------------------------------------------------------!
  !   Crash guard: cap the per-cell Stefan source so its equivalent      !
  !   dilatation Courant number cannot exceed MASS_TRANSFER_CFL.  The    !
  !   macro (q_1-q_0)/L source is otherwise unbounded, and a             !
  !   freshly-flipped near-wall front cell can drive an unclamped        !
  !   within-SIMPLE feedback - m_dot raises the pressure volume          !
  !   source, which sharpens the near-wall velocity, which sharpens      !
  !   the gradient, which raises m_dot - that detonates inside a         !
  !   single time step.  Bounding m_dot at the source node caps the      !
  !   feedback gain.                                                     !
  !                                                                      !
  !   This is a STABILITY guard, not a physics fix.  While it bites,     !
  !   the clamped cells violate m_dot*L = (q_1-q_0) because the energy   !
  !   Robin term is not clamped, so a non-zero clamp count marks a run   !
  !   as diagnostic-only.  The cap scales as vol/dt, so the CFL value    !
  !   is resolution dependent.                                           !
  !----------------------------------------------------------------------!
  call Vof % Get_Vapour_And_Liquid_Phase(vapour, liquid)
  c_exp = abs(1.0/Vof % phase_dens(vapour) - 1.0/Vof % phase_dens(liquid))
  cfl   = Flow % mass_transfer_cfl
  n_clamped = 0
  if(c_exp > TINY .and. Flow % dt > 0.0 .and. cfl > 0.0) then
    do c = Cells_In_Domain_And_Buffers()
      if(Front % elem_in_cell(c) > 0) then
        mdot_cap = cfl * Grid % vol(c) / (Flow % dt * c_exp)
        if(abs(Vof % m_dot(c)) > mdot_cap) then
          Vof % m_dot(c) = sign(mdot_cap, Vof % m_dot(c))
          n_clamped = n_clamped + 1
        end if
      end if
    end do
  end if
  call Global % Sum_Int(n_clamped)
  if(First_Proc() .and. n_clamped > 0 .and. Iter % Current() .eq. 1) then
    print '(a,i0,a,1pe9.2)',                                       &
      ' # Mass transfer: m_dot CLAMPED (over-feed guard) in ',     &
      n_clamped, ' cell(s); CFL=', cfl
  end if

  ! Sum
  sum_m_dot = 0.0
  sum_nab_t = 0.0
  elem = 0
  do c = Cells_In_Domain()
    e  = Front % elem_in_cell(c)
    if(e > 0) then
      sum_m_dot = sum_m_dot + Vof % m_dot(c)
      sum_nab_t = sum_nab_t + norm2((/t % x(c), t % y(c), t % z(c)/))
      elem = elem + 1
    end if
  end do
  call Global % Sum_Real(sum_m_dot)
  call Global % Sum_Real(sum_nab_t)
  call Global % Sum_Int(elem)

  if(DEBUG) then
    fname = "m_dot-"
    write(fname(7:11), '(i5.5)') Time % Curr_Dt()
    call Grid % Save_Debug_Vtu(fname,                            &
                               scalar_cell = Vof % m_dot,        &
                               scalar_name = "vof_m_dot")
  end if

  call Work % Disconnect_Real_Cell(elem_sx, elem_sy, elem_sz)
  call Work % Disconnect_Int_Cell(elem_used)

  end subroutine
