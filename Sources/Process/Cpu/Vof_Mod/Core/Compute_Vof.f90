!==============================================================================!
  subroutine Compute_Vof(Vof, Sol, dt)
!------------------------------------------------------------------------------!
!   Solves Volume Fraction equation using UPWIND ADVECTION and CICSAM          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type),   target :: Vof
  type(Solver_Type), target :: Sol
  real                      :: dt
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: Flow
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: fun
  type(Face_Type),   pointer :: v_flux
  type(Var_Type),    pointer :: u, v, w
  type(Matrix_Type), pointer :: A
  character(SL)              :: fname
  real, contiguous,  pointer :: b(:)
  real, contiguous,  pointer :: beta_f(:)
  real, contiguous,  pointer :: c_d(:)
  real                       :: phi_c1, phi_c2, fs
  real                       :: u_orig_f, v_orig_f, w_orig_f
  real                       :: courant_max
  integer                    :: i_sub, n_sub, wrong_vf,  &
                                n_wrong_vf0, n_wrong_vf1
  integer                    :: s, c, c1, c2, fu, corr, reg, s_in_reg
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!==============================================================================!

  call Profiler % Start('Compute_Vof (without solvers)')

  call User_Mod_Beginning_Of_Compute_Vof(Vof, Sol)

  ! Take aliases
  Flow   => Vof % pnt_flow
  Grid   => Flow % pnt_grid
  v_flux => Flow % v_flux
  fun    => Vof % fun
  beta_f => Vof % beta_f
  c_d    => Vof % c_d
  A      => Sol % Nat % A
  b      => Sol % Nat % b % val
  u      => Flow % u
  v      => Flow % v
  w      => Flow % w

  ! Velocity flux extrapolation over the interface is disabled:
  ! call Vof % Extrapolate_Normal_To_Front(Flow, v_flux % n , towards=0)
  ! Discretize and Predict_Beta now use Flow % v_flux directly.
  Vof % v_flux1 % n(1:) = v_flux % n(1:)

  if(.false. .and.                                                &
     (Flow % mass_transfer_model .eq. TEMPERATURE_GRADIENTS)) then
    ! Initialize ...
    Vof % v_flux1 % n(1:) = v_flux % n(1:)  ! initialize extrapolated flux
    Vof % u_0 % n(1:) = u % n(1:)
    Vof % v_0 % n(1:) = v % n(1:)
    Vof % w_0 % n(1:) = w % n(1:)

    ! ... and extrapolate to front
    ! call Vof % Extrapolate_Normal_To_Front(Flow, Vof % v_flux1 % n ,  &
    !                                        towards=0)

    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % u_0 % n , towards=0)
    ! DEBUG: Check for NaN after u extrapolation
    do c = 1, Grid % n_cells
      if(Vof % u_0 % n(c) /= Vof % u_0 % n(c)) then
        if(First_Proc()) print '(a,i0,a,i0)',                        &
          ' # [NAN TRACE] u_0 NaN after extrapolation at cell ', c,  &
          ' step ', Time % Curr_Dt()
        exit
      end if
    end do
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % v_0 % n , towards=0)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % w_0 % n , towards=0)

    ! On the boundary perform the extrapolation
    do reg = Boundary_Regions()
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)
        Vof % u_0 % n(c2) = u % n(c2)
        Vof % v_0 % n(c2) = v % n(c2)
        Vof % w_0 % n(c2) = w % n(c2)
        Vof % u_0 % n(c1) = u % n(c1)
        Vof % v_0 % n(c1) = v % n(c1)
        Vof % w_0 % n(c1) = w % n(c1)
      end do  ! faces
    end do    ! region

    do s = Faces_In_Domain_And_At_Buffers()
      c1 = Grid % faces_c(1, s)
      c2 = Grid % faces_c(2, s)
      fs = Grid % f(s)
      phi_c1 = fun % n(c1)
      phi_c2 = fun % n(c2)
      ! Assertions disabled - VOF can reach 0.5 during ...
      ! ... bubble growth
      ! Assert(phi_c1 .ne. 0.5)
      ! Assert(phi_c2 .ne. 0.5)
      if (c2>0) then
        ! Compute v_flux1 from extrapolated velocities
        if (Vof % u_0 % n(c2) .ne. u % n(c2) .or.  &
            Vof % u_0 % n(c1) .ne. u % n(c1) .or.  &
            Vof % v_0 % n(c2) .ne. v % n(c2) .or.  &
            Vof % v_0 % n(c1) .ne. v % n(c1) .or.  &
            Vof % w_0 % n(c2) .ne. w % n(c2) .or.  &
            Vof % w_0 % n(c1) .ne. w % n(c1)) then
          Vof % u_f % n(s) = fs * Vof % u_0 % n(c1)  &
                           + (1.0 - fs) * Vof % u_0 % n(c2)
          Vof % v_f % n(s) = fs * Vof % v_0 % n(c1)  &
                           + (1.0 - fs) * Vof % v_0 % n(c2)
          Vof % w_f % n(s) = fs * Vof % w_0 % n(c1)  &
                           + (1.0 - fs) * Vof % w_0 % n(c2)
          ! Original face velocity (same interpolation as Rhie-Chow)
          u_orig_f = fs * u % n(c1) + (1.0 - fs) * u % n(c2)
          v_orig_f = fs * v % n(c1) + (1.0 - fs) * v % n(c2)
          w_orig_f = fs * w % n(c1) + (1.0 - fs) * w % n(c2)

          ! Correction-based flux preserving Rhie-Chow
          Vof % v_flux1 % n(s) = v_flux % n(s)              &
            + (Vof % u_f % n(s) - u_orig_f) * Grid % sx(s)  &
            + (Vof % v_f % n(s) - v_orig_f) * Grid % sy(s)  &
            + (Vof % w_f % n(s) - w_orig_f) * Grid % sz(s)
        end if
      else
        continue
      end if
    end do

    ! Set boundary face fluxes
    do reg = Boundary_Regions()
      do s = Faces_In_Region(reg)
        Vof % v_flux1 % n(s)  = v_flux % n(s)
      end do  ! faces
    end do    ! region
  else
    Vof % v_flux1 % n = v_flux % n
  end if

  if(fun % adv_scheme .eq. CICSAM .or. &
     fun % adv_scheme .eq. STACS) then

    if(fun % adv_scheme .eq. CICSAM) then
      ! Compute Courant number close to the interface:
      call Vof % Max_Courant_Number(dt, 1, courant_max)

      n_sub = min(max(ceiling(courant_max / Vof % courant_max_param), 1),  &
                  Vof % n_sub_param)

      ! Warning if Courant number is exceeded
      if(n_sub > 1) then
        if(First_Proc()) then
          call File % Append_For_Writing_Ascii('alert-dt-vof.dat', fu)
          write(fu,*) 'Courant Number exceded at iteration: ', Time % Curr_Dt()
          write(fu,*) 'Co_max = ', courant_max
          write(fu,*) 'Try reducing time step'
          close(fu)
        end if
      end if
    else
      n_sub = 1
    end if
  else
    ! Old volume fraction:
    fun % o(:) = fun % n(:)
  end if

  if(DEBUG) then
    fname = "ureal-"
    write(fname(8:12), '(i5.5)') Time % Curr_Dt()
    call Grid % Save_Debug_Vtu(fname,                &
                               scalar_cell = u % n,  &
                               scalar_name = "v_flux1")
    fname = "u0real-"
    write(fname(8:12), '(i5.5)') Time % Curr_Dt()
    call Grid % Save_Debug_Vtu(fname,                        &
                               scalar_cell = Vof % u_0 % n,  &
                               scalar_name = "v_flux")
  end if

  if(fun % adv_scheme .eq. UPWIND) then

    !------------------------------------------------!
    !   Discretize the system for the vof function   !
    !------------------------------------------------!
    call Vof % Discretize(A, b, dt)

    ! Solve System
    call Vof % Solve_System(Sol, b)

    call Grid % Exchange_Cells_Real(fun % n)

    !-----------------------------!
    !   Correct Volume Fraction   !
    !-----------------------------!
    do c = Cells_In_Domain()
      fun % n(c) = max(min(fun % n(c),1.0),0.0)
    end do

  else if(fun % adv_scheme .eq. CICSAM .or. &
          fun % adv_scheme .eq. STACS) then

    do i_sub = 1, n_sub

      ! Courant number full domain:
      call Vof % Max_Courant_Number(dt / real(n_sub), 0, courant_max)

      !---------------------------!
      !   Predict beta at faces   !
      !---------------------------!

      do reg = Boundary_Regions()

        if(Grid % region % type(reg) .ne. INFLOW) then
          do s_in_reg = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s_in_reg)
            c2 = Grid % faces_c(2,s_in_reg)

            fun % n(c2) = fun % n(c1)
          end do  ! faces
        end if    ! inflow
      end do      ! region

      ! Old volume fraction:
      fun % o(:) = fun % n(:)

      ! Compute gradient:
      call Flow % Grad_Variable(fun)

      call Vof % Predict_Beta()

      do corr = 1, Vof % corr_num_max

        !------------------------------------------------!
        !   Discretize the system for the vof function   !
        !------------------------------------------------!
        call Vof % Discretize(A, b, dt / real(n_sub))

        ! Solve System
        call Vof % Solve_System(Sol, b)

        do reg = Boundary_Regions()

          if(Grid % region % type(reg) .ne. INFLOW) then
            do s_in_reg = Faces_In_Region(reg)
              c1 = Grid % faces_c(1,s_in_reg)
              c2 = Grid % faces_c(2,s_in_reg)

              fun % n(c2) = fun % n(c1)
            end do  ! faces
          end if    ! inflow
        end do      ! region

        n_wrong_vf0 = 0
        n_wrong_vf1 = 0
        wrong_vf = 0

        ! Determine if 0 <= fun <= 1.0
        do c = Cells_In_Domain()
          if(fun % n(c) < -fun % tol) then
            n_wrong_vf0 = n_wrong_vf0 + 1
          end if

          if(fun % n(c) - 1.0 > fun % tol) then
            n_wrong_vf1 = n_wrong_vf1 + 1
          end if
        end do

        call Grid % Exchange_Cells_Real(fun % n)

        !---------------------------!
        !   Correct beta at faces   !
        !---------------------------!

        if(n_wrong_vf0 > 0 .or. n_wrong_vf1 > 0) then
          wrong_vf = 1
        end if

        call Global % Sum_Int(wrong_vf)

        if(wrong_vf == 0) then
          goto 1
        else
          call Vof % Correct_Beta()
        end if

      end do
1     continue

      !------------------------------------------------------------------------!
      !   Isotropic-beta diagnostic: did the imposed constant hold?            !
      !                                                                        !
      !   A constant beta_f bypasses the scheme's own CBC limiting, so the     !
      !   solution can leave [0,1].  Correct_Beta() then repairs beta_f face   !
      !   by face, which quietly restores the orientation dependence the       !
      !   constant was meant to remove.  corr = 1 means the first sweep was    !
      !   already bounded and the compression really was uniform; a larger     !
      !   value bounds how far the run departed from the intended              !
      !   intervention, and corr_num_max means it never converged.             !
      !                                                                        !
      !   Printed for the scheme's own beta as well, not only for the          !
      !   imposed constant: a sweep count is uninterpretable without the       !
      !   baseline.  Exhaustion may simply be how this case always runs.       !
      !   With corr_num_max = 1 the loop runs once and any Correct_Beta edit   !
      !   is discarded (Predict_Beta overwrites beta_f next step), so the      !
      !   scheme is PURE -- but alpha is then free to leave [0,1] and the      !
      !   violation has to be reported or the trade is invisible.              !
      !------------------------------------------------------------------------!
      call Global % Sum_Int(n_wrong_vf0)
      call Global % Sum_Int(n_wrong_vf1)
      if(First_Proc()) then
        if(Vof % beta_iso >= 0.0) then
          print '(a,f6.4,a,i0,a,i0,a,i0,a,i0)',                     &
            ' # VOF_CORR: iso beta=', Vof % beta_iso,               &
            '  sweeps=', corr, ' of max ', Vof % corr_num_max,      &
            '  unbounded lo=', n_wrong_vf0, ' hi=', n_wrong_vf1
        else
          print '(a,i0,a,i0,a,i0,a,i0)',                            &
            ' # VOF_CORR: scheme beta    sweeps=', corr,            &
            ' of max ', Vof % corr_num_max,                         &
            '  unbounded lo=', n_wrong_vf0, ' hi=', n_wrong_vf1
        end if
      end if

      !------------------------!
      !   Correct boundaries   !
      !------------------------!
      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .ne. INFLOW) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            if(fun % n(c2) < FEMTO) then
              fun % n(c2) = 0.0
            end if
            if(fun % n(c2) - 1.0 >= FEMTO) then
              fun % n(c2) = 1.0
            end if
          end do  ! faces
        end if    ! inflow
      end do      ! region

      !--------------------------------------!
      !   Correct Interior Volume Fraction   !
      !--------------------------------------!
      do c = Cells_In_Domain()
        if(fun % n(c) < FEMTO) then
          fun % n(c) = 0.0
        end if
        if(fun % n(c) - 1.0 >= FEMTO) then
          fun % n(c) = 1.0
        end if
      end do

      call Grid % Exchange_Cells_Real(fun % n)
    end do

  end if

  !------------------------------!
  !   Volume fraction at faces   !
  !------------------------------!

  if(fun % adv_scheme .eq. UPWIND) then

    ! At boundaries
    do reg = Boundary_Regions()

      if(Grid % region % type(reg) .ne. INFLOW) then
        do s_in_reg = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s_in_reg)
          c2 = Grid % faces_c(2,s_in_reg)

          fun % n(c2) = fun % n(c1)
        end do  ! faces
      end if    ! inflow
    end do      ! region

  else if(fun % adv_scheme .eq. CICSAM .or. fun % adv_scheme .eq. STACS) then

    ! At boundaries
    do reg = Boundary_Regions()

      if(Grid % region % type(reg) .ne. INFLOW) then
        do s_in_reg = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s_in_reg)
          c2 = Grid % faces_c(2,s_in_reg)

          fun % n(c2) = fun % n(c1)
        end do  ! faces
      end if    ! inflow
    end do      ! region

  end if

  call Flow % Grad_Variable(fun)

  call User_Mod_End_Of_Compute_Vof(Vof, Sol)

  call Profiler % Stop('Compute_Vof (without solvers)')

  end subroutine
