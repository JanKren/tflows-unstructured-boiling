!==============================================================================!
  subroutine Calculate_Grad_Matrix_With_Front(Vof)
!------------------------------------------------------------------------------!
!   Calculates gradient matrix for cases with front tracking (boiling)         !
!                                                                              !
!   When the front-modified stencil produces a singular gradient matrix,       !
!   falls back to the standard gradient matrix for that cell.                  !
!                                                                              !
!   (Closely related to Field_Mod/Gradients/Calculate_Grad_Matrix)             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Front_Type), pointer :: Front
  integer                   :: c, c1, c2, s, cn, i_ext
  real                      :: dx_c1, dy_c1, dz_c1, dx_c2, dy_c2, dz_c2
  real                      :: jac, g_inv(6)
  real                      :: w_c1, w_c2, w_e, p_w, d_c1, d_c2
  real                      :: dx_e, dy_e, dz_e, d_e, f_c, f_n, w_f
  real, contiguous, pointer :: g1(:), g2(:), g3(:), g4(:), g5(:), g6(:)
  real, contiguous, pointer :: g1_orig(:), g2_orig(:), g3_orig(:)
  real, contiguous, pointer :: g4_orig(:), g5_orig(:), g6_orig(:)
  integer                   :: fallback_count
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!==============================================================================!

  if(DEBUG) then
    call Work % Connect_Real_Cell(g1, g2, g3, g4, g5, g6)
  end if

  ! Allocate storage for original gradient matrix as fallback
  call Work % Connect_Real_Cell(g1_orig, g2_orig, g3_orig,  &
                                g4_orig, g5_orig, g6_orig)

  ! Take alias
  Grid  => Vof % pnt_grid
  Flow  => Vof % pnt_flow
  Front => Vof % Front

  ! Least-squares weight exponent; zero reproduces the unweighted stencil
  p_w = Vof % ext_weight_p

  !----------------------------------------------------------!
  !   Compute and store the standard gradient matrix first   !
  !----------------------------------------------------------!
  call Flow % Calculate_Grad_Matrix()
  do c = 1, Grid % n_cells
    g1_orig(c) = Flow % grad_c2c(1,c)
    g2_orig(c) = Flow % grad_c2c(2,c)
    g3_orig(c) = Flow % grad_c2c(3,c)
    g4_orig(c) = Flow % grad_c2c(4,c)
    g5_orig(c) = Flow % grad_c2c(5,c)
    g6_orig(c) = Flow % grad_c2c(6,c)
  end do

  !--------------------------------------------!
  !   Initialize gradient matrices for cells   !
  !--------------------------------------------!
  do c = Cells_In_Domain_And_Buffers()
    Flow % grad_c2c(1,c) = 0.0
    Flow % grad_c2c(2,c) = 0.0
    Flow % grad_c2c(3,c) = 0.0
    Flow % grad_c2c(4,c) = 0.0
    Flow % grad_c2c(5,c) = 0.0
    Flow % grad_c2c(6,c) = 0.0
  end do

  !----------------------------------------------------------------------!
  !   Compute the gradient matrix for all cells browsing through faces   !
  !----------------------------------------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    dx_c1 = Grid % dx(s)
    dy_c1 = Grid % dy(s)
    dz_c1 = Grid % dz(s)
    dx_c2 = Grid % dx(s)
    dy_c2 = Grid % dy(s)
    dz_c2 = Grid % dz(s)

    ! If face is at the front, reduce the extents of the stencil
    if(Front % intersects_face(s)) then
      dx_c1 = Front % xs(s) - Grid % xc(c1)
      dy_c1 = Front % ys(s) - Grid % yc(c1)
      dz_c1 = Front % zs(s) - Grid % zc(c1)
      dx_c2 = Grid % xc(c2) - Front % xs(s)
      dy_c2 = Grid % yc(c2) - Front % ys(s)
      dz_c2 = Grid % zc(c2) - Front % zs(s)
    end if

    ! Least-squares weights w = |r|^-p.  With p = 0 (the default) both
    ! weights are exactly one and the sums below are unchanged.
    w_c1 = 1.0
    w_c2 = 1.0
    if(p_w > 0.0) then
      d_c1 = sqrt(dx_c1**2 + dy_c1**2 + dz_c1**2)
      d_c2 = sqrt(dx_c2**2 + dy_c2**2 + dz_c2**2)
      if(d_c1 > FEMTO) w_c1 = d_c1 ** (-p_w)
      if(d_c2 > FEMTO) w_c2 = d_c2 ** (-p_w)
    end if

    Flow % grad_c2c(1,c1) = Flow % grad_c2c(1,c1) + w_c1*dx_c1*dx_c1
    Flow % grad_c2c(2,c1) = Flow % grad_c2c(2,c1) + w_c1*dy_c1*dy_c1
    Flow % grad_c2c(3,c1) = Flow % grad_c2c(3,c1) + w_c1*dz_c1*dz_c1
    Flow % grad_c2c(4,c1) = Flow % grad_c2c(4,c1) + w_c1*dx_c1*dy_c1
    Flow % grad_c2c(5,c1) = Flow % grad_c2c(5,c1) + w_c1*dx_c1*dz_c1
    Flow % grad_c2c(6,c1) = Flow % grad_c2c(6,c1) + w_c1*dy_c1*dz_c1
    if(c2 > 0) then  ! this is enough even for parallel
      Flow % grad_c2c(1,c2) = Flow % grad_c2c(1,c2) + w_c2*dx_c2*dx_c2
      Flow % grad_c2c(2,c2) = Flow % grad_c2c(2,c2) + w_c2*dy_c2*dy_c2
      Flow % grad_c2c(3,c2) = Flow % grad_c2c(3,c2) + w_c2*dz_c2*dz_c2
      Flow % grad_c2c(4,c2) = Flow % grad_c2c(4,c2) + w_c2*dx_c2*dy_c2
      Flow % grad_c2c(5,c2) = Flow % grad_c2c(5,c2) + w_c2*dx_c2*dz_c2
      Flow % grad_c2c(6,c2) = Flow % grad_c2c(6,c2) + w_c2*dy_c2*dz_c2
    end if

  end do

  !-------------------------------------------------------------------!
  !   Extended stencil: add the node-connected (edge and corner)      !
  !   neighbours on top of the face contributions above.              !
  !                                                                   !
  !   Each pair is visited once from each side, so only cell c is     !
  !   updated here.  The interface shortening follows the same rule   !
  !   as Mark_Cells_And_Faces: the connection is cut where the        !
  !   volume fraction crosses 0.5, and the far value becomes T_sat.   !
  !-------------------------------------------------------------------!
  if(Vof % ext_stencil) then

    do c = Cells_In_Domain()
      f_c = Vof % fun % n(c)

      do i_ext = 1, Vof % ext_n(c)
        cn = Vof % ext_c(i_ext, c)
        if(cn < 1) cycle

        dx_e = Grid % xc(cn) - Grid % xc(c)
        dy_e = Grid % yc(cn) - Grid % yc(c)
        dz_e = Grid % zc(cn) - Grid % zc(c)

        ! Cut the connection at the interface, exactly as for faces
        f_n = Vof % fun % n(cn)
        if( (f_c - 0.5) * (f_n - 0.5) < 0.0 ) then
          if(abs(f_c - f_n) > FEMTO) then
            w_f  = abs(f_c - 0.5) / abs(f_c - f_n)
            dx_e = w_f * dx_e
            dy_e = w_f * dy_e
            dz_e = w_f * dz_e
          end if
        end if

        w_e = 1.0
        d_e = sqrt(dx_e**2 + dy_e**2 + dz_e**2)
        if(d_e < FEMTO) cycle
        if(p_w > 0.0) w_e = d_e ** (-p_w)

        Flow % grad_c2c(1,c) = Flow % grad_c2c(1,c) + w_e*dx_e*dx_e
        Flow % grad_c2c(2,c) = Flow % grad_c2c(2,c) + w_e*dy_e*dy_e
        Flow % grad_c2c(3,c) = Flow % grad_c2c(3,c) + w_e*dz_e*dz_e
        Flow % grad_c2c(4,c) = Flow % grad_c2c(4,c) + w_e*dx_e*dy_e
        Flow % grad_c2c(5,c) = Flow % grad_c2c(5,c) + w_e*dx_e*dz_e
        Flow % grad_c2c(6,c) = Flow % grad_c2c(6,c) + w_e*dy_e*dz_e
      end do
    end do

  end if

  !--------------------------------!
  !   Find the inverse of matrix   !
  !--------------------------------!
  fallback_count = 0
  do c = Cells_In_Domain()
    jac = Flow % grad_c2c(1,c) * Flow % grad_c2c(2,c) * Flow % grad_c2c(3,c)  &
        - Flow % grad_c2c(1,c) * Flow % grad_c2c(6,c) * Flow % grad_c2c(6,c)  &
        - Flow % grad_c2c(4,c) * Flow % grad_c2c(4,c) * Flow % grad_c2c(3,c)  &
        + Flow % grad_c2c(4,c) * Flow % grad_c2c(5,c) * Flow % grad_c2c(6,c)  &
        + Flow % grad_c2c(4,c) * Flow % grad_c2c(5,c) * Flow % grad_c2c(6,c)  &
        - Flow % grad_c2c(5,c) * Flow % grad_c2c(5,c) * Flow % grad_c2c(2,c)

    ! Fall back to standard matrix if front-modified one is singular
    if(jac <= 0.0) then
      Flow % grad_c2c(1,c) = g1_orig(c)
      Flow % grad_c2c(2,c) = g2_orig(c)
      Flow % grad_c2c(3,c) = g3_orig(c)
      Flow % grad_c2c(4,c) = g4_orig(c)
      Flow % grad_c2c(5,c) = g5_orig(c)
      Flow % grad_c2c(6,c) = g6_orig(c)
      fallback_count = fallback_count + 1
      cycle
    end if

    g_inv(1) = +(  Flow % grad_c2c(2,c) * Flow % grad_c2c(3,c)  &
                 - Flow % grad_c2c(6,c) * Flow % grad_c2c(6,c) ) / jac
    g_inv(2) = +(  Flow % grad_c2c(1,c) * Flow % grad_c2c(3,c)  &
                 - Flow % grad_c2c(5,c) * Flow % grad_c2c(5,c) ) / jac
    g_inv(3) = +(  Flow % grad_c2c(1,c) * Flow % grad_c2c(2,c)  &
                 - Flow % grad_c2c(4,c) * Flow % grad_c2c(4,c) ) / jac
    g_inv(4) = -(  Flow % grad_c2c(4,c) * Flow % grad_c2c(3,c)  &
                 - Flow % grad_c2c(5,c) * Flow % grad_c2c(6,c) ) / jac
    g_inv(5) = +(  Flow % grad_c2c(4,c) * Flow % grad_c2c(6,c)  &
                 - Flow % grad_c2c(5,c) * Flow % grad_c2c(2,c) ) / jac
    g_inv(6) = -(  Flow % grad_c2c(1,c) * Flow % grad_c2c(6,c)  &
                 - Flow % grad_c2c(4,c) * Flow % grad_c2c(5,c) ) / jac

    Flow % grad_c2c(1,c) = g_inv(1)
    Flow % grad_c2c(2,c) = g_inv(2)
    Flow % grad_c2c(3,c) = g_inv(3)
    Flow % grad_c2c(4,c) = g_inv(4)
    Flow % grad_c2c(5,c) = g_inv(5)
    Flow % grad_c2c(6,c) = g_inv(6)
  end do

  ! Report fallback usage across all processors
  call Global % Sum_Int(fallback_count)
  if(First_Proc() .and. fallback_count > 0) then
    print *, '# GRAD MATRIX: Fallback to standard for',  &
             fallback_count, 'cells (singular front stencil)'
  end if

  if(DEBUG) then
    do c = Cells_In_Domain_And_Buffers()
      g1(c) = Flow % grad_c2c(1,c)
      g2(c) = Flow % grad_c2c(2,c)
      g3(c) = Flow % grad_c2c(3,c)
      g4(c) = Flow % grad_c2c(4,c)
      g5(c) = Flow % grad_c2c(5,c)
      g6(c) = Flow % grad_c2c(6,c)
    end do
    call Grid % Save_Debug_Vtu("grad_matrix_with_front",                  &
                               tensor_cell = (/g1, g2, g3, g4, g5, g6/),  &
                               tensor_comp = 6,                           &
                               tensor_name = "grad_matrix_with_front")

    call Work % Disconnect_Real_Cell(g1, g2, g3, g4, g5, g6)
  end if

  ! Clean up fallback storage
  call Work % Disconnect_Real_Cell(g1_orig, g2_orig, g3_orig,  &
                                   g4_orig, g5_orig, g6_orig)

  end subroutine
