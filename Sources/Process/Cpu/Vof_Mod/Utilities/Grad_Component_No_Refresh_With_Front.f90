!==============================================================================!
  subroutine Grad_Component_No_Refresh_With_Front(Vof, phi, i, phii, phif)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi by a least squares method,     !
!   without refreshing the buffers.                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  real,    intent(in)     :: phi (-Vof % pnt_grid % n_bnd_cells:  &
                                   Vof % pnt_grid % n_cells)
  integer, intent(in)     :: i
  real,    intent(out)    :: phii(-Vof % pnt_grid % n_bnd_cells:  &
                                   Vof % pnt_grid % n_cells)
  real,    intent(in)     :: phif  ! phi at the front
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Front_Type), pointer :: Front
  integer                   :: s, c1, c2, c, cn, i_ext
  real                      :: dphi1, dphi2
  real                      :: dx_c1, dy_c1, dz_c1, dx_c2, dy_c2, dz_c2
  real                      :: w_c1, w_c2, p_w, d_c1, d_c2
  real                      :: dx_e, dy_e, dz_e, d_e, w_e, f_c, f_n, w_f
  real                      :: dphi_e
  integer                   :: n_nan
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAX_NAN_REPORTS = 5
  integer, dimension(3,3), parameter :: MAP = reshape((/ 1, 4, 5,  &
                                                         4, 2, 6,  &
                                                         5, 6, 3 /), shape(MAP))
!==============================================================================!

  ! Take alias
  Grid  => Vof % pnt_grid
  Flow  => Vof % pnt_flow
  Front => Vof % Front

  ! Least-squares weight exponent; zero reproduces the unweighted stencil
  p_w = Vof % ext_weight_p

  ! Initialize gradients
  phii(1:Grid % n_cells) = 0.

  ! Cap on the NaN reports below.  Once a field goes non-finite EVERY face
  ! touching it reports, so on a 618k-cell case this printed tens of
  ! millions of lines per step: two runs produced 131 GB and 201 GB logs,
  ! and one of them filled the filesystem to zero bytes and took a
  ! neighbouring job with it.  A handful of faces names the location just as
  ! well as all of them, and the count that follows says how bad it is.
  n_nan = 0

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! DEBUG: Check input phi values for NaN
    if(phi(c1) /= phi(c1) .or. phi(c2) /= phi(c2)) then
      n_nan = n_nan + 1
      if(First_Proc() .and. n_nan <= MAX_NAN_REPORTS) then
        print '(a,i0,a,i0,a,i0)', &
              ' # [GRAD NAN] Input phi NaN at face ', s, &
              ' c1=', c1, ' c2=', c2
      end if
      cycle  ! skip this face
    end if

    ! Take differences as if front doesn't exist
    dphi1 = phi(c2)-phi(c1)
    dphi2 = phi(c2)-phi(c1)
    dx_c1 = Grid % dx(s)
    dy_c1 = Grid % dy(s)
    dz_c1 = Grid % dz(s)
    dx_c2 = Grid % dx(s)
    dy_c2 = Grid % dy(s)
    dz_c2 = Grid % dz(s)

    ! If face is at the front, reduce the extents of the stencil
    if(Front % intersects_face(s)) then
      dphi1 = phif - phi(c1)
      dphi2 = phi(c2) - phif
      dx_c1 = Vof % Front % xs(s) - Grid % xc(c1)
      dy_c1 = Vof % Front % ys(s) - Grid % yc(c1)
      dz_c1 = Vof % Front % zs(s) - Grid % zc(c1)
      dx_c2 = Grid % xc(c2) - Vof % Front % xs(s)
      dy_c2 = Grid % yc(c2) - Vof % Front % ys(s)
      dz_c2 = Grid % zc(c2) - Vof % Front % zs(s)

      ! DEBUG: Check for NaN in front intersection coordinates
      if(dx_c1 /= dx_c1 .or. dy_c1 /= dy_c1 .or. dz_c1 /= dz_c1 .or.  &
         dx_c2 /= dx_c2 .or. dy_c2 /= dy_c2 .or. dz_c2 /= dz_c2) then
        n_nan = n_nan + 1
        if(First_Proc() .and. n_nan <= MAX_NAN_REPORTS) then
          print '(a,i0,a,i0,a,i0)', &
                ' # [GRAD NAN] Front intersection NaN at face ', s, &
                ' c1=', c1, ' c2=', c2
        end if
        cycle  ! skip this face to avoid propagating NaN
      end if
    end if

    ! Least-squares weights, matching Calculate_Grad_Matrix_With_Front
    w_c1 = 1.0
    w_c2 = 1.0
    if(p_w > 0.0) then
      d_c1 = sqrt(dx_c1**2 + dy_c1**2 + dz_c1**2)
      d_c2 = sqrt(dx_c2**2 + dy_c2**2 + dz_c2**2)
      if(d_c1 > FEMTO) w_c1 = d_c1 ** (-p_w)
      if(d_c2 > FEMTO) w_c2 = d_c2 ** (-p_w)
    end if

    ! On the boundaries
    if(c2 < 0) then
      if(Grid % Bnd_Cond_Type(c2) .eq. SYMMETRY) then
        dphi1 = 0.
      end if

      phii(c1) = phii(c1)                                              &
               + w_c1 * dphi1 * (  Flow % grad_c2c(MAP(i,1),c1) * dx_c1   &
                                 + Flow % grad_c2c(MAP(i,2),c1) * dy_c1   &
                                 + Flow % grad_c2c(MAP(i,3),c1) * dz_c1)
    end if

    ! Inside the domain
    if(c2 > 0) then
      phii(c1) = phii(c1)                                              &
               + w_c1 * dphi1 * (  Flow % grad_c2c(MAP(i,1),c1) * dx_c1   &
                                 + Flow % grad_c2c(MAP(i,2),c1) * dy_c1   &
                                 + Flow % grad_c2c(MAP(i,3),c1) * dz_c1)

      phii(c2) = phii(c2)                                              &
               + w_c2 * dphi2 * (  Flow % grad_c2c(MAP(i,1),c2) * dx_c2   &
                                 + Flow % grad_c2c(MAP(i,2),c2) * dy_c2   &
                                 + Flow % grad_c2c(MAP(i,3),c2) * dz_c2)
    end if
  end do

  !----------------------------------------------------------------!
  !   Extended stencil: the same edge and corner neighbours that   !
  !   Calculate_Grad_Matrix_With_Front added to the matrix.        !
  !----------------------------------------------------------------!
  if(Vof % ext_stencil) then

    do c = Cells_In_Domain()
      f_c = Vof % fun % n(c)

      do i_ext = 1, Vof % ext_n(c)
        cn = Vof % ext_c(i_ext, c)
        if(cn < 1) cycle

        dx_e = Grid % xc(cn) - Grid % xc(c)
        dy_e = Grid % yc(cn) - Grid % yc(c)
        dz_e = Grid % zc(cn) - Grid % zc(c)

        dphi_e = phi(cn) - phi(c)

        ! Cut at the interface and impose the front value there
        f_n = Vof % fun % n(cn)
        if( (f_c - 0.5) * (f_n - 0.5) < 0.0 ) then
          if(abs(f_c - f_n) > FEMTO) then
            w_f    = abs(f_c - 0.5) / abs(f_c - f_n)
            dx_e   = w_f * dx_e
            dy_e   = w_f * dy_e
            dz_e   = w_f * dz_e
            dphi_e = phif - phi(c)
          end if
        end if

        w_e = 1.0
        d_e = sqrt(dx_e**2 + dy_e**2 + dz_e**2)
        if(d_e < FEMTO) cycle
        if(p_w > 0.0) w_e = d_e ** (-p_w)

        phii(c) = phii(c)                                            &
                + w_e * dphi_e * (  Flow % grad_c2c(MAP(i,1),c) * dx_e  &
                                  + Flow % grad_c2c(MAP(i,2),c) * dy_e  &
                                  + Flow % grad_c2c(MAP(i,3),c) * dz_e)
      end do
    end do

  end if

  ! Say how bad it was.  The cap above prints the first few faces, which
  ! names the location; this names the extent, and the difference between
  ! five and five million is the difference between a stray cell and a
  ! field that has already gone.
  if(First_Proc() .and. n_nan > MAX_NAN_REPORTS) then
    print '(a,i0,a)', ' # [GRAD NAN] ... and ', n_nan - MAX_NAN_REPORTS,  &
          ' further NaN faces this call (reports capped)'
  end if

  end subroutine
