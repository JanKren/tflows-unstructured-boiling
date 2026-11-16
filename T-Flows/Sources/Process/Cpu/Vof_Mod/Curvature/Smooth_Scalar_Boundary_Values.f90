!==============================================================================!
  subroutine Smooth_Scalar_Boundary_Values(Vof, Grid, smooth_var)
!------------------------------------------------------------------------------!
!   Values of the smoothed colour in the boundary cells, for Smooth_Scalar.    !
!                                                                              !
!   Away from walls a boundary cell takes the value of its cell.  At a wall    !
!   that value makes the smoothed colour meet the wall at 90 degrees, what-    !
!   ever the contact angle: each pass of the smoothing draws the colour of     !
!   the wall layers towards 90 degrees, their normals and curvature with it,   !
!   and after the twelve passes of MAX_SMOOTHING_CYCLES_CURVATURE_VOF 12 a     !
!   cap at 45 degrees spreads.  At a wall the value is extrapolated from the   !
!   cell with its gradient turned to the contact angle, as Brackbill turns     !
!   the normal in Curvature_Csf: magnitude and tangential direction kept,      !
!   the angle with the wall normal set to theta.  The gradient is taken with   !
!   the wall values of the call before, so the passes of the smoothing also    !
!   bring the two together.                                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  type(Grid_Type)         :: Grid
  real                    :: smooth_var(-Grid % n_bnd_cells : Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  integer                   :: nb, nc, reg, s, c1, c2
  real                      :: nwallx, nwally, nwallz, dotprod, theta
  real                      :: ntanx, ntany, ntanz, mag_grad, mag_ntan
  real                      :: gwallx, gwally, gwallz
  real, contiguous, pointer :: grad_x(:), grad_y(:), grad_z(:)
  logical                   :: walls
!==============================================================================!

  Flow => Vof % pnt_flow

  nb = Grid % n_bnd_cells
  nc = Grid % n_cells

  ! Regions are the same on all ranks, and so is the decision to go
  ! into the exchanges of the gradient below
  walls = .false.
  do reg = Boundary_Regions()
    if(     Grid % region % type(reg) .eq. WALL     &
       .or. Grid % region % type(reg) .eq. WALLFL) walls = .true.
  end do

  !-----------------------------------------!
  !   Other boundaries: value of the cell   !
  !-----------------------------------------!
  do reg = Boundary_Regions()
    if(     Grid % region % type(reg) .eq. WALL     &
       .or. Grid % region % type(reg) .eq. WALLFL) cycle
    do s = Faces_In_Region(reg)
      ! Need to check with Bojan: at some runs I get s = -1, which ...
      ! ... causes a crash (hence the s > 0 guard)
      if(s>0) then
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)
        smooth_var(c2) = smooth_var(c1)
      end if
    end do
  end do

  if(.not. walls) return

  !----------------------------------------------!
  !   Walls: extrapolated at the contact angle   !
  !----------------------------------------------!
  call Work % Connect_Real_Cell(grad_x, grad_y, grad_z)

  call Flow % Grad_Component(Grid, smooth_var(-nb:nc), 1, grad_x(-nb:nc))
  call Flow % Grad_Component(Grid, smooth_var(-nb:nc), 2, grad_y(-nb:nc))
  call Flow % Grad_Component(Grid, smooth_var(-nb:nc), 3, grad_z(-nb:nc))

  ! Next to a symmetry plane a cell sees its mirror image, as in a full domain
  call Vof % Mirror_Gradient(Grid, smooth_var(-nb:nc), grad_x(-nb:nc),  &
                             grad_y(-nb:nc), grad_z(-nb:nc))

  do reg = Boundary_Regions()
    if(      Grid % region % type(reg) .ne. WALL     &
       .and. Grid % region % type(reg) .ne. WALLFL) cycle
    do s = Faces_In_Region(reg)
      if(s>0) then
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        mag_grad = sqrt(grad_x(c1)**2 + grad_y(c1)**2 + grad_z(c1)**2)

        ! Away from the interface: the value of the cell, as before
        if(mag_grad < FEMTO) then
          smooth_var(c2) = smooth_var(c1)
          cycle
        end if

        ! Unit wall normal vector directed into wall
        nwallx = Grid % sx(s) / Grid % s(s)
        nwally = Grid % sy(s) / Grid % s(s)
        nwallz = Grid % sz(s) / Grid % s(s)

        ! Direction in the wall normal to the contact line (Curvature_Csf)
        dotprod = grad_x(c1) * nwallx + grad_y(c1) * nwally  &
                + grad_z(c1) * nwallz
        ntanx = grad_x(c1) - dotprod * nwallx
        ntany = grad_y(c1) - dotprod * nwally
        ntanz = grad_z(c1) - dotprod * nwallz
        mag_ntan = sqrt(ntanx**2 + ntany**2 + ntanz**2 + FEMTO)
        ntanx = ntanx / mag_ntan
        ntany = ntany / mag_ntan
        ntanz = ntanz / mag_ntan

        ! The gradient at the contact angle
        theta  = Vof % fun % q(c2) * PI / 180.0
        gwallx = mag_grad * (nwallx * cos(theta) + ntanx * sin(theta))
        gwally = mag_grad * (nwally * cos(theta) + ntany * sin(theta))
        gwallz = mag_grad * (nwallz * cos(theta) + ntanz * sin(theta))

        ! Value at the face centre; a colour, so bounded, which an
        ! extrapolation need not be
        smooth_var(c2) = smooth_var(c1) + gwallx * Grid % dx(s)  &
                                        + gwally * Grid % dy(s)  &
                                        + gwallz * Grid % dz(s)
        smooth_var(c2) = max(0.0, min(1.0, smooth_var(c2)))
      end if
    end do
  end do

  call Work % Disconnect_Real_Cell(grad_x, grad_y, grad_z)

  end subroutine
