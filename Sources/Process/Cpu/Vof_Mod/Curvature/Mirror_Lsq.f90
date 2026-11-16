!==============================================================================!
  subroutine Mirror_Lsq(Vof, Grid, qx, qy, qz, reflect)
!------------------------------------------------------------------------------!
!   Least-squares gradients of the fields qx, qy and qz in the cells with a    !
!   symmetry face, with the mirror image of the cell as the neighbour behind   !
!   the face (see Mirror_Setup).  The gradients, nine numbers per cell, are    !
!   left in Vof % mirror_acc:                                                  !
!                                                                              !
!     (1:3) gradient of qx,  (4:6) gradient of qy,  (7:9) gradient of qz.      !
!                                                                              !
!   reflect = .false.:  qx, qy and qz are scalars; the image has the value of  !
!                       the cell itself.                                       !
!   reflect = .true.:   (qx, qy, qz) is a vector; the image has it reflected   !
!                       in the plane, n - 2 (n . s) s, with s the plane's unit !
!                       normal.                                                !
!                                                                              !
!   Everything but the symmetry faces is summed as in Grad_Component, so       !
!   that a cell next to a wall or to another cell is treated like any other.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  type(Grid_Type)         :: Grid
  real                    :: qx(-Grid % n_bnd_cells : Grid % n_cells)
  real                    :: qy(-Grid % n_bnd_cells : Grid % n_cells)
  real                    :: qz(-Grid % n_bnd_cells : Grid % n_cells)
  logical                 :: reflect
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c1, c2, i, k, reg
  real    :: da, db, dc, dx, dy, dz, sx, sy, sz, dot, g(6), r(3)
  logical :: is_plane
!==============================================================================!

  Vof % mirror_acc(:,:) = 0.0

  !---------------------------------------------------------------!
  !   Faces on the boundaries, update the cell next to them only   !
  !---------------------------------------------------------------!
  do reg = Boundary_Regions()
    is_plane = (Grid % region % type(reg) .eq. SYMMETRY)

    ! A scalar is its own image, which adds nothing to the right-hand side
    if(is_plane .and. .not. reflect) cycle

    do s = Faces_In_Region(reg)
      if(s < 1) cycle
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      i  = Vof % mirror_idx(c1)
      if(i .eq. 0) cycle

      if(is_plane) then
        sx  = Grid % sx(s) / Grid % s(s)
        sy  = Grid % sy(s) / Grid % s(s)
        sz  = Grid % sz(s) / Grid % s(s)
        dot = qx(c1) * sx + qy(c1) * sy + qz(c1) * sz
        da  = -2.0 * dot * sx
        db  = -2.0 * dot * sy
        dc  = -2.0 * dot * sz
        dx  = 2.0 * Grid % dx(s)
        dy  = 2.0 * Grid % dy(s)
        dz  = 2.0 * Grid % dz(s)
      else
        da  = qx(c2) - qx(c1)
        db  = qy(c2) - qy(c1)
        dc  = qz(c2) - qz(c1)
        dx  = Grid % dx(s)
        dy  = Grid % dy(s)
        dz  = Grid % dz(s)
      end if

      Vof % mirror_acc(1,i) = Vof % mirror_acc(1,i) + dx * da
      Vof % mirror_acc(2,i) = Vof % mirror_acc(2,i) + dy * da
      Vof % mirror_acc(3,i) = Vof % mirror_acc(3,i) + dz * da
      Vof % mirror_acc(4,i) = Vof % mirror_acc(4,i) + dx * db
      Vof % mirror_acc(5,i) = Vof % mirror_acc(5,i) + dy * db
      Vof % mirror_acc(6,i) = Vof % mirror_acc(6,i) + dz * db
      Vof % mirror_acc(7,i) = Vof % mirror_acc(7,i) + dx * dc
      Vof % mirror_acc(8,i) = Vof % mirror_acc(8,i) + dy * dc
      Vof % mirror_acc(9,i) = Vof % mirror_acc(9,i) + dz * dc
    end do
  end do

  !---------------------------!
  !   Faces inside the domain   !
  !---------------------------!
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    da = qx(c2) - qx(c1)
    db = qy(c2) - qy(c1)
    dc = qz(c2) - qz(c1)
    dx = Grid % dx(s)
    dy = Grid % dy(s)
    dz = Grid % dz(s)

    do k = 1, 2
      if(k .eq. 1) i = Vof % mirror_idx(c1)
      if(k .eq. 2) i = Vof % mirror_idx(c2)
      if(i .eq. 0) cycle
      Vof % mirror_acc(1,i) = Vof % mirror_acc(1,i) + dx * da
      Vof % mirror_acc(2,i) = Vof % mirror_acc(2,i) + dy * da
      Vof % mirror_acc(3,i) = Vof % mirror_acc(3,i) + dz * da
      Vof % mirror_acc(4,i) = Vof % mirror_acc(4,i) + dx * db
      Vof % mirror_acc(5,i) = Vof % mirror_acc(5,i) + dy * db
      Vof % mirror_acc(6,i) = Vof % mirror_acc(6,i) + dz * db
      Vof % mirror_acc(7,i) = Vof % mirror_acc(7,i) + dx * dc
      Vof % mirror_acc(8,i) = Vof % mirror_acc(8,i) + dy * dc
      Vof % mirror_acc(9,i) = Vof % mirror_acc(9,i) + dz * dc
    end do
  end do

  !-----------------------------------------!
  !   Multiply the sums with the inverse     !
  !-----------------------------------------!
  do i = 1, Vof % n_mirror
    g(1:6) = Vof % mirror_inv(1:6,i)
    do k = 0, 6, 3
      r(1:3) = Vof % mirror_acc(k+1:k+3,i)
      Vof % mirror_acc(k+1,i) = g(1) * r(1) + g(4) * r(2) + g(5) * r(3)
      Vof % mirror_acc(k+2,i) = g(4) * r(1) + g(2) * r(2) + g(6) * r(3)
      Vof % mirror_acc(k+3,i) = g(5) * r(1) + g(6) * r(2) + g(3) * r(3)
    end do
  end do

  end subroutine
