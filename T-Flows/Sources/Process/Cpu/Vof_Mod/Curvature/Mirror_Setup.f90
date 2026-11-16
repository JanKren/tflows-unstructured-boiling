!==============================================================================!
  subroutine Mirror_Setup(Vof)
!------------------------------------------------------------------------------!
!   Prepares the gradients of the cells that have a symmetry face.             !
!                                                                              !
!   The least-squares gradient of the Field module counts a symmetry face as   !
!   a neighbour at the centre of the boundary cell, half way to the plane,     !
!   and gives that neighbour the value of the cell itself.  What a full        !
!   domain has on the other side of the plane is the mirror image of the       !
!   cell: twice as far, with the same colour and with the interface normal     !
!   reflected in the plane.  With the neighbour half way the gradient normal   !
!   to the plane is pulled towards zero, and the curvature next to a          !
!   symmetry plane comes out 20 to 40 % low.                                   !
!                                                                              !
!   The cells with a symmetry face get a matrix of their own here: the one     !
!   of the Field module, with each symmetry face counted at twice the          !
!   distance.  Mirror_Gradient and Mirror_Divergence use it.                   !
!                                                                              !
!   It is collective: whether there is a symmetry plane on any rank is         !
!   decided for all of them, because the routines that use it exchange.        !
!   Call it again after the grid has changed (mirror_ready = .false.).         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c, c1, c2, i, s, reg, k
  real                     :: w, dx, dy, dz, jac, a(6)
  logical                  :: has_plane
!==============================================================================!

  Grid => Vof % pnt_grid

  if(allocated(Vof % mirror_cell)) deallocate(Vof % mirror_cell)
  if(allocated(Vof % mirror_idx))  deallocate(Vof % mirror_idx)
  if(allocated(Vof % mirror_inv))  deallocate(Vof % mirror_inv)
  if(allocated(Vof % mirror_acc))  deallocate(Vof % mirror_acc)

  allocate(Vof % mirror_idx(Grid % n_cells))
  Vof % mirror_idx(:) = 0

  !----------------------------------------------!
  !   Find the cells which have a symmetry face   !
  !----------------------------------------------!
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .ne. SYMMETRY) cycle
    do s = Faces_In_Region(reg)
      if(s < 1) cycle
      Vof % mirror_idx(Grid % faces_c(1,s)) = -1
    end do
  end do

  Vof % n_mirror = 0
  do c = Cells_In_Domain()
    if(Vof % mirror_idx(c) .eq. -1) then
      Vof % n_mirror = Vof % n_mirror + 1
      Vof % mirror_idx(c) = Vof % n_mirror
    end if
  end do
  where(Vof % mirror_idx < 0) Vof % mirror_idx = 0

  allocate(Vof % mirror_cell(Vof % n_mirror))
  allocate(Vof % mirror_inv (6, Vof % n_mirror))
  allocate(Vof % mirror_acc (9, Vof % n_mirror))
  do c = Cells_In_Domain()
    if(Vof % mirror_idx(c) > 0) Vof % mirror_cell(Vof % mirror_idx(c)) = c
  end do

  !-------------------------------------------------------------------!
  !   Matrix of these cells, as in Calculate_Grad_Matrix, but with    !
  !   the neighbour behind a symmetry face at twice the distance      !
  !-------------------------------------------------------------------!
  Vof % mirror_inv(:,:) = 0.0

  do reg = Boundary_Regions()
    w = 1.0
    if(Grid % region % type(reg) .eq. SYMMETRY) w = 4.0  ! (2 d) (2 d)
    do s = Faces_In_Region(reg)
      if(s < 1) cycle
      i = Vof % mirror_idx(Grid % faces_c(1,s))
      if(i .eq. 0) cycle
      dx = Grid % dx(s)
      dy = Grid % dy(s)
      dz = Grid % dz(s)
      Vof % mirror_inv(1,i) = Vof % mirror_inv(1,i) + w * dx * dx
      Vof % mirror_inv(2,i) = Vof % mirror_inv(2,i) + w * dy * dy
      Vof % mirror_inv(3,i) = Vof % mirror_inv(3,i) + w * dz * dz
      Vof % mirror_inv(4,i) = Vof % mirror_inv(4,i) + w * dx * dy
      Vof % mirror_inv(5,i) = Vof % mirror_inv(5,i) + w * dx * dz
      Vof % mirror_inv(6,i) = Vof % mirror_inv(6,i) + w * dy * dz
    end do
  end do

  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    dx = Grid % dx(s)
    dy = Grid % dy(s)
    dz = Grid % dz(s)
    do k = 1, 2
      if(k .eq. 1) i = Vof % mirror_idx(c1)
      if(k .eq. 2) i = Vof % mirror_idx(c2)
      if(i .eq. 0) cycle
      Vof % mirror_inv(1,i) = Vof % mirror_inv(1,i) + dx * dx
      Vof % mirror_inv(2,i) = Vof % mirror_inv(2,i) + dy * dy
      Vof % mirror_inv(3,i) = Vof % mirror_inv(3,i) + dz * dz
      Vof % mirror_inv(4,i) = Vof % mirror_inv(4,i) + dx * dy
      Vof % mirror_inv(5,i) = Vof % mirror_inv(5,i) + dx * dz
      Vof % mirror_inv(6,i) = Vof % mirror_inv(6,i) + dy * dz
    end do
  end do

  !--------------------------------!
  !   Find the inverse of matrix   !
  !--------------------------------!
  do i = 1, Vof % n_mirror
    a(1:6) = Vof % mirror_inv(1:6,i)
    jac = a(1) * a(2) * a(3)         &
        - a(1) * a(6) * a(6)         &
        - a(4) * a(4) * a(3)         &
        + a(4) * a(5) * a(6) * 2.0   &
        - a(5) * a(5) * a(2)
    Assert(jac > 0.0)
    Vof % mirror_inv(1,i) = +(a(2) * a(3) - a(6) * a(6)) / jac
    Vof % mirror_inv(2,i) = +(a(1) * a(3) - a(5) * a(5)) / jac
    Vof % mirror_inv(3,i) = +(a(1) * a(2) - a(4) * a(4)) / jac
    Vof % mirror_inv(4,i) = -(a(4) * a(3) - a(5) * a(6)) / jac
    Vof % mirror_inv(5,i) = +(a(4) * a(6) - a(5) * a(2)) / jac
    Vof % mirror_inv(6,i) = -(a(1) * a(6) - a(4) * a(5)) / jac
  end do

  !-----------------------------------------------!
  !   Is there a symmetry plane on any processor?   !
  !-----------------------------------------------!
  has_plane = (Vof % n_mirror > 0)
  call Global % Lor_Log(has_plane)
  Vof % mirror_active = has_plane
  Vof % mirror_ready  = .true.

  end subroutine
