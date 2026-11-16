!==============================================================================!
  subroutine Mirror_Gradient(Vof, Grid, phi, phi_x, phi_y, phi_z)
!------------------------------------------------------------------------------!
!   Gradient of a scalar in the cells next to a symmetry plane, with the       !
!   mirror image of the cell as the neighbour behind the plane.  Call it       !
!   after the least-squares gradient of the Field module, which it corrects.   !
!   The scalar has to be even in the plane, like a colour function.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  type(Grid_Type)         :: Grid
  real                    :: phi  (-Grid % n_bnd_cells : Grid % n_cells)
  real                    :: phi_x(-Grid % n_bnd_cells : Grid % n_cells)
  real                    :: phi_y(-Grid % n_bnd_cells : Grid % n_cells)
  real                    :: phi_z(-Grid % n_bnd_cells : Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i
!==============================================================================!

  if(.not. Vof % mirror_ready) call Vof % Mirror_Setup()

  ! The same on all ranks, so that they all go on to the exchanges
  if(.not. Vof % mirror_active) return

  call Vof % Mirror_Lsq(Grid, phi, phi, phi, .false.)

  do i = 1, Vof % n_mirror
    c = Vof % mirror_cell(i)
    phi_x(c) = Vof % mirror_acc(1,i)
    phi_y(c) = Vof % mirror_acc(2,i)
    phi_z(c) = Vof % mirror_acc(3,i)
  end do

  call Grid % Exchange_Cells_Real(phi_x)
  call Grid % Exchange_Cells_Real(phi_y)
  call Grid % Exchange_Cells_Real(phi_z)

  end subroutine
