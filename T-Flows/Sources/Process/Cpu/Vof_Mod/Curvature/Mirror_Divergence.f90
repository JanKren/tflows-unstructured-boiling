!==============================================================================!
  subroutine Mirror_Divergence(Vof, Grid, nx, ny, nz, div_x, div_y, div_z)
!------------------------------------------------------------------------------!
!   Derivatives d nx / dx, d ny / dy and d nz / dz of a vector field in the    !
!   cells next to a symmetry plane, with the mirror image of the cell,         !
!   reflected in the plane, as the neighbour behind it.  Call it after the     !
!   least-squares gradients of the Field module, which it corrects.            !
!   Their sum is the divergence.                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  type(Grid_Type)         :: Grid
  real                    :: nx   (-Grid % n_bnd_cells : Grid % n_cells)
  real                    :: ny   (-Grid % n_bnd_cells : Grid % n_cells)
  real                    :: nz   (-Grid % n_bnd_cells : Grid % n_cells)
  real                    :: div_x(-Grid % n_bnd_cells : Grid % n_cells)
  real                    :: div_y(-Grid % n_bnd_cells : Grid % n_cells)
  real                    :: div_z(-Grid % n_bnd_cells : Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, i
!==============================================================================!

  if(.not. Vof % mirror_ready) call Vof % Mirror_Setup()

  ! The same on all ranks, so that they all go on to the exchanges
  if(.not. Vof % mirror_active) return

  call Vof % Mirror_Lsq(Grid, nx, ny, nz, .true.)

  do i = 1, Vof % n_mirror
    c = Vof % mirror_cell(i)
    div_x(c) = Vof % mirror_acc(1,i)  ! d nx / dx
    div_y(c) = Vof % mirror_acc(5,i)  ! d ny / dy
    div_z(c) = Vof % mirror_acc(9,i)  ! d nz / dz
  end do

  call Grid % Exchange_Cells_Real(div_x)
  call Grid % Exchange_Cells_Real(div_y)
  call Grid % Exchange_Cells_Real(div_z)

  end subroutine
