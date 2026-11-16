!==============================================================================!
  function Interpolate_Cell_To_Point(Exporter, Grid, field, x, y, z) result(value)
!------------------------------------------------------------------------------!
!   Interpolate a cell-centered field to a point (placeholder)               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surface_Exporter_Type), intent(inout) :: Exporter
  type(Grid_Type),             intent(in)    :: Grid
  real,                        intent(in)    :: field(:)
  real,                        intent(in)    :: x, y, z
  real                                        :: value
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, nearest_cell
  real    :: min_distance, distance
!==============================================================================!

  min_distance = 1.0e30
  nearest_cell = 1

  ! Find the nearest cell center
  do c = 1, Grid % n_cells
    distance = sqrt((Grid % xc(c) - x)**2 + (Grid % yc(c) - y)**2 + (Grid % zc(c) - z)**2)
    if (distance < min_distance) then
      min_distance = distance
      nearest_cell = c
    end if
  end do

  ! Return the field value of the nearest cell
  value = field(nearest_cell)

  end function