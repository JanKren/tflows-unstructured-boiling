!==============================================================================!
  function Point_Plane_Distance(Exporter, cut, x, y, z) result(distance)
!------------------------------------------------------------------------------!
!   Calculate signed distance from a point to the cutting plane              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surface_Exporter_Type), intent(inout) :: Exporter
  type(Surface_Cut_Type),       intent(in)    :: cut
  real,                        intent(in)    :: x, y, z
  real                                        :: distance
!==============================================================================!

  distance = (x - cut % origin(1)) * cut % normal(1) + &
             (y - cut % origin(2)) * cut % normal(2) + &
             (z - cut % origin(3)) * cut % normal(3)

  end function