!==============================================================================!
  subroutine Add_Cut_Plane(Exporter, name, origin, normal, tolerance, &
                           extraction_method, search_radius)
!------------------------------------------------------------------------------!
!   Add a new cutting plane definition to the surface exporter               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surface_Exporter_Type), intent(inout) :: Exporter
  character(*),                 intent(in)    :: name
  real,                        intent(in)    :: origin(3)
  real,                        intent(in)    :: normal(3)
  real,                        intent(in)    :: tolerance
  integer,                     intent(in), optional :: extraction_method
  real,                        intent(in), optional :: search_radius
!-----------------------------------[Locals]-----------------------------------!
  type(Surface_Cut_Type), allocatable :: temp_cuts(:)
  real                                :: norm_magnitude
  integer                             :: i
!==============================================================================!

  ! Normalize the normal vector
  norm_magnitude = sqrt(normal(1)**2 + normal(2)**2 + normal(3)**2)

  if (norm_magnitude < 1.0e-12) then
    print *, 'ERROR: Normal vector magnitude is too small!'
    return
  end if

  ! Resize cuts array
  if (allocated(Exporter % cuts)) then
    allocate(temp_cuts(Exporter % n_cuts))
    temp_cuts(1:Exporter % n_cuts) = Exporter % cuts(1:Exporter % n_cuts)
    deallocate(Exporter % cuts)
    allocate(Exporter % cuts(Exporter % n_cuts + 1))
    Exporter % cuts(1:Exporter % n_cuts) = temp_cuts(1:Exporter % n_cuts)
    deallocate(temp_cuts)
  else
    allocate(Exporter % cuts(1))
  end if

  Exporter % n_cuts = Exporter % n_cuts + 1
  i = Exporter % n_cuts

  ! Set up the new cut
  Exporter % cuts(i) % name      = trim(name)
  Exporter % cuts(i) % origin(1) = origin(1)
  Exporter % cuts(i) % origin(2) = origin(2)
  Exporter % cuts(i) % origin(3) = origin(3)
  Exporter % cuts(i) % normal(1) = normal(1) / norm_magnitude
  Exporter % cuts(i) % normal(2) = normal(2) / norm_magnitude
  Exporter % cuts(i) % normal(3) = normal(3) / norm_magnitude
  Exporter % cuts(i) % tolerance = tolerance
  Exporter % cuts(i) % active    = .true.

  ! Set extraction parameters with defaults
  if (present(extraction_method)) then
    Exporter % cuts(i) % extraction_method = extraction_method
  else
    Exporter % cuts(i) % extraction_method = 1  ! Default: cell centers
  end if

  if (present(search_radius)) then
    Exporter % cuts(i) % search_radius = search_radius
  else
    Exporter % cuts(i) % search_radius = tolerance * 10.0  ! Default: 10x tolerance
  end if

  print *, 'Added surface cut: ', trim(name)
  print *, '  Origin: ', origin
  print *, '  Normal: ', Exporter % cuts(i) % normal
  print *, '  Tolerance: ', tolerance

  end subroutine