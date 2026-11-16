!==============================================================================!
  subroutine Extract_Surface(Exporter, cut_index, Grid, Flow, Vof)
!------------------------------------------------------------------------------!
!   Extract surface points and triangles for a specific cutting plane        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surface_Exporter_Type), intent(inout) :: Exporter
  integer,                     intent(in)    :: cut_index
  type(Grid_Type),             intent(in)    :: Grid
  type(Field_Type),            intent(in), optional :: Flow
  type(Vof_Type),              intent(in), optional :: Vof
!-----------------------------------[Locals]-----------------------------------!
  integer                         :: s, c1, c2
  real                           :: dist1, dist2, t_intersect, face_dist
  real                           :: x_int, y_int, z_int
  real                           :: dist_c1, dist_c2, weight_c1, weight_c2
  integer                        :: n_intersections
  real, allocatable              :: temp_points(:,:)
  integer, allocatable           :: temp_triangles(:,:)
  real, allocatable              :: temp_fields(:,:)
  integer                        :: max_points, max_triangles, max_fields
  real                           :: u_int, v_int, w_int, t_int, vof_int, p_int
  integer                        :: field_count, pt_index, method1_points, method2_points
!==============================================================================!

  if (cut_index < 1 .or. cut_index > Exporter % n_cuts) return

  if (.not. Exporter % cuts(cut_index) % active) return

  ! Initialize temporary storage
  max_points = Grid % n_faces * 2  ! Conservative estimate
  max_triangles = Grid % n_faces

  ! Calculate maximum number of fields
  max_fields = 0
  if (present(Flow)) max_fields = max_fields + 4  ! u, v, w, pressure
  if (present(Flow) .and. allocated(Flow % t % n)) max_fields = max_fields + 1  ! temperature
  if (present(Vof) .and. allocated(Vof % fun % n)) max_fields = max_fields + 1  ! VOF
  if (present(Vof) .and. allocated(Vof % m_dot)) max_fields = max_fields + 1  ! mass transfer rate

  allocate(temp_points(3, max_points))
  allocate(temp_triangles(3, max_triangles))
  if (max_fields > 0) allocate(temp_fields(max_fields, max_points))

  n_intersections = 0
  Exporter % cuts(cut_index) % n_triangles = 0

  method1_points = 0
  method2_points = 0

  if (First_Proc() .and. Exporter % first_extraction) then
    print *, 'Surface_Exporter: tolerance = ', Exporter % cuts(cut_index) % tolerance
    print *, 'Surface_Exporter: search_radius = ', Exporter % cuts(cut_index) % search_radius
    print *, 'Surface_Exporter: extraction_method = ', Exporter % cuts(cut_index) % extraction_method
  end if

  ! ParaView-style cell-based intersection method
  if (Exporter % cuts(cut_index) % extraction_method == 1 .or. &
      Exporter % cuts(cut_index) % extraction_method == 3) then
    ! Method 1: Process each cell to find plane intersections
    do c1 = 1, Grid % n_cells
      ! Calculate distance from cell center to cutting plane
      dist1 = Exporter % Point_Plane_Distance(Exporter % cuts(cut_index), Grid % xc(c1), Grid % yc(c1), Grid % zc(c1))

      ! Check if cell is close enough to the plane to potentially intersect
      if (abs(dist1) <= Exporter % cuts(cut_index) % tolerance) then
        ! For cells very close to the plane, add the cell center projected onto the plane
        x_int = Grid % xc(c1) - dist1 * Exporter % cuts(cut_index) % normal(1)
        y_int = Grid % yc(c1) - dist1 * Exporter % cuts(cut_index) % normal(2)
        z_int = Grid % zc(c1) - dist1 * Exporter % cuts(cut_index) % normal(3)

        n_intersections = n_intersections + 1
        method1_points = method1_points + 1
        if (n_intersections <= max_points) then
          temp_points(1, n_intersections) = x_int
          temp_points(2, n_intersections) = y_int
          temp_points(3, n_intersections) = z_int

          ! Use field values directly from the cell (ParaView-style)
          if (max_fields > 0) then
            field_count = 0

            if (present(Flow)) then
              field_count = field_count + 1
              temp_fields(field_count, n_intersections) = Flow % u % n(c1)
              field_count = field_count + 1
              temp_fields(field_count, n_intersections) = Flow % v % n(c1)
              field_count = field_count + 1
              temp_fields(field_count, n_intersections) = Flow % w % n(c1)
              field_count = field_count + 1
              temp_fields(field_count, n_intersections) = Flow % p % n(c1)

              ! Temperature (if available)
              if (allocated(Flow % t % n)) then
                field_count = field_count + 1
                temp_fields(field_count, n_intersections) = Flow % t % n(c1)
              end if
            end if

            ! VOF (if available)
            if (present(Vof) .and. allocated(Vof % fun % n)) then
              field_count = field_count + 1
              temp_fields(field_count, n_intersections) = Vof % fun % n(c1)
            end if

            ! Mass transfer rate (if available)
            if (present(Vof) .and. allocated(Vof % m_dot)) then
              field_count = field_count + 1
              temp_fields(field_count, n_intersections) = Vof % m_dot(c1)
            end if
          end if
        end if
      end if
    end do
  end if

  ! Method 2: Denser cell sampling (ParaView-style with larger search radius)
  if (Exporter % cuts(cut_index) % extraction_method == 2 .or. &
      Exporter % cuts(cut_index) % extraction_method == 3) then

    ! Method 2 uses a larger search radius to get more cells
    do c1 = 1, Grid % n_cells
      ! Calculate distance from cell center to cutting plane
      dist1 = Exporter % Point_Plane_Distance(Exporter % cuts(cut_index), Grid % xc(c1), Grid % yc(c1), Grid % zc(c1))

      ! Use search_radius for Method 2 (larger tolerance for denser sampling)
      if (abs(dist1) <= Exporter % cuts(cut_index) % search_radius) then
        ! Project cell center onto the cutting plane
        x_int = Grid % xc(c1) - dist1 * Exporter % cuts(cut_index) % normal(1)
        y_int = Grid % yc(c1) - dist1 * Exporter % cuts(cut_index) % normal(2)
        z_int = Grid % zc(c1) - dist1 * Exporter % cuts(cut_index) % normal(3)

        n_intersections = n_intersections + 1
        method2_points = method2_points + 1
        if (n_intersections <= max_points) then
          temp_points(1, n_intersections) = x_int
          temp_points(2, n_intersections) = y_int
          temp_points(3, n_intersections) = z_int

          ! Use field values directly from the cell (ParaView-style)
          if (max_fields > 0) then
            field_count = 0

            if (present(Flow)) then
              field_count = field_count + 1
              temp_fields(field_count, n_intersections) = Flow % u % n(c1)
              field_count = field_count + 1
              temp_fields(field_count, n_intersections) = Flow % v % n(c1)
              field_count = field_count + 1
              temp_fields(field_count, n_intersections) = Flow % w % n(c1)
              field_count = field_count + 1
              temp_fields(field_count, n_intersections) = Flow % p % n(c1)

              ! Temperature (if available)
              if (allocated(Flow % t % n)) then
                field_count = field_count + 1
                temp_fields(field_count, n_intersections) = Flow % t % n(c1)
              end if
            end if

            ! VOF (if available)
            if (present(Vof) .and. allocated(Vof % fun % n)) then
              field_count = field_count + 1
              temp_fields(field_count, n_intersections) = Vof % fun % n(c1)
            end if

            ! Mass transfer rate (if available)
            if (present(Vof) .and. allocated(Vof % m_dot)) then
              field_count = field_count + 1
              temp_fields(field_count, n_intersections) = Vof % m_dot(c1)
            end if
          end if
        end if
      end if
    end do
  end if

  ! Store the extracted points
  Exporter % cuts(cut_index) % n_points = n_intersections
  if (Exporter % cuts(cut_index) % n_points > 0) then
    if (allocated(Exporter % cuts(cut_index) % points)) deallocate(Exporter % cuts(cut_index) % points)
    allocate(Exporter % cuts(cut_index) % points(3, Exporter % cuts(cut_index) % n_points))
    Exporter % cuts(cut_index) % points(:, 1:Exporter % cuts(cut_index) % n_points) = temp_points(:, 1:Exporter % cuts(cut_index) % n_points)

    ! Count available fields and allocate field data storage
    field_count = 0
    if (present(Flow)) field_count = field_count + 4  ! u, v, w, pressure
    if (present(Flow) .and. allocated(Flow % t % n)) field_count = field_count + 1  ! temperature
    if (present(Vof) .and. allocated(Vof % fun % n)) field_count = field_count + 1  ! VOF
    if (present(Vof) .and. allocated(Vof % m_dot)) field_count = field_count + 1  ! mass transfer rate

    if (field_count > 0) then
      if (allocated(Exporter % cuts(cut_index) % field_data)) deallocate(Exporter % cuts(cut_index) % field_data)
      if (allocated(Exporter % cuts(cut_index) % field_names)) deallocate(Exporter % cuts(cut_index) % field_names)
      allocate(Exporter % cuts(cut_index) % field_data(field_count, Exporter % cuts(cut_index) % n_points))
      allocate(Exporter % cuts(cut_index) % field_names(field_count))
      Exporter % cuts(cut_index) % n_fields = field_count

      ! Set field names
      field_count = 0
      if (present(Flow)) then
        field_count = field_count + 1
        Exporter % cuts(cut_index) % field_names(field_count) = 'Velocity_U'
        field_count = field_count + 1
        Exporter % cuts(cut_index) % field_names(field_count) = 'Velocity_V'
        field_count = field_count + 1
        Exporter % cuts(cut_index) % field_names(field_count) = 'Velocity_W'
        field_count = field_count + 1
        Exporter % cuts(cut_index) % field_names(field_count) = 'Pressure'
        if (allocated(Flow % t % n)) then
          field_count = field_count + 1
          Exporter % cuts(cut_index) % field_names(field_count) = 'Temperature'
        end if
      end if
      if (present(Vof) .and. allocated(Vof % fun % n)) then
        field_count = field_count + 1
        Exporter % cuts(cut_index) % field_names(field_count) = 'VOF_Sharp'
      end if
      if (present(Vof) .and. allocated(Vof % m_dot)) then
        field_count = field_count + 1
        Exporter % cuts(cut_index) % field_names(field_count) = 'Mass_Transfer_Rate'
      end if

      ! Copy field data from temporary storage
      if (max_fields > 0) then
        Exporter % cuts(cut_index) % field_data(:, 1:Exporter % cuts(cut_index) % n_points) = &
          temp_fields(:, 1:Exporter % cuts(cut_index) % n_points)
      end if
    end if
  end if

  deallocate(temp_points)
  deallocate(temp_triangles)
  if (max_fields > 0) deallocate(temp_fields)

  if (First_Proc() .and. Exporter % first_extraction) then
    print *, 'Surface_Exporter: Method 1 contributed ', method1_points, ' points'
    print *, 'Surface_Exporter: Method 2 contributed ', method2_points, ' points'
    Exporter % first_extraction = .false.
  end if
  ! Only print extraction summary occasionally or when points > 0
  if (Exporter % cuts(cut_index) % n_points > 0) then
    print *, 'Extracted ', Exporter % cuts(cut_index) % n_points, ' points for surface: ', trim(Exporter % cuts(cut_index) % name)
  end if

  end subroutine
