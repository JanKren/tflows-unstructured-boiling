!==============================================================================!
  subroutine Vof_Initialization_Ellipsoid(Vof)
!------------------------------------------------------------------------------!
!   Initialize dependent variables.  (It is a bit of a mess still)             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  integer                   :: c, n, fu, i_nod
  integer                   :: ee, n_ellipses
  integer                   :: i, j, k, n_in
  real                      :: radius_1, radius_2, radius_3
  real                      :: cent_x, cent_y, cent_z
  real                      :: x_min, x_max, y_min, y_max, z_min, z_max
  real                      :: xs, ys, zs
  integer, parameter        :: N_SUB = 32  ! samples per cell edge
  real, contiguous, pointer :: prelim_vof(:), min_dist(:), max_dist(:)
  real, contiguous, pointer :: dist_node(:)
!==============================================================================!

  call Work % Connect_Real_Cell(prelim_vof, min_dist, max_dist)
  call Work % Connect_Real_Node(dist_node)

  ! First take aliases
  Grid => Vof % pnt_grid

  ! Open file to read Ellipsoid parameters:
  call File % Open_For_Reading_Ascii('ellipsoid_parameters.ini', fu)

  call File % Read_Line(fu)
  read(line % tokens(1), *) n_ellipses

  do ee = 1, n_ellipses

    ! Initialize working arrays
    prelim_vof(:) = 0.0

    call File % Read_Line(fu)
    read(line % tokens(1), *) radius_1
    read(line % tokens(2), *) radius_2
    read(line % tokens(3), *) radius_3

    call File % Read_Line(fu)
    read(line % tokens(1), *) cent_x
    read(line % tokens(2), *) cent_y
    read(line % tokens(3), *) cent_z

    ! Normalized distance from ellipsoid center in nodes
    dist_node (:) = 0.0
    do n = 1, Grid % n_nodes
      dist_node(n) = sqrt(  ((Grid % xn(n) - cent_x) / radius_1)**2   &
                          + ((Grid % yn(n) - cent_y) / radius_2)**2   &
                          + ((Grid % zn(n) - cent_z) / radius_3)**2)
    end do

    ! Minimum and maximum normalized distance in cells.  Polyhedral cells
    ! store their node count as a negative number; without abs() the loop
    ! never ran for them and every polyhedral cell counted as inside.
    min_dist(:) = +HUGE
    max_dist(:) = -HUGE
    do c = Cells_In_Domain_And_Buffers()

      ! A polyhedral cell is flagged by a NEGATIVE node count, so the count
      ! must be taken in absolute value.  Without the abs() the loop runs zero
      ! times on every polyhedral cell, min_dist and max_dist keep their
      ! sentinels, and the test below reads -HUGE <= 1.0 as "cell entirely
      ! inside the bubble" -- which silently fills the whole domain with
      ! vapour instead of raising anything.
      do i_nod = 1, abs(Grid % cells_n_nodes(c))
        n = Grid % cells_n(i_nod, c)

        min_dist(c)= min(dist_node(n), min_dist(c))
        max_dist(c)= max(dist_node(n), max_dist(c))
      end do

      ! Refuse to classify a cell whose nodes were never visited
      Assert(max_dist(c) > -HUGE*0.5)
    end do

    ! Volume fraction of cut hexahedra by sampling N_SUB**3 points of the
    ! cell's bounding box.  The linear estimate (1 - d_min) / (d_max - d_min)
    ! is exact only for a plane along the cell diagonal; on a sphere it left
    ! the first reconstruction with 33-38 % excess interface area, which
    ! the phase-change source (applied over that area) turned into a start-up
    ! over-delivery.  The bounding box is the cell only on the axis-aligned
    ! hexahedra of the uniform meshes; other cells keep the linear estimate.
    do c = Cells_In_Domain_And_Buffers()

      ! Since surface is at 1.0 this checks if cell crosses the surface
      if (min_dist(c) < 1.0 .and. max_dist(c) > 1.0 .and.  &
          Grid % cells_n_nodes(c) .ne. 8) then
        prelim_vof(c) = (1.0 - min_dist(c))  &
                      / (max_dist(c)-min_dist(c))

      else if (min_dist(c) < 1.0 .and. max_dist(c) > 1.0) then
        x_min = +HUGE;  x_max = -HUGE
        y_min = +HUGE;  y_max = -HUGE
        z_min = +HUGE;  z_max = -HUGE
        do i_nod = 1, Grid % cells_n_nodes(c)
          n = Grid % cells_n(i_nod, c)
          x_min = min(x_min, Grid % xn(n));  x_max = max(x_max, Grid % xn(n))
          y_min = min(y_min, Grid % yn(n));  y_max = max(y_max, Grid % yn(n))
          z_min = min(z_min, Grid % zn(n));  z_max = max(z_max, Grid % zn(n))
        end do
        n_in = 0
        do k = 1, N_SUB
          zs = z_min + (real(k) - 0.5) * (z_max - z_min) / real(N_SUB)
          do j = 1, N_SUB
            ys = y_min + (real(j) - 0.5) * (y_max - y_min) / real(N_SUB)
            do i = 1, N_SUB
              xs = x_min + (real(i) - 0.5) * (x_max - x_min) / real(N_SUB)
              if(  ((xs - cent_x) / radius_1)**2                  &
                 + ((ys - cent_y) / radius_2)**2                  &
                 + ((zs - cent_z) / radius_3)**2 < 1.0) n_in = n_in + 1
            end do
          end do
        end do
        prelim_vof(c) = real(n_in) / real(N_SUB**3)

      else if (max_dist(c) <= 1.0) then
        prelim_vof(c) = 1.0
      end if

    end do

    ! Precision (what on earth?)
    do c = Cells_In_Domain_And_Buffers()
      Vof % fun % n(c) = max(prelim_vof(c), Vof % fun % n(c))
    end do

  end do

  close(fu)

  ! Reverse
  do c = Cells_In_Domain_And_Buffers()
    Vof % fun % n(c) = 1.0 - Vof % fun % n(c)
  end do

  call Work % Disconnect_Real_Cell(prelim_vof, min_dist, max_dist)
  call Work % Disconnect_Real_Node(dist_node)

  end subroutine
