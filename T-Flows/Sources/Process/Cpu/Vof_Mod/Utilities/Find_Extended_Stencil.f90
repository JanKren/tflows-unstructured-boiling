!==============================================================================!
  subroutine Find_Extended_Stencil(Vof)
!------------------------------------------------------------------------------!
!   Builds, for every cell, the list of node-connected neighbours which are    !
!   NOT already face neighbours.  On a hexahedral mesh these are the 12 edge   !
!   and 8 corner neighbours, so the front-modified least-squares stencil can   !
!   be widened from 6 to 26 points.                                            !
!                                                                              !
!   The list is a pure mesh property, so it is built once at start-up.         !
!                                                                              !
!   Face neighbours are deliberately excluded: the gradient routines keep      !
!   their original face loop untouched and merely add these extra terms on     !
!   top of it.  With EXTENDED_GRADIENT_STENCIL off, nothing here is used and   !
!   the gradient is bit-identical to the face-only formulation.                !
!                                                                              !
!   Periodicity: Find_Nodes_Cells merges node cell lists across periodic       !
!   boundaries, but the coordinates of those partners are not shifted, so a    !
!   raw centre-to-centre vector would span the whole domain.  Such pairs are   !
!   rejected by a distance guard (EXT_MAX_DIST cell sizes).  Face neighbours   !
!   are unaffected because they keep using Grid % dx(s), which is periodic-    !
!   aware.                                                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!------------------------------[Local parameters]------------------------------!
  real, parameter :: EXT_MAX_DIST = 2.2  !! reject neighbours farther than this
                                         !! many local cell sizes.  A hexahedral
                                         !! corner neighbour sits at sqrt(3) =
                                         !! 1.73, so 2.2 keeps every genuine
                                         !! neighbour with margin.  It must stay
                                         !! well below the smallest periodic
                                         !! period: a slab three cells thick has
                                         !! its unshifted periodic partner at
                                         !! exactly 3 cell sizes, which a
                                         !! threshold of 3.0 would NOT reject.
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c, cn, i_nod, i_cel, i_fac, n, run, cnt, max_ext
  integer                  :: nc, nb, s, c1, c2, max_fn
  real                     :: dx, dy, dz, dist, cell_size, guard
  integer, allocatable     :: seen(:), n_fn(:), fn(:,:)
!==============================================================================!

  Grid => Vof % pnt_grid

  nc = Grid % n_cells
  nb = Grid % n_bnd_cells

  !---------------------------------------------------------------!
  !   Nothing to do unless the extended stencil was switched on   !
  !---------------------------------------------------------------!
  if(.not. Vof % ext_stencil) return

  !-------------------------------------------------------------------!
  !   Build the face adjacency locally.                               !
  !                                                                   !
  !   Grid % cells_c cannot be used here: it is filled by             !
  !   Find_Cells_Cells, which needs Find_Cells_Faces, and Main_Pro    !
  !   does not call those.  The array is empty, so the face           !
  !   neighbours are collected here in one sweep over the faces       !
  !   instead.                                                        !
  !-------------------------------------------------------------------!
  allocate(n_fn(1:nc))
  n_fn(:) = 0

  do run = 1, 2
    n_fn(:) = 0
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c1 > 0 .and. c1 <= nc .and. c2 > 0 .and. c2 <= nc) then
        n_fn(c1) = n_fn(c1) + 1
        n_fn(c2) = n_fn(c2) + 1
        if(run .eq. 2) then
          fn(n_fn(c1), c1) = c2
          fn(n_fn(c2), c2) = c1
        end if
      end if
    end do
    if(run .eq. 1) then
      max_fn = 1
      if(nc > 0) max_fn = max(maxval(n_fn(1:nc)), 1)
      allocate(fn(1:max_fn, 1:nc))
      fn(:,:) = 0
    end if
  end do

  !--------------------------------------------------------------!
  !   "seen" marks cells already collected for the current one   !
  !   so the same neighbour is not stored twice                  !
  !--------------------------------------------------------------!
  allocate(seen(-nb:nc))
  seen(:) = 0

  allocate(Vof % ext_n(1:nc))
  Vof % ext_n(:) = 0

  !-----------------------------------------------------------!
  !   Run 1: count the extra neighbours to size the storage   !
  !   Run 2: store them                                       !
  !-----------------------------------------------------------!
  do run = 1, 2

    ! Only cells inside the domain: the gradient routines browse with
    ! Cells_In_Domain() and buffer gradients are overwritten by the
    ! exchange in Grad_Variable_With_Front, so buffer lists are never
    ! read.  Building them would also inflate the count, because a
    ! buffer cell's faces are not all stored on this rank, so the face
    ! sweep above under-counts its face neighbours and too few of its
    ! node-connected cells get excluded (20 becomes 23-24 in parallel,
    ! while a serial run correctly reports 20).
    do c = Cells_In_Domain()

      cell_size = Grid % vol(c) ** (1.0 / 3.0)
      guard     = EXT_MAX_DIST * cell_size

      ! Mark this cell and all of its face neighbours as already covered
      seen(c) = c
      do i_fac = 1, n_fn(c)
        cn = fn(i_fac, c)
        if(cn .gt. 0) seen(cn) = c
      end do

      cnt = 0

      ! Browse all cells sharing a node with this cell
      do i_nod = 1, abs(Grid % cells_n_nodes(c))
        n = Grid % cells_n(i_nod, c)

        do i_cel = 1, Grid % nodes_n_cells(n)
          cn = Grid % nodes_c(i_cel, n)

          ! Skip boundary cells and anything already taken
          if(cn .le. 0)      cycle
          if(cn .gt. nc)     cycle
          if(seen(cn) .eq. c) cycle

          ! Reject unshifted periodic partners (and anything absurd)
          dx = Grid % xc(cn) - Grid % xc(c)
          dy = Grid % yc(cn) - Grid % yc(c)
          dz = Grid % zc(cn) - Grid % zc(c)
          dist = sqrt(dx*dx + dy*dy + dz*dz)
          if(dist .gt. guard) cycle

          ! Genuine edge or corner neighbour
          seen(cn) = c
          cnt = cnt + 1
          if(run .eq. 2) Vof % ext_c(cnt, c) = cn
        end do
      end do

      Vof % ext_n(c) = cnt

    end do  ! through cells

    ! Allocate storage between the two runs
    if(run .eq. 1) then
      max_ext = maxval(Vof % ext_n(1:nc))
      call Global % Max_Int(max_ext)
      max_ext = max(max_ext, 1)
      allocate(Vof % ext_c(1:max_ext, 1:nc))
      Vof % ext_c(:,:) = 0
      Vof % ext_n(:) = 0
      seen(:) = 0
    end if

  end do  ! runs

  deallocate(seen)
  deallocate(n_fn)
  deallocate(fn)

  !---------------------------!
  !   Report what was built   !
  !---------------------------!
  cnt = 0
  if(nc .gt. 0) cnt = maxval(Vof % ext_n(1:nc))
  call Global % Max_Int(cnt)
  if(First_Proc()) then
    print '(a,i0,a)', ' # EXTENDED STENCIL: up to ', cnt,  &
                      ' extra (edge/corner) neighbours per cell'
    print '(a,f5.2,a)', ' # EXTENDED STENCIL: LSQ weight exponent p = ',  &
                        Vof % ext_weight_p, ' (w = |r|^-p)'
  end if

  end subroutine
