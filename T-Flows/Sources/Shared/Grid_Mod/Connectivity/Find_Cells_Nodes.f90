!==============================================================================!
  subroutine Find_Cells_Nodes(Grid)
!------------------------------------------------------------------------------!
!>  Rebuilds cells_n and cells_n_nodes for every cell as the union of the nodes
!>  of the faces which bound it.
!------------------------------------------------------------------------------!
!   For any closed polyhedron whatsoever, the node list of the cell IS the     !
!   union of the node lists of its faces.  That makes this routine a complete  !
!   and shape-agnostic way to recover cells_n, which is exactly what is wanted !
!   after refinement: the alternative, patching cells_n cell by cell as faces  !
!   are split, has to get the bookkeeping right separately for the cell that   !
!   was refined and for its unrefined neighbours, and is easy to get subtly    !
!   wrong.                                                                     !
!                                                                              !
!   The node ORDER that comes out is arbitrary, so the cells are flagged       !
!   polyhedral (cells_n_nodes negative) and the output paths use the general   !
!   face-based route.  For a grid that has just been refined that is not a     !
!   loss of information but a statement of fact - the cells along a refinement !
!   front genuinely are polyhedra, whatever they started life as.              !
!                                                                              !
!   Boundary cells are left alone.  They are one per boundary face and         !
!   Sort_Faces_By_Region has already given them their face's nodes.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), intent(inout) :: Grid  !! computational grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: c, s, i_fac, i_nod, n, cnt, max_cn
  integer, allocatable :: seen(:)
!==============================================================================!

  call Profiler % Start('Find_Cells_Nodes')

  ! "seen" marks the cell for which a node has already been collected, so a
  ! node shared by several faces of the same cell is stored only once.  The
  ! stamp is the cell number itself, which saves clearing the array per cell.
  allocate(seen(Grid % n_nodes));  seen(:) = 0

  !-----------------------------------------------------!
  !   First pass: how wide does cells_n have to be?     !
  !-----------------------------------------------------!
  max_cn = 0
  do c = 1, Grid % n_cells
    cnt = 0
    do i_fac = 1, abs(Grid % cells_n_faces(c))
      s = Grid % cells_f(i_fac, c)
      do i_nod = 1, Grid % faces_n_nodes(s)
        n = Grid % faces_n(i_nod, s)
        if(seen(n) .ne. c) then
          seen(n) = c
          cnt = cnt + 1
        end if
      end do
    end do
    max_cn = max(max_cn, cnt)
  end do

  call Enlarge % Matrix_Int(Grid % cells_n,               &
                            i=(/1, max(max_cn, 1)/),      &
                            j=(/-Grid % n_bnd_cells, Grid % n_cells/))

  !--------------------------------------!
  !   Second pass: store them for real   !
  !--------------------------------------!
  seen(:) = 0

  do c = 1, Grid % n_cells
    cnt = 0
    do i_fac = 1, abs(Grid % cells_n_faces(c))
      s = Grid % cells_f(i_fac, c)
      do i_nod = 1, Grid % faces_n_nodes(s)
        n = Grid % faces_n(i_nod, s)
        if(seen(n) .ne. c) then
          seen(n) = c
          cnt = cnt + 1
          Grid % cells_n(cnt, c) = n
        end if
      end do
    end do

    ! Negative means polyhedral, which is what an arbitrary node order is
    Grid % cells_n_nodes(c) = -cnt
  end do

  Grid % polyhedral = .true.

  deallocate(seen)

  call Profiler % Stop('Find_Cells_Nodes')

  end subroutine
