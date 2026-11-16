!==============================================================================!
  subroutine Find_Cells_Cells(Grid)
!------------------------------------------------------------------------------!
!>  This subroutine searches for cells neighbours (cells_c, based on the
!>  information already stored in the faces surrounding each cell (cells_f)
!------------------------------------------------------------------------------!
!   Process                                                                    !
!                                                                              !
!   * Initialization:                                                          !
!     - Resets the count of cells' neighbours to zero for all cells.           !
!   * Cell Mapping:                                                            !
!     - Browse through all cells (boundary and inside).                        !
!     - For each cell, browse through its faces (i_fac)                        !
!     - Use data stored in matrix faces_c (cells surrounding each face to      !
!       fill up the cells_c matrices.                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type) :: Grid  !! grid under consideration
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, c1, c2, s, i_fac, i_cel
!==============================================================================!

  ! This function is supposed to be used from Processor only
  Assert(PROGRAM_NAME(1:6) .eq. 'Proces')

  call Profiler % Start('Find_Cells_Cells')

  ! Reset number of neighbours for each cell
  Grid % cells_n_cells(:) = 0
  !Grid % cells_c = Grid % cells_f
  !------------------------------!
  !   Browse through all cells   !
  !------------------------------!
  do c = -Grid % n_bnd_cells, Grid % n_cells
    Grid % cells_n_cells(c) = Grid % cells_n_faces(c)
    !write(6,*) 'Cells-n-cells:',  Grid % cells_n_cells(c)
    if(Grid % cells_n_cells(c) .eq. 0) then
      Grid % cells_n_cells(c) = 1
    endif
    call Enlarge % Matrix_Int(Grid % cells_c, i=(/1,Grid % cells_n_cells(c)/))

    ! Browse through faces of the cell c
    do i_fac = 1, Grid % cells_n_faces(c)
      i_cel = i_fac
      s  = Grid % cells_f(i_fac, c)  ! take true face number
      c1 = Grid % faces_c(1, s)      ! face's c1
      c2 = Grid % faces_c(2, s)      ! face's c2
      if(c1 .eq. c) Grid % cells_c(i_cel, c) = c2
      if(c2 .eq. c) Grid % cells_c(i_cel, c) = c1
    end do
  end do

  call Profiler % Stop('Find_Cells_Cells')

  end subroutine
