!==============================================================================!
  subroutine Allocate_Work(Work, Grid)
!------------------------------------------------------------------------------!
!>  This subroutine is responsible for allocating memory for various working
!>  arrays within the Work object.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * The subroutine first records the required number of each type of array   !
!     in the Work object.                                                      !
!   * It then proceeds to call specific allocation subroutines (such as        !
!     Allocate_Real_Cell, Allocate_Real_Face, etc.) for each array type,       !
!     passing the Grid object and the number of arrays needed.                 !
!   * After the allocation process, the subroutine updates the Work object's   !
!     allocated status to true, signifying that the working arrays are now set !
!     up and ready for use.                                                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Work_Type) :: Work      !! parent; the singleton Work object
  type(Grid_Type)  :: Grid(:)   !! grids on which the Work will be used
!==============================================================================!

  !-------------------------------------------------!
  !   Fetch the maximum number of cells, boundary   !
  !   cells, faces and nodes across all the grids   !
  !-------------------------------------------------!

  ! Get the maximum number of cells and boundary cells
  Work % max_nc = maxval(Grid(1:size(Grid)) % n_cells)
  Work % max_nb = maxval(Grid(1:size(Grid)) % n_bnd_cells)

  ! Get the maximum number of faces
  Work % max_nf = maxval(Grid(1:size(Grid)) % n_faces)

  ! Get the maximum number of nodes
  Work % max_nn = maxval(Grid(1:size(Grid)) % n_nodes)

  !------
  !
  !------

  !--------------------------------------------------------------------!
  !   Made re-callable, so a refinement sweep can resize the work      !
  !   space.  Dropping the containers drops the arrays inside them     !
  !   too, which is the point: those were sized from the old grid and  !
  !   are shorter than the loops that would index them.                !
  !--------------------------------------------------------------------!
  if(allocated(Work % i_cell)) deallocate(Work % i_cell)
  if(allocated(Work % i_face)) deallocate(Work % i_face)
  if(allocated(Work % i_node)) deallocate(Work % i_node)
  if(allocated(Work % r_cell)) deallocate(Work % r_cell)
  if(allocated(Work % r_face)) deallocate(Work % r_face)
  if(allocated(Work % r_node)) deallocate(Work % r_node)

  ! Allocate the basic container's work space
  allocate(Work % i_cell(MAX_WORK_ARRAYS));  Work % last_i_cell = 0
  allocate(Work % i_face(MAX_WORK_ARRAYS));  Work % last_i_face = 0
  allocate(Work % i_node(MAX_WORK_ARRAYS));  Work % last_i_node = 0
  allocate(Work % r_cell(MAX_WORK_ARRAYS));  Work % last_r_cell = 0
  allocate(Work % r_face(MAX_WORK_ARRAYS));  Work % last_r_face = 0
  allocate(Work % r_node(MAX_WORK_ARRAYS));  Work % last_r_node = 0

  end subroutine
