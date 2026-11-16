!==============================================================================!
  subroutine Allocate_Nodes(Grid, nn, margin)
!------------------------------------------------------------------------------!
!>  Allocates memory for node-based data (arrays and matrices), for geometrical
!>  (xn, yn, zn ...) and connectivity data (new_n, old_n, Comm % node_glo).
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)              :: Grid    !! grid being expanded
  integer, intent(in)           :: nn      !! number of nodes in the grid
  integer, intent(in), optional :: margin  !! margin for allocation
!-----------------------------------[Locals]-----------------------------------!
  integer :: n, nn_m, n_old
!==============================================================================!

  ! Generator has growing number of faces, don't set them here
  if(PROGRAM_NAME(1:6) .ne. 'Genera') then
    Grid % n_nodes = nn
  end if

  ! If node-based arrays are allocated and they
  ! are bigger than requested, get out of here
  if(allocated(Grid % xn)) then
    if(nn .le. ubound(Grid % xn, 1)) then
      return
    end if
  end if

  ! Process the margin if specified
  nn_m = nn
  if(present(margin)) then
    Assert(margin .ge. 0)
    nn_m = nn + margin
  end if

  ! I can't figure out how to print something meaningful in parallel here
  if(Sequential_Run()) then
    print '(a,i9)', ' # Expanding memory for nodes to size: ', nn_m
  end if

  ! Allocate memory for node coordinates
  call Enlarge % Array_Real(Grid % xn, i=(/1,nn_m/))
  call Enlarge % Array_Real(Grid % yn, i=(/1,nn_m/))
  call Enlarge % Array_Real(Grid % zn, i=(/1,nn_m/))

  !----------------------------------------------------------------------!
  !   Global node numbers.  Only the NEW tail is given the identity.     !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - +
  !   The identity is the right answer for a grid which is being built here    !
  !   (Convert, Generate) and for a sequential run, where local and global     !
  !   numbering are the same thing.  It is emphatically NOT the right answer   !
  !   for a sub-domain: Load_Cfn reads the real numbering from the .cfn file   !
  !   right after calling this routine, and a later call - refinement enlarges !
  !   the node arrays on every sweep - used to overwrite the whole array and   !
  !   throw that numbering away.                                              !
  !                                                                            !
  !   The damage was silent.  Subdivide_Cells orders a marked cell's nodes by  !
  !   global number so that two ranks holding the same cell number its         !
  !   children identically; with node_glo reset to the identity that sort ran  !
  !   on LOCAL numbers, which differ between ranks, and the agreement the       !
  !   whole parallel refinement design rests on was down to luck.               !
  !                                                                            !
  !   Enlarge % Array_Int copies the existing contents across, so filling only  !
  !   the tail leaves every program that never had anything but the identity    !
  !   with exactly the array it had before.                                     !
  !----------------------------------------------------------------------------!
  n_old = 0
  if(allocated(Grid % Comm % node_glo)) then
    n_old = ubound(Grid % Comm % node_glo, 1)
  end if

  call Enlarge % Array_Int(Grid % Comm % node_glo, i=(/1,nn_m/))
  do n = n_old + 1, nn_m
    Grid % Comm % node_glo(n) = n
  end do

  ! Allocate new and old numbers (this is so often used, maybe is better here)
  if(PROGRAM_NAME(1:6) .ne. 'Proces') then
    call Enlarge % Array_Int(Grid % new_n, i=(/1,nn_m/))
    call Enlarge % Array_Int(Grid % old_n, i=(/1,nn_m/))
  end if

  if(PROGRAM_NAME(1:6) .eq. 'Genera') then
    ! This variable used to be in Gen_Mod for a long time and
    ! is used only in Generate to take care of periodicity
    call Enlarge % Matrix_Int(Grid % twin_n, i=(/1,nn_m/), j=(/0,8/))

    ! This is used only in initial stages of Generate, inside Domain_Mod really
    call Enlarge % Array_Log(Grid % node_positioned, i=(/1,nn_m/))
  end if

  ! Nodal coordinates of homogeneous planes
  Assert(Grid % n_x_planes .ge. 0)
  Assert(Grid % n_y_planes .ge. 0)
  Assert(Grid % n_z_planes .ge. 0)
  if(Grid % n_x_planes .gt. 0) then
    call Enlarge % Array_Real(Grid % x_coord_plane, i=(/1,Grid % n_x_planes/))
  end if
  if(Grid % n_y_planes .gt. 0) then
    call Enlarge % Array_Real(Grid % y_coord_plane, i=(/1,Grid % n_y_planes/))
  end if
  if(Grid % n_z_planes .gt. 0) then
    call Enlarge % Array_Real(Grid % z_coord_plane, i=(/1,Grid % n_z_planes/))
  end if

  end subroutine
