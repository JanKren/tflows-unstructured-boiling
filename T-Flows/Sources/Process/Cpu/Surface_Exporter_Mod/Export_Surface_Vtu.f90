!==============================================================================!
  subroutine Export_Surface_Vtu(Exporter, cut_index)
!------------------------------------------------------------------------------!
!   Export extracted surface to VTU format (parallel-aware)                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surface_Exporter_Type), intent(inout) :: Exporter
  integer,                     intent(in)    :: cut_index
!-----------------------------------[Locals]-----------------------------------!
  character(SL)                  :: pvtu_name, vtu_name, piece_name
  integer                        :: fu_pvtu, fu_vtu, i, field_count, proc
!==============================================================================!

  if (cut_index < 1 .or. cut_index > Exporter % n_cuts) return

  if (.not. Exporter % cuts(cut_index) % active) return

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  !--------------------------------!
  !   Create PVTU and VTU files    !
  !--------------------------------!

  ! Set file names following T-Flows pattern
  call File % Set_Name(pvtu_name,                    &
                       time_step = Time % Curr_Dt(), &
                       appendix  = '-' // trim(Exporter % cuts(cut_index) % name), &
                       extension = '.pvtu',           &
                       domain    = 1)

  call File % Set_Name(vtu_name,                               &
                       processor = (/This_Proc(), N_Procs()/), &
                       time_step = Time % Curr_Dt(),           &
                       appendix  = '-' // trim(Exporter % cuts(cut_index) % name), &
                       extension = '.vtu',                     &
                       domain    = 1)


  !--------------------------!
  !   Write PVTU file        !
  !--------------------------!
  ! Only first processor writes the PVTU file (for parallel runs)
  if(First_Proc()) then
    call File % Open_For_Writing_Ascii(pvtu_name, fu_pvtu)

    ! Write PVTU header
    write(fu_pvtu,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(fu_pvtu,'(a,a)') IN_0, '<VTKFile type="PUnstructuredGrid" version="1.0"' //  &
                                 ' byte_order="LittleEndian">'
    write(fu_pvtu,'(a,a)') IN_1, '<PUnstructuredGrid GhostLevel="0">'

    ! Write field information
    write(fu_pvtu,'(a,a)') IN_2, '<PPointData>'
    write(fu_pvtu,'(a,a)') IN_3, '<DataArray type="Int32" Name="PointIndex"/>'

    ! Add field data arrays - use standard field names regardless of local point count
    ! This ensures PVTU has consistent field declarations across all processors
    write(fu_pvtu,'(a,a)') IN_3, '<DataArray type="Float64" Name="Velocity_U"/>'
    write(fu_pvtu,'(a,a)') IN_3, '<DataArray type="Float64" Name="Velocity_V"/>'
    write(fu_pvtu,'(a,a)') IN_3, '<DataArray type="Float64" Name="Velocity_W"/>'
    write(fu_pvtu,'(a,a)') IN_3, '<DataArray type="Float64" Name="Pressure"/>'
    write(fu_pvtu,'(a,a)') IN_3, '<DataArray type="Float64" Name="Temperature"/>'
    write(fu_pvtu,'(a,a)') IN_3, '<DataArray type="Float64" Name="VOF_Sharp"/>'
    write(fu_pvtu,'(a,a)') IN_3, '<DataArray type="Float64" Name="Mass_Transfer_Rate"/>'

    write(fu_pvtu,'(a,a)') IN_2, '</PPointData>'
    write(fu_pvtu,'(a,a)') IN_2, '<PPoints>'
    write(fu_pvtu,'(a,a)') IN_3, '<DataArray type="Float64" NumberOfComponents="3"/>'
    write(fu_pvtu,'(a,a)') IN_2, '</PPoints>'

    ! Write piece information for each processor
    do proc = 1, N_Procs()
      call File % Set_Name(piece_name,                      &
                           processor = (/proc, N_Procs()/), &
                           time_step = Time % Curr_Dt(),    &
                           appendix  = '-' // trim(Exporter % cuts(cut_index) % name), &
                           extension = '.vtu',              &
                           domain    = 1)
      write(fu_pvtu,'(a,a)') IN_2, '<Piece Source="' // trim(piece_name) // '"/>'
    end do

    write(fu_pvtu,'(a,a)') IN_1, '</PUnstructuredGrid>'
    write(fu_pvtu,'(a,a)') IN_0, '</VTKFile>'
    close(fu_pvtu)
  end if

  !--------------------------!
  !   Write VTU file         !
  !--------------------------!
  ! Each processor writes its own VTU file (even if no points)
  call File % Open_For_Writing_Ascii(vtu_name, fu_vtu)

  ! Write VTU header
  write(fu_vtu,'(a,a)') IN_0, '<?xml version="1.0"?>'
  write(fu_vtu,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid"'  //  &
                              ' version="0.1"'                    //  &
                              ' byte_order="LittleEndian">'
  write(fu_vtu,'(a,a)') IN_1, '<UnstructuredGrid>'

  write(fu_vtu,'(a,a,i0,a,i0,a)')   &
                IN_2, '<Piece NumberOfPoints="', Exporter % cuts(cut_index) % n_points,  &
                           '" NumberOfCells="', Exporter % cuts(cut_index) % n_points, '">'

  ! Write points (even if zero points)
  write(fu_vtu,'(a,a)') IN_3, '<Points>'
  write(fu_vtu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                              ' NumberOfComponents'       //  &
                              '="3" format="ascii">'
  if (Exporter % cuts(cut_index) % n_points > 0) then
    do i = 1, Exporter % cuts(cut_index) % n_points
      write(fu_vtu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')            &
                    IN_5, Exporter % cuts(cut_index) % points(1,i), &
                          Exporter % cuts(cut_index) % points(2,i), &
                          Exporter % cuts(cut_index) % points(3,i)
    end do
  end if
  write(fu_vtu,'(a,a)') IN_4, '</DataArray>'
  write(fu_vtu,'(a,a)') IN_3, '</Points>'

  ! Write point data
  write(fu_vtu,'(a,a)') IN_3, '<PointData>'

  ! Write PointIndex
  write(fu_vtu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                              ' Name="PointIndex" '      //  &
                              'format="ascii">'
  if (Exporter % cuts(cut_index) % n_points > 0) then
    do i = 1, Exporter % cuts(cut_index) % n_points
      write(fu_vtu,'(a,i9)') IN_5, i
    end do
  end if
  write(fu_vtu,'(a,a)') IN_4, '</DataArray>'

  ! Write field data if available
  if (Exporter % cuts(cut_index) % n_fields > 0 .and. Exporter % cuts(cut_index) % n_points > 0) then
    do field_count = 1, Exporter % cuts(cut_index) % n_fields
      write(fu_vtu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                                  ' Name="' // trim(Exporter % cuts(cut_index) % field_names(field_count)) // '"' //  &
                                  ' format="ascii">'
      do i = 1, Exporter % cuts(cut_index) % n_points
        write(fu_vtu,'(a,1pe16.6e4)') IN_5, Exporter % cuts(cut_index) % field_data(field_count, i)
      end do
      write(fu_vtu,'(a,a)') IN_4, '</DataArray>'
    end do
  end if

  write(fu_vtu,'(a,a)') IN_3, '</PointData>'

  ! Write cells (for point cloud, use vertices)
  write(fu_vtu,'(a,a)') IN_3, '<Cells>'
  write(fu_vtu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                              ' Name="connectivity"'    //  &
                              ' format="ascii">'
  if (Exporter % cuts(cut_index) % n_points > 0) then
    ! Create individual vertices for point cloud
    do i = 1, Exporter % cuts(cut_index) % n_points
      write(fu_vtu,'(a,i9)') IN_5, i-1  ! VTK uses 0-based indexing
    end do
  end if
  write(fu_vtu,'(a,a)') IN_4, '</DataArray>'

  ! Cell offsets
  write(fu_vtu,'(a,a)') IN_4, '<DataArray type='//intp//' Name="offsets"' //  &
                              ' format="ascii">'
  if (Exporter % cuts(cut_index) % n_points > 0) then
    do i = 1, Exporter % cuts(cut_index) % n_points
      write(fu_vtu,'(a,i9)') IN_5, i
    end do
  end if
  write(fu_vtu,'(a,a)') IN_4, '</DataArray>'

  ! Cell types (VTK_VERTEX = 1)
  write(fu_vtu,'(a,a)') IN_4, '<DataArray type='//intp//' Name="types"' //  &
                              ' format="ascii">'
  if (Exporter % cuts(cut_index) % n_points > 0) then
    do i = 1, Exporter % cuts(cut_index) % n_points
      write(fu_vtu,'(a,i9)') IN_5, VTK_VERTEX
    end do
  end if
  write(fu_vtu,'(a,a)') IN_4, '</DataArray>'
  write(fu_vtu,'(a,a)') IN_3, '</Cells>'

  write(fu_vtu,'(a,a)') IN_3, '<CellData>'
  write(fu_vtu,'(a,a)') IN_3, '</CellData>'

  write(fu_vtu,'(a,a)') IN_2, '</Piece>'
  write(fu_vtu,'(a,a)') IN_1, '</UnstructuredGrid>'
  write(fu_vtu,'(a,a)') IN_0, '</VTKFile>'
  close(fu_vtu)

  if (First_Proc()) then
    print *, 'Exported surface to: ', trim(pvtu_name), ' (and processor-specific VTU files)'
  end if

  end subroutine

