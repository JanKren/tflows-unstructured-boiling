!==============================================================================!
  subroutine Save_Vtu_Scalar_Int(Results,                &
                                 var_name, plot_inside,  &
                                 val, fs, fp,            &
                                 data_offset, sweep)
!------------------------------------------------------------------------------!
!>  The Save_Vtu_Scalar_Int subroutine is designed to write a single intefer
!>  scalar variable defined over cells into a .vtu file as a single record.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * Preparatory work: Writes the supplied integer record directly      !
!     and sets the plotting precision.                                         !
!   * Data handling: Based on the 'sweep' stage, it either writes the header   !
!     data for the variable or the actual scalar values to the file.           !
!   * Data writing: Scalar values are written into the file in an appended     !
!     format, taking into account whether the data is for internal or boundary !
!     cells.                                                                   !
!   * Data offset management: Updates the data offset after writing, ensuring  !
!     accurate placement of subsequent data in the file.                       !
!------------------------------------------------------------------------------!
!   Workflow                                                                   !
!                                                                              !
!   * Input section: The supplied cells are streamed without a cell buffer.    !
!   * Header writing: In the first sweep, the subroutine writes the XML header !
!     for the scalar variable, including its name and data type.               !
!   * Data writing: In the second sweep, it writes the actual scalar values.   !
!     This includes processing values for either internal or boundary cells    !
!     based on the 'plot_inside' flag.                                         !
!   * Closing operations: Finally, it updates the appended-data offset.        !
!                                                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Results_Type) :: Results      !! parent class
  character(len=*)    :: var_name     !! name of the variabl
  logical             :: plot_inside  !! true to plots inside,
                                      !! false to plot on the boundary
  integer             :: val(:)       !! variable's values
  integer             :: fs, fp       !! file unit sequential and parallel
  integer             :: data_offset  !! data offset in the .vtu file
  integer             :: sweep        !! is it the first or second sweep
!-----------------------------------[Locals]-----------------------------------!
  integer(SP)                  :: data_size
  integer                      :: c2, c_f, c_l
  character(SL)                :: str1
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Results)
!==============================================================================!


  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  data_size = 0

  c_f = lbound(val, 1)
  c_l = ubound(val, 1)

  ! Header
  if(sweep .eq. 1) then
    if(Parallel_Run() .and. First_Proc()) then
      write(fs) IN_4                                  //  &
                '<PDataArray type='//intp             //  &
                ' Name="' // trim(var_name) // '"/>'  // LF
    end if

    write(str1, '(i0.0)') data_offset
    if(data_offset .eq. 0) write(str1, '(i1)') data_offset
    write(fp) IN_4 // '<DataArray type='//intp           //  &
                      ' Name="' // trim(var_name) // '"' //  &
                      ' format="appended"'               //  &
                      ' offset="' // trim(str1) //'">'   // LF
    write(fp) IN_4 // '</DataArray>' // LF
  end if

  ! Data
  if(sweep .eq. 2) then

    call Profiler % Start('Save_Vtu_Results (scalar int)')

    ! The input already contains exactly the requested cells. Boundary
    ! records can be longer than the positive range of a cell work array
    ! on thin meshes or small MPI partitions; copying to buffer(1:n) then
    ! overwrites the heap. Stream the input section directly, preserving
    ! the existing payload order and precision (including strided inputs).
    data_size = int((c_l-c_f+1) * IP, SP)
    write(fp) data_size
    write(fp) val(c_f:c_l)

    call Profiler % Stop('Save_Vtu_Results (scalar int)')

  end if

  ! Update data_offset
  if(sweep .eq. 1) then
    if(plot_inside) then
      data_offset = data_offset + (c_l-c_f+1) * IP
    else
      do c2 = c_f, c_l
        data_offset = data_offset + IP
      end do
    end if
    data_offset = data_offset + SP
  end if


  end subroutine
