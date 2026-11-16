!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Vof(Vof, Sol)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Vof function.                !
!                                                                              !
!   The paper's hook (one-dimensionality: the axis cell's value is copied to   !
!   its transverse line), unchanged, followed by a record of the mass transfer !
!   the volume-fraction source of this step used.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Field_Type),  pointer :: Flow
  type(Var_Type),    pointer :: fun
  integer                    :: c, d, fu
  real                       :: x_ref, f_ref, mdot_vof, vol_liq
!==============================================================================!

  ! Take aliases
  Flow => Vof  % pnt_flow
  Grid => Flow % pnt_grid
  fun  => Vof % fun

  do c = Cells_In_Domain_And_Buffers()

    ! Find relevant cell
    if( Math % Approx_Real(Grid % yc(c), 0.0) .and.  &
        Math % Approx_Real(Grid % zc(c), 0.0) ) then
      x_ref = Grid % xc(c)
      f_ref = fun % n(c)

      ! Browse through all other cells and homogenize the values
      do d = 1, Grid % n_cells
        if(d .ne. c) then
          if(Math % Approx_Real(Grid % xc(d), x_ref)) then
            fun % n(d) = f_ref
          end if
        end if
      end do
    end if
  end do

  ! Mass transfer the source of this VOF step was built from, written to
  ! caf-vof.dat (this hook is compiled into Vof_Mod and cannot see User_Mod)
  mdot_vof = 0.0
  vol_liq  = 0.0
  do c = Cells_In_Domain()
    mdot_vof = mdot_vof + Vof % m_dot(c)
    vol_liq  = vol_liq  + fun % n(c) * Grid % vol(c)
  end do
  call Global % Sum_Real(mdot_vof)
  call Global % Sum_Real(vol_liq)
  if(First_Proc()) then
    call File % Append_For_Writing_Ascii('caf-vof.dat', fu)
    write(fu,'(i8,3(1x,es22.14e3))') Time % Curr_Dt(), Time % Get_Time(),  &
                                     mdot_vof, vol_liq
    close(fu)
  end if

  end subroutine
