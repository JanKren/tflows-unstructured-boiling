!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Vof(Vof, Sol)
!------------------------------------------------------------------------------!
!   Mean-preserving variant of the Stefan one-dimensionality enforcement.       !
!                                                                              !
!   The shipped hook holds the problem one-dimensional by copying the value of  !
!   the axis cell (yc = zc = 0) onto every cell sharing its x.  That is a       !
!   projection, but not a conservative one: one cell's alpha replaces the whole !
!   line, so the line's volume integral of alpha -- its mass -- changes.        !
!   Measured against the global budget it costs about 1e-3 in M_s/M_0.          !
!                                                                              !
!   This variant enforces exactly the same one-dimensionality, but assigns the  !
!   VOLUME-WEIGHTED MEAN of the line instead of one member of it.  The mean is  !
!   the orthogonal projection onto x-only fields and leaves sum(alpha V)        !
!   unchanged by construction, so it is conservative.                           !
!                                                                              !
!   The two differ only in that choice, which is what makes the pair a          !
!   controlled measurement of what the enforcement itself costs.                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Field_Type),  pointer :: Flow
  type(Var_Type),    pointer :: fun
  integer                    :: c, d
  real                       :: x_ref, f_sum, v_sum, f_avg
!==============================================================================!

  Flow => Vof  % pnt_flow
  Grid => Flow % pnt_grid
  fun  => Vof % fun

  do c = Cells_In_Domain_And_Buffers()
    if( Math % Approx_Real(Grid % yc(c), 0.0) .and.  &
        Math % Approx_Real(Grid % zc(c), 0.0) ) then
      x_ref = Grid % xc(c)

      f_sum = 0.0
      v_sum = 0.0
      do d = 1, Grid % n_cells
        if(Math % Approx_Real(Grid % xc(d), x_ref)) then
          f_sum = f_sum + fun % n(d) * Grid % vol(d)
          v_sum = v_sum + Grid % vol(d)
        end if
      end do

      if(v_sum > TINY) then
        f_avg = f_sum / v_sum
        do d = 1, Grid % n_cells
          if(Math % Approx_Real(Grid % xc(d), x_ref)) then
            fun % n(d) = f_avg
          end if
        end do
      end if

    end if
  end do

  end subroutine
