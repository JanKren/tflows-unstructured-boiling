!==============================================================================!
  subroutine User_Mod_Beginning_Of_Iteration(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   Uniform-temperature check (caf_mode 1): impose the interfacial mass flux   !
!   on the cells that hold a front element.  The pressure source of this       !
!   iteration then makes div u non-zero there, while the mass transfer         !
!   estimated from the (uniform) temperature stays zero.                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c, e
  real                     :: area
!==============================================================================!

  Grid => Flow % pnt_grid

  if(caf_mode .ne. 1) return

  do c = Cells_In_Domain_And_Buffers()
    e = Vof % Front % elem_in_cell(c)
    if(e .ne. 0) then
      area = sqrt(  Vof % Front % Elem(e) % sx**2  &
                  + Vof % Front % Elem(e) % sy**2  &
                  + Vof % Front % Elem(e) % sz**2)
      Vof % m_dot(c) = caf_flux * area
    end if
  end do

  end subroutine
