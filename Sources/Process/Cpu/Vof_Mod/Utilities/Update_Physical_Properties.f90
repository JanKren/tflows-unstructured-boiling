!==============================================================================!
  subroutine Update_Physical_Properties(Vof)
!------------------------------------------------------------------------------!
!   Update physical properties based on volume fraction variable               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  type(Var_Type),   pointer :: col
  integer                   :: c, c1, c2, s, reg
  real                      :: f
!==============================================================================!

  ! Take aliases
  Flow => Vof  % pnt_flow
  Grid => Flow % pnt_grid
  col  => Vof % fun
  ! col  => Vof % smooth

  ! Density and viscosity in cells
  do c = Cells_In_Domain_And_Buffers()
    Flow % density(c)      = col % n(c)         * Vof % phase_dens(1)  &
                           + (1.0 - col % n(c)) * Vof % phase_dens(0)
    Flow % viscosity(c)    = col % n(c)         * Vof % phase_visc(1)  &
                           + (1.0 - col % n(c)) * Vof % phase_visc(0)
    if(Flow % heat_transfer) then
      Flow % capacity(c)     = col % n(c)         * Vof % phase_capa(1)  &
                             + (1.0 - col % n(c)) * Vof % phase_capa(0)
      Flow % conductivity(c) = col % n(c)         * Vof % phase_cond(1)  &
                             + (1.0 - col % n(c)) * Vof % phase_cond(0)
    end if
  end do

  !---------------------------------------------------------------!
  !   Properties in boundary cells.  Face coefficients are        !
  !   harmonic means with them (Turb_Mod/Face_Vis.f90 and         !
  !   Face_Cond_And_Stress.f90), so a boundary cell left at its   !
  !   initial zero removes the wall shear and wall conduction.    !
  !   Inflows carry the incoming mixture; every other boundary    !
  !   takes the material in the cell next to it, which is what    !
  !   a single-phase run sees with constant properties.           !
  !---------------------------------------------------------------!
  do reg = Boundary_Regions()
    do s = Faces_In_Region(reg)
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(Grid % region % type(reg) .eq. INFLOW) then
        f = col % n(c2)
      else
        f = col % n(c1)
      end if
      Flow % density(c2)   = f         * Vof % phase_dens(1)  &
                           + (1.0 - f) * Vof % phase_dens(0)
      Flow % viscosity(c2) = f         * Vof % phase_visc(1)  &
                           + (1.0 - f) * Vof % phase_visc(0)
      if(Flow % heat_transfer) then
        Flow % capacity(c2)     = f         * Vof % phase_capa(1)  &
                                + (1.0 - f) * Vof % phase_capa(0)
        Flow % conductivity(c2) = f         * Vof % phase_cond(1)  &
                                + (1.0 - f) * Vof % phase_cond(0)
      end if
    end do
  end do

  call Grid % Exchange_Cells_Real(Flow % density)
  call Grid % Exchange_Cells_Real(Flow % viscosity)
  if(Flow % heat_transfer) then
    call Grid % Exchange_Cells_Real(Flow % capacity)
    call Grid % Exchange_Cells_Real(Flow % conductivity)
  end if

  ! Surface tension is a property too, the same constant in every cell
  Vof % surf_tens(-Grid % n_bnd_cells:Grid % n_cells) = Vof % surface_tension

  end subroutine
