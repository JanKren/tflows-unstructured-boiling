#include "./Vof_Initialization_Ellipsoid.f90"

!==============================================================================!
  subroutine User_Mod_Initialize_Variables(Flow, Turb, Vof, Swarm, Sol)
!------------------------------------------------------------------------------!
!   Seed the three-dimensional static droplet analytically.                     !
!                                                                              !
!   The quasi-two-dimensional study seeded from an STL, which does not scale    !
!   here: the cost of the distance calculation is O(cells x facets), and a      !
!   sphere resolved finely enough for the 128 level needs thousands of facets.  !
!   Measured, that initialisation reached 1-2% in forty minutes.  The analytic  !
!   ellipsoid is exact, independent of any faceting, and immediate.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: Swarm
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: fun
  integer                  :: c
!==============================================================================!

  Grid => Flow % pnt_grid
  fun  => Vof % fun

  do c = Cells_In_Domain_And_Buffers()
    fun % n(c) = 0.0
  end do

  call Vof_Initialization_Ellipsoid(Vof)
  call Grid % Exchange_Cells_Real(fun % n)

  fun % o (:) = fun % n(:)
  fun % oo(:) = fun % n(:)

  end subroutine
