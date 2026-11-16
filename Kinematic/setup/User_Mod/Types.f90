!==============================================================================!
!   Kinematic advection tests (CAF review, comment 2.8): user variables        !
!==============================================================================!

  integer :: kin_mode = 1       ! 1 expansion, 2 translation, 3 deformation
  integer :: kin_init = 0       ! 0 linear estimate (production), 1 sub-sampled
  real    :: kin_c(3) = 0.0     ! sphere centre
  real    :: kin_rad  = 0.0     ! sphere radius
  integer :: kin_cyl  = 0       ! 1: a cylinder along z (two-dimensional disc)
  real    :: kin_p(4) = 0.0     ! mode parameters
  real    :: kin_vol0 = 0.0     ! initial volume of the sphere phase
  real    :: kin_w0   = 0.0     ! initial volume of cells with 0.01 < a < 0.99
  real, allocatable :: kin_fun0(:)  ! initial colour function
