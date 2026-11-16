!==============================================================================!
!   Module-level data shared by the hooks of the CAF review instrumentation.   !
!                                                                              !
!   caf_mode  0 : the case as it is, with the per-step record                  !
!             1 : uniform-temperature check; the temperature starts at T_sat   !
!                 everywhere and an interfacial mass flux caf_flux [kg/m^2/s]  !
!                 is imposed on the pressure source, so that div u is not zero !
!                 while nothing but the convective operator can move T         !
!   Read once from caf_mode.ini (two numbers: mode, flux); absent means 0.     !
!------------------------------------------------------------------------------!
  integer :: caf_mode   = -1       ! -1: not read yet
  real    :: caf_flux   =  0.0     ! imposed mass flux in mode 1 [kg/m^2/s]
  real    :: caf_mdot_vof = 0.0    ! sum of m_dot the VOF source used [kg/s]
  real    :: caf_dev_max  = 0.0    ! running max |T - T_sat|, mode 1 [K]
