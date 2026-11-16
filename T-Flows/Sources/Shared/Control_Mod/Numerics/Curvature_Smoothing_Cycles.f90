!==============================================================================!
  subroutine Curvature_Smoothing_Cycles(Control, val, verbose)
!------------------------------------------------------------------------------!
!>  How many times the curvature field itself is smoothed, after it has been
!>  computed.  This is NOT the same as MAX_SMOOTHING_CYCLES_CURVATURE_VOF,
!>  which smooths the COLOUR field the curvature is derived from.
!>
!>  The two act in opposite directions and only one of them helps.  Smoothing
!>  the colour field more widens the interface band; the band then spans a
!>  range of radii over which the true curvature genuinely varies, so the
!>  SCATTER of kappa goes UP even as its bias goes down.  Smoothing kappa
!>  itself attacks the scatter directly.
!>
!>  That distinction matters because parasitic currents are driven by the
!>  scatter of kappa and are indifferent to a uniform bias in it - a constant
!>  error only rescales a Laplace jump that the pressure still balances.
!>
!>  Default 1 reproduces the historical behaviour exactly.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Control_Type)  :: Control  !! parent class
  integer, intent(out) :: val      !! number of curvature smoothing passes
  logical, optional    :: verbose  !! controls output verbosity
!==============================================================================!

  call Control % Read_Int_Item('CURVATURE_SMOOTHING_CYCLES', 1, val, verbose)

  end subroutine
