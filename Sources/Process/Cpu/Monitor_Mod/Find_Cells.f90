!==============================================================================!
  subroutine Find_Cells(Monitor)
!------------------------------------------------------------------------------!
!>  Works out which cell each monitoring point falls in.
!------------------------------------------------------------------------------!
!   A monitoring point is a place, but what gets stored is a cell number, and a !
!   cell number only means anything for one particular grid.  Refinement        !
!   renumbers the cells underneath the monitor, so the number it is holding     !
!   afterwards belongs to whatever cell has since taken that slot.               !
!                                                                              !
!   That is a quiet failure of the worst kind: nothing crashes and the file     !
!   keeps filling, but the numbers in it are from somewhere else.  On the       !
!   hydrostatic cavity the reported pressure jumped from -0.275 to 0.375 at the !
!   refinement step, which is the correct answer for a point two thirds of the  !
!   way down the box and no answer at all for the one that was asked about.     !
!                                                                              !
!   This is split out of Initialize so it can be called again on its own.       !
!   Initialize cannot: it allocates, opens files and writes headers, and none   !
!   of that wants doing twice.                                                  !
!                                                                              !
!   The reduction at the end is what settles ownership in parallel.  Every rank !
!   finds its own nearest cell; only the one that is nearest overall keeps the  !
!   point, and the others set their cell to zero so they stay quiet.  It has to !
!   be redone after refinement too, since a rank that gains finer cells there   !
!   may now hold a cell closer to the point than the rank that had it before.   !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Monitor_Type) :: Monitor  !! parent class of Monitor_Type
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: c, m
  real                     :: curr_dist, min_dist_all, min_dist
!==============================================================================!

  if(Monitor % n_points .eq. 0) return

  Grid => Monitor % pnt_grid

  do m = 1, Monitor % n_points

    min_dist = HUGE
    Monitor % cell(m) = 0

    do c = Cells_In_Domain_And_Buffers()
      curr_dist = Math % Distance( Monitor % x(m),  &
                                   Monitor % y(m),  &
                                   Monitor % z(m),  &
                                   Grid % xc(c),    &
                                   Grid % yc(c),    &
                                   Grid % zc(c))
      if(curr_dist < min_dist) then
        Monitor % cell(m) = c
        min_dist          = curr_dist
      end if
    end do

    ! Somebody else may be closer; if so, this rank lets the point go
    if(Parallel_Run()) then
      min_dist_all = min_dist
      call Global % Min_Real(min_dist_all)
      if(abs(min_dist_all - min_dist) >= TINY) Monitor % cell(m) = 0
    end if

  end do

  end subroutine
