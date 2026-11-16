!==============================================================================!
  subroutine Calculate_Bnd_Cell_Centers(Grid)
!------------------------------------------------------------------------------!
!>  Places every boundary cell centre at the foot of the perpendicular dropped
!>  from its inside neighbour's centre onto the plane of the boundary face.
!------------------------------------------------------------------------------!
!   This reproduces what Convert does when it builds a grid from scratch, and  !
!   it is needed after refinement because Calculate_Cell_Centers only covers   !
!   cells 1 to n_cells - boundary cells carry negative indices and are left    !
!   untouched by it.  A freshly created boundary cell would otherwise sit at   !
!   the origin, which quietly wrecks every boundary condition that reads a     !
!   distance or a direction.                                                   !
!                                                                              !
!   The foot of the perpendicular, rather than the face centre, is the         !
!   convention the rest of the code assumes: it makes the cell-to-cell vector  !
!   across a boundary face normal to that face, so the wall-normal gradient    !
!   is a plain difference over a plain distance and needs no skewness          !
!   correction.  Convert offers the face-centre variant behind an option, but  !
!   the perpendicular foot is the default and is what is used here.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), intent(inout) :: Grid  !! computational grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c1, c2
  real    :: tot_surf, t
!==============================================================================!

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)

    if(c2 .lt. 0) then

      tot_surf = sqrt(  Grid % sx(s) * Grid % sx(s)     &
                      + Grid % sy(s) * Grid % sy(s)     &
                      + Grid % sz(s) * Grid % sz(s) )

      ! A zero area face would mean the subdivision produced a degenerate
      ! sub-face, which is worth stopping for rather than dividing by
      if(tot_surf .lt. TINY) then
        call Message % Error(72,                                            &
           'A boundary face with zero surface area was found while '//      &
           'placing boundary cell centres.  This normally means a '//       &
           'refinement sweep produced a degenerate sub-face.  Exiting!',    &
           file=__FILE__, line=__LINE__)
      end if

      t = (  Grid % sx(s) * (Grid % xf(s) - Grid % xc(c1))     &
           + Grid % sy(s) * (Grid % yf(s) - Grid % yc(c1))     &
           + Grid % sz(s) * (Grid % zf(s) - Grid % zc(c1)) ) / tot_surf

      Grid % xc(c2) = Grid % xc(c1) + Grid % sx(s) * t / tot_surf
      Grid % yc(c2) = Grid % yc(c1) + Grid % sy(s) * t / tot_surf
      Grid % zc(c2) = Grid % zc(c1) + Grid % sz(s) * t / tot_surf

    end if
  end do

  end subroutine
