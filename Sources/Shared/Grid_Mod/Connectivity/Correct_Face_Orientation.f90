!==============================================================================!
  subroutine Correct_Face_Orientation(Grid)
!------------------------------------------------------------------------------!
!>  Reverses the node order of any face whose surface vector points the wrong
!>  way, so that afterwards every inside face has its normal pointing from c1
!>  to c2 and every boundary face has its normal pointing out of the domain.
!------------------------------------------------------------------------------!
!   Calculate_Face_Geometry insists on this and stops the run if it is          !
!   violated, and rightly so: the sign of the convective flux through a face    !
!   is read straight off that dot product, so a single face wound the other way !
!   feeds mass into the domain where it should take it out.                     !
!                                                                              !
!   Getting the winding right at the moment a face is created is possible but   !
!   fragile.  A sub-face cut from the corner of a parent face inherits its      !
!   parent's winding, which is fine, but a face created between two children of !
!   the same cell has no parent to inherit from, and which of its two           !
!   neighbours ends up as c1 is decided later still, when the cells are         !
!   renumbered and the smaller number is put first.  Repairing the orientation  !
!   once, after the geometry is known, is both shorter and harder to get wrong  !
!   than threading the correct winding through all of that.                     !
!                                                                              !
!   Reversing the stored node order is what actually flips the normal, since    !
!   Calculate_Face_Surfaces derives the surface vector from that order; the     !
!   surface vector is negated here as well so the two stay consistent without   !
!   the caller having to recompute it.                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type), intent(inout) :: Grid  !! computational grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c1, c2, k, i, tmp, n_flip
  real    :: dot
!==============================================================================!

  n_flip = 0

  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)

    if(c2 .gt. 0) then
      ! Inside face: the normal should follow the cell to cell connection
      dot =   Grid % sx(s) * (Grid % xc(c2) - Grid % xc(c1))     &
            + Grid % sy(s) * (Grid % yc(c2) - Grid % yc(c1))     &
            + Grid % sz(s) * (Grid % zc(c2) - Grid % zc(c1))
    else
      ! Boundary face: the normal should point away from the inside cell
      dot =   Grid % sx(s) * (Grid % xf(s) - Grid % xc(c1))     &
            + Grid % sy(s) * (Grid % yf(s) - Grid % yc(c1))     &
            + Grid % sz(s) * (Grid % zf(s) - Grid % zc(c1))
    end if

    if(dot .lt. 0.0) then

      k = Grid % faces_n_nodes(s)

      ! Reverse the node order, which is what turns the normal around
      do i = 1, k / 2
        tmp = Grid % faces_n(i, s)
        Grid % faces_n(i, s) = Grid % faces_n(k - i + 1, s)
        Grid % faces_n(k - i + 1, s) = tmp
      end do

      ! ... and keep the already computed surface vector in step with it
      Grid % sx(s) = -Grid % sx(s)
      Grid % sy(s) = -Grid % sy(s)
      Grid % sz(s) = -Grid % sz(s)

      n_flip = n_flip + 1
    end if

  end do

  if(First_Proc() .and. n_flip .gt. 0) then
    print '(a,i9,a)', ' # Reoriented ', n_flip, ' face(s)'
  end if

  end subroutine
