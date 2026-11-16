!==============================================================================!
  subroutine Smooth_Scalar(Vof, Grid, var, smooth_var, n_conv)
!------------------------------------------------------------------------------!
!   Smoothes scalar using a Laplacian smoother.  The values of the boundary    !
!   cells, the contact angle at walls, are set in Smooth_Scalar_Boundary_      !
!   Values; it is only called for the colour function.                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  type(Grid_Type)         :: Grid
  integer                 :: n_conv
  real                    :: var(-Grid % n_bnd_cells: Grid % n_cells)
  real                    :: smooth_var(-Grid % n_bnd_cells    &
                                       : Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer                   :: nb, nc, reg
  integer                   :: s, c, c1, c2, c_iter
  real                      :: fs, vol_face
  real, contiguous, pointer :: sum_vol_area(:), sum_area(:)
!==============================================================================!

  call Work % Connect_Real_Cell(sum_vol_area, sum_area)

  ! Take aliases
  nb = Grid % n_bnd_cells
  nc = Grid % n_cells

  ! Copy the values from phi % n to local variable
  smooth_var(:) = var(:)

  ! Boundary cells start from the value of their cell: the walls among
  ! them keep the values of a pass for the gradient of the next
  do reg = Boundary_Regions()
    do s = Faces_In_Region(reg)
      if(s>0) then
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)
        smooth_var(c2) = smooth_var(c1)
      end if
    end do
  end do

  do c_iter = 1, n_conv

    sum_vol_area(-nb:nc) = 0.0
    sum_area    (-nb:nc) = 0.0

    !-------------------------------!
    !   Extrapolate to boundaries   !
    !-------------------------------!
    call Vof % Smooth_Scalar_Boundary_Values(Grid, smooth_var)

    ! At boundaries
    do reg = Boundary_Regions()
      do s = Faces_In_Region(reg)
        if(s>0) then
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)

          sum_vol_area(c1) = sum_vol_area(c1) + smooth_var(c2) * Grid % s(s)
          sum_area(c1) = sum_area(c1) + Grid % s(s)
        end if
      end do
    end do

    do s = Faces_In_Domain_And_At_Buffers()
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      fs = Grid % f(s)
      vol_face         = fs * smooth_var(c1) + (1.0-fs) * smooth_var(c2)
      sum_vol_area(c1) = sum_vol_area(c1) + vol_face * Grid % s(s)
      sum_vol_area(c2) = sum_vol_area(c2) + vol_face * Grid % s(s)
      sum_area(c1)     = sum_area(c1) + Grid % s(s)
      sum_area(c2)     = sum_area(c2) + Grid % s(s)
    end do

    call Grid % Exchange_Cells_Real(sum_vol_area(-nb:nc))
    call Grid % Exchange_Cells_Real(sum_area    (-nb:nc))

    do c = Cells_In_Domain()
      ! No clipping here, deliberately.  A volume-weighted face average of a
      ! bounded colour function is itself bounded, so a clip would only ever
      ! act when var already overshoots [0,1] - masking that instead of
      ! fixing it.  (An unreachable max(min(...,1.0),0.0) used to sit above
      ! this line and never took effect.)  Bounding smooth_var is a numerics
      ! change, not a cleanup: it feeds both surface-tension terms.
      smooth_var(c) = sum_vol_area(c) / sum_area(c)
    end do
    call Grid % Exchange_Cells_Real(smooth_var)

  end do

  ! At boundaries
  call Vof % Smooth_Scalar_Boundary_Values(Grid, smooth_var)
  call Grid % Exchange_Cells_Real(smooth_var(-nb:nc))

  call Work % Disconnect_Real_Cell(sum_vol_area, sum_area)

  end subroutine
