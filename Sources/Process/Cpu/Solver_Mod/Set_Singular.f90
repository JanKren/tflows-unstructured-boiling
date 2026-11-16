!==============================================================================!
  subroutine Set_Singular(Sol, phi)
!------------------------------------------------------------------------------!
!>  The Set_Singular subroutine in the Solver_Mod module plays an important
!>  role in preparing the solver for handling singular matrices. This
!>  functionality is especially relevant when dealing with matrices that might
!>  be deficient in rank or lack full invertibility, which is typically the
!>  case for the pressure correction solution.  Just like her sister "Run",
!>  the subroutine shows how encapsulation of decision making on wheather to
!>  use native or PETSc solvers works in practice in T-Flows.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Solver_Type) :: Sol  !! parent class
  type(Var_Type)     :: phi  !! unknown variable
!-----------------------------------[Locals]-----------------------------------!
  type(Matrix_Type), pointer :: A
  type(Grid_Type),   pointer :: Grid
  integer                    :: i, ni, reference_mode, reference_gid, n_reference
!==============================================================================!

  ! Take aliases
  A => phi % pnt_matrix
  Grid => A % pnt_grid

  ! Prepare solver for the singular matrix
  if(Sol % solvers == PETSC) then
    call C_Petsc_Mat_Set_Null_Space(phi % Pet % A)
  else
    ni = A % pnt_grid % n_cells - A % pnt_grid % Comm % n_buff_cells
    call Control % Read_Int_Item('NATIVE_PRESSURE_REFERENCE_MODE',0,reference_mode,.false.)
    if(reference_mode==0) then
      do i = 1, Grid % n_cells
        A % val(A % dia(i)) = A % val(A % dia(i)) * (1.0 + MICRO)
      end do
    else if(reference_mode==1) then
      ! For a connected Neumann pressure system with compatible RHS, adding
      ! a diagonal at one global cell fixes only the pressure gauge: summing
      ! the equations gives p(reference)=0. Perturbing every diagonal instead
      ! leaves artificial local divergence proportional to the correction.
      reference_gid=not(ishft(1,bit_size(0)-1))
      do i=1,ni
        reference_gid=min(reference_gid,Grid % Comm % cell_glo(i))
      end do
      call Global % Min_Int(reference_gid)
      n_reference=0
      do i=1,Grid % n_cells
        if(Grid % Comm % cell_glo(i)/=reference_gid) cycle
        if(A % val(A % dia(i))<=0.) call Message % Error(80, &
          'Native pressure reference requires a positive diagonal.', &
          file=__FILE__,line=__LINE__,one_proc=.true.)
        A % val(A % dia(i))=2.*A % val(A % dia(i))
        if(i<=ni) n_reference=n_reference+1
      end do
      call Global % Sum_Int(n_reference)
      if(n_reference/=1) call Message % Error(80, &
        'Native pressure reference requires exactly one global owned cell.', &
        file=__FILE__,line=__LINE__,one_proc=.true.)
    else
      call Message % Error(80,'Unsupported native pressure reference mode.', &
                           file=__FILE__,line=__LINE__,one_proc=.true.)
    end if
  end if

  end subroutine
