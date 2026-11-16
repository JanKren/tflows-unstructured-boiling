!==============================================================================!
  subroutine User_Mod_Beginning_Of_Compute_Vof(Vof, Sol)
!------------------------------------------------------------------------------!
!   Kinematic advection tests (CAF review, comment 2.8).                       !
!                                                                              !
!   Hand-sets the face volume fluxes just before the volume fraction is        !
!   advected.  With zero outer iterations (MIN = MAX_SIMPLE_ITERATIONS = 0)    !
!   momentum, pressure and energy are never solved, so nothing else touches    !
!   Flow % v_flux.                                                             !
!                                                                              !
!   Every flux field is discretely solenoidal BY CONSTRUCTION, which is what   !
!   makes hand-setting legitimate:                                             !
!                                                                              !
!   mode 1  radial expansion: the exact flux of a point source through the     !
!           face, Q / (4 pi) times the solid angle the face subtends at the    !
!           source.  The solid angles of a closed polyhedron which does not    !
!           contain the source sum to zero.  (Variants inner = 1, 2, 3 scale   !
!           this flux in the vapour and are deliberately NOT solenoidal in     !
!           the interface band; see where they are applied.)                   !
!   mode 2  uniform translation: U . S, and the face vectors of a closed       !
!           cell sum to zero.                                                  !
!   mode 3  three-dimensional deformation (LeVeque): the circulation of the    !
!           vector potential around the face edges, and each edge of a         !
!           closed cell is traversed once in each direction.                   !
!                                                                              !
!   The measured defect is written to kin-flux.dat at every step.              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: MAXN = 256   ! largest number of nodes in a face
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Field_Type), pointer :: Flow
  logical, save             :: first = .true.
  integer, save             :: mode, inner, init
  real,    save             :: cen(3), rad, par(4)
  real,    save             :: div_max      ! max |sum geo| / vol, regular cells
  real,    save             :: src_sum      ! sum of geo over the source cells
  integer, save             :: n_src        ! number of source cells
  real, allocatable, save   :: geo(:)       ! geometric flux factor of a face
  real, allocatable, save   :: rf(:)        ! distance of a face from centre
  real, allocatable         :: div(:), sgeo(:), bgeo(:)
  integer                   :: s, c, c1, c2, i, nn, fu, reg, n_open
  integer                   :: nod(MAXN)
  real                      :: xyz(3, MAXN)
  real                      :: t0, t1, tm, dt, tv, r0, r1, rm, q, fac, defect
  real                      :: beta, d_l, wgt
!==============================================================================!

  Unused(Sol)

  ! Take aliases
  Flow => Vof  % pnt_flow
  Grid => Flow % pnt_grid

  !-----------------------------------------------------------------------!
  !                                                                       !
  !   First call: read the parameters and form the geometric factors      !
  !                                                                       !
  !-----------------------------------------------------------------------!
  if(first) then
    first = .false.

    open(newunit=fu, file='kinematic.in', action='read', status='old')
    read(fu, *) mode
    read(fu, *) cen(1), cen(2), cen(3), rad
    read(fu, *) par(1), par(2), par(3), par(4)
    read(fu, *) inner
    read(fu, *) init
    close(fu)

    allocate(geo(Grid % n_faces));  geo(:) = 0.0
    allocate(rf (Grid % n_faces));  rf (:) = 0.0

    do s = 1, Grid % n_faces
      nn = abs(Grid % faces_n_nodes(s))
      if(nn < 3) cycle
      Assert(nn <= MAXN)
      do i = 1, nn
        xyz(1, i) = Grid % xn(Grid % faces_n(i, s))
        xyz(2, i) = Grid % yn(Grid % faces_n(i, s))
        xyz(3, i) = Grid % zn(Grid % faces_n(i, s))
      end do

      ! Nodes in rotational order, counter-clockwise about the face vector.
      ! (Taking the node loop as stored in the mesh instead, reversed where
      ! needed, was tried on the polyhedral 75 mesh and gives the same
      ! closure statistics.)
      call Kin_Sort_Face_Nodes(nn, xyz,  &
                               (/Grid % sx(s), Grid % sy(s), Grid % sz(s)/))

      rf(s) = norm2((/sum(xyz(1,1:nn)), sum(xyz(2,1:nn)),  &
                      sum(xyz(3,1:nn))/) / real(nn) - cen(1:3))

      if(mode .eq. 1) then
        geo(s) = Kin_Solid_Angle(nn, xyz, cen) / (4.0 * PI)
      else if(mode .eq. 2) then
        geo(s) = par(1) * Grid % sx(s)  &
               + par(2) * Grid % sy(s)  &
               + par(3) * Grid % sz(s)
      else if(mode .eq. 3) then
        geo(s) = Kin_Circulation(nn, xyz)
      end if
    end do

    !-----------------------------------------------!
    !   Discrete divergence of the imposed fluxes   !
    !-----------------------------------------------!
    allocate(div (-Grid % n_bnd_cells:Grid % n_cells));  div (:) = 0.0
    allocate(sgeo(-Grid % n_bnd_cells:Grid % n_cells));  sgeo(:) = TINY
    do reg = Boundary_Regions()
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        div (c1) = div (c1) + geo(s)
        sgeo(c1) = sgeo(c1) + abs(geo(s))
      end do
    end do
    do s = Faces_In_Domain_And_At_Buffers()
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      div (c1) = div (c1) + geo(s)
      div (c2) = div (c2) - geo(s)
      sgeo(c1) = sgeo(c1) + abs(geo(s))
      sgeo(c2) = sgeo(c2) + abs(geo(s))
    end do

    div_max = 0.0
    src_sum = 0.0
    n_src   = 0
    n_open  = 0

    ! The faces of a cell at a corner of the domain in a polyhedral dual do
    ! not close (the corner itself is cut off).  In a uniform flow such a
    ! cell is a source of volume, and of spurious phase with it.  The flux
    ! through the missing piece is put on the boundary faces of the cell, in
    ! proportion to their own flux, which closes it.
    if(mode .eq. 2) then
      allocate(bgeo(-Grid % n_bnd_cells:Grid % n_cells));  bgeo(:) = 0.0
      do reg = Boundary_Regions()
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          bgeo(c1) = bgeo(c1) + abs(geo(s))
        end do
      end do
      do reg = Boundary_Regions()
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          if(abs(div(c1)) > 1.0e-9 * sgeo(c1) .and. bgeo(c1) > 0.0) then
            geo(s) = geo(s) - div(c1) * abs(geo(s)) / bgeo(c1)
          end if
        end do
      end do
      do c = Cells_In_Domain()
        if(abs(div(c)) > 1.0e-9 * sgeo(c)) then
          n_open = n_open + 1
          print '(a,i3,a,3es12.4,a,es12.4)',                              &
            ' # KIN: closed the cell on rank ', This_Proc(), ' at ',      &
            Grid % xc(c), Grid % yc(c), Grid % zc(c),                     &
            '  |sum geo| / sum |geo| was ', abs(div(c)) / sgeo(c)
          if(bgeo(c) > 0.0) div(c) = 0.0
        end if
      end do
      deallocate(bgeo)
    end if

    do c = Cells_In_Domain()
      ! A cell holding the point source (or a share of it) carries an
      ! O(1) fraction of the unit source.  The cells at the corners of the
      ! domain in a polyhedral dual are not closed by their faces and carry
      ! about 1e-6 of it.  Both are listed; everything else is round-off
      if(mode .eq. 1 .and. abs(div(c)) > 1.0e-9) then
        n_src   = n_src + 1
        src_sum = src_sum + div(c)
        print '(a,i3,a,3es12.4,a,es22.14)',                            &
          ' # KIN: source cell on rank ', This_Proc(), ' at ',         &
          Grid % xc(c), Grid % yc(c), Grid % zc(c), '  share ', div(c)
      else
        div_max = max(div_max, abs(div(c)) / Grid % vol(c))
      end if
    end do
    call Global % Sum_Int(n_open)
    if(First_Proc() .and. mode .eq. 2)  &
      print '(a,i6)', ' # KIN: corner cells closed: ', n_open
    call Global % Max_Real(div_max)
    call Global % Sum_Real(src_sum)
    call Global % Sum_Int (n_src)
    deallocate(div, sgeo)

    if(First_Proc()) then
      print '(a,i2,a,i2,a,i2)',  ' # KIN: mode = ', mode,  &
                                 '  inner = ', inner, '  init = ', init
      print '(a,es12.4,a,i4,a,es22.14)',                              &
        ' # KIN: max |sum geo| / vol in regular cells = ', div_max,   &
        '   source cells: ', n_src, '   their share of the source: ', src_sum
      open(newunit=fu, file='kin-flux.dat', action='write', status='replace')
      write(fu,'(a)') '# step  time  scale(Q or prefactor)  R_exact  '  //  &
                      'V_exact  defect = scale*dt*max|div|/vol  n_src  src_sum'
      close(fu)
    end if
  end if

  !------------------------------------------!
  !                                          !
  !   Every call: impose the face fluxes     !
  !                                          !
  !------------------------------------------!
  dt = Flow % dt
  t1 = Time % Get_Time()   ! time has been advanced before the VOF step
  t0 = t1 - dt
  tm = 0.5 * (t0 + t1)
  r1 = 0.0

  ! Radial expansion: the volume flow which makes the exact sphere
  ! follow the Scriven radius R_S(t) = 2 beta sqrt(D_l (t + t_v)) exactly
  if(mode .eq. 1) then
    beta = par(1)
    d_l  = par(2)
    tv   = rad**2 / (4.0 * beta**2 * d_l)
    r0   = 2.0 * beta * sqrt(d_l * (t0 + tv))
    r1   = 2.0 * beta * sqrt(d_l * (t1 + tv))
    rm   = 2.0 * beta * sqrt(d_l * (tm + tv))
    q    = 4.0 * PI / 3.0 * (r1**3 - r0**3) / dt
    fac  = q
    do s = 1, Grid % n_faces
      Flow % v_flux % n(s) = q * geo(s)

      ! inner = 0: the liquid velocity everywhere, which is solenoidal
      !            (the point source sits in pure vapour)
      ! inner = 1: no flux inside the exact sphere R_S(t)
      ! inner = 2: the one-fluid (volume-averaged) velocity with the vapour
      !            at rest, u = vof u_l, the face value of vof being the
      !            mean of the two cells.  Its divergence, u_l . grad(vof),
      !            lies in the interface band, as in the coupled solver
      ! inner = 3: as 2, with a step at vof = 1/2 in place of the mean
      if(inner .eq. 1) then
        if(rf(s) < rm) Flow % v_flux % n(s) = 0.0
      else if(inner .eq. 2 .or. inner .eq. 3) then
        c1  = Grid % faces_c(1,s)
        c2  = Grid % faces_c(2,s)
        wgt = 0.5 * (  min(max(Vof % fun % n(c1), 0.0), 1.0)  &
                     + min(max(Vof % fun % n(c2), 0.0), 1.0))
        if(inner .eq. 3) wgt = merge(1.0, 0.0, wgt > 0.5)
        Flow % v_flux % n(s) = wgt * Flow % v_flux % n(s)
      end if
    end do

  ! Uniform translation, optionally reversed at t = par(4)
  else if(mode .eq. 2) then
    fac = 1.0
    if(par(4) > 0.0 .and. tm > par(4)) fac = -1.0
    do s = 1, Grid % n_faces
      Flow % v_flux % n(s) = fac * geo(s)
    end do

  ! Deformation, reversing with the cosine of period par(1)
  else if(mode .eq. 3) then
    fac = cos(2.0 * PI * tm / par(1))
    do s = 1, Grid % n_faces
      Flow % v_flux % n(s) = fac * geo(s)
    end do
  end if

  defect = abs(fac) * dt * div_max

  if(First_Proc()) then
    call File % Append_For_Writing_Ascii('kin-flux.dat', fu)
    write(fu,'(i8,5(1x,es18.10),i6,1x,es22.14)')                     &
      Time % Curr_Dt(), t1, fac, r1, 4.0*PI/3.0*r1**3, defect, n_src, src_sum
    close(fu)
  end if

  end subroutine

!==============================================================================!
  subroutine Kin_Sort_Face_Nodes(nn, xyz, sv)
!------------------------------------------------------------------------------!
!   Puts the nodes of a face in rotational order, counter-clockwise when       !
!   seen against the face vector sv, by their angle about the node average.    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(in)    :: nn
  real,    intent(inout) :: xyz(3, nn)
  real,    intent(in)    :: sv(3)
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j
  real    :: cen(3), nv(3), e1(3), e2(3), d(3), ang(nn), a, tmp(3)
!==============================================================================!

  do i = 1, 3
    cen(i) = sum(xyz(i, 1:nn)) / real(nn)
  end do
  nv(:) = sv(:) / norm2(sv)

  ! In-plane basis
  e1(:) = xyz(:, 1) - cen(:)
  e1(:) = e1(:) - dot_product(e1, nv) * nv(:)
  e1(:) = e1(:) / norm2(e1)
  e2(1) = nv(2) * e1(3) - nv(3) * e1(2)
  e2(2) = nv(3) * e1(1) - nv(1) * e1(3)
  e2(3) = nv(1) * e1(2) - nv(2) * e1(1)

  do i = 1, nn
    d(:)   = xyz(:, i) - cen(:)
    ang(i) = atan2(dot_product(d, e2), dot_product(d, e1))
    if(ang(i) < -1.0e-12) ang(i) = ang(i) + 2.0 * PI
  end do

  ! Insertion sort on the angle
  do i = 2, nn
    a      = ang(i)
    tmp(:) = xyz(:, i)
    j = i - 1
    do while(j >= 1)
      if(ang(j) <= a) exit
      ang(j+1)    = ang(j)
      xyz(:, j+1) = xyz(:, j)
      j = j - 1
    end do
    ang(j+1)    = a
    xyz(:, j+1) = tmp(:)
  end do

  end subroutine

!==============================================================================!
  real function Kin_Solid_Angle(nn, xyz, p)
!------------------------------------------------------------------------------!
!   Signed solid angle which a face subtends at point p, positive when the     !
!   face vector points away from p.  The face is fanned into triangles about   !
!   the average of its nodes, and each triangle is evaluated with the formula  !
!   of Van Oosterom and Strackee (1983).                                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(in) :: nn
  real,    intent(in) :: xyz(3, nn)   ! nodes, counter-clockwise about S
  real,    intent(in) :: p(3)
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j
  real    :: c(3), a(3), b(3), la, lb, lc, num, den, len
!==============================================================================!

  do i = 1, 3
    c(i) = sum(xyz(i, 1:nn)) / real(nn) - p(i)
  end do
  lc = norm2(c)

  ! Size of the face, to recognise a source lying on it
  len = 0.0
  do i = 1, nn
    len = max(len, norm2(xyz(:, i) - p(:) - c(:)))
  end do

  Kin_Solid_Angle = 0.0
  do i = 1, nn
    j = mod(i, nn) + 1
    a(:) = xyz(:, i) - p(:)
    b(:) = xyz(:, j) - p(:)
    la = norm2(a)
    lb = norm2(b)

    num = c(1) * (a(2)*b(3) - a(3)*b(2))  &
        + c(2) * (a(3)*b(1) - a(1)*b(3))  &
        + c(3) * (a(1)*b(2) - a(2)*b(1))
    den = lc*la*lb + dot_product(c, a) * lb  &
                   + dot_product(c, b) * la  &
                   + dot_product(a, b) * lc

    ! The source lies in the plane of the triangle (on a node, on an edge
    ! or in the face): the flux through it vanishes by symmetry
    if(min(la, lb, lc) < 1.0e-9 * len) cycle
    if(abs(num) < 1.0e-12 * lc * la * lb) cycle

    Kin_Solid_Angle = Kin_Solid_Angle + 2.0 * atan2(num, den)
  end do

  end function

!==============================================================================!
  real function Kin_Circulation(nn, xyz)
!------------------------------------------------------------------------------!
!   Circulation around a face of the vector potential of the deformation       !
!   field of LeVeque (1996),                                                   !
!                                                                              !
!     u =  2 sin^2(pi x) sin(2 pi y) sin(2 pi z)                               !
!     v = -  sin(2 pi x) sin^2(pi y) sin(2 pi z)                               !
!     w = -  sin(2 pi x) sin(2 pi y) sin^2(pi z)                               !
!                                                                              !
!   which is the curl of                                                       !
!                                                                              !
!     A = (0, -sin^2(pi x) sin(2 pi y) sin^2(pi z),                            !
!              sin^2(pi x) sin^2(pi y) sin(2 pi z)) / pi                       !
!                                                                              !
!   Each edge is integrated from its lexicographically smaller end, with a     !
!   symmetric Gauss rule, so both faces sharing it get the same number.        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer, intent(in) :: nn
  real,    intent(in) :: xyz(3, nn)   ! nodes, counter-clockwise about S
!-----------------------------------[Locals]-----------------------------------!
  integer, parameter :: NG = 5
  real,    parameter :: XG(NG) = (/ -0.9061798459386640, -0.5384693101056831, &
                                     0.0,                                     &
                                     0.5384693101056831,  0.9061798459386640 /)
  real,    parameter :: WG(NG) = (/  0.2369268850561891,  0.4786286704993665, &
                                     0.5688888888888889,                      &
                                     0.4786286704993665,  0.2369268850561891 /)
  integer :: i, j, g, k
  real    :: a(3), b(3), x(3), d(3), pot(3), sgn, edge
!==============================================================================!

  Kin_Circulation = 0.0

  do i = 1, nn
    j = mod(i, nn) + 1
    a(:) = xyz(:, i)
    b(:) = xyz(:, j)

    ! Canonical direction
    sgn = 1.0
    do k = 1, 3
      if(a(k) < b(k)) exit
      if(a(k) > b(k)) then
        a(:) = xyz(:, j)
        b(:) = xyz(:, i)
        sgn  = -1.0
        exit
      end if
    end do

    d(:) = b(:) - a(:)
    edge = 0.0
    do g = 1, NG
      x(:)   = 0.5 * (a(:) + b(:)) + 0.5 * XG(g) * d(:)
      pot(1) = 0.0
      pot(2) = - sin(PI*x(1))**2 * sin(2.0*PI*x(2)) * sin(PI*x(3))**2 / PI
      pot(3) =   sin(PI*x(1))**2 * sin(PI*x(2))**2 * sin(2.0*PI*x(3)) / PI
      edge   = edge + 0.5 * WG(g) * dot_product(pot, d)
    end do

    Kin_Circulation = Kin_Circulation + sgn * edge
  end do

  end function
