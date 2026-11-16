#include "../../Shared/Assert.h90"
#include "../../Shared/Browse.h90"
#include "../../Shared/Unused.h90"
#include "../../Shared/Macros.h90"

!==============================================================================!
  module Vof_Mod
!------------------------------------------------------------------------------!
!   Definition of variables used for all VOF multiphase modelling              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Stl_Mod
  use Surf_Mod  ! inherited from Front_Mod
  use Turb_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------!
  !               !
  !   Vof class   !
  !               !
  !---------------!
  type Vof_Type

    ! Stores the name of the STL file for initialization
    character(SL) :: name_stl = ''
    logical       :: init_stl = .false.   ! is it initialized from STL?
    logical       :: oncevof = .true.     ! is it initialized from STL?

    type(Grid_Type),   pointer :: pnt_grid    ! its grid
    type(Field_Type),  pointer :: pnt_flow    ! its flow
    type(Matrix_Type), pointer :: pnt_matrix  ! its matrix
    type(Front_Type)           :: Front       ! pointer to front
    type(Surf_Type)            :: Surf        ! pointer to surface

    ! Volume fraction (colour function) and its smooth variant
    type(Var_Type) :: fun
    type(Var_Type) :: smooth

    ! Additional variables for calculation of vof function (fun)
    real, allocatable :: beta_f(:)
    real, allocatable :: c_d(:)

    ! Surface curvature
    real, allocatable :: curv(:)

    ! Surface normals
    real, allocatable :: nx(:), ny(:), nz(:)

    ! Surface tension force
    real, allocatable :: surf_fx(:), surf_fy(:), surf_fz(:)  ! [N/m^3]

    real, allocatable :: fract_tmp(:)

    ! Physical properties in case of multiphase flow
    real :: phase_visc(0:1), phase_dens(0:1)
    real :: phase_capa(0:1), phase_cond(0:1)
    real :: surface_tension

    ! Cell-centred surface tension [N/m], the value the solver uses.  It
    ! holds the constant surface_tension in every cell (filled in
    ! Update_Physical_Properties).
    real, allocatable :: surf_tens(:)

    ! Phase change (called mass transfer to be consistent
    ! with heat transfer in the rest of the code)
    real    :: m_d, m_ini, m_s, m_s_acc

    ! For phase change
    real :: t_sat, latent_heat  ! [K, J/kg]

    ! Neighbouring coefficients due to presence of saturated front
    real, allocatable :: a12(:)
    real, allocatable :: a21(:)
    real, allocatable :: m_dot(:)       ! [kg/s]
    real, allocatable :: q_0(:)         ! [kg/s]
    real, allocatable :: q_1(:)         ! [kg/s]
    real, allocatable :: m_dot0(:)      ! [kg/s]

    type(Var_Type) :: t_0, t_1, u_0, v_0, w_0
    type(Face_Type) :: v_flux1  ! extrapolated flux for mass transfer
    type(Face_Type) :: u_f, v_f, w_f

    ! User define parameters for vof function (fun)
    real    :: courant_max_param
    integer :: n_sub_param, corr_num_max
    integer :: n_smooth_for_curv_csf

    ! Passes of Smooth_Curvature applied to kappa itself.  Distinct from
    ! n_smooth_for_curv_csf, which smooths the colour field; that one
    ! lowers the bias of kappa but RAISES its scatter, and scatter is
    ! what drives parasitic currents.  Default 1 = historical behaviour.
    integer :: n_smooth_curv = 1

    ! Gradients next to symmetry planes, where a cell has its own mirror
    ! image for a neighbour (see Mirror_Setup).  Prepared on first use.
    logical              :: mirror_ready  = .false.
    logical              :: mirror_active = .false.  ! plane on any rank
    integer              :: n_mirror      = 0        ! cells with a plane
    integer, allocatable :: mirror_cell(:)           ! these cells
    integer, allocatable :: mirror_idx (:)           ! place in the list, or 0
    real,    allocatable :: mirror_inv (:,:)         ! their inverse matrices
    real,    allocatable :: mirror_acc (:,:)         ! their gradients

    ! Averaging
    integer, allocatable :: avg_cells(:,:)

    ! Mesh the front or surface
    logical :: track_front    ! simple but ugly
    logical :: track_surface  ! complicated but beautiful

    ! Extended (node-connected) stencil for the front-modified gradient.
    ! With ext_stencil off and ext_weight_p zero, the front-modified
    ! gradient is bit-identical to the original face-only formulation.
    ! See Vof_Mod/Utilities/Find_Extended_Stencil.f90
    logical :: ext_stencil   = .false.  ! add edge/corner neighbours

    ! Diagnostic: replace the computed interfacial mass flux by a spatially
    ! UNIFORM one carrying the same total, so that every point of the
    ! interface receives identical mass per unit area.  Any deformation that
    ! survives is then attributable to the transport of alpha alone, not to
    ! angular variation of the mass transfer rate.  Default off.
    logical :: uniform_m_dot = .false.

    ! Diagnostic: override the CICSAM/STACS blend weight gamma_f, which is the
    ! only explicitly orientation-dependent term in the limiter --- cos^n of the
    ! angle between grad(alpha) and the cell connection, and on a Cartesian mesh
    ! there are only six connection directions.  Setting it to 1 gives pure
    ! Hyper-C (maximally compressive, and bounded by construction since
    ! alfa_cbc = min(1, alfa_d_til/Co) IS the CBC limit); 0 gives pure
    ! ULTIMATE-QUICKEST for CICSAM, or pure STOIC for STACS.  The two bracket
    ! the blend.  Note this does not remove every orientation dependence: the
    ! donor Courant number entering alfa_cbc still varies with face direction.
    ! Negative (the default) leaves the computed blend untouched.
    real    :: blend_gamma   = -1.0

    ! Diagnostic: impose a spatially ISOTROPIC compression factor.  CICSAM and
    ! STACS blend a compressive and a bounded face value with a weight
    ! gamma_f = cos^n(angle between grad(alpha) and the cell connection), and
    ! the donor Courant number entering alfa_cbc varies with face orientation
    ! too.  On a Cartesian mesh there are only six connection directions, so
    ! both make beta_f modulate with cubic symmetry however exact the normal
    ! is.  Overwriting beta_f by one constant on every interface face removes
    ! that modulation while keeping the compression strength, so the interface
    ! stays sharp and the bubble grows at the same rate -- unlike an equal
    ! density ablation, which removes the growth along with the advection.
    ! Negative (the default) leaves the computed value untouched.
    !
    ! WARNING.  A constant bypasses the scheme's own CBC limiting, so alpha can
    ! leave [0,1]; Correct_Beta() then repairs beta_f face by face and quietly
    ! reinstates the orientation dependence the constant was meant to remove.
    ! Imposing the mean of a distribution pointwise is NOT reproducing its mean
    ! behaviour: CICSAM varies beta_f precisely to satisfy Courant boundedness,
    ! so faces needing beta_f ~ 0 get the mean instead and overshoot.  At
    ! beta_iso = 0.3096, the measured <beta_f> on this problem, the corrector
    ! exhausted on 99.7% of steps against 5% for the scheme's own beta.  Read
    ! the '# VOF_CORR' line Compute_Vof prints before trusting any such run.
    real    :: beta_iso      = -1.0

    real    :: ext_weight_p  = 0.0      ! LSQ weight w = |r|^-p
    integer, allocatable :: ext_n(:)    ! number of extra neighbours per cell
    integer, allocatable :: ext_c(:,:)  ! the extra neighbours themselves

    contains

      !----------------------------------------!
      !   Procedures to advance vof function   !
      !----------------------------------------!
      procedure          :: Check_Controls
      procedure          :: Create_Vof
      procedure          :: Main_Vof
      procedure, private :: Compute_Vof
      procedure, private :: Discretize
      procedure, private :: Correct_Beta
      procedure          :: Initialize_From_Stl
      procedure, private :: Max_Courant_Number
      procedure, private :: Predict_Beta
      procedure, private :: Solve_System

      !------------------------------------------!
      !   Procedures helping to prepare smooth   !
      !   (convoluted) variant of vof function   !
      !   for eventual estimation of curvature   !
      !------------------------------------------!
      procedure, private :: Curvature_Csf
      procedure, private :: Smooth_Curvature
      procedure, private :: Smooth_Scalar
      procedure, private :: Smooth_Scalar_Boundary_Values
      procedure, private :: Mirror_Setup
      procedure, private :: Mirror_Lsq
      procedure, private :: Mirror_Gradient
      procedure, private :: Mirror_Divergence
      procedure          :: Smooth_Vof_And_Compute_Surface_Normals

      !----------------------------------------------!
      !   Procedures to be called by other modules   !
      !----------------------------------------------!
      procedure :: Calculate_Grad_Matrix_With_Front
      procedure :: Extrapolate_Normal_To_Front
      procedure :: Find_Extended_Stencil
      procedure :: Get_Vapour_And_Liquid_Phase
      procedure :: Cap_Dens
      procedure :: Grad_Component_No_Refresh_With_Front
      procedure :: Grad_Variable_With_Front
      procedure :: Mass_Transfer_Added_Volume
      procedure :: Mass_Transfer_Estimate
      procedure :: Mass_Transfer_Pressure_Source
      procedure :: Mass_Transfer_Vof_Source
      procedure :: Surface_Tension_Force
      procedure :: Update_Physical_Properties

  end type

  ! type(Vof_Type) :: Vof

  contains

    !----------------------------------------!
    !   Procedures to advance vof function   !
    !----------------------------------------!
#   include "Vof_Mod/Core/Check_Controls.f90"
#   include "Vof_Mod/Core/Create_Vof.f90"
#   include "Vof_Mod/Core/Main_Vof.f90"
#   include "Vof_Mod/Core/Compute_Vof.f90"
#   include "Vof_Mod/Core/Discretize.f90"
#   include "Vof_Mod/Core/Correct_Beta.f90"
#   include "Vof_Mod/Core/Initialize_From_Stl.f90"
#   include "Vof_Mod/Core/Max_Courant_Number.f90"
#   include "Vof_Mod/Core/Predict_Beta.f90"
#   include "Vof_Mod/Core/Solve_System.f90"

    !------------------------------------------!
    !   Procedures helping to prepare smooth   !
    !   (convoluted) variant of vof function   !
    !   for eventual estimation of curvature   !
    !------------------------------------------!
#   include "Vof_Mod/Curvature/Curvature_Csf.f90"
#   include "Vof_Mod/Curvature/Smooth_Curvature.f90"
#   include "Vof_Mod/Curvature/Smooth_Vof_And_Compute_Surface_Normals.f90"
#   include "Vof_Mod/Curvature/Smooth_Scalar.f90"
#   include "Vof_Mod/Curvature/Smooth_Scalar_Boundary_Values.f90"
#   include "Vof_Mod/Curvature/Mirror_Setup.f90"
#   include "Vof_Mod/Curvature/Mirror_Lsq.f90"
#   include "Vof_Mod/Curvature/Mirror_Gradient.f90"
#   include "Vof_Mod/Curvature/Mirror_Divergence.f90"

    !----------------------------------------------!
    !   Procedures to be called by other modules   !
    !----------------------------------------------!
#   include "Vof_Mod/Utilities/Calculate_Grad_Matrix_With_Front.f90"
#   include "Vof_Mod/Utilities/Extrapolate_Normal_To_Front.f90"
#   include "Vof_Mod/Utilities/Find_Extended_Stencil.f90"
#   include "Vof_Mod/Utilities/Get_Gas_And_Liquid_Phase.f90"
#   include "Vof_Mod/Utilities/Cap_Dens.f90"
#   include "Vof_Mod/Utilities/Grad_Component_No_Refresh_With_Front.f90"
#   include "Vof_Mod/Utilities/Grad_Variable_With_Front.f90"
#   include "Vof_Mod/Utilities/Mass_Transfer_Added_Volume.f90"
#   include "Vof_Mod/Utilities/Mass_Transfer_Estimate.f90"
#   include "Vof_Mod/Utilities/Mass_Transfer_Pressure_Source.f90"
#   include "Vof_Mod/Utilities/Mass_Transfer_Vof_Source.f90"
#   include "Vof_Mod/Utilities/Surface_Tension_Force.f90"
#   include "Vof_Mod/Utilities/Update_Physical_Properties.f90"

    !------------------------------------!
    !   User functions for this module   !
    !------------------------------------!
#   include "User_Mod/Beginning_Of_Compute_Vof.f90"
#   include "User_Mod/End_Of_Compute_Vof.f90"

  end module
