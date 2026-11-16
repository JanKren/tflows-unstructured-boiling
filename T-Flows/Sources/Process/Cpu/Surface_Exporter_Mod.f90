#include "../../Shared/Assert.h90"
#include "../../Shared/Browse.h90"
#include "../../Shared/Unused.h90"
#include "../../Shared/Macros.h90"

!==============================================================================!
  module Surface_Exporter_Mod
!------------------------------------------------------------------------------!
!   Surface extraction and export functionality for T-Flows                   !
!------------------------------------------------------------------------------!
!   This module provides functionality to extract and export specific         !
!   surfaces from the computational domain based on cutting planes defined   !
!   by origin point and normal vector, similar to OpenFOAM's surface        !
!   extraction capability.                                                    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Field_Mod
  use File_Mod
  use Grid_Mod
  use Time_Mod
  use Vof_Mod
  use Vtk_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------------------!
  !                            !
  !   Surface Cut Type         !
  !                            !
  !----------------------------!
  type Surface_Cut_Type

    ! Plane definition
    real :: origin(3)                ! origin point coordinates (x0, y0, z0)
    real :: normal(3)                ! normal vector components (nx, ny, nz)
    real :: tolerance               ! tolerance for point-plane distance

    ! Extraction parameters
    integer :: extraction_method    ! 1=cell centers, 2=face centers, 3=both
    real    :: search_radius        ! radius around plane to search for points

    ! Output settings
    character(SL) :: name           ! name for output files
    logical       :: active         ! whether this cut is active

    ! Extracted surface data
    integer :: n_points    = 0      ! number of points on extracted surface
    integer :: n_triangles = 0      ! number of triangles on extracted surface

    ! Point coordinates
    real, allocatable :: points(:,:)     ! point coordinates (3, n_points)

    ! Triangle connectivity
    integer, allocatable :: triangles(:,:)  ! triangle node indices

    ! Field data on surface
    real, allocatable :: field_data(:,:)    ! field values at points
    integer           :: n_fields = 0       ! number of fields stored
    character(SL), allocatable :: field_names(:)  ! names of stored fields

  end type

  !------------------------------!
  !                              !
  !   Surface Exporter Type      !
  !                              !
  !------------------------------!
  type Surface_Exporter_Type

    integer :: n_cuts = 0                          ! number of defined cuts
    type(Surface_Cut_Type), allocatable :: cuts(:) ! array of surface cuts

    ! Save frequency control
    integer :: save_interval = 1                   ! surface save interval
    integer :: last_saved_step = -1                ! last saved step
    logical :: first_extraction = .true.           ! flag for verbose output

    contains
      procedure :: Add_Cut_Plane
      procedure :: Extract_All_Surfaces
      procedure :: Extract_Surface
      procedure :: Export_All_Surfaces
      procedure :: Export_Surface_Vtu
      procedure :: Set_Save_Interval
      procedure :: Should_Save_Now

      procedure, private :: Find_Face_Plane_Intersection
      procedure, private :: Interpolate_Cell_To_Point
      procedure, private :: Point_Plane_Distance

  end type

  ! Global instance
  type(Surface_Exporter_Type) :: Surface_Exporter

  contains

#   include "Surface_Exporter_Mod/Add_Cut_Plane.f90"
#   include "Surface_Exporter_Mod/Export_All_Surfaces.f90"
#   include "Surface_Exporter_Mod/Export_Surface_Vtu.f90"
#   include "Surface_Exporter_Mod/Extract_All_Surfaces.f90"
#   include "Surface_Exporter_Mod/Extract_Surface.f90"
#   include "Surface_Exporter_Mod/Find_Face_Plane_Intersection.f90"
#   include "Surface_Exporter_Mod/Interpolate_Cell_To_Point.f90"
#   include "Surface_Exporter_Mod/Point_Plane_Distance.f90"
#   include "Surface_Exporter_Mod/Set_Save_Interval.f90"
#   include "Surface_Exporter_Mod/Should_Save_Now.f90"

  end module