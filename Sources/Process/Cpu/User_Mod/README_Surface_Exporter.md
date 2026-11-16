# T-Flows Surface Exporter

## Overview

The Surface Exporter provides functionality similar to OpenFOAM's surface extraction capability, allowing users to extract and export specific surfaces from the computational domain based on cutting planes defined by an origin point and normal vector.

## Files

- `Surface_Exporter_Mod.f90` - Main surface exporter module
- `Surface_Exporter_Example.f90` - User-modifiable example showing how to define cutting planes
- `Save_Results.f90` - Integration point in User_Mod

## How to Use

### 1. Enable Surface Export

In `Save_Results.f90`, uncomment the line:
```fortran
call User_Mod_Surface_Exporter_Example(Flow, Turb, Vof, Swarm, domain)
```

### 2. Define Cutting Planes

Edit `Surface_Exporter_Example.f90` and uncomment/modify the cutting plane definitions:

```fortran
! Example: Vertical plane at x = 0.5
origin(1) = 0.5   ! x-coordinate where plane passes through
origin(2) = 0.0   ! y-coordinate (any value for YZ plane)
origin(3) = 0.0   ! z-coordinate (any value for YZ plane)
normal(1) = 1.0   ! normal pointing in x-direction
normal(2) = 0.0
normal(3) = 0.0
call Surface_Exporter % Add_Cut_Plane('x_plane_05', origin, normal, 1.0e-6)
```

### 3. Cutting Plane Parameters

- **Name**: Unique identifier for the cutting plane (used in output filenames)
- **Origin**: Point (x, y, z) that the plane passes through
- **Normal**: Vector (nx, ny, nz) perpendicular to the plane (automatically normalized)
- **Tolerance**: Distance tolerance for point-plane intersection (typically 1.0e-6)

### 4. Common Cutting Planes

**Coordinate Planes:**
- **X = constant (YZ plane)**: `normal = [1.0, 0.0, 0.0]`
- **Y = constant (XZ plane)**: `normal = [0.0, 1.0, 0.0]`
- **Z = constant (XY plane)**: `normal = [0.0, 0.0, 1.0]`

**Diagonal Planes:**
- **45° in XY**: `normal = [1.0, 1.0, 0.0]`
- **45° in XZ**: `normal = [1.0, 0.0, 1.0]`
- **45° in YZ**: `normal = [0.0, 1.0, 1.0]`

### 5. Output

The extracted surfaces are saved as VTU files with names like:
- `results-00001234-x_plane_05.vtu`
- `results-00001234-y_plane_025.vtu`

These files can be opened and visualized in ParaView or VisIt.

## Example Use Cases

1. **Flow visualization**: Extract velocity profiles on specific planes
2. **Temperature analysis**: Show temperature distribution on cutting planes
3. **Pressure gradients**: Visualize pressure fields on selected surfaces
4. **Interface tracking**: Extract VOF interfaces at specific locations

## Future Enhancements

- Field interpolation onto extracted surfaces
- Support for curved surfaces (spheres, cylinders)
- Integration with existing T-Flows variable output system
- Parallel efficiency improvements

## Notes

- Surface extraction runs every time results are saved (based on RESULTS_SAVE_INTERVAL)
- Only the first processor prints status messages to avoid cluttered output
- The cutting plane algorithm currently uses a simple cell-center based approach
- For best results, ensure cutting planes intersect the computational domain
