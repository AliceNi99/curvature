# Curvature Computation on 3D Meshes

This project provides a C++ implementation for computing curvature on 3D meshes using different weighting schemes. It supports input files in STL, OBJ, and OFF formats and outputs in OBJ and PLY formats. 
OBJ output files will contain per-vertex curvature information encoded as RGB `float`colors ranging in [0,1]. 
PLY per-vertex colors are saved as `uchar` RGB in [0,255]. Additionally, PLY files will contain per-vertex curvature information as a scalar property
## Authors:
- Alice Nicoletta
  
    <img src="https://github.com/user-attachments/assets/18c6f562-099b-4ed0-be91-f566f00690ce" alt="Curvature" width="200">

# Requirements 
    - g++ 15.2+ (MinGW-w64 for Windows)
    - CMake 3.15+
    - Ninja 1.13.1+
    - Eigen 3.3+ (either installed system-wide or included in the `external` folder). A submodule is provided already.

These requirements are quite conservative; older versions should also (maybe?) work fine.

# How to use
To run the code in this project, follow these steps:
1. **Clone the Repository**: 
    If you haven't already, clone the repository to your local machine using:
    ```
    git clone https://github.com/AliceNi99/curvature.git mean_curvature
    cd mean_curvature
    git submodule update --init --recursive
    ```
2. **Build the Project**:
   - Open a terminal and navigate to the root directory of the project.
   - For release build: 
   ```
   cmake --preset mingw-release ..
   cmake --build build/mingw-release
   ``` 
   - For debug build:
   ```
    cmake --preset mingw-debug ..
    cmake --build build/mingw-debug
   ```
    The name of the executable will be `curvature`, but feel free to change it from the CMakeLists.txt file. You can also modify the build presets in the `CMakePresets.json` file to suit your environment.
    I suggest making a bash/batch script to automate the build process (just in case you need to rebuild often).
3. **Run the Executable**: 
    After building, you can run the generated executable (located in the `build/mingw-release` (or `build/mingw-debug`) directory) with the appropriate command-line arguments:
   - `path/input_file`: Path and name of the input file, can support STL, OBJ, or OFF formats.
   - `path/output_file`: Path and name of the output file, can support OBJ and PLY format. <br>
   - `option`: optional; curvature can be calculated as `uniform`, `area` or `angle` (default is `area`). `area` or `angle` affect only the per-vertex normal, thus the curvature final sign and both will compute the curve through the cotangent formula.
    
    Example on Linux and Mac OS: 
    ```
    ./build/mingw-release/curvature ./mesh/bunny.obj ./mesh/bunny_curvature.ply
    ```
    Example on Windows:
    ```
    .\build\mingw-release\curvature.exe ./mesh/bunny.obj ./mesh/bunny_curvature.ply
    ```
4. **View the Output**: 
    The curvature values can be visualized as colors on the mesh, with different colors representing different curvature values.
    Red is for high negative curvature (concave), blue indicates high positive curvature (convex), and green indicates low or zero curvature.
    You can view the output files using any 3D model viewer that supports OBJ or PLY formats, such as MeshLab or Blender. Unfortunately 3D viewers like `Windows 3D Viewer` or `Paint 3D` do not support per-vertex colors in OBJ or PLY files.
    
# Limitations
- Only mean curvature is computed, not Gaussian curvature.
- The input mesh must be a 2-Manifold 3D triangular mesh.
- The program does not handle non-triangular meshes.
- The program assumes the input mesh is well-formed and does not contain degenerate faces.
- The curvature computation methods are based on discrete approximations and may not be accurate for all mesh geometries.
- The program takes a mesh and returns a mesh, it does not provide visualization capabilities.

# Future Work
- Implement additional curvature computation methods.
- Optimize performance for large meshes.
- Add support for more input/output file formats.
