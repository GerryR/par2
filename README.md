# PAR²

Parallel Random Walk Particle Tracking Method for Solute Transport in Porous Media

## Build

The following softwares and libraries must be installed:

* cmake (version 3.8 or higher)
* CUDA Toolkit (tested with version 9.0)
* yaml-cpp library
* spdlog library (included in the source code)

### Compile on Linux

1. Make sure to have a valid c++ compiler (e.g., gcc)
2. Create a build directory:

        mkdir Build  
        cd Build

2. Create makefile:

        cmake -DCMAKE_BUILD_TYPE=Release ..

3. Compile:

        make

### Compile on Windows

1. Make sure to have Visual Studio 2017 installed
2. Create a build directory:

        mkdir Build  
        cd Build

2. Create MSVC solution:

        cmake -G "Visual Studio 15 2017 Win64" -T v140 -DCMAKE_BUILD_TYPE=Release -DYAML_ROOT=C:/path/to/yaml-cpp ..

3. Compile using the Developer Command Prompt for VS:

        devenv /Build Release par2.sln

## Citations

Rizzo, C. B., Nakano, A., and de Barros, F. P. J. [PAR²: Parallel Random Walk Particle Tracking Method for Solute Transport in Porous Media](https://doi.org/10.1016/j.cpc.2019.01.013) Computer Physics Communications

## Directory Structure
```
.
├── CMakeLists.txt
├── maingpu.cu (Entry point for PAR2)
├── README.md
|
├── CMakeModules
│   └── FindYamlCpp.cmake (CMake module to find the yaml-cpp library)
|
├── Example
│   ├── example.ftl (MODFLOW output containing the velocity field)
│   ├── example.yaml (Example of YAML parameter file)
│   └── output (Folder where the output is written)
|
├── Geometry
│   ├── CartesianGrid.cuh (Data stuctures and functions of a Cartesian grid)
│   ├── CellField.cuh (A cellfield is a field that is defined at the center of every cells of the grid)
│   ├── CornerField.cuh (A cornerfield is a field that is defined at the each corner of every cells of the grid)
│   ├── FaceField.cuh (A facefield is a variable that is defined at the center of of each cell interface of the grid)
│   ├── Interpolation.cuh (Interpolation methods)
│   ├── Point.cuh (Operations for points)
│   └── Vector.cuh (Operations for vectors)
|
├── Particles
│   ├── MoveParticle.cuh (Kernel used to update particle position)
│   ├── PParticles.cu (Implementation file for PParticles class)
│   └── PParticles.cuh (Header file for PParticles class)
|
└── Utilities
    ├── Parameters.h (Parse the YAML configuration file and create the simulation parameters)
    └── spdlog (Folder containing the spdlog library for logging)
```
