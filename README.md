# PAR²
PAR² is a Lagrangian solute transport simulator using a parallelized Random Walk Particle Tracking (RWPT) method.

## Getting Started
You can download the latest release of PAR² for Windows or Linux executable from [here](https://github.com/GerryR/par2/releases). Make sure your computer is equipped with an NVIDIA GPU and NVIDIA drivers are updated. You can control the simulation parameters through a YAML configuration file. Look inside the Examples folder to get started.

## Build
The following software and libraries must be installed:

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
