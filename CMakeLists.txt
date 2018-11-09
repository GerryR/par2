cmake_minimum_required(VERSION 3.8)
project(par2 LANGUAGES CXX CUDA)

# LINUX
# Create makefile:
# > cmake -DCMAKE_BUILD_TYPE=Release ..
# Compile:
# > make
#
# WINDOWS
# Create MSVC solution:
# > cmake -G "Visual Studio 15 2017 Win64" -T v140 -DCMAKE_BUILD_TYPE=Release -DYAML_ROOT=C:/path/to/yaml-cpp ..
# Compile using the Developer Command Prompt for VS:
# > devenv /Build Release par2.sln

set(CUDA_PROPAGATE_HOST_FLAGS OFF)

# Add custom cmake module for YamlCpp
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules
                      ${CMAKE_MODULE_PATH})

# Search YamlCpp library
set(YAMLCPP_USE_STATIC_LIBS ON)
find_package(YamlCpp REQUIRED)
include_directories(${YAMLCPP_INCLUDE_DIR})

# In MSVC, use static compilation to have a standalone executable
if(MSVC)
    foreach(flag_var
            CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
            CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO
            CMAKE_CUDA_FLAGS CMAKE_CUDA_FLAGS_DEBUG CMAKE_CUDA_FLAGS_RELEASE)
        if(${flag_var} MATCHES "/MD")
            string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
            set(${flag_var} "${${flag_var}}" CACHE INTERNAL "")
        endif(${flag_var} MATCHES "/MD")
        if(${flag_var} MATCHES "-MD")
            string(REGEX REPLACE "-MD" "-MT" ${flag_var} "${${flag_var}}")
            set(${flag_var} "${${flag_var}}" CACHE INTERNAL "")
        endif(${flag_var} MATCHES "-MD")
    endforeach(flag_var)
endif()

# Compiler definitions
add_definitions(-DPAR2_FLOAT=float)
add_definitions(-arch=sm_30
 -gencode=arch=compute_30,code=sm_30
 -gencode=arch=compute_50,code=sm_50
 -gencode=arch=compute_52,code=sm_52
 -gencode=arch=compute_60,code=sm_60
 -gencode=arch=compute_61,code=sm_61
 -gencode=arch=compute_62,code=sm_62
 -gencode=arch=compute_70,code=sm_70
 -gencode=arch=compute_70,code=compute_70)

# Create executable
add_executable(par2 maingpu.cu)

# Enable C++11 features
target_compile_features(par2 PUBLIC cxx_std_11)
set_property(TARGET par2 PROPERTY CUDA_STANDARD 11)

# Link libraries to executable
target_link_libraries(par2 ${CUDA_curand_LIBRARY} ${CUDA_LIBRARIES} ${YAMLCPP_LIBRARY})
