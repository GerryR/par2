# https://github.com/Astron/Astron/blob/master/cmake/modules/FindYamlCpp.cmake
#
# What follows is the Modified BSD License.
# Copyright (c) 2014, the Astron Authors. All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# Neither the name of the copyright holder nor the names of its authors or contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE AUTHORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Locate yaml-cpp
#
# This module defines
#  YAMLCPP_FOUND, if false, do not try to link to yaml-cpp
#  YAMLCPP_LIBNAME, name of yaml library
#  YAMLCPP_LIBRARY, where to find yaml-cpp
#  YAMLCPP_LIBRARY_RELEASE, where to find Release or RelWithDebInfo yaml-cpp
#  YAMLCPP_LIBRARY_DEBUG, where to find Debug yaml-cpp
#  YAMLCPP_INCLUDE_DIR, where to find yaml.h
#  YAMLCPP_LIBRARY_DIR, the directories to find YAMLCPP_LIBRARY
#
# By default, the dynamic libraries of yaml-cpp will be found. To find the static ones instead,
# you must set the YAMLCPP_USE_STATIC_LIBS variable to TRUE before calling find_package(YamlCpp ...)

# attempt to find static library first if this is set
if(YAMLCPP_USE_STATIC_LIBS)
    set(YAMLCPP_STATIC libyaml-cpp.a)
    set(YAMLCPP_STATIC_DEBUG libyaml-cpp-dbg.a)
endif()

if(${MSVC})    ### Set Yaml libary name for Windows
  set(YAMLCPP_LIBNAME "libyaml-cppmt" CACHE STRING "Name of YAML library")
  set(YAMLCPP_LIBNAME optimized ${YAMLCPP_LIBNAME} debug ${YAMLCPP_LIBNAME}d)
else()                      ### Set Yaml libary name for Unix, Linux, OS X, etc
  set(YAMLCPP_LIBNAME "yaml-cpp" CACHE STRING "Name of YAML library")
endif()

# find the yaml-cpp include directory
find_path(YAMLCPP_INCLUDE_DIR
  NAMES yaml-cpp/yaml.h
  PATH_SUFFIXES include
  PATHS
    ${PROJECT_SOURCE_DIR}/dependencies/yaml-cpp-0.5.1/include
    ~/Library/Frameworks/yaml-cpp/include/
    /Library/Frameworks/yaml-cpp/include/
    /usr/local/include/
    /usr/include/
    /sw/yaml-cpp/         # Fink
    /opt/local/yaml-cpp/  # DarwinPorts
    /opt/csw/yaml-cpp/    # Blastwave
    /opt/yaml-cpp/
	${YAML_ROOT}/include/)

# find the release yaml-cpp library
find_library(YAMLCPP_LIBRARY_RELEASE
  NAMES ${YAMLCPP_STATIC} yaml-cpp libyaml-cppmt.lib
  PATH_SUFFIXES lib64 lib Release RelWithDebInfo
  PATHS
    ${PROJECT_SOURCE_DIR}/dependencies/yaml-cpp-0.5.1/
    ${PROJECT_SOURCE_DIR}/dependencies/yaml-cpp-0.5.1/build
    ~/Library/Frameworks
    /Library/Frameworks
    /usr/local
    /usr
    /sw
    /opt/local
    /opt/csw
    /opt
	${YAML_ROOT}/bin/Release)

# find the debug yaml-cpp library
find_library(YAMLCPP_LIBRARY_DEBUG
  NAMES ${YAMLCPP_STATIC_DEBUG} yaml-cpp-dbg libyaml-cppmtd.lib
  PATH_SUFFIXES lib64 lib Debug
  PATHS
    ${PROJECT_SOURCE_DIR}/dependencies/yaml-cpp-0.5.1/
    ${PROJECT_SOURCE_DIR}/dependencies/yaml-cpp-0.5.1/build
    ~/Library/Frameworks
    /Library/Frameworks
    /usr/local
    /usr
    /sw
    /opt/local
    /opt/csw
    /opt
	${YAML_ROOT}/bin/Debug)

# set library vars
set(YAMLCPP_LIBRARY ${YAMLCPP_LIBRARY_RELEASE})
if(CMAKE_BUILD_TYPE MATCHES Debug AND EXISTS ${YAMLCPP_LIBRARY_DEBUG})
  set(YAMLCPP_LIBRARY ${YAMLCPP_LIBRARY_DEBUG})
endif()

get_filename_component(YAMLCPP_LIBRARY_RELEASE_DIR ${YAMLCPP_LIBRARY_RELEASE} PATH)
get_filename_component(YAMLCPP_LIBRARY_DEBUG_DIR ${YAMLCPP_LIBRARY_DEBUG} PATH)
set(YAMLCPP_LIBRARY_DIR ${YAMLCPP_LIBRARY_RELEASE_DIR} ${YAMLCPP_LIBRARY_DEBUG_DIR})

# handle the QUIETLY and REQUIRED arguments and set YAMLCPP_FOUND to TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(YamlCpp DEFAULT_MSG
  YAMLCPP_INCLUDE_DIR
  YAMLCPP_LIBRARY
  YAMLCPP_LIBRARY_DIR)
mark_as_advanced(
  YAMLCPP_INCLUDE_DIR
  YAMLCPP_LIBRARY_DIR
  YAMLCPP_LIBRARY
  YAMLCPP_LIBRARY_RELEASE
  YAMLCPP_LIBRARY_RELEASE_DIR
  YAMLCPP_LIBRARY_DEBUG
  YAMLCPP_LIBRARY_DEBUG_DIR)
