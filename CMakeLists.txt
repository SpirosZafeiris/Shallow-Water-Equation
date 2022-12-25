# Set Compiler 
cmake_host_system_information(RESULT host QUERY HOSTNAME)

set(ENV{FC} mpiifort)

cmake_minimum_required(VERSION 3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(Shallow_Water_Equation  VERSION 1.0
      	   DESCRIPTION "D'Alembert Solver"
                         LANGUAGES Fortran)

# SET compiled module path
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/modules/)

# SET BUILD OPTIONS
set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release")
endif()

# Source Dirs
add_subdirectory(src)

# SET COMPILER FOR THE TWO BUILD TYPES
set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -qopenmp") #-qmkl=parallel 
set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -pg -qmkl -C -init=snan -traceback -fpe0 -fpp -diag-enable=all -warn all")

# DEFINE EXECUTABLES
add_executable (Swe ${OBJS})
