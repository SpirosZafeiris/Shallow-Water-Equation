# Set Compiler 
cmake_host_system_information(RESULT host QUERY HOSTNAME)

#set(ENV{FC} ifort)



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
set(LIBS "")

option(CGNS "Include CGNS library" OFF)
if (CGNS)
   set(CGNS_LIB_DIR ${CMAKE_SOURCE_DIR}/lib/cgns/lib)
   set(CGNS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/lib/cgns/include)
   find_library(LIBCGNS NAMES libcgns.so HINTS ${CMAKE_SOURCE_DIR}/lib/cgns/lib)
   string(APPEND CMAKE_Fortran_FLAGS " -I ${CGNS_INCLUDE_DIR} \
                       -L ${CGNS_LIB_DIR} ")
                      list(APPEND LIBS ${LIBCGNS})
endif()




# Source Dirs
add_subdirectory(src)

# Set some default compilation settings for Fortran compiler
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
   set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
   set(CMAKE_Fortran_FLAGS_DEBUG "-fbounds-check -fbacktrace -fdefault-real-8 -ffree-line-length-none")
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
   set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
   set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -pg -qmkl -C -init=snan -traceback -fpe0 -fpp -diag-enable=all -warn all")
endif()

# DEFINE EXECUTABLES
add_executable (shallow_water_equation.out ${OBJS} ${OBJS_CGNS})

target_link_libraries(shallow_water_equation.out PRIVATE ${LIBS})
