# OpenCMISS (Open Continuum Mechanics, Imaging, Signal processing and System identification)
# is a mathematical modelling environment that enables the application of finite element
# analysis techniques to a variety of complex bioengineering problems.
# 
# The OpenCMISS project website can be found at http://www.opencmiss.org
# For more information see http://www.opencmiss.org/documentation

# Include the OpenCMISS initalisation & macro definition file.
include(./OpenCMISS.cmake)

# Project setup
# -------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(OpenCMISS-Example VERSION 1.0 LANGUAGES C)
enable_language(Fortran OPTIONAL)
enable_testing()

# Find OpenCMISS
# --------------
#
# This call tries to find a matching OpenCMISS (SDK) installation for your (default) Toolchain and MPI choice. 
# Currently available OpenCMISS package components: Iron Zinc
# Requiring Iron will also include the iron_c bindings target if available.
#
# This function creates a link target 'opencmiss' that contains references to the requested OpenCMISS components.
# For Iron, it also adds the necessary MPI information - no further find_package(MPI ..) and setting of extra include
# paths needed!               
find_package(OpenCMISS 1.1 REQUIRED COMPONENTS Iron CONFIG)

# CMake application code
# ----------------------
#add_subdirectory(C)
if (CMAKE_Fortran_COMPILER)
    add_subdirectory(Fortran)
endif()
