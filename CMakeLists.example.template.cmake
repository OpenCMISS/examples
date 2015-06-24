# This is the main build file for OpenCMISS examples using the CMake build system
#
# If standard procedure has been followed building OpenCMISS using CMake, all you
# need to do is set OPENCMISS_INSTALL_DIR to the <OPENCMISS_ROOT>/install/[release|debug|...] directory.
# The script will do the rest.
#
# Otherwise, if the FindOpenCMISS.cmake module is located elsewhere on your system
# (it is placed inside the OpenCMISS installation folder by default), you need to additionally add that path to
# the CMAKE_MODULE_PATH variable.

# Get the build context from the OpenCMISS installation
if (NOT OPENCMISS_INSTALL_DIR)
    set(OPENCMISS_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../install)
    get_filename_component(OPENCMISS_INSTALL_DIR ${OPENCMISS_INSTALL_DIR} ABSOLUTE)
endif()
set(CMAKE_MODULE_PATH ${OPENCMISS_INSTALL_DIR})
find_package(OpenCMISS REQUIRED)

#################### Actual example code ####################
cmake_minimum_required(VERSION ${OPENCMISS_CMAKE_MIN_VERSION} FATAL_ERROR)
project(OpenCMISS-Example VERSION 1.0 LANGUAGES Fortran C CXX)

OPENCMISS_IMPORT()

# Get sources in /src
file(GLOB SRC src/*.f90 src/*.c ../input/*.f90)
set(EXAMPLE_TARGET run)
# Add example executable
add_executable(${EXAMPLE_TARGET} ${SRC})
# Link to opencmiss - contains forward refs to all other necessary libs
target_link_libraries(${EXAMPLE_TARGET} PRIVATE opencmiss)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
if (WIN32)
    target_compile_definitions(${EXAMPLE_TARGET} PRIVATE NOMPIMOD)
endif()
install(TARGETS ${EXAMPLE_TARGET} DESTINATION .)