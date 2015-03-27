cmake_minimum_required(VERSION 3.0)

set(OPENCMISS_ROOT /home/dwirtz/software/opencmiss)
set(ARCH_DIR /x86_64_linux/gnu-4.8/openmpi/static/release)
# Have CMake find the package components
list(APPEND CMAKE_PREFIX_PATH ${OPENCMISS_ROOT}/install/${ARCH_DIR})

project(OpenCMISS-Example LANGUAGES C Fortran VERSION 1.0)
SET(CMAKE_NO_SYSTEM_FROM_IMPORTED YES)

find_package(Iron REQUIRED) #${IRON_VERSION}
find_package(LIBCELLML REQUIRED)
find_package(CELLML-API REQUIRED)
find_package(PETSC REQUIRED)
find_package(SUNDIALS REQUIRED)
find_package(HYPRE REQUIRED)
find_package(MUMPS REQUIRED)
find_package(SCALAPACK REQUIRED)
find_package(PASTIX REQUIRED)
find_package(SUITESPARSE REQUIRED)
find_package(PTSCOTCH REQUIRED)
find_package(PARMETIS REQUIRED)
find_package(SUPERLU REQUIRED)
find_package(SUPERLU_DIST REQUIRED)

get_filename_component(EXAMPLE_PATH ${CMAKE_CURRENT_BINARY_DIR} PATH)
get_filename_component(EXAMPLE_NAME ${CMAKE_CURRENT_BINARY_DIR} NAME)
string(TOLOWER ${EXAMPLE_NAME} EXAMPLE_NAME)

file(GLOB SRC src/*.f90 src/*.c ../input/*.f90)
set(TARGET_NAME ${EXAMPLE_NAME})
add_executable(${TARGET_NAME} ${SRC})
target_link_libraries(${TARGET_NAME} PUBLIC iron)
set(CMAKE_Fortran_FLAGS "-cpp")
 