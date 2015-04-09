cmake_minimum_required(VERSION 3.0)
# Get lowercase name
string(TOLOWER "@EXAMPLE_NAME@" EXAMPLE_NAME)

project(OpenCMISSExample-${EXAMPLE_NAME} LANGUAGES C CXX Fortran VERSION 1.0)
SET(CMAKE_NO_SYSTEM_FROM_IMPORTED YES)

# Have CMake find the package components
list(APPEND CMAKE_PREFIX_PATH @OPENCMISS_PREFIX_PATH@) # var defined in ~CMakeLists.main.template.cmake:90

# Grab config files for packages
find_package(Iron REQUIRED)
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

# Get sources in /src
file(GLOB SRC src/*.f90 src/*.c ../input/*.f90)
# Add example executable
add_executable(${EXAMPLE_NAME} ${SRC})
# Link to iron - contains forward refs to all other necessary libs
target_link_libraries(${EXAMPLE_NAME} PUBLIC iron)
set(CMAKE_Fortran_FLAGS "-cpp")