
# Add example executable
add_executable(ReactionDiffusionNoSource1D  src/ReactionDiffusionNoSource1DExample.f90)

# Turn on Fortran preprocessing (#include directives)
if (MSVC)
    set(PROC_OPT "/fpp")
else()
    set(PROC_OPT "-cpp")
endif()

# in debug mode, add flag that make Floating point exceptions fatal
set(OPT "")
if (FALSE)
    if(CMAKE_BUILD_TYPE MATCHES DEBUG)
        set(OPT "-ffpe-trap=invalid,zero,underflow,overflow,denormal")		# GCC: make floating point exceptions fatal
    endif()
endif()

target_compile_options(ReactionDiffusionNoSource1D PRIVATE ${PROC_OPT} ${OPT} "-fbacktrace")

# Link to opencmiss
# -----------------
#
# This simply uses the opencmiss target created by OC_INIT() above.
#
# Alternatively, you can also directly invoke 'find_package(Iron <IRON_VERSION>)' to explicitly
# require that version. But then you have to deal with setting up the correct MPI include directives and library paths.
#
# If required, add any other required link libraries (cellml, petsc ..) here, too.
# For example, if you needed PetSC functionality, issue
#
#     find_package(PETSC <PETSC_VERSION> REQUIRED)
#     target_link_libraries(${EXAMPLE_BINARY} PRIVATE petsc)
#
# All the OpenCMISS dependencies provide a target you can link against corresponding to the (lowercase) component name.
target_link_libraries(ReactionDiffusionNoSource1D PRIVATE opencmiss)

# Add a simple test that runs the executable
add_test(NAME test_ReactionDiffusionNoSource1D COMMAND ReactionDiffusionNoSource1D)
add_opencmiss_environment(test_ReactionDiffusionNoSource1D)

###################
# Developer notice!
#
# If you write Fortran code and use MPI, you need to use the following MPI directives:
#
# #ifndef NOMPIMOD
#   USE MPI
# #endif
# [...]
#   IMPLICIT NONE
# [...]
# #ifdef NOMPIMOD
#     #include "mpif.h"
# #endif
#
# Reasoning: In some cases like Windows/MPICH2 there sometimes is no mpi.mod file. In order to yet make
# the example work the build system adds the definition 'NOMPIMOD', which can be checked and acted to accordingly.
#
