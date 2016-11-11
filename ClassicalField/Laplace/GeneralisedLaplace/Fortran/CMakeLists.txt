# Add example executable
add_executable(example1.x src/FortranExample.f90)

# Turn on Fortran preprocessing (#include directives)
if (MSVC)
    set(PROC_OPT "/fpp")
else()
    set(PROC_OPT "-cpp")
endif()
target_compile_options(example1.x PRIVATE ${PROC_OPT})

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
target_link_libraries(example1.x PRIVATE opencmiss)

# Add a simple test that runs the executable
add_test(NAME Laplace_Fortran COMMAND example1.x)
add_opencmiss_environment(Laplace_Fortran)

add_executable(example2.x src/FortranExample_2D_testproblem.f90)
target_compile_options(example2.x PRIVATE ${PROC_OPT})
target_link_libraries(example2.x PRIVATE opencmiss)

add_executable(example3.x src/FortranExample_3D_testproblem.f90)
target_compile_options(example3.x PRIVATE ${PROC_OPT})
target_link_libraries(example3.x PRIVATE opencmiss)

add_executable(example4.x src/FortranExample_3D_tests.f90)
target_compile_options(example4.x PRIVATE ${PROC_OPT})
target_link_libraries(example4.x PRIVATE opencmiss)

