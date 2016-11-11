# Add example executable
add_executable(laplace_c src/CExample.c)

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
target_link_libraries(laplace_c PRIVATE opencmiss)

# Add a simple test that runs the executable
add_test(NAME Laplace_C COMMAND laplace_c)
add_opencmiss_environment(Laplace_C)

