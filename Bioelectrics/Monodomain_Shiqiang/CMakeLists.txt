# Add example executable
add_executable(monodomain_shiqiang.x src/Monodomain_ShiqiangExample.f90)

# Turn on Fortran preprocessing (#include directives)
if (MSVC)
    set(PROC_OPT "/fpp")
else()
    set(PROC_OPT "-cpp")
endif()
target_compile_options(monodomain_shiqiang.x PRIVATE ${PROC_OPT})

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
target_link_libraries(monodomain_shiqiang.x PRIVATE opencmiss)

# Add a simple test that runs the executable
add_test(NAME MonoDomain_Shiqiang COMMAND monodomain_shiqiang.x)
add_opencmiss_environment(MonoDomain_Shiqiang)

