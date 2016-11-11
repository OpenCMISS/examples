# Turn on Fortran preprocessing (#include directives)
if (MSVC)
    set(PROC_OPT "/fpp")
else()
    set(PROC_OPT "-cpp")
endif()

# Build the Fortran example
add_executable(laplace_fortran Fortran/src/FortranExample.f90)
target_compile_options(laplace_fortran PRIVATE ${PROC_OPT})
target_link_libraries(laplace_fortran PRIVATE opencmiss)

# Build the C example
add_executable(laplace_c C/src/CExample.c)
target_link_libraries(laplace_c PRIVATE opencmiss)

#add_test(NAME Laplace_Fortran COMMAND laplace_fortran)
#add_opencmiss_environment(Laplace_Fortran)

