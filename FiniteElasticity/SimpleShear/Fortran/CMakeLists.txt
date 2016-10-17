# Turn on Fortran preprocessing (#include directives)
if (MSVC)
    set(PROC_OPT "/fpp")
else()
    set(PROC_OPT "-cpp")
endif()

# Add first example executable
add_executable(shear_ex1.x src/FortranExample.f90)
target_compile_options(shear_ex1.x PRIVATE ${PROC_OPT})
target_link_libraries(shear_ex1.x PRIVATE opencmiss)

# Add second example executable
add_executable(shear_ex2.x src/FortranExample_incomp.f90)
target_compile_options(shear_ex2.x PRIVATE ${PROC_OPT})
target_link_libraries(shear_ex2.x PRIVATE opencmiss)

# Add third example executable
add_executable(shear_ex3.x src/FortranExample_nearly.f90)
target_compile_options(shear_ex3.x PRIVATE ${PROC_OPT})
target_link_libraries(shear_ex3.x PRIVATE opencmiss)

