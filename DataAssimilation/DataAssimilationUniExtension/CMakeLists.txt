# Add example executable
add_executable(data_assim_uni1.x src/data_assimilation_uniextension.cpp)
target_link_libraries(data_assim_uni1.x PRIVATE opencmiss)

# Add example executable
add_executable(data_assim_uni2.x src/data_assimilation_uniextension_OpenCMISS_test.cpp)
target_link_libraries(data_assim_uni2.x PRIVATE opencmiss)

# Add example executable
add_executable(data_assim_uni3.x src/UniAxialExtensionExample.f90)
target_link_libraries(data_assim_uni3.x PRIVATE opencmiss)
