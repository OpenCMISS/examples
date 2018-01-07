# CMake generated Testfile for 
# Source directory: /store/software/opencmiss/OpenCMISS-examples/FiniteElasticity/RightBiceps/Fortran
# Build directory: /store/software/opencmiss/OpenCMISS-examples/FiniteElasticity/RightBiceps/build_debug/Fortran
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_RightBiceps "/store/software/opencmiss/OpenCMISS-examples/FiniteElasticity/RightBiceps/build_debug/Fortran/RightBiceps")
set_tests_properties(test_RightBiceps PROPERTIES  ENVIRONMENT "LD_LIBRARY_PATH=/store/software/opencmiss/iron_maierbn/install/x86_64_linux/gnu-5.4-F5.4/openmpi_release/debug/bin::/home/maierbn/perl5/perlbrew/build/perl-5.10.0:/store/software/opendihu/dependencies/petsc/install/lib:/store/software/simbody/install:/store/software/opensim/opensim-core/install")
