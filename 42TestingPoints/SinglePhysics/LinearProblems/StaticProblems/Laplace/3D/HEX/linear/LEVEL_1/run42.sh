#! /bin/bash
$OPENCMISS_ROOT/cm/examples/ClassicalField/Laplace/42Master/bin/x86_64-linux/mpich2/gnu/42MasterExample -3D -hex -linearbasis -nx 10 -ny 10 -nz 10
mv Laplace* output/
