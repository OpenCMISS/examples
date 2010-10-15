#! /bin/bash
$OPENCMISS_ROOT/cm/examples/ClassicalField/Laplace/42Master/bin/x86_64-linux/mpich2/gnu/42MasterExample -3D -hermite -nx 4 -ny 4 -nz 4
mv Laplace* output/
