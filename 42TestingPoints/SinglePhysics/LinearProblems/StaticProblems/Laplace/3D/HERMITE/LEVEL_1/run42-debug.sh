#! /bin/bash
$OPENCMISS_ROOT/cm/examples/ClassicalField/Laplace/42Master/bin/x86_64-linux/mpich2/gnu_4.4/42MasterExample-debug -3D -hermite -nx 8 -ny 8 -nz 8
mv Laplace* output/
