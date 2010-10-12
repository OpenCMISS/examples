#! /bin/bash
$OPENCMISS_ROOT/cm/examples/ClassicalField/Laplace/42Master/bin/x86_64-linux/mpich2/gnu/42MasterExample-debug -2D -quad -linearbasis -nx 10 -ny 10
mv Laplace* output/
mv Diagnostics.diag output/
