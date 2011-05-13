#! /bin/bash
$OPENCMISS_ROOT/cm/examples/ClassicalField/Laplace/42Master/bin/x86_64-linux/mpich2/gnu_4.4/42MasterExample -2D -tri -quadraticbasis -nx 10 -ny 10
mv Laplace* output/
