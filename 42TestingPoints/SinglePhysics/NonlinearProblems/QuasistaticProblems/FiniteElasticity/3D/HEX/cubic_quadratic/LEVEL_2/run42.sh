#!/bin/bash
$OPENCMISS_ROOT/cm/examples/FiniteElasticity/testingPoints/bin/x86_64-linux/mpich2/gnu/testingPointsExample  -DIM=3D -ELEM=HEX -BASIS_1=cubic -BASIS_2=quadratic -LEVEL=2
mv *.exnode *.exelem output/