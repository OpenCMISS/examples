#!/bin/bash
$OPENCMISS_ROOT/cm/examples/FiniteElasticity/testingPoints/bin/x86_64-linux/mpich2/gnu/testingPointsExample  -DIM=3D -ELEM=TET -BASIS=cubic -LEVEL=1
mv *.exnode *.exelem output/
