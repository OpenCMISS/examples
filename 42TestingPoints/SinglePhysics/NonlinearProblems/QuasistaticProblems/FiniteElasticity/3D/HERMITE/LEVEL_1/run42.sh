#!/bin/bash
$OPENCMISS_ROOT/cm/examples/FiniteElasticity/testingPoints/bin/x86_64-linux/mpich2/gnu_4.4/testingPointsExample  -DIM=3D -ELEM=HERMITE -BASIS=hermite -LEVEL=1
mv *.exnode *.exelem output/
