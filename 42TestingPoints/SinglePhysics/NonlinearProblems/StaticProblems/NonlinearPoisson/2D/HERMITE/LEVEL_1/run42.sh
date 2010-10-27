#!/bin/bash
$OPENCMISS_ROOT/cm/examples/ClassicalField/NonlinearPoisson/AnalyticNonlinearPoisson/bin/x86_64-linux/mpich2/gnu/AnalyticNonlinearPoissonExample 3 3 0 4
mv *.exnode output/
mv *.exelem output/
