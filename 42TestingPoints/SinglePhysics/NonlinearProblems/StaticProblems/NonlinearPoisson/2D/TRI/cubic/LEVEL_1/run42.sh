#!/bin/bash
$OPENCMISS_ROOT/cm/examples/ClassicalField/NonlinearPoisson/AnalyticNonlinearPoisson/bin/x86_64-linux/mpich2/gnu/AnalyticNonlinearPoissonExample 5 5 0 9
mv *.exnode output/
mv *.exelem output/
