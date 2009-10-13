# -*- coding: utf-8 -*-
import os, sys

cwd = os.getcwd();
logDir = cwd + "/../../build/logs";
rootUrl = "http://autotest.bioeng.auckland.ac.nz/opencmiss-build/"

if not os.path.isdir(logDir):
  os.mkdir(cwd + "/../../build")
  os.mkdir(logDir);
compiler = sys.argv[1];
f = open(logDir+'/failedBuilds',"w")

def buildExample(path) :
   global compiler,logDir;
   newDir = logDir
   for folder in path.split('/') :
     newDir = newDir + '/' + folder
     if not os.path.isdir(newDir):
       os.mkdir(newDir)
   os.chdir(path)
   if os.path.exists(newDir + "/build-" + compiler) :
     os.remove(newDir + "/build-" + compiler)
   err=os.system("make COMPILER=" + compiler + " > " + newDir + "/build-" + compiler +" 2>&1")
   if err==0 :
     print "Building %s: <a class='success' href='%slogs_x86_64-linux/%s/build-%s'>success</a><br>" %(path,rootUrl,path,compiler)
   else :
     f.write(path+'\n')
     print "Building %s: <a class='fail' href='%slogs_x86_64-linux/%s/build-%s'>failed</a><br>" %(path,rootUrl,path,compiler)
   os.chdir(cwd)
   return;

buildExample("ClassicalField/AnalyticLaplace")
buildExample("ClassicalField/Diffusion")
buildExample("ClassicalField/Helmholtz")
buildExample("ClassicalField/Laplace")
buildExample("ClassicalField/NonlinearPoisson")

buildExample("Bioelectrics/Monodomain")

buildExample("FluidMechanics/Stokes/HEX_CHANNEL/SteadyState")
buildExample("FluidMechanics/Stokes/HEX_CHANNEL/Transient")
buildExample("FluidMechanics/Stokes/SINGLE_ELEMENTS/SteadyState")
buildExample("FluidMechanics/Stokes/SINGLE_ELEMENTS/Transient")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/hex/Cubic/DirectSolver/SteadyState")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/hex/Cubic/DirectSolver/Transient")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/hex/Cubic/IterativeSolver/SteadyState")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/hex/Cubic/IterativeSolver/Transient")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/hex/Quadratic/DirectSolver/SteadyState")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/hex/Quadratic/DirectSolver/Transient")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/hex/Quadratic/IterativeSolver/SteadyState")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/hex/Quadratic/IterativeSolver/Transient")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/quad/Cubic/DirectSolver/SteadyState")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/quad/Cubic/DirectSolver/Transient")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/quad/Cubic/IterativeSolver/SteadyState")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/quad/Cubic/IterativeSolver/Transient")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/quad/Quadratic/DirectSolver/SteadyState")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/quad/Quadratic/DirectSolver/Transient")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/quad/Quadratic/IterativeSolver/SteadyState")
buildExample("FluidMechanics/Stokes/TESTING_ELEMENTS/quad/Quadratic/IterativeSolver/Transient")

buildExample("FluidMechanics/NavierStokes/HEX_CHANNEL/SteadyState")
buildExample("FluidMechanics/NavierStokes/HEX_CHANNEL/Transient")
buildExample("FluidMechanics/NavierStokes/SINGLE_ELEMENTS/SteadyState")
buildExample("FluidMechanics/NavierStokes/SINGLE_ELEMENTS/Transient")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/hex/Cubic/DirectSolver/SteadyState")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/hex/Cubic/DirectSolver/Transient")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/hex/Cubic/IterativeSolver/SteadyState")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/hex/Cubic/IterativeSolver/Transient")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/hex/Quadratic/DirectSolver/SteadyState")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/hex/Quadratic/DirectSolver/Transient")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/hex/Quadratic/IterativeSolver/SteadyState")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/hex/Quadratic/IterativeSolver/Transient")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/quad/Cubic/DirectSolver/SteadyState")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/quad/Cubic/DirectSolver/Transient")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/quad/Cubic/IterativeSolver/SteadyState")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/quad/Cubic/IterativeSolver/Transient")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/quad/Quadratic/DirectSolver/SteadyState")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/quad/Quadratic/DirectSolver/Transient")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/quad/Quadratic/IterativeSolver/SteadyState")
buildExample("FluidMechanics/NavierStokes/TESTING_ELEMENTS/quad/Quadratic/IterativeSolver/Transient")

buildExample("FluidMechanics/Darcy/ConvergenceStudy")
buildExample("FluidMechanics/Darcy/FiveSpotProblem")
buildExample("FluidMechanics/Darcy/VenousCompartment")

buildExample("FiniteElasticity/UniAxialExtension")
buildExample("FiniteElasticity/TwoElementTriLinear")
buildExample("FiniteElasticity/MixedBoundaryConditions")
buildExample("FiniteElasticity/TriCubicAxialExtension")

buildExample("LinearElasticity/2DPlaneStressLagrangeBasis")
buildExample("LinearElasticity/2DPlaneStressLagrangeBasisAnalytic")
buildExample("LinearElasticity/3DLagrangeBasis")
buildExample("LinearElasticity/3DLagrangeBasisAnisotropicFibreField")
buildExample("LinearElasticity/3DCubicHermiteBasis")

# TODO Group them
buildExample("LagrangeSimplexMesh")
buildExample("cellml")
buildExample("define-geometry-and-export")
buildExample("MoreComplexMesh")
buildExample("simple-field-manipulation-direct-access")
buildExample("SimplexMesh")
buildExample("TwoRegions")

f.close()
