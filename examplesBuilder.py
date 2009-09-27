import os, sys

cwd = os.getcwd();
logDir = cwd + "/../../build/logs";
rootUrl = "http://autotest.bioeng.auckland.ac.nz/opencmiss-build/"

if not os.path.isdir(logDir):
  os.mkdir(cwd + "/../../build")
  os.mkdir(logDir);
compiler = sys.argv[1];


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
     print "Building %s: <a class='fail' href='%slogs_x86_64-linux/%s/build-%s'>failed</a><br>" %(path,rootUrl,path,compiler)
   os.chdir(cwd)
   return;

buildExample("ClassicalField/AnalyticLaplace")
buildExample("ClassicalField/Diffusion")
buildExample("ClassicalField/Helmholtz")
buildExample("ClassicalField/Laplace")
buildExample("ClassicalField/NonlinearPoisson")

buildExample("Bioelectrics/Monodomain")
  
buildExample("FluidMechanics/Stokes/HexChannel")
buildExample("FluidMechanics/Stokes/HexPipe")
buildExample("FluidMechanics/Stokes/SingleElement")
buildExample("FluidMechanics/Stokes/VesselPipe")

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


