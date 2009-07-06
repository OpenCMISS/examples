import os, sys

cwd = os.getcwd();
success = 1;
logDir = cwd + "/../../build/logs";
rootUrl = "https://autotest.bioeng.auckland.ac.nz/opencmiss-admin/"

if not os.path.isdir(logDir):
  os.mkdir(cwd + "/../../build")
  os.mkdir(logDir);
compiler = sys.argv[1];


def buildExample(path) :
   global success,compiler,logDir;
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
     print "Building " + path+': success'
   else :
     success=0
     print "Building " + path+': fail'
   os.chdir(cwd)
   return;

buildExample("ClassicalField/AnalyticLaplace")
buildExample("ClassicalField/Diffusion")
buildExample("ClassicalField/Helmholtz")
buildExample("ClassicalField/Laplace")
buildExample("ClassicalField/NonlinearPoisson")

buildExample("Bioelectrics/Monodomain")
  
buildExample("FluidMechanics/StokesFlow/HexChannel")
buildExample("FluidMechanics/StokesFlow/HexPipe")
buildExample("FluidMechanics/StokesFlow/SingleElement")
buildExample("FluidMechanics/StokesFlow/VesselPipe")
buildExample("FluidMechanics/Darcy")

buildExample("FiniteElasticity/UniAxialExtension")
buildExample("FiniteElasticity/TwoElementTriLinear")
buildExample("FiniteElasticity/MixedBoundaryConditions")
buildExample("FiniteElasticity/TriCubicAxialExtension")

# TODO Group them
buildExample("LagrangeSimplexMesh")
buildExample("LinearElasticity")
buildExample("cellml")
buildExample("define-geometry-and-export")
buildExample("MoreComplexMesh")
buildExample("simple-field-manipulation-direct-access")
buildExample("SimplexMesh")
buildExample("TwoRegions")

print "See %slogs_x86_64-linux-%s for detail" % (rootUrl, compiler) 
if success==0 :
  raise RuntimeError


