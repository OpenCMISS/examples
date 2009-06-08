import os, sys

cwd = os.getcwd();
success = 1;

if os.path.isdir("logs"):
  os.system("rm -r logs")
os.mkdir("logs");
logFile = open("logs/log.txt","w")
compiler = sys.argv[1];
failedExample = '';

def buildExample(path) :
  global failedExample, success, compiler;
  os.chdir(path)
  if compiler == 'gnu' :
    os.system("make COMPILER=gnu > " + cwd + "/logs/" + path.replace('/', '_'))
  elif compiler == 'intel' :
    err=os.system("make > " + cwd + "/logs/" + path.replace('/', '_'))
    if err==0 :
      logFile.write(path.replace('/', '_')+'=success\n') 
    else :
      success=0
      logFile.write(path.replace('/', '_')+'=fail\n') 
      failedExample += path.replace('/', ' - ') + ' '
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

# TODO Group them
buildExample("LagrangeSimplexMesh")
buildExample("LinearElasticity")
buildExample("cellml")
buildExample("Darcy")
buildExample("define-geometry-and-export")
buildExample("FiniteElasticity")
buildExample("MoreComplexMesh")
buildExample("simple-field-manipulation-direct-access")
buildExample("SimplexMesh")
buildExample("TwoRegions")

logFile.close()
if success==0 :
  raise RuntimeError, 'Failed in %s' % (failedExample)


