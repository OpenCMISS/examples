import os, sys

cwd = os.getcwd();
success = 1;
logDir = cwd + "/../../build/logs";
rootUrl = "http://autotest.bioeng.auckland.ac.nz:8888/"

if os.path.isdir(logDir):
  os.system("rm -r " + logDir)
  os.rmdir(cwd + "/../../build")
os.mkdir(cwd + "/../../build")
os.mkdir(logDir);
logFile = open(logDir + "/log.txt","w")
compiler = sys.argv[1];
failedExample = '';
successedExample = ''

def buildExample(path) :
  global failedExample,successedExample,success,compiler,logDir;
  os.chdir(path)
  err=0
  if compiler == 'gnu' :
    err=os.system("make COMPILER=gnu > " + logDir + "/" + path.replace('/', '_')+ " 2>&1")
  elif compiler == 'intel' :
    err=os.system("make > " + logDir + "/" + path.replace('/', '_') + " 2>&1")
  if err==0 :
    logFile.write(path.replace('/', '_')+'=success\n')
    print path.replace('/', '_')+': success'
    successedExample += path.replace('/', ' - ') + '\n'
  else :
    success=0
    logFile.write(path.replace('/', '_')+'=fail\n') 
    print path.replace('/', '_')+': fail\n'
    failedExample += path.replace('/', ' - ') + '\n'
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

buildExample("FiniteElasticity/UniAxialExtension")
buildExample("FiniteElasticity/TriCubicAxialExtension")

# TODO Group them
buildExample("LagrangeSimplexMesh")
buildExample("LinearElasticity")
buildExample("cellml")
buildExample("Darcy")
buildExample("define-geometry-and-export")
buildExample("MoreComplexMesh")
buildExample("simple-field-manipulation-direct-access")
buildExample("SimplexMesh")
buildExample("TwoRegions")

logFile.close()
print "See %slogs_x86_64-linux-%s for detail" % (rootUrl, compiler) 
print "Success Examples: \n%s" % (successedExample) 
if success==0 :
  raise RuntimeError, 'Failed in \n%s' % (failedExample)


