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
compiler = sys.argv[1];

def buildExample(path) :
  global success,compiler,logDir;
  os.chdir(path)
  err=0
  if compiler == 'gnu' :
    err=os.system("make COMPILER=gnu > " + logDir + "/build-" + path.replace('/', '_')+ " 2>&1")
  elif compiler == 'intel' :
    err=os.system("make > " + logDir + "/" + path.replace('/', '_') + " 2>&1")
  if err==0 :
    print "Building " + path.replace('/', '_')+': success'
  else :
    success=0
    print "Building " + path.replace('/', '_')+': fail'
  os.chdir(cwd)
  return;
  
def testExample(path, nodes) :
  global success,compiler,logDir;
  os.putenv('PATH', os.environ['PATH']+':'+cwd+'/../../../opencmissextras/cm/external/x86_64-linux-debug-'+compiler+'/bin')
  os.chdir(path)
  err=os.system('mpiexec -n '+nodes+' bin/x86_64-linux/'+path.rpartition('/')[2]+'Example-debug-'+compiler+' 10 10 0 '+nodes+" > " + logDir + "/test-" + path.replace('/', '_')+ " 2>&1")
  if err==0 :
    print "Testing" + path.replace('/', '_')+': success'
  else :
    success=0
    print "Testing" + path.replace('/', '_')+': fail'
  os.chdir(cwd)
  return;

buildExample("ClassicalField/AnalyticLaplace")
buildExample("ClassicalField/Diffusion")
buildExample("ClassicalField/Helmholtz")
buildExample("ClassicalField/Laplace")
testExample("ClassicalField/Laplace", '2')
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


