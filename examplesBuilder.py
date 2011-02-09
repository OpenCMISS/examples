# -*- coding: utf-8 -*-
import os, sys, commands
from time import strftime
import socket

cwd = os.getcwd();
logDir = cwd + "/../../build/logs";
rootUrl = "http://autotest.bioeng.auckland.ac.nz/opencmiss-build/"
hostname = socket.gethostname()

def getVersion() :
  global compiler;
  operating_system = commands.getoutput('uname') 
  if operating_system == 'Linux':
    if compiler == 'gnu':
      version_info = commands.getoutput('gfortran -v')
      if version_info.find('gcc version 4.5') != -1 :
        return '_4.5'
      elif version_info.find('gcc version 4.4') != -1 :
        return '_4.4'
  return ''

if not os.path.isdir(logDir):
  os.mkdir(cwd + "/../../build")
  os.mkdir(logDir);
compiler = sys.argv[1];
compiler_version = getVersion()
f = open(logDir+'/build.log',"a")

def insertTag(tag,outfile) :
  outfile.write(tag+"\n")
  

def buildExample(path) :
   global compiler,compiler_version,logDir,libSuccess;
   if (libSuccess!=-1) :
     newDir = logDir
     for folder in path.split('/') :
       newDir = newDir + '/' + folder
       if not os.path.isdir(newDir):
         os.mkdir(newDir)
     os.chdir(path)
     if os.path.exists(newDir + "/temp") :
       os.remove(newDir + "/temp")
     err=os.system("make COMPILER=" + compiler + " > " + newDir + "/temp" +" 2>&1")
     htmlWrapper(newDir + "/temp",newDir + "/build-" + compiler + compiler_version)
     if os.path.exists(newDir + "/history-build-" + compiler + compiler_version) :
       history = open(newDir + "/history-build-" + compiler + compiler_version,"a")
     else :
       history = open(newDir + "/history-build-" + compiler + compiler_version,"w")
       history.write("<pre>Completed Time\t\tStatus\tHostname\n")
     if err==0 :
       f.write(compiler+compiler_version+'_'+path+'_build|success|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
       print "<a class='success' href='%slogs_x86_64-linux/%s/history-build-%s'>Pass</a>: building <a href='%slogs_x86_64-linux/%s/build-%s'>%s</a><br>" %(rootUrl,path,compiler+compiler_version,rootUrl,path,compiler+compiler_version,path)
       history.write(strftime("%Y-%m-%d %H:%M:%S")+'\tsuccess\t'+hostname+'\n')
     else :
       f.write(compiler+compiler_version+'_'+path+'_build|fail|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
       print "<a class='fail' href='%slogs_x86_64-linux/%s/history-build-%s'>Fail</a>: building <a href='%slogs_x86_64-linux/%s/build-%s'>%s</a><br>" %(rootUrl,path,compiler+compiler_version,rootUrl,path,compiler+compiler_version,path)
       history.write(strftime("%Y-%m-%d %H:%M:%S")+'\tfailed\t'+hostname+'\n')
     os.chdir(cwd)
   else :
     f.write(compiler+compiler_version+'_'+path+'_build|fail|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
     print "<a class='fail' href='%slogs_x86_64-linux/%s/history-build-%s'>Fail</a>: building %s failed due to library build failure<br>" %(rootUrl,path,compiler+compiler_version,path)
   return;

def buildLibrary() :
   global compiler,compiler_version,logDir,hostname;
   newDir = logDir
   os.chdir('..')
   if os.path.exists(newDir + "/temp") :
     os.remove(newDir + "/temp")
   err=os.system("make USEFIELDML=true COMPILER=" + compiler + " > " + newDir + "/temp" +" 2>&1")
   htmlWrapper(newDir + "/temp",newDir + "/build-" + compiler+compiler_version)
   if os.path.exists(newDir + "/history-" + compiler+compiler_version) :
     history = open(newDir + "/history-" + compiler+compiler_version,"a")
   else :
     history = open(newDir + "/history-" + compiler+compiler_version,"w")
     history.write("<pre>Completed Time\t\tStatus\tHostname\n")
   if err==0 :
     f.write(compiler+compiler_version+'_library_build|success|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
     print "<a class='success' href='%slogs_x86_64-linux/history-%s'>Pass</a>: building <a href='%slogs_x86_64-linux/build-%s'>OpenCMISS Library</a><br>" %(rootUrl,compiler+compiler_version,rootUrl,compiler+compiler_version)
     history.write(strftime("%Y-%m-%d %H:%M:%S")+'\tsuccess\t'+hostname+'\n')
     os.chdir(cwd)
     return 0;
   else :
     f.write(compiler+compiler_version+'_library_build|fail|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
     print "<a class='fail' href='%slogs_x86_64-linux/history-%s'>Fail</a>: building <a href='%slogs_x86_64-linux/build-%s'>OpenCMISS Library</a><br>" %(rootUrl,compiler+compiler_version,rootUrl,compiler+compiler_version)
     history.write(strftime("%Y-%m-%d %H:%M:%S")+'\tfailed\t'+hostname+'\n')
     os.chdir(cwd)
     return -1;

def htmlWrapper(infile, outfile) :
  stext1 = "<";
  stext2 = ">";
  rtext1 = "&lt;";
  rtext2 = "&gt;";
  input = sys.stdin
  output = sys.stdout
  input = open(infile)
  output = open(outfile, 'w')
  insertTag("<pre>",output)
  for s in input.xreadlines():
    output.write(s.replace(stext1, rtext1).replace(stext2,rtext2))
  insertTag("</pre>",output)
  output.close()
   
   
libSuccess = buildLibrary()

if (compiler!="intel") :
  ##buildExample("FluidMechanics/Stokes/42Master & FluidMechanics/NavierStokes/42Master")
  buildExample("ClassicalField/Laplace/42Master")
  buildExample("FluidMechanics/NavierStokes/42Master")
  buildExample("FluidMechanics/Stokes/42Master")
  
  #buildExample("ClassicalField/AdvectionDiffusion")
  buildExample("ClassicalField/AdvectionDiffusion/AdvectionDiffusionIO")
  buildExample("ClassicalField/AdvectionDiffusion/StaticAdvectionDiffusion")
  buildExample("ClassicalField/AdvectionDiffusion/StaticAdvectionDiffusion_FieldML")
  buildExample("ClassicalField/Diffusion/Diffusion")
  buildExample("ClassicalField/Diffusion/LinearConvergenceTest")
  buildExample("ClassicalField/Diffusion/QuadraticConvergenceTest")
  buildExample("ClassicalField/Diffusion/CubicConvergenceTest") 
  buildExample("ClassicalField/Diffusion/DiffusionConstantSource")
  #buildExample("ClassicalField/Helmholtz")
  buildExample("ClassicalField/Poisson/AnalyticNonlinearPoisson")
  buildExample("ClassicalField/Poisson/42Master")
  buildExample("ClassicalField/Laplace/AnalyticLaplace")
  buildExample("ClassicalField/Laplace/Laplace")
  buildExample("ClassicalField/Laplace/NumberLaplace")

  #buildExample("Bioelectrics/Monodomain")

  buildExample("FluidMechanics/Stokes/RoutineCheck/ALE")
  buildExample("FluidMechanics/Stokes/RoutineCheck/Static")
  buildExample("FluidMechanics/Stokes/RoutineCheck/Dynamic")

  buildExample("FluidMechanics/NavierStokes/RoutineCheck/ALE")
  buildExample("FluidMechanics/NavierStokes/RoutineCheck/Static")
  buildExample("FluidMechanics/NavierStokes/RoutineCheck/Dynamic")
  buildExample("FluidMechanics/NavierStokes/Static_FieldML")
  buildExample("FluidMechanics/NavierStokes/1DTransient")

  buildExample("FluidMechanics/Darcy/Analytic")
  buildExample("FluidMechanics/Darcy/Static")
  buildExample("FluidMechanics/Darcy/QuasistaticMaterial")
  #buildExample("FluidMechanics/Darcy/VenousCompartment")

  buildExample("FiniteElasticity/UniAxialExtension")
  #buildExample("FiniteElasticity/TwoElementTriLinear")
  #buildExample("FiniteElasticity/MixedBoundaryConditions")
  buildExample("FiniteElasticity/TriCubicAxialExtension")
  buildExample("FiniteElasticity/CylinderInflation")
  buildExample("FiniteElasticity/LargeUniAxialExtension")
  buildExample("FiniteElasticity/SimplexElements/LargeQuadraticTet")
  buildExample("FiniteElasticity/testingPoints")

  buildExample("LinearElasticity/2DExtensionPlaneStressLagrangeBasis")
  buildExample("LinearElasticity/3DExtensionLagrangeBasis")
  buildExample("LinearElasticity/3DExtensionHermiteBasis")
  buildExample("LinearElasticity/Analytic/Extension")

  buildExample("MultiPhysics/Poroelasticity/FiniteElasticityDarcy/SameRegionSameMesh")
  buildExample("MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDarcySharedVars")
  buildExample("MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDarcySharedVarsTimeBCs")

  # TODO Group them
  #buildExample("LagrangeSimplexMesh")
  #buildExample("cellml")
  #buildExample("define-geometry-and-export")
  #buildExample("MoreComplexMesh")
  #buildExample("simple-field-manipulation-direct-access")
  #buildExample("SimplexMesh")
  #buildExample("TwoRegions")

f.close()

