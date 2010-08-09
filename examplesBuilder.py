# -*- coding: utf-8 -*-
import os, sys
from time import strftime

cwd = os.getcwd();
logDir = cwd + "/../../build/logs";
rootUrl = "http://autotest.bioeng.auckland.ac.nz/opencmiss-build/"

if not os.path.isdir(logDir):
  os.mkdir(cwd + "/../../build")
  os.mkdir(logDir);
compiler = sys.argv[1];
f = open(logDir+'/build.log',"a")

def insertTag(tag,filePath,isOpen) :
  if (isOpen) :
    logFile = open(filePath,"w")
  else :
    logFile = open(filePath,"a")
  logFile.write(tag+"\n")
  logFile.close()
  

def buildExample(path) :
   global compiler,logDir,libSuccess;
   if (libSuccess!=-1) :
     newDir = logDir
     for folder in path.split('/') :
       newDir = newDir + '/' + folder
       if not os.path.isdir(newDir):
         os.mkdir(newDir)
     os.chdir(path)
     if os.path.exists(newDir + "/build-" + compiler) :
       os.remove(newDir + "/build-" + compiler)
     insertTag("<pre>",newDir + "/build-" + compiler,True)
     err=os.system("make COMPILER=" + compiler + " >> " + newDir + "/build-" + compiler +" 2>&1")
     insertTag("</pre>",newDir + "/build-" + compiler,False)
     if err==0 :
       f.write(compiler+'_'+path+'_build|success|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
       print "Building %s: <a class='success' href='%slogs_x86_64-linux/%s/build-%s'>success</a><br>" %(path,rootUrl,path,compiler)
     else :
       f.write(compiler+'_'+path+'_build|fail|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
       print "Building %s: <a class='fail' href='%slogs_x86_64-linux/%s/build-%s'>failed</a><br>" %(path,rootUrl,path,compiler)
     os.chdir(cwd)
   else :
     f.write(compiler+'_'+path+'_build|fail|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
     print "Building %s: <a class='fail'>failed</a> due to library build failure<br>" %(path)
   return;

def buildLibrary() :
   global compiler,logDir;
   newDir = logDir
   os.chdir('..')
   if os.path.exists(newDir + "/build-" + compiler) :
     os.remove(newDir + "/build-" + compiler)
   insertTag("<pre>",newDir + "/build-" + compiler,True)
   err=os.system("make USEFIELDML=true COMPILER=" + compiler + " >> " + newDir + "/build-" + compiler +" 2>&1")
   insertTag("</pre>",newDir + "/build-" + compiler,False)
   if err==0 :
     f.write(compiler+'_library_build|success|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
     print "Building OpenCMISS Library: <a class='success' href='%slogs_x86_64-linux/build-%s'>success</a><br>" %(rootUrl,compiler)
     os.chdir(cwd)
     return 0;
   else :
     f.write(compiler+'_library_build|fail|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
     print "Building OpenCMISS Library: <a class='fail' href='%slogs_x86_64-linux/build-%s'>failed</a><br>" %(rootUrl,compiler)
     os.chdir(cwd)
     return -1;
   
   
libSuccess = buildLibrary()

if (compiler!="intel") :
  #buildExample("ClassicalField/AdvectionDiffusion")
  buildExample("ClassicalField/AdvectionDiffusion/AdvectionDiffusionIO")
  buildExample("ClassicalField/AdvectionDiffusion/StaticAdvectionDiffusion")
  buildExample("ClassicalField/AdvectionDiffusion/StaticAdvectionDiffusion_FieldML")
  buildExample("ClassicalField/Diffusion")
  buildExample("ClassicalField/Diffusion/LinearConvergenceTest")
  buildExample("ClassicalField/Diffusion/QuadraticConvergenceTest")
  buildExample("ClassicalField/Diffusion/CubicConvergenceTest") 
  buildExample("ClassicalField/DiffusionConstantSource")
  #buildExample("ClassicalField/Helmholtz")
  buildExample("ClassicalField/NonlinearPoisson")
  buildExample("ClassicalField/Laplace/AnalyticLaplace")
  buildExample("ClassicalField/Laplace/Laplace")
  buildExample("ClassicalField/Laplace/NumberLaplace")

  #buildExample("Bioelectrics/Monodomain")

  buildExample("FluidMechanics/Stokes/ALE")
  buildExample("FluidMechanics/Stokes/Static")
  buildExample("FluidMechanics/Stokes/Dynamic")

  buildExample("FluidMechanics/NavierStokes/ALE")
  buildExample("FluidMechanics/NavierStokes/Static")
  buildExample("FluidMechanics/NavierStokes/Dynamic")
  buildExample("FluidMechanics/NavierStokes/Static_FieldML")


  buildExample("FluidMechanics/Darcy/Analytic")
  buildExample("FluidMechanics/Darcy/Static")
  buildExample("FluidMechanics/Darcy/QuasistaticMaterial")
  #buildExample("FluidMechanics/Darcy/VenousCompartment")

  #buildExample("FiniteElasticity/UniAxialExtension")
  #buildExample("FiniteElasticity/TwoElementTriLinear")
  #buildExample("FiniteElasticity/MixedBoundaryConditions")
  #buildExample("FiniteElasticity/TriCubicAxialExtension")
  buildExample("FiniteElasticity/SimplexElements/LargeQuadraticTet")

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
