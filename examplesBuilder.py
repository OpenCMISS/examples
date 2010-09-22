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

def insertTag(tag,outfile) :
  outfile.write(tag+"\n")
  

def buildExample(path) :
   global compiler,logDir,libSuccess;
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
     htmlWrapper(newDir + "/temp",newDir + "/build-" + compiler)
     if err==0 :
       f.write(compiler+'_'+path+'_build|success|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
       print "<a class='success' href='%slogs_x86_64-linux/%s/build-%s'>Pass</a>: building %s<br>" %(rootUrl,path,compiler,path)
     else :
       f.write(compiler+'_'+path+'_build|fail|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
       print "<a class='fail' href='%slogs_x86_64-linux/%s/build-%s'>Fail</a>: building %s<br>" %(rootUrl,path,compiler,path)
     os.chdir(cwd)
   else :
     f.write(compiler+'_'+path+'_build|fail|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
     print "<a class='fail'>Fail</a>: building %s failed due to library build failure<br>" %(path)
   return;

def buildLibrary() :
   global compiler,logDir;
   newDir = logDir
   os.chdir('..')
   if os.path.exists(newDir + "/temp") :
     os.remove(newDir + "/temp")
   err=os.system("make USEFIELDML=true COMPILER=" + compiler + " > " + newDir + "/temp" +" 2>&1")
   htmlWrapper(newDir + "/temp",newDir + "/build-" + compiler)
   if err==0 :
     f.write(compiler+'_library_build|success|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
     print "<a class='success' href='%slogs_x86_64-linux/build-%s'>Pass</a>: building OpenCMISS Library<br>" %(rootUrl,compiler)
     os.chdir(cwd)
     return 0;
   else :
     f.write(compiler+'_library_build|fail|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
     print "<a class='fail' href='%slogs_x86_64-linux/build-%s'>Fail</a>: building OpenCMISS Library<br>" %(rootUrl,compiler)
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

  buildExample("FluidMechanics/Stokes/RoutineCheck/ALE")
  buildExample("FluidMechanics/Stokes/RoutineCheck/Static")
  buildExample("FluidMechanics/Stokes/RoutineCheck/Dynamic")

  buildExample("FluidMechanics/NavierStokes/RoutineCheck/ALE")
  buildExample("FluidMechanics/NavierStokes/RoutineCheck/Static")
  buildExample("FluidMechanics/NavierStokes/RoutineCheck/Dynamic")
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
