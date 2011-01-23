# -*- coding: utf-8 -*-
import os, sys, subprocess
from time import strftime
import socket

cwd = os.getcwd();
logDir = cwd + "/../../build/logs";
ndiff = cwd + "/../utils/ndiff"
rootUrl = "http://autotest.bioeng.auckland.ac.nz/opencmiss-build/"
hostname = socket.gethostname()

if not os.path.isdir(logDir):
  os.mkdir(cwd + "/../../build")
  os.mkdir(logDir);
compiler = sys.argv[1];
os.putenv('HOME', '/home/autotest')
mpidir = cwd+'/../../../opencmissextras/cm/external/x86_64-linux-debug/mpich2/'+compiler+'/bin'
os.system('python ' + mpidir + '/mpd.py &')
f = open(logDir+'/build.log',"r")
successbuilds = f.read()
f.close()
f = open(logDir+'/build.log',"a")

def insertTag(tag,outfile) :
  outfile.write(tag+"\n")

def testExample(id, path, nodes, input=None, args=None, ndiffDir=None,outputDir=None,master42=None) :
   global compiler,logDir,successbuilds,f;
   if master42!=None :
     index = successbuilds.find(compiler+'_'+master42)
     path='42TestingPoints/'+path
   else :
     index = successbuilds.find(compiler+'_'+path)
   index = successbuilds.find('|',index)
   if (successbuilds.startswith('success',index+1)) :
     newDir = logDir
     for folder in path.split('/') :
       newDir = newDir + '/' + folder
       if not os.path.isdir(newDir):
         os.mkdir(newDir)
     os.chdir(path)
     if os.path.exists(newDir + "/temp") :
       os.remove(newDir + "/temp")
     if master42!=None:
       execPath='./run42-debug.sh'
     else :
       execPath='bin/x86_64-linux/mpich2/'+compiler+'/'+path.rpartition('/')[2]+'Example-debug'
     if nodes == '1' :
       if input != None :
         inputPipe = subprocess.Popen(["echo", input], stdout=subprocess.PIPE)
         f1 = open(newDir + "/temp","w")
         execArgs = [execPath]
         if args != None :
           execArgs.extend(args.split(' '))
         execCommand = subprocess.Popen(args=execArgs, stdin=inputPipe.stdout, stdout=f1,stderr=f1)
         f1.close()
         err = os.waitpid(execCommand.pid, 0)[1]
       elif args==None :
         err=os.system(execPath +" > " + newDir + "/temp" + " 2>&1")
       else :
         err=os.system(execPath + ' ' + args +" >> " + newDir + "/temp" + " 2>&1")
     else :
       if input != None :
         inputPipe = subprocess.Popen(["echo", input], stdout=subprocess.PIPE)
         f1 = open(newDir + "/temp","w")
         execArgs = ["mpiexec","-n",nodes,execPath]
         if args != None :
           execArgs.extend(args.split(' '))
         execCommand = subprocess.Popen(args=execArgs, stdin=inputPipe.stdout, stdout=f1,stderr=subprocess.PIPE)
         f1.close()
         err = os.waitpid(execCommand.pid, 0)[1]
       elif args==None :
         err=os.system('python ' + mpidir + '/mpiexec.py -n ' + nodes + ' ' + execPath +" > " + newDir + "/temp" + " 2>&1")
       else :
         err=os.system('python ' + mpidir + '/mpiexec.py -n ' + nodes + " " + execPath + ' ' + args+" > " + newDir + "/temp" + " 2>&1")
     if err==0 and ndiffDir != None and outputDir != None :
     	for outputFile in os.listdir(ndiffDir) :
		if outputFile!='.svn' :
       			error=os.system(ndiff + ' --tolerance=1e-10 ' + ndiffDir +'/'+outputFile + ' ' + outputDir +'/'+outputFile + ' >> '  + newDir + "/temp" + " 2>&1")
			if err==0 :
				err=error
     if not os.path.exists(execPath) :
       err=-1
     htmlWrapper(newDir + "/temp",newDir + "/test" + id + "-" + compiler)
     outputfile = open(newDir + "/test" + id + "-" + compiler, 'r')
     if "ERROR:" in outputfile.read() :
       err=-1
     if os.path.exists(newDir + "/history-test"+id+"-" + compiler) :
       history = open(newDir + "/history-test"+id+"-" + compiler,"a")
     else :
       history = open(newDir + "/history-test"+id+"-" + compiler,"w")
       history.write("<pre>Completed Time\t\tStatus\tHostname\n")
     if err==0 :
       f.write(compiler+'_'+path+'_test|success|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
       print "<a class='success' href='%slogs_x86_64-linux/%s/history-test%s-%s'>Pass</a>: testing <a href='%slogs_x86_64-linux/%s/test%s-%s'>%s%s</a><br>" %(rootUrl,path,id,compiler,rootUrl,path,id,compiler,path,id)
       history.write(strftime("%Y-%m-%d %H:%M:%S")+'\tsuccess\t'+hostname+'\n')
     else :
       f.write(compiler+'_'+path+'_test|fail|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
       print "<a class='fail' href='%slogs_x86_64-linux/%s/history-test%s-%s'>Fail</a>: testing <a href='%slogs_x86_64-linux/%s/test%s-%s'>%s%s</a> <br>" %(rootUrl,path,id,compiler,rootUrl,path,id,compiler,path,id)
       history.write(strftime("%Y-%m-%d %H:%M:%S")+'\tfailed\t'+hostname+'\n')
     os.chdir(cwd)
   else :
     f.write(compiler+'_'+path+'_test|fail|'+ strftime("%Y-%m-%d %H:%M:%S")+'\n')
     print "<a class='fail' href='%slogs_x86_64-linux/%s/history-test%s-%s'>Fail</a>: testing %s%s failed due to build failure<br>" %(rootUrl,path,id,compiler,path,id)
   return;

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

# ##################################################################################
# START: 42 TESTING MATRIX
# ##################################################################################

# Laplace 1D
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/1D/linear/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/1D/quadratic/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/1D/cubic/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/1D/hermite/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
# Laplace 2D
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/2D/QUAD/linear/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/2D/QUAD/quadratic/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/2D/QUAD/cubic/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/2D/TRI/linear/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/2D/TRI/quadratic/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
#testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/2D/TRI/cubic/LEVEL_1",master42="ClassicalField/Laplace/42Master",
#nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/2D/HERMITE/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
# Laplace 3D
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/3D/HEX/linear/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/3D/HEX/quadratic/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/3D/HEX/cubic/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/3D/TET/linear/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/3D/TET/quadratic/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
#testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/3D/TET/cubic/LEVEL_1",master42="ClassicalField/Laplace/42Master",
#nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/3D/HERMITE/LEVEL_1",master42="ClassicalField/Laplace/42Master",
nodes='1',ndiffDir='expected_results',outputDir='output')
# Analytic Laplace 2D
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/2D/QUAD/linear/LEVEL_2",master42="ClassicalField/Laplace/AnalyticLaplace",
nodes='1',ndiffDir='expected_results',outputDir='output')
testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/2D/TRI/linear/LEVEL_2",master42="ClassicalField/Laplace/AnalyticLaplace",
nodes='1',ndiffDir='expected_results',outputDir='output')
#TODO: Hermite is currently not working
#testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Laplace/2D/HERMITE/linear/LEVEL_2",master42="ClassicalField/Laplace/AnalyticLaplace",
#nodes='1',ndiffDir='expected_results',outputDir='output')

#Linear Poisson problem, constant source
#1D
for interp in ['linear','quadratic','cubic','hermite']:
    testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Poisson/1D/"+interp+"/LEVEL_1",master42="ClassicalField/Poisson/42Master",
    nodes='1',ndiffDir='expected_results',outputDir='output')
#2D,3D
for dims in ['2D','3D']:
    if dims == '2D': elements = ['QUAD','TRI']
    elif dims == '3D': elements = ['HEX','TET']
    for els in elements:
        for interp in ['linear','quadratic']:
            testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Poisson/"+dims+"/"+els+"/"+interp+"/LEVEL_1",master42="ClassicalField/Poisson/42Master",
            nodes='1',ndiffDir='expected_results',outputDir='output')
        if els in ['QUAD','HEX']: #cubic broken for tet/tri
            interp = 'cubic'
            testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Poisson/"+dims+"/"+els+"/"+interp+"/LEVEL_1",master42="ClassicalField/Poisson/42Master",
            nodes='1',ndiffDir='expected_results',outputDir='output')
    els = 'HERMITE'
    testExample(id='1',path="SinglePhysics/LinearProblems/StaticProblems/Poisson/"+dims+"/"+els+"/LEVEL_1",master42="ClassicalField/Poisson/42Master",
    nodes='1',ndiffDir='expected_results',outputDir='output')

#Nonlinear Poisson problem, exponential source
#2D only and no Hermite for analytic solution
dims = '2D'
for els in ['QUAD','TRI']:
    for interp in ['linear','quadratic']:
        testExample(id='1',path="SinglePhysics/NonlinearProblems/StaticProblems/NonlinearPoisson/"+dims+"/"+els+"/"+interp+"/LEVEL_1",master42="ClassicalField/Poisson/AnalyticNonlinearPoisson",
        nodes='1',ndiffDir='expected_results',outputDir='output')
    if els == 'QUAD': #cubic broken for tet
        interp = 'cubic'
        testExample(id='1',path="SinglePhysics/NonlinearProblems/StaticProblems/NonlinearPoisson/"+dims+"/"+els+"/"+interp+"/LEVEL_1",master42="ClassicalField/Poisson/AnalyticNonlinearPoisson",
        nodes='1',ndiffDir='expected_results',outputDir='output')

#Linear Stokes problem, static/dynamic, 2D/3D, HEX/TET, 221/331/332
for times in ['StaticProblems','DynamicProblems']:
  for dims in ['2D','3D']:
    if dims == '2D': elements = ['QUAD','TRI']
    elif dims == '3D': elements = ['HEX','TET']
    for els in elements:
        for interp in ['QuadraticVelocityLinearPressure','CubicVelocityLinearPressure','CubicVelocityQuadraticPressure']:
          for level in ['LEVEL_1','LEVEL_2']:
            testExample(id='1',path="SinglePhysics/LinearProblems/"+times+"/Stokes/"+dims+"/"+els+"/"+interp+"/"+level+"",master42="FluidMechanics/Stokes/42Master",nodes='1',ndiffDir='expected_results',outputDir='output')     

#Nonlinear Navier-Stokes problem, static/dynamic, 2D/3D, HEX/TET, 221/331/332
for times in ['StaticProblems','DynamicProblems']:
  for dims in ['2D','3D']:
    if dims == '2D': elements = ['QUAD','TRI']
    elif dims == '3D': elements = ['HEX','TET']
    for els in elements:
        for interp in ['QuadraticVelocityLinearPressure','CubicVelocityLinearPressure','CubicVelocityQuadraticPressure']:
          for level in ['LEVEL_1','LEVEL_2']:
            testExample(id='1',path="SinglePhysics/NonlinearProblems/"+times+"/NavierStokes/"+dims+"/"+els+"/"+interp+"/"+level+"",master42="FluidMechanics/NavierStokes/42Master",nodes='1',ndiffDir='expected_results',outputDir='output')     
            
#Nonlinear FiniteElasticity, Quasistatic, 3D, HEX, 21/31,32
bases=['cubic','quadratic','linear']
for els in ['HEX']:
  for interp1 in range(2):
    for interp2 in range(interp1+1,3):
      for level in ['LEVEL_1','LEVEL_2','LEVEL_3']:
        path_prefix="SinglePhysics/NonlinearProblems/QuasistaticProblems/FiniteElasticity/3D/"
        if level=='LEVEL_3':
          testExample(id='1',path=path_prefix+els+"/"+bases[interp1]+"_"+bases[interp2]+"/"+level+"",master42="FiniteElasticity/testingPoints",nodes='2',ndiffDir='expected_results',outputDir='output')
        else:
          testExample(id='1',path=path_prefix+els+"/"+bases[interp1]+"_"+bases[interp2]+"/"+level+"",master42="FiniteElasticity/testingPoints",nodes='1',ndiffDir='expected_results',outputDir='output')



# ##################################################################################
# END: 42 TESTING MATRIX
# ##################################################################################

#testExample(id='1', path="ClassicalField/AdvectionDiffusion", nodes='1')
testExample(id='1', path="ClassicalField/AdvectionDiffusion/AdvectionDiffusionIO", nodes='1')
testExample(id='1', path="ClassicalField/AdvectionDiffusion/StaticAdvectionDiffusion", nodes='1')
testExample(id='1', path="ClassicalField/AdvectionDiffusion/StaticAdvectionDiffusion_FieldML", nodes='1')
testExample(id='1', path="ClassicalField/Diffusion/Diffusion", nodes='1')
testExample(id='1', path="ClassicalField/Diffusion/LinearConvergenceTest", nodes='1', ndiffDir='expected_files', outputDir='.')
testExample(id='1', path="ClassicalField/Diffusion/QuadraticConvergenceTest", nodes='1', ndiffDir='expected_files', outputDir='.')
testExample(id='1', path="ClassicalField/Diffusion/CubicConvergenceTest", nodes='1', ndiffDir='expected_files', outputDir='.')
testExample(id='1', path="ClassicalField/Diffusion/DiffusionConstantSource", nodes='1')
#testExample(id='1', path="ClassicalField/Helmholtz", nodes='1',input='4\n4\n0\n1')
testExample(id='1', path="ClassicalField/NonlinearPoisson/AnalyticNonlinearPoisson",nodes='1',args='4 4 0 1')
testExample(id='2', path="ClassicalField/Laplace/NumberLaplace", nodes='1',args='4 4 0 1')

#testExample(id='1',path="Bioelectrics/Monodomain",nodes='1',input='4\n4\n0\n1')
  
testExample(id='1',path="FluidMechanics/Stokes/RoutineCheck/ALE",nodes='1',input='\n')
testExample(id='1',path="FluidMechanics/Stokes/RoutineCheck/Static",nodes='1',input='\n')
testExample(id='1',path="FluidMechanics/Stokes/RoutineCheck/Dynamic",nodes='1',input='\n')
testExample(id='1',path="FluidMechanics/NavierStokes/RoutineCheck/ALE",nodes='1',input='\n')
testExample(id='1',path="FluidMechanics/NavierStokes/RoutineCheck/Static",nodes='1',input='\n')
testExample(id='1',path="FluidMechanics/NavierStokes/Static_FieldML",nodes='1',input='\n')
testExample(id='1',path="FluidMechanics/NavierStokes/RoutineCheck/Dynamic",nodes='1',input='\n')
testExample(id='1',path="FluidMechanics/NavierStokes/1DTransient",nodes='1',input='\n')


testExample(id='1',path="FluidMechanics/Darcy/Analytic",nodes='1',input='\n',args='1 test1')
#testExample(id='1',path="FluidMechanics/Darcy/Analytic",nodes='1',input='\n',args='3 test3')
#testExample(id='1',path="FluidMechanics/Darcy/Analytic",nodes='1',input='\n',args='1 test1',ndiffDir="expected_results",outputDir=".")
testExample(id='2',path="FluidMechanics/Darcy/Analytic",nodes='1',input='\n',args='3 test3',ndiffDir="expected_results",outputDir=".")
testExample(id='1',path="FluidMechanics/Darcy/Static",nodes='1',input='\n')
testExample(id='1',path="FluidMechanics/Darcy/QuasistaticMaterial",nodes='1',input='\n',ndiffDir="expected_results",outputDir="output")
#testExample(id='1',path="FluidMechanics/Darcy/QuasistaticMaterial",nodes='1',input='\n')
#testExample(id='1',path="FluidMechanics/Darcy/VenousCompartment",nodes='1',input='4\n4\n0\n1')

testExample(id='1',path="FiniteElasticity/UniAxialExtension",nodes='1',input='\n',ndiffDir="expected_results",outputDir=".")
#testExample(id='1',path="FiniteElasticity/TwoElementTriLinear",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="FiniteElasticity/MixedBoundaryConditions",nodes='1',input='4\n4\n0\n1')
testExample(id='1',path="FiniteElasticity/TriCubicAxialExtension",nodes='1',input='\n',ndiffDir="expected_results",outputDir=".")
testExample(id='1',path="FiniteElasticity/LargeUniAxialExtension",nodes='1',input='\n',ndiffDir="expected_results",outputDir=".")
testExample(id='1',path="FiniteElasticity/CylinderInflation",nodes='1',input='\n',ndiffDir="expected_results/1proc",outputDir="./outputs")
testExample(id='2',path="FiniteElasticity/CylinderInflation",nodes='2',ndiffDir="expected_results/2proc",outputDir="./outputs")
testExample(id='1',path="FiniteElasticity/SimplexElements/LargeQuadraticTet",nodes='1',input='\n')

#testExample(id='1',path="LinearElasticity/2DAnalytic1",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="LinearElasticity/2DPlaneStressLagrangeBasis",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="LinearElasticity/3DAnalytic1",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="LinearElasticity/3DLagrangeBasis",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="LinearElasticity/3DCubicHermiteBasis",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="LinearElasticity/3DLagrangeBasisAnisotropicFibreField",nodes='1',input='4\n4\n0\n1')

testExample(id='1',path="MultiPhysics/Poroelasticity/FiniteElasticityDarcy/SameRegionSameMesh",nodes='1',input='\n')
testExample(id='1',path="MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDarcySharedVars",nodes='1',input='\n')
testExample(id='1',path="MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDarcySharedVarsTimeBCs",nodes='1',input='\n')


# TODO Group them
#testExample(id='1',path="LagrangeSimplexMesh",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="cellml",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="define-geometry-and-export",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="MoreComplexMesh",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="simple-field-manipulation-direct-access",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="SimplexMesh",nodes='1',input='4\n4\n0\n1')
#testExample(id='1',path="TwoRegions",nodes='1',input='4\n4\n0\n1')

f.close()

os.system('python ' + mpidir + '/mpdallexit.py')

