import os, subprocess,sys
sys.path.append(os.environ['OPENCMISS_ROOT']+"/cm/examples")
from noseTester import *

masterPath=os.environ['OPENCMISS_ROOT']+"/cm/examples/FluidMechanics/Stokes/42Master"
system = os.uname()[0].lower()
arch = os.uname()[4]
testArgs = "-density 1.0 -viscosity 1.0 -velocity 1.0 0.0 0.0 -dynamic FALSE -starttime 0.0 -stoptime 0.0 -timeincrement 0.0 -ALE FALSE -analytic FALSE -directsolver TRUE"
ndiffDir='expected_results'
outputDir='output'

def test_build():
  if system=="linux" :
    yield check_build, system, "gnu", masterPath
    #yield NoseTester.check_build, "intel"

def test_run():
  if system=="linux" :
    yield check_run, system, arch, "gnu", masterPath, testArgs

def test_output():
  if system=="linux" :
    yield check_output, ndiffDir, outputDir


