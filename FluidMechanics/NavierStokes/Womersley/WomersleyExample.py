#!/usr/bin/env python

#> \file
#> \author David Ladd
#> \brief This is an OpenCMISS script to solve Navier-Stokes oscillatory flow through a cylinder and validate the solution against Womersley's analytic solution for the velocity profile. 
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is openCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s):
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. if you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. if you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>

#> \example FluidMechanics/NavierStokes/Womersley/WomersleyExample.py
## Python OpenCMISS script to solve Navier-Stokes oscillatory flow through a cylinder and validate the solution against Womersley's analytic solution for the velocity profile.
## \par Latest Builds:
#<


# Add Python bindings directory to PATH
import sys, os
sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

import numpy
import gzip
import pylab
import time
import re
import math
import matplotlib.pyplot as plt
import contextlib
import WomersleyAnalytic

# Intialise OpenCMISS
from opencmiss import CMISS


@contextlib.contextmanager
def ChangeDirectory(path):
    """A context manager which changes the working directory to the given
    path, and then changes it back to its previous value on exit.
    """
    prev_cwd = os.getcwd()
    os.chdir(path)
    try:
        yield path
    finally:
        os.chdir(prev_cwd)

# Get the computational nodes information
numberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet()
computationalNodeNumber = CMISS.ComputationalNodeNumberGet()

# -----------------------------------------------
#  Get the mesh information from FieldML data
# -----------------------------------------------

meshName = 'hexCylinder14'            
inputDir = './input/hexCylinder14/'
length = 10.
radius = 0.5
axialComponent = 1
fieldmlInput = inputDir + meshName + '.xml'
print('FieldML input file: ' + fieldmlInput)

#Wall boundary nodes
filename=inputDir + 'bc/wallNodes.dat'
try:
    with open(filename):
        f = open(filename,"r")
        numberOfWallNodes=int(f.readline())
        wallNodes=map(int,(re.split(',',f.read())))
        f.close()
except IOError:
   print ('Could not open Wall boundary node file: ' + filename)
#Inlet boundary nodes
filename=inputDir + 'bc/inletVelocityNodes.dat'
try:
    with open(filename):
        f = open(filename,"r")
        numberOfInletNodes=int(f.readline())
        inletNodes=map(int,(re.split(',',f.read())))
        f.close()
except IOError:
   print ('Could not open Inlet boundary node file: ' + filename)
#Outlet boundary nodes
filename=inputDir + 'bc/outletNodes.dat'
try:
    with open(filename):
        f = open(filename,"r")
        numberOfOutletNodes=int(f.readline())
        outletNodes=map(int,(re.split(',',f.read())))
        f.close()
except IOError:
   print ('Could not open Outlet boundary node file: ' + filename)

print("Inlet Nodes:")
print(inletNodes)

print("outlet Nodes:")
print(outletNodes)

# -----------------------------------------------
#  Set up problem
# -----------------------------------------------

(coordinateSystemUserNumber,
 regionUserNumber,
 linearBasisUserNumber,
 quadraticBasisUserNumber,
 generatedMeshUserNumber,
 meshUserNumber,
 decompositionUserNumber,
 geometricFieldUserNumber,
 equationsSetFieldUserNumber,
 dependentFieldUserNumber,
 materialsFieldUserNumber,
 analyticFieldUserNumber,
 equationsSetUserNumber,
 problemUserNumber) = range(1,15)

#Initialise fieldML IO
fieldmlInfo=CMISS.FieldMLIO()
fieldmlInfo.InputCreateFromFile(fieldmlInput)

# Creation a RC coordinate system
coordinateSystem = CMISS.CoordinateSystem()
fieldmlInfo.InputCoordinateSystemCreateStart("CylinderMesh.coordinates",coordinateSystem,coordinateSystemUserNumber)
#fieldmlInfo.InputCoordinateSystemCreateStart("CylinderMesh.coordinates",coordinateSystemUserNumber)
coordinateSystem.CreateFinish()
numberOfDimensions = coordinateSystem.DimensionGet()

# Create a region
region = CMISS.Region()
region.CreateStart(regionUserNumber,CMISS.WorldRegion)
region.label = "Cylinder"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create nodes
nodes=CMISS.Nodes()
fieldmlInfo.InputNodesCreateStart("CylinderMesh.nodes.argument",region,nodes)
nodes.CreateFinish()
numberOfNodes = nodes.numberOfNodes
print("number of nodes: " + str(numberOfNodes))

# Create bases
basisNumberQuadratic = 1
basisNumberLinear = 2
meshComponentQuadratic = 1
meshComponentLinear = 2
gaussQuadrature = [3,3,3]
fieldmlInfo.InputBasisCreateStartNum("CylinderMesh.triquadratic_lagrange",basisNumberQuadratic)
CMISS.Basis_QuadratureNumberOfGaussXiSetNum(basisNumberQuadratic,gaussQuadrature)
CMISS.Basis_QuadratureLocalFaceGaussEvaluateSetNum(basisNumberQuadratic,True)
CMISS.Basis_CreateFinishNum(basisNumberQuadratic)
fieldmlInfo.InputBasisCreateStartNum("CylinderMesh.trilinear_lagrange",basisNumberLinear)
CMISS.Basis_QuadratureNumberOfGaussXiSetNum(basisNumberLinear,gaussQuadrature)
CMISS.Basis_QuadratureLocalFaceGaussEvaluateSetNum(basisNumberLinear,True)
CMISS.Basis_CreateFinishNum(basisNumberLinear)

# Create Mesh
numberOfMeshComponents=2
meshComponentQuadratic=1
meshComponentLinear=2
mesh = CMISS.Mesh()
fieldmlInfo.InputMeshCreateStart("CylinderMesh.mesh.argument",mesh,meshUserNumber,region)
mesh.NumberOfComponentsSet(numberOfMeshComponents)
fieldmlInfo.InputCreateMeshComponent(mesh,meshComponentQuadratic,"CylinderMesh.template.triquadratic")
fieldmlInfo.InputCreateMeshComponent(mesh,meshComponentLinear,"CylinderMesh.template.trilinear")
mesh.CreateFinish()
numberOfElements = mesh.numberOfElements
print("number of elements: " + str(numberOfElements))

# Create a decomposition for the mesh
decomposition = CMISS.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = CMISS.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CalculateFacesSet(True)
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = CMISS.Field()
fieldmlInfo.InputFieldCreateStart(region,decomposition,geometricFieldUserNumber,
                                  geometricField,CMISS.FieldVariableTypes.U,
                                  "CylinderMesh.coordinates")
geometricField.CreateFinish()
fieldmlInfo.InputFieldParametersUpdate(geometricField,"CylinderMesh.node.coordinates",
                                       CMISS.FieldVariableTypes.U,
                                       CMISS.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateStart(CMISS.FieldVariableTypes.U,
                                       CMISS.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(CMISS.FieldVariableTypes.U,
                                       CMISS.FieldParameterSetTypes.VALUES)
fieldmlInfo.Finalise()

# Export mesh geometry
fields = CMISS.Fields()
fields.CreateRegion(region)
fields.NodesExport("Geometry","FORTRAN")
fields.ElementsExport("Geometry","FORTRAN")
fields.Finalise()
print("Exported Geometric Mesh")

# Create standard Navier-Stokes equations set
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        CMISS.EquationsSetClasses.FLUID_MECHANICS,
        CMISS.EquationsSetTypes.NAVIER_STOKES_EQUATION,
        CMISS.EquationsSetSubtypes.TRANSIENT_SUPG_NAVIER_STOKES,
        equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()


def solveProblem(transient,viscosity,density,offset,amplitude,period):
    """ Sets up the problem and solve with the provided parameter values


          Oscillatory flow through a rigid cylinder

                                     u=0
                  ------------------------------------------- R = 0.5
                                             \
                                              \  
        p = offset + A*sin(2*pi*(t/period))   | u(r,t)        p = 0
                                              /
                                             /
                  ------------------------------------------- L = 10
                                     u=0
    """

    period = math.pi/2.0
    angularFrequency = 2.0*math.pi/period
    Womersley = 0.5*math.sqrt(angularFrequency*density/viscosity)
    if computationalNodeNumber == 0:
        print("-----------------------------------------------")
        print("Setting up problem for Womersley number: " + str(Womersley))
        print("-----------------------------------------------")

    # Create dependent field
    dependentField = CMISS.Field()
    equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
    # velocity
    for component in range(1,4):
        dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,component,meshComponentQuadratic)        
        dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN,component,meshComponentQuadratic) 
    # pressure
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,4,meshComponentLinear)        
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN,4,meshComponentLinear) 
    dependentField.DOFOrderTypeSet(CMISS.FieldVariableTypes.U,CMISS.FieldDOFOrderTypes.SEPARATED)
    dependentField.DOFOrderTypeSet(CMISS.FieldVariableTypes.DELUDELN,CMISS.FieldDOFOrderTypes.SEPARATED)
    equationsSet.DependentCreateFinish()
    # Initialise dependent field to 0
    for component in range(1,5):
        dependentField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,component,0.0)
        dependentField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.DELUDELN,CMISS.FieldParameterSetTypes.VALUES,component,0.0)

    # Initialise dependent field to analytic values
    initialiseAnalytic = False
    if initialiseAnalytic:
        for node in range(1,numberOfNodes+1):
            sumPositionSq = 0.
            nodeNumber = nodes.UserNumberGet(node)
            nodeDomain=decomposition.NodeDomainGet(nodeNumber,meshComponentQuadratic)
            if (nodeDomain == computationalNodeNumber):
                for component in range(1,4):
                    if component != axialComponent:
                        value=geometricField.ParameterSetGetNodeDP(CMISS.FieldVariableTypes.U,
                                                                      CMISS.FieldParameterSetTypes.VALUES,
                                                                      1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,component)
                        sumPositionSq += value**2
                radialNodePosition=math.sqrt(sumPositionSq)
                for component in range(1,4):
                    if component == axialComponent + 1:
                        value = WomersleyAnalytic.WomersleyAxialVelocity(transient[0],offset,amplitude,radius,
                                                                         radialNodePosition,period,viscosity,
                                                                         Womersley,length)
                    else:
                        value = 0.0
                    dependentField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,
                                                            1,1,nodeNumber,component,value)

    # Create materials field
    materialsField = CMISS.Field()
    equationsSet.MaterialsCreateStart(materialsFieldUserNumber,materialsField)
    equationsSet.MaterialsCreateFinish()
    # Initialise materials field parameters
    materialsField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,viscosity)
    materialsField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,2,density)

    # Create analytic field (allows for time-dependent calculation of sinusoid-based pressure waveform during solve)
    analytic = False
    if analytic:
        analyticField = CMISS.Field()
        equationsSet.AnalyticCreateStart(CMISS.NavierStokesAnalyticFunctionTypes.Sinusoid,analyticFieldUserNumber,analyticField)
        equationsSet.AnalyticCreateFinish()
        # Initialise analytic field parameters: (1-4) Dependent params, 5 amplitude, 6 offset, 7 period
        analyticField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,0.0)
        analyticField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,2,0.0)
        analyticField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,3,0.0)
        analyticField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,4,1.0)
        analyticField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,5,amplitude)
        analyticField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,6,offset)
        analyticField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,7,period)

    # Create equations
    equations = CMISS.Equations()
    equationsSet.EquationsCreateStart(equations)
    equations.sparsityType = CMISS.EquationsSparsityTypes.SPARSE
    equations.outputType = CMISS.EquationsOutputTypes.NONE
    equationsSet.EquationsCreateFinish()

    # Create Navier-Stokes problem
    problem = CMISS.Problem()
    problem.CreateStart(problemUserNumber)
    problem.SpecificationSet(CMISS.ProblemClasses.FLUID_MECHANICS,
                             CMISS.ProblemTypes.NAVIER_STOKES_EQUATION,
                             CMISS.ProblemSubTypes.TRANSIENT_SUPG_NAVIER_STOKES)
    problem.CreateFinish()

    # Create control loops
    problem.ControlLoopCreateStart()
    controlLoop = CMISS.ControlLoop()
    problem.ControlLoopGet([CMISS.ControlLoopIdentifiers.NODE],controlLoop)
    controlLoop.TimesSet(transient[0],transient[1],transient[2])
    controlLoop.TimeOutputSet(transient[3])
    problem.ControlLoopCreateFinish()

    # Create problem solver
    dynamicSolver = CMISS.Solver()
    problem.SolversCreateStart()
    problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],1,dynamicSolver)
    dynamicSolver.outputType = CMISS.SolverOutputTypes.NONE
    dynamicSolver.dynamicTheta = [1.0]
    nonlinearSolver = CMISS.Solver()
    dynamicSolver.DynamicNonlinearSolverGet(nonlinearSolver)
    nonlinearSolver.newtonJacobianCalculationType = CMISS.JacobianCalculationTypes.EQUATIONS
    nonlinearSolver.outputType = CMISS.SolverOutputTypes.PROGRESS
    nonlinearSolver.newtonAbsoluteTolerance = 1.0E-7
    nonlinearSolver.newtonRelativeTolerance = 1.0E-7
    nonlinearSolver.newtonSolutionTolerance = 1.0E-7
    nonlinearSolver.newtonMaximumFunctionEvaluations = 10000
    linearSolver = CMISS.Solver()
    nonlinearSolver.NewtonLinearSolverGet(linearSolver)
    linearSolver.outputType = CMISS.SolverOutputTypes.NONE
    linearSolver.linearType = CMISS.LinearSolverTypes.DIRECT
    linearSolver.libraryType = CMISS.SolverLibraries.MUMPS
    problem.SolversCreateFinish()

    # Create solver equations and add equations set to solver equations
    solver = CMISS.Solver()
    solverEquations = CMISS.SolverEquations()
    problem.SolverEquationsCreateStart()
    problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],1,solver)
    solver.SolverEquationsGet(solverEquations)
    solverEquations.sparsityType = CMISS.SolverEquationsSparsityTypes.SPARSE
    equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
    problem.SolverEquationsCreateFinish()

    # Create boundary conditions
    boundaryConditions = CMISS.BoundaryConditions()
    solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
    # Wall boundary nodes u = 0 (no-slip)
    value=0.0
    for nodeNumber in wallNodes:
        nodeDomain=decomposition.NodeDomainGet(nodeNumber,meshComponentQuadratic)
        if (nodeDomain == computationalNodeNumber):
            for component in range(numberOfDimensions):
                componentId = component + 1
                boundaryConditions.SetNode(dependentField,CMISS.FieldVariableTypes.U,
                                           1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                           nodeNumber,componentId,CMISS.BoundaryConditionsTypes.FIXED,value)


    # Note: inlet/outlet nodes are pressure-based so only defined on linear nodes
    # outlet boundary nodes p = 0 
    value=0.0
    for nodeNumber in outletNodes:
        nodeDomain=decomposition.NodeDomainGet(nodeNumber,meshComponentLinear)
        if (nodeDomain == computationalNodeNumber):
            boundaryConditions.SetNode(dependentField,CMISS.FieldVariableTypes.U,
                                       1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                       nodeNumber,4,CMISS.BoundaryConditionsTypes.FIXED,value)
    # inlet boundary nodes p = f(t) - will be updated in pre-solve
    # value = 1.0
    # for nodeNumber in inletNodes:
    #     nodeDomain=decomposition.NodeDomainGet(nodeNumber,meshComponentQuadratic)
    #     if (nodeDomain == computationalNodeNumber):
    #         boundaryConditions.SetNode(dependentField,CMISS.FieldVariableTypes.U,
    #                                    1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
    #                                    nodeNumber,2,CMISS.BoundaryConditionsTypes.FIXED_INLET,value)
    # solverEquations.BoundaryConditionsCreateFinish()

    value = 1.0
    for nodeNumber in inletNodes:
        nodeDomain=decomposition.NodeDomainGet(nodeNumber,meshComponentQuadratic)
        if (nodeDomain == computationalNodeNumber):
            for component in range (1,4):
                if component == 2:
                    value = 1.0
                else:
                    value = 0.0
                boundaryConditions.SetNode(dependentField,CMISS.FieldVariableTypes.U,
                                           1,CMISS.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                           nodeNumber,component,CMISS.BoundaryConditionsTypes.FIXED_INLET,value)
    solverEquations.BoundaryConditionsCreateFinish()

    # Solve the problem
    print("solving problem...")
    problem.Solve()

    # Clear fields so can run in batch mode on this region
    materialsField.Destroy()
    analyticField.Destroy()
    equationsSet.Destroy()
    problem.Destroy()



#==========================================================
# P r o b l e m     C o n t r o l
#==========================================================

# Problem defaults
offset = 0.0
density = 1.0
amplitude = 1.0
period = math.pi/2

# transient parameters: startTime,stopTime,timeIncrement,outputFrequency
#transient = [3.0/4.0*period,7.0/4.0*period+0.00001,period/40,1]
#transient = [0.0,period+0.00001,period/10,1]
transient = [0.0,2.0001,1.0,1]

WomersleyNumbers = [1.0]

for Wo in WomersleyNumbers:
    # determine Wo using viscosity:density ratio (fixed angular frequency and radius)
    viscosity = density/(Wo**2.0)
    # make a new output directory if necessary
    outputDirectory = "./output/Wo" + str(Wo) + meshName + "/"
    try:
        os.makedirs(outputDirectory)
    except OSError, e:
        if e.errno != 17:
            raise   
    # change to new directory and solve problem (note will return to original directory on exit)
    with ChangeDirectory(outputDirectory):
        solveProblem(transient,viscosity,density,offset,amplitude,period)

CMISS.Finalise()





