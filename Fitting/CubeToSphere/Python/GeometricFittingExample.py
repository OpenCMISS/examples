#!/usr/bin/env python

#> \file
#> \author David Ladd
#> \brief This is an example to use linear fitting to fit a generated cube mesh surface to a sphere.
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
#> The Original Code is OpenCMISS
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

#> \example /Fitting/GeometricFitting/GeometricFittingExample.py
## Example script to fit a generated cube mesh to a sphere using OpenCMISS calls in python.
## \par Latest Builds:
#<

# Add Python bindings directory to PATH
import sys, os
sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

import exfile
import numpy
from numpy import linalg
import math
import random

# Intialise OpenCMISS
from opencmiss import iron

def writeExdataFile(filename,dataPointLocations,offset):
    "Writes data points to an exdata file"

    numberOfDimensions = dataPointLocations[1].shape[0]
    try:
        f = open(filename,"w")    
        header = '''Group name: DataPoints
 #Fields=1
 1) data_coordinates, coordinate, rectangular cartesian, #Components='''+str(numberOfDimensions)+'''
  1.  Value index=1, #Derivatives=0, #Versions=1
  2.  Value index=2, #Derivatives=0, #Versions=1
'''
        if numberOfDimensions == 3:
            header+= '''  3.  Value index=3, #Derivatives=0, #Versions=1
'''
        f.write(header)

        numberOfDataPoints = len(dataPointLocations)
        for i in range(numberOfDataPoints):
            line = " Node: " + str(offset+i+1) + '\n'
            f.write(line)
            for j in range (numberOfDimensions):
                line = ' ' + str(dataPointLocations[i,j]) + '\n'
                f.write(line)
        f.close()
            
    except IOError:
        print ('Could not open file: ' + filename)


#=================================================================
# Control Panel
#=================================================================

#-----------------------------------------------------------------------------------------------
# ***NOTE: running this in debug mode with FPE detection can incur large performance overheads.
#          Suggest using an optimized version of the library if possible.
#-----------------------------------------------------------------------------------------------

# Set cube dimensions
numberOfDimensions = 3
length = 1.25
meshDimensions = []
meshOrigin = []
for d in range(numberOfDimensions):
    meshDimensions.append(length)
    meshOrigin.append(-length/2.0)

# Set cube generated mesh resolution/type
meshResolution = [3]*numberOfDimensions
hermite = False

# Set data point resolution (will be randomly placed on surface of a sphere)
numberOfDataPoints = 1000
radius = 1.0
origin = [0.,0.,0.]

# fix interior nodes so that fitting only applies to surface
fixInterior = True

# iteratively fit the cube to sphere- default 1 for automated testing
numberOfIterations = 1

# If start iteration > 1, read in geometry from a previous fit iteration
iteration = 1
if iteration > 1:
    exfileMesh = True
    exnode = exfile.Exnode("DeformedGeometry" + str(iteration-1) + ".part0.exnode")
    exelem = exfile.Exelem("UndeformedGeometry.part0.exelem")
else:
    exfileMesh = False

# Set Sobolev smoothing parameters
tau = 0.5
kappa = 0.1

numberOfGaussXi = numberOfDimensions
zeroTolerance = 0.00001

#=================================================================

(coordinateSystemUserNumber,
    regionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    meshUserNumber,
    decompositionUserNumber,
    geometricFieldUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    independentFieldUserNumber,
    dataPointFieldUserNumber,
    materialFieldUserNumber,
    analyticFieldUserNumber,
    dependentDataFieldUserNumber,
    dataProjectionUserNumber,
    equationsSetUserNumber,
    problemUserNumber) = range(1,18)

# Get the computational nodes information
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Create a RC coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = numberOfDimensions
coordinateSystem.CreateFinish()

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "FittingRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

#=================================================================
# Mesh
#=================================================================
# Calc number of elements
numberOfElements = 1
for dimension in range(numberOfDimensions):
    numberOfElements = numberOfElements*meshResolution[dimension]

# Create a lagrange basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = numberOfDimensions
if hermite:
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*numberOfDimensions
else:
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfDimensions

basis.quadratureNumberOfGaussXi = [numberOfGaussXi]*numberOfDimensions
basis.CreateFinish()

if (exfileMesh):
    # Read previous mesh
    mesh = iron.Mesh()
    mesh.CreateStart(meshUserNumber, region, numberOfDimensions)
    mesh.NumberOfComponentsSet(1)
    mesh.NumberOfElementsSet(exelem.num_elements)
    # Define nodes for the mesh
    nodes = iron.Nodes()
    nodes.CreateStart(region, exnode.num_nodes)
    nodes.CreateFinish()
    # Define elements for the mesh
    elements = iron.MeshElements()
    meshComponentNumber = 1
    elements.CreateStart(mesh, meshComponentNumber, basis)
    for elem in exelem.elements:
        elements.NodesSet(elem.number, elem.nodes)
    elements.CreateFinish()
    mesh.CreateFinish()
else:
    # Create a generated mesh
    generatedMesh = iron.GeneratedMesh()
    generatedMesh.CreateStart(generatedMeshUserNumber,region)
    generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
    generatedMesh.basis = [basis]
    generatedMesh.extent = meshDimensions
    generatedMesh.origin = meshOrigin

    generatedMesh.numberOfElements = meshResolution
    mesh = iron.Mesh()
    generatedMesh.CreateFinish(meshUserNumber,mesh)

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

#=================================================================
# Geometric Field
#=================================================================

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.meshDecomposition = decomposition
for dimension in range(numberOfDimensions):
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,dimension+1,1)
geometricField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
geometricField.CreateFinish()

# Get nodes
nodes = iron.Nodes()
region.NodesGet(nodes)
numberOfNodes = nodes.numberOfNodes

# Get or calculate geometric parameters
if (exfileMesh):
    # Read the geometric field from the exnode file
    geometricField.ParameterSetUpdateStart(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    for node_num in range(1, exnode.num_nodes + 1):
        version = 1
        derivative = 1
        for component in range(1, numberOfDimensions + 1):
            component_name = ["x", "y", "z"][component - 1]
            value = exnode.node_value("Coordinate", component_name, node_num, derivative)
            geometricField.ParameterSetUpdateNode(
                    iron.FieldVariableTypes.U,
                    iron.FieldParameterSetTypes.VALUES,
                    version, derivative, node_num, component, value)
    geometricField.ParameterSetUpdateFinish(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
else:
    # Create undeformed geometry from the generated mesh
    generatedMesh.GeometricParametersCalculate(geometricField)
    # Export undeformed mesh geometry
    print("Writing undeformed geometry")
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("UndeformedGeometry","FORTRAN")
    fields.ElementsExport("UndeformedGeometry","FORTRAN")
    fields.Finalise()

#=================================================================
# Data Points
#=================================================================

# Create the data points
dataPoints = iron.DataPoints()
dataPoints.CreateStart(region,numberOfDataPoints)

localNumberOfDataPoints = 0
dataPointLocations = numpy.zeros((numberOfDataPoints,numberOfDimensions))
print("Number of data points: " + str(numberOfDataPoints))

# Calculate data point locations of points on a sphere
random.seed(1)
for i in range(numberOfDataPoints):
    if numberOfDimensions == 3:
        theta = 2.0*math.pi*random.uniform(0.0,radius)
        phi = math.acos(2.0*random.uniform(0.0,radius)-radius)
        x = math.cos(theta)*math.sin(phi)
        y = math.sin(theta)*math.sin(phi)
        z = math.cos(phi)
        dataPointLocations[i,:] = [x,y,z]
    if numberOfDimensions == 2:
        theta = random.uniform(0.0,2.0*math.pi)
        x = radius*math.cos(theta)
        y = radius*math.sin(theta)
        dataPointLocations[i,:] = [x,y]

# Set up CMISS data points with geometric values
for dataPoint in range(numberOfDataPoints):
    dataPointId = dataPoint + 1
    dataList = dataPointLocations[dataPoint,:]
    dataPoints.ValuesSet(dataPointId,dataList)

dataPoints.CreateFinish()

# write data points to exdata file for CMGUI
offset = 0
writeExdataFile("DataPoints.part"+str(computationalNodeNumber)+".exdata",dataPointLocations,offset)

#=================================================================
# Data Projection on Geometric Field
#=================================================================

print("Projecting data points onto geometric field")
# Set up data projection
dataProjection = iron.DataProjection()
dataProjection.CreateStart(dataProjectionUserNumber,dataPoints,mesh)
dataProjection.projectionType = iron.DataProjectionProjectionTypes.ALL_ELEMENTS
dataProjection.CreateFinish()

# Evaluate data projection based on geometric field
dataProjection.DataPointsProjectionEvaluate(geometricField)
# Create mesh topology for data projection
mesh.TopologyDataPointsCalculateProjection(dataProjection)
# Create decomposition topology for data projection
decomposition.TopologyDataProjectionCalculate()
print("Projection complete")

# # MPI DEBUG: Check data point info
# numberOfLocalDataPoints = 0
# localDataPointLocations = []
# localDataPoints = []
# elementDataPoints = numpy.zeros((numberOfElements,numberOfDataPoints,numberOfDimensions))
# for element in range(1,numberOfElements+1):
#     elementDomain = decomposition.ElementDomainGet(element)
#     if (elementDomain == computationalNodeNumber):
#         numberOfProjectedDataPoints = decomposition.TopologyNumberOfElementDataPointsGet(element)
#         for dataPoint in range(1,numberOfProjectedDataPoints+1):
#             dataList = numpy.zeros((numberOfDimensions))
#             dataPointNumber = decomposition.TopologyElementDataPointUserNumberGet(element,dataPoint)
#             dataList = dataPoints.ValuesGet(dataPointNumber,numberOfDimensions)
#             numberOfLocalDataPoints += 1
#             localDataPointLocations.append(dataList)
#             localDataPoints.append(dataPointNumber)
#             elementDataPoints[element-1,dataPoint-1,:] = dataList
# localDataPointLocations = numpy.array((localDataPointLocations))
# xi = numpy.zeros((numberOfDimensions))
# for d in range(numberOfLocalDataPoints):
#     localDataPointLocations[d,:] = localDataPointLocations[d]
#     print('Data point: ' + str(localDataPoints[d]))
#     elementNumber = dataProjection.ResultElementNumberGet(localDataPoints[d])
#     print('    ElementNumber: ' + str(elementNumber))
#     xi = dataProjection.ResultXiGet(localDataPoints[d],numberOfDimensions)
#     print('    Xi: ' + str(xi))
#     distance = dataProjection.ResultDistanceGet(localDataPoints[d])
#     print('    Distance: ' + str(distance))
 
#=================================================================
# Equations Set
#=================================================================

# Create vector fitting equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.FITTING,
                             iron.EquationsSetTypes.DATA_FITTING_EQUATION,
                             iron.EquationsSetSubtypes.DATA_POINT_VECTOR_STATIC_FITTING]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

#=================================================================
# Dependent Field
#=================================================================

# Create dependent field (will be deformed fitted values based on data point locations)
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
dependentField.ScalingTypeSet(iron.FieldScalingTypes.UNIT)
equationsSet.DependentCreateFinish()
# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)

# Initialise dependent field to undeformed geometric field
for component in range (1,numberOfDimensions+1):
    geometricField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                            component, dependentField, iron.FieldVariableTypes.U,
                                                            iron.FieldParameterSetTypes.VALUES, component)

#=================================================================
# Independent Field
#=================================================================

# Create data point field (independent field, with vector values stored at the data points)
independentField = iron.Field()
equationsSet.IndependentCreateStart(independentFieldUserNumber,independentField)
independentField.VariableLabelSet(iron.FieldVariableTypes.U,"data point vector")
independentField.VariableLabelSet(iron.FieldVariableTypes.V,"data point weight")
independentField.DataProjectionSet(dataProjection)
equationsSet.IndependentCreateFinish()
# Initialise data point vector field to 0
independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)
# Initialise data point weight field to 1
independentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,1,1.0)

# loop over each element's data points and set independent field values to data point locations on surface of the sphere
for element in range(numberOfElements):
    elementId = element + 1
    elementDomain = decomposition.ElementDomainGet(elementId)
    if (elementDomain == computationalNodeNumber):
        numberOfProjectedDataPoints = decomposition.TopologyNumberOfElementDataPointsGet(elementId)
        for dataPoint in range(numberOfProjectedDataPoints):
            dataPointId = dataPoint + 1
            dataPointNumber = decomposition.TopologyElementDataPointUserNumberGet(elementId,dataPointId)
            dataList = dataPoints.ValuesGet(dataPointNumber,numberOfDimensions)
            # set data point field values
            for component in range(numberOfDimensions):
                componentId = component + 1
                dataPointNumberIndex = dataPointNumber - 1
                value = dataList[component]
                independentField.ParameterSetUpdateElementDataPointDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,elementId,dataPointId,componentId,value)

#=================================================================
# Material Field
#=================================================================

# Create material field (Sobolev parameters)
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
materialField.VariableLabelSet(iron.FieldVariableTypes.U,"Smoothing Parameters")
equationsSet.MaterialsCreateFinish()

# Set kappa and tau - Sobolev smoothing parameters
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,tau)
materialField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,kappa)

#=================================================================
# Equations
#=================================================================

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.FULL
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

#=================================================================
# Problem setup
#=================================================================

# Create fitting problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.FITTING,
                        iron.ProblemTypes.DATA_FITTING,
                        iron.ProblemSubtypes.DATA_POINT_VECTOR_STATIC_FITTING]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.outputType = iron.SolverOutputTypes.NONE # NONE / MATRIX
solver.linearType = iron.LinearSolverTypes.ITERATIVE
solver.LibraryTypeSet(iron.SolverLibraries.UMFPACK) # UMFPACK/SUPERLU
solver.linearIterativeAbsoluteTolerance = 1.0E-10
solver.linearIterativeRelativeTolerance = 1.0E-05
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.FULL
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

#=================================================================
# Boundary Conditions
#=================================================================

# Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

version = 1
meshComponent = decomposition.MeshComponentGet()
# Fix the interior nodes- use to only apply fit to surface nodes
if (fixInterior):
    # first find which nodes are non-surface nodes
    for node in range(numberOfNodes):
        nodeId = node + 1
        nodeDomain = decomposition.NodeDomainGet(nodeId,meshComponent)
        if (nodeDomain == computationalNodeNumber):
            geometricValue = numpy.zeros((numberOfDimensions))
            for component in range(numberOfDimensions):
                componentId = component+1
                geometricValue[component]=geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                                               iron.FieldParameterSetTypes.VALUES,
                                                                               1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,
                                                                               nodeId,componentId)

            derivList = []
            deriv = 0
            if hermite:
                if numberOfDimensions ==3:
                    if (abs(geometricValue[0]) < (abs(meshOrigin[0]) - zeroTolerance)) and (abs(geometricValue[1]) < (abs(meshOrigin[1]) - zeroTolerance)) and (abs(geometricValue[2]) < (abs(meshOrigin[2]) - zeroTolerance)):
                        # Interior nodes
                        derivList = [1,2,3,4,5,6,7,8]
                    # Radial nodes
                    elif abs(geometricValue[0]) < zeroTolerance and abs(geometricValue[1]) < zeroTolerance:
                        deriv = 5
                    elif abs(geometricValue[1]) < zeroTolerance and abs(geometricValue[2]) < zeroTolerance:
                        deriv = 2
                    elif abs(geometricValue[0]) < zeroTolerance and abs(geometricValue[2]) < zeroTolerance:
                        deriv = 3

                    if deriv > 0 and deriv not in derivList:
                        derivList.append(deriv)

                elif numberOfDimensions ==2:
                    if (abs(geometricValue[0]) < (abs(meshOrigin[0]) - zeroTolerance)) and (abs(geometricValue[1]) < (abs(meshOrigin[1]) - zeroTolerance)):
                        # Interior nodes
                        derivList = [1,2,3,4]
                    # Radial nodes
                    elif abs(geometricValue[1]) < zeroTolerance:
                        deriv = 2
                    elif abs(geometricValue[0]) < zeroTolerance:
                        deriv = 3

                    if deriv > 0 and deriv not in derivList:
                        derivList.append(deriv)

            elif numberOfDimensions == 3:
                if (abs(geometricValue[0]) < (abs(meshOrigin[0]) - zeroTolerance)) and (abs(geometricValue[1]) < (abs(meshOrigin[1]) - zeroTolerance)) and (abs(geometricValue[2]) < (abs(meshOrigin[2]) - zeroTolerance)):
                    derivList = [1]
            elif numberOfDimensions == 2:
                if (abs(geometricValue[0]) < (abs(meshOrigin[0]) - zeroTolerance)) and (abs(geometricValue[1]) < (abs(meshOrigin[1]) - zeroTolerance)):
                    derivList = [1]

            for globalDeriv in derivList: 
                for component in range(1,numberOfDimensions+1):
                    value=geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                               iron.FieldParameterSetTypes.VALUES,
                                                               version,globalDeriv,nodeId,component)
                    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,
                                               version,globalDeriv,nodeId,component,
                                               iron.BoundaryConditionsTypes.FIXED,value)

solverEquations.BoundaryConditionsCreateFinish()


#=================================================================
# S o l v e    a n d    E x p o r t    D a t a
#=================================================================
for iteration in range (1,numberOfIterations+1):

    # Solve the problem
    print("Solving fitting problem, iteration: " + str(iteration))
    problem.Solve()

    # Copy dependent field to geometric 
    for component in range(1,numberOfDimensions+1):
        dependentField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,component,geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,component)

    # Export fields
    print("Writing deformed geometry")
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("DeformedGeometry" + str(iteration),"FORTRAN")
    fields.ElementsExport("DeformedGeometry" + str(iteration),"FORTRAN")
    fields.Finalise()

#-----------------------------------------------------------------

iron.Finalise()
