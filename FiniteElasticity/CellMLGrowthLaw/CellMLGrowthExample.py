#!/usr/bin/env python

#> \file
#> \author Chris Bradley
#> \brief This is an example script to solve a finite elasticity problem with a growth and constituation law in CellML
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

#> Main script
# Add Python bindings directory to PATH
import sys, os

sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

# Intialise OpenCMISS
from opencmiss import CMISS

# Set problem parameters
height = 1.0
width = 1.0
length = 1.0

NumberOfGaussXi = 2

startTime = 0.0
stopTime = 10.0
timeIncrement = 1.0

coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
generatedMeshUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
geometricFieldUserNumber = 1
fibreFieldUserNumber = 2
materialFieldUserNumber = 3
dependentFieldUserNumber = 4
equationsSetUserNumber = 1
equationsSetFieldUserNumber = 5
growthCellMLUserNumber = 1
growthCellMLModelsFieldUserNumber = 6
growthCellMLStateFieldUserNumber = 7
growthCellMLParametersFieldUserNumber = 8
constituativeCellMLUserNumber = 2
constituativeCellMLModelsFieldUserNumber = 9
constituativeCellMLParametersFieldUserNumber = 10
constituativeCellMLIntermediateFieldUserNumber = 11
problemUserNumber = 1

InterpolationType = 1

# Get the number of computational nodes and this computational node number
numberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet()
computationalNodeNumber = CMISS.ComputationalNodeNumberGet()

# Create a 3D rectangular cartesian coordinate system
coordinateSystem = CMISS.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = CMISS.Region()
region.CreateStart(regionUserNumber,CMISS.WorldRegion)
region.LabelSet("Region")
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Define basis
basis = CMISS.Basis()
basis.CreateStart(basisUserNumber)
if InterpolationType in (1,2,3,4):
    basis.type = CMISS.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = 3
basis.interpolationXi = [CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
if(NumberOfGaussXi>0):
    basis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*3
basis.CreateFinish()

# Start the creation of a manually generated mesh in the region
mesh = CMISS.Mesh()
mesh.CreateStart(meshUserNumber,region,3)
mesh.NumberOfComponentsSet(1)
mesh.NumberOfElementsSet(1)

#Define nodes for the mesh
nodes = CMISS.Nodes()
nodes.CreateStart(region,8)
nodes.CreateFinish()

elements = CMISS.MeshElements()
elements.CreateStart(mesh,1,basis)
elements.NodesSet(1,[1,2,3,4,5,6,7,8])
elements.CreateFinish()

mesh.CreateFinish() 

# Create a decomposition for the mesh
decomposition = CMISS.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = CMISS.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = CMISS.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(CMISS.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(CMISS.FieldVariableTypes.U,"Geometry")
geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,3,1)
if InterpolationType == 4:
    geometricField.fieldScalingType = CMISS.FieldScalingTypes.ARITHMETIC_MEAN
geometricField.CreateFinish()

# Update the geometric field parameters manually
geometricField.ParameterSetUpdateStart(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)
# node 1
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,1,1,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,1,2,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,1,3,0.0)
# node 2
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,2,1,height)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,2,2,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,2,3,0.0)
# node 3
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,3,1,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,3,2,width)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,3,3,0.0)
# node 4
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,4,1,height)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,4,2,width)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,4,3,0.0)
# node 5
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,5,1,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,5,2,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,5,3,length)
# node 6
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,6,1,height)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,6,2,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,6,3,length)
# node 7
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,7,1,0.0)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,7,2,width)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,7,3,length)
# node 8
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,8,1,height)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,8,2,width)
geometricField.ParameterSetUpdateNodeDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,1,8,3,length)
geometricField.ParameterSetUpdateFinish(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)

# Create a fibre field and attach it to the geometric field
fibreField = CMISS.Field()
fibreField.CreateStart(fibreFieldUserNumber,region)
fibreField.TypeSet(CMISS.FieldTypes.FIBRE)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(CMISS.FieldVariableTypes.U,"Fibre")
if InterpolationType == 4:
    fibreField.fieldScalingType = CMISS.FieldScalingTypes.ARITHMETIC_MEAN
fibreField.CreateFinish()

# Create the material field
materialField = CMISS.Field()
materialField.CreateStart(materialFieldUserNumber,region)
materialField.TypeSet(CMISS.FieldTypes.MATERIAL)
materialField.MeshDecompositionSet(decomposition)
materialField.GeometricFieldSet(geometricField)
materialField.NumberOfVariablesSet(1)
materialField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U,2)
materialField.VariableLabelSet(CMISS.FieldVariableTypes.U,"Material")
materialField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,1,1)
materialField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,2,1)
if InterpolationType == 4:
    materialField.fieldScalingType = CMISS.FieldScalingTypes.ARITHMETIC_MEAN
materialField.CreateFinish()

# Set Mooney-Rivlin constants c10 and c01 respectively.
CMISS.Field.ComponentValuesInitialiseDP(
    materialField,CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,2.0)
CMISS.Field.ComponentValuesInitialiseDP(
    materialField,CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,2,6.0)

# Create the dependent field
dependentField = CMISS.Field()
dependentField.CreateStart(dependentFieldUserNumber,region)
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.U,"Dependent")
dependentField.TypeSet(CMISS.FieldTypes.GEOMETRIC_GENERAL)  
dependentField.MeshDecompositionSet(decomposition)
dependentField.GeometricFieldSet(geometricField) 
dependentField.DependentTypeSet(CMISS.FieldDependentTypes.DEPENDENT) 
dependentField.NumberOfVariablesSet(5)
dependentField.VariableTypesSet([CMISS.FieldVariableTypes.U,CMISS.FieldVariableTypes.DELUDELN,CMISS.FieldVariableTypes.U1,CMISS.FieldVariableTypes.U2,CMISS.FieldVariableTypes.U3])
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U,3)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.DELUDELN,3)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U1,6)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U2,6)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U3,3)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U1,1,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U1,2,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U1,3,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U1,4,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U1,5,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U1,6,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U2,1,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U2,2,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U2,3,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U2,4,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U2,5,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U2,6,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U3,1,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U3,2,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U3,3,CMISS.FieldInterpolationTypes.GAUSS_POINT_BASED)
if InterpolationType == 4:
    dependentField.fieldScalingType = CMISS.FieldScalingTypes.ARITHMETIC_MEAN
dependentField.CreateFinish()

# Initialise dependent field from undeformed geometry
CMISS.Field.ParametersToFieldParametersComponentCopy(
    geometricField,CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,
    dependentField,CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1)
CMISS.Field.ParametersToFieldParametersComponentCopy(
    geometricField,CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,2,
    dependentField,CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,2)
CMISS.Field.ParametersToFieldParametersComponentCopy(
    geometricField,CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,3,
    dependentField,CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,3)

# Create the equations_set
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equationsSetUserNumber,region,fibreField,
    CMISS.EquationsSetClasses.ELASTICITY,
    CMISS.EquationsSetTypes.FINITE_ELASTICITY,
    CMISS.EquationsSetSubtypes.CONSTIT_AND_GROWTH_LAW_IN_CELLML,
    equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
equationsSet.MaterialsCreateFinish()

equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
equationsSet.DependentCreateFinish()

# Create the CellML environment for the growth law
growthCellML = CMISS.CellML()
growthCellML.CreateStart(growthCellMLUserNumber,region)
growthCellMLIdx = growthCellML.ModelImport("simplegrowth.cellml")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/S11")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/S22")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/S33")
#growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/bff")
#growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/bss")
#growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/bnn")
growthCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps
growthCellML.FieldMapsCreateStart()
growthCellML.CreateFieldToCellMLMap(dependentField,CMISS.FieldVariableTypes.U2,1,CMISS.FieldParameterSetTypes.VALUES,
    growthCellMLIdx,"Main/S11",CMISS.FieldParameterSetTypes.VALUES)
growthCellML.CreateFieldToCellMLMap(dependentField,CMISS.FieldVariableTypes.U2,4,CMISS.FieldParameterSetTypes.VALUES,
    growthCellMLIdx,"Main/S22",CMISS.FieldParameterSetTypes.VALUES)
growthCellML.CreateFieldToCellMLMap(dependentField,CMISS.FieldVariableTypes.U2,6,CMISS.FieldParameterSetTypes.VALUES,
    growthCellMLIdx,"Main/S33",CMISS.FieldParameterSetTypes.VALUES)
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda1",CMISS.FieldParameterSetTypes.VALUES,
    dependentField,CMISS.FieldVariableTypes.U3,1,CMISS.FieldParameterSetTypes.VALUES)
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda2",CMISS.FieldParameterSetTypes.VALUES,
    dependentField,CMISS.FieldVariableTypes.U3,2,CMISS.FieldParameterSetTypes.VALUES)
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda3",CMISS.FieldParameterSetTypes.VALUES,
    dependentField,CMISS.FieldVariableTypes.U3,3,CMISS.FieldParameterSetTypes.VALUES)
growthCellML.FieldMapsCreateFinish()

# Create the CELL models field
growthCellMLModelsField = CMISS.Field()
growthCellML.ModelsFieldCreateStart(growthCellMLModelsFieldUserNumber,growthCellMLModelsField)
growthCellML.ModelsFieldCreateFinish()

# Create the CELL parameters field
growthCellMLParametersField = CMISS.Field()
growthCellML.ParametersFieldCreateStart(growthCellMLParametersFieldUserNumber,growthCellMLParametersField)
growthCellML.ParametersFieldCreateFinish()

# Create the CELL state field
growthCellMLStateField = CMISS.Field()
growthCellML.StateFieldCreateStart(growthCellMLStateFieldUserNumber,growthCellMLStateField)
growthCellML.StateFieldCreateFinish()

# Create the CellML environment for the consitutative law
constituativeCellML = CMISS.CellML()
constituativeCellML.CreateStart(constituativeCellMLUserNumber,region)
constituativeCellMLIdx = constituativeCellML.ModelImport("mooneyrivlin.cellml")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/E11")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/E12")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/E13")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/E22")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/E23")
constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/E33")
#constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/c1")
#constituativeCellML.VariableSetAsKnown(constituativeCellMLIdx,"equations/c2")
constituativeCellML.VariableSetAsWanted(constituativeCellMLIdx,"equations/Tdev11")
constituativeCellML.VariableSetAsWanted(constituativeCellMLIdx,"equations/Tdev12")
constituativeCellML.VariableSetAsWanted(constituativeCellMLIdx,"equations/Tdev13")
constituativeCellML.VariableSetAsWanted(constituativeCellMLIdx,"equations/Tdev22")
constituativeCellML.VariableSetAsWanted(constituativeCellMLIdx,"equations/Tdev23")
constituativeCellML.VariableSetAsWanted(constituativeCellMLIdx,"equations/Tdev33")
constituativeCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps
constituativeCellML.FieldMapsCreateStart()
constituativeCellML.CreateFieldToCellMLMap(dependentField,CMISS.FieldVariableTypes.U1,1,CMISS.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/E11",CMISS.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,CMISS.FieldVariableTypes.U1,2,CMISS.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/E12",CMISS.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,CMISS.FieldVariableTypes.U1,3,CMISS.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/E13",CMISS.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,CMISS.FieldVariableTypes.U1,4,CMISS.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/E22",CMISS.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,CMISS.FieldVariableTypes.U1,5,CMISS.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/E23",CMISS.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,CMISS.FieldVariableTypes.U1,6,CMISS.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/E33",CMISS.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev11",CMISS.FieldParameterSetTypes.VALUES,
    dependentField,CMISS.FieldVariableTypes.U2,1,CMISS.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev12",CMISS.FieldParameterSetTypes.VALUES,
    dependentField,CMISS.FieldVariableTypes.U2,2,CMISS.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev13",CMISS.FieldParameterSetTypes.VALUES,
    dependentField,CMISS.FieldVariableTypes.U2,3,CMISS.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev22",CMISS.FieldParameterSetTypes.VALUES,
    dependentField,CMISS.FieldVariableTypes.U2,4,CMISS.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev23",CMISS.FieldParameterSetTypes.VALUES,
    dependentField,CMISS.FieldVariableTypes.U2,5,CMISS.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev33",CMISS.FieldParameterSetTypes.VALUES,
    dependentField,CMISS.FieldVariableTypes.U2,6,CMISS.FieldParameterSetTypes.VALUES)
constituativeCellML.FieldMapsCreateFinish()

# Create the CELL models field
constituativeCellMLModelsField = CMISS.Field()
constituativeCellML.ModelsFieldCreateStart(constituativeCellMLModelsFieldUserNumber,constituativeCellMLModelsField)
constituativeCellML.ModelsFieldCreateFinish()

# Create the CELL parameters field
constituativeCellMLParametersField = CMISS.Field()
constituativeCellML.ParametersFieldCreateStart(constituativeCellMLParametersFieldUserNumber,constituativeCellMLParametersField)
constituativeCellML.ParametersFieldCreateFinish()

# Create the CELL intermediate field
constituativeCellMLIntermediateField = CMISS.Field()
constituativeCellML.IntermediateFieldCreateStart(constituativeCellMLIntermediateFieldUserNumber,constituativeCellMLIntermediateField)
constituativeCellML.IntermediateFieldCreateFinish()

# Create equations
equations = CMISS.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = CMISS.EquationsSparsityTypes.SPARSE
equations.outputType = CMISS.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Define the problem
problem = CMISS.Problem()
problem.CreateStart(problemUserNumber)
problem.SpecificationSet(CMISS.ProblemClasses.ELASTICITY,
        CMISS.ProblemTypes.FINITE_ELASTICITY,
        CMISS.ProblemSubTypes.FINITE_ELASTICITY_WITH_GROWTH_CELLML)
problem.CreateFinish()

# Create control loops
timeLoop = CMISS.ControlLoop()
problem.ControlLoopCreateStart()
problem.ControlLoopGet([CMISS.ControlLoopIdentifiers.NODE],timeLoop)
timeLoop.TimesSet(startTime,stopTime,timeIncrement)
problem.ControlLoopCreateFinish()

# Create problem solvers
odeIntegrationSolver = CMISS.Solver()
nonlinearSolver = CMISS.Solver()
linearSolver = CMISS.Solver()
cellMLEvaluationSolver = CMISS.Solver()
problem.SolversCreateStart()
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],1,odeIntegrationSolver)
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],2,nonlinearSolver)
nonlinearSolver.outputType = CMISS.SolverOutputTypes.PROGRESS
nonlinearSolver.NewtonJacobianCalculationTypeSet(CMISS.JacobianCalculationTypes.FD)
nonlinearSolver.NewtonCellMLSolverGet(cellMLEvaluationSolver)
nonlinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.linearType = CMISS.LinearSolverTypes.DIRECT
problem.SolversCreateFinish()

# Create nonlinear equations and add equations set to solver equations
nonlinearEquations = CMISS.SolverEquations()
problem.SolverEquationsCreateStart()
nonlinearSolver.SolverEquationsGet(nonlinearEquations)
nonlinearEquations.sparsityType = CMISS.SolverEquationsSparsityTypes.SPARSE
nonlinearEquationsSetIndex = nonlinearEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Create CellML equations and add growth and constituative equations to the solvers
growthEquations = CMISS.CellMLEquations()
constituativeEquations = CMISS.CellMLEquations()
problem.CellMLEquationsCreateStart()
odeIntegrationSolver.CellMLEquationsGet(growthEquations)
growthEquationsIndex = growthEquations.CellMLAdd(growthCellML)
cellMLEvaluationSolver.CellMLEquationsGet(constituativeEquations)
constituativeEquationsIndex = constituativeEquations.CellMLAdd(constituativeCellML)
problem.CellMLEquationsCreateFinish()

# Prescribe boundary conditions (absolute nodal parameters)
boundaryConditions = CMISS.BoundaryConditions()
nonlinearEquations.BoundaryConditionsCreateStart(boundaryConditions)

#Set x=0 nodes to no x displacment in x. Set x=width nodes to 10% x displacement
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,1,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,3,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,5,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,7,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)

boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,2,1,CMISS.BoundaryConditionsTypes.FIXED,0.1*width)
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,4,1,CMISS.BoundaryConditionsTypes.FIXED,0.1*width)
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,6,1,CMISS.BoundaryConditionsTypes.FIXED,0.1*width)
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,8,1,CMISS.BoundaryConditionsTypes.FIXED,0.1*width)

# Set y=0 nodes to no y displacement
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,1,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,2,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,5,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,6,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)

# Set z=0 nodes to no y displacement
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,1,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,2,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,3,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,4,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)

nonlinearEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Export results
fields = CMISS.Fields()
fields.CreateRegion(region)
fields.NodesExport("CellMLGrowth","FORTRAN")
fields.ElementsExport("CellMLGrowth","FORTRAN")
fields.Finalise()

