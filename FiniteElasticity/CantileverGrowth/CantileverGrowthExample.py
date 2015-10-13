#!/usr/bin/env python

#> \file
#> \author Chris Bradley
#> \brief This is an example script to solve a finite elasticity cantilever problem with a growth and constituative law in CellML. The growth occurs just at the bottom of the cantilever in order to cause upward bending. 
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

#> Main script
# Add Python bindings directory to PATH
import sys, os

sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

# Intialise OpenCMISS
from opencmiss import iron

# Set the physical size of the cantilever
length = 5.0
width = 1.0
height = 1.0

# Set the number of elements in the cantilever
numberXElements = 5
numberYElements = 1
numberZElements = 1

# Set the number of Gauss points to use in each direction
numberOfGaussXi = 3

# Set the growth rates
xGrowthRate = 0.02
yGrowthRate = 0.0
zGrowthRate = 0.0

# Set the similation times.
startTime = 0.0
stopTime = 10.0
timeIncrement = 1.0

# Set the user numbers
coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
generatedMeshUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
geometricFieldUserNumber = 1
fibreFieldUserNumber = 2
dependentFieldUserNumber = 3
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

# Get the number of computational nodes and this computational node number
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Create a 3D rectangular cartesian coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.LabelSet("Region")
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Define a basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = 3
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
basis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
basis.CreateFinish()

# Start the creation of a generated mesh in the region
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [length,width,height]
generatedMesh.numberOfElements = [numberXElements,numberYElements,numberZElements]
# Finish the creation of a generated mesh in the region
mesh = iron.Mesh()
generatedMesh.CreateFinish(meshUserNumber,mesh)

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
geometricField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
geometricField.CreateFinish()

# Update the geometric field parameters from generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create a fibre field and attach it to the geometric field
fibreField = iron.Field()
fibreField.CreateStart(fibreFieldUserNumber,region)
fibreField.TypeSet(iron.FieldTypes.FIBRE)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(iron.FieldVariableTypes.U,"Fibre")
fibreField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
fibreField.CreateFinish()

# Create the dependent field
dependentField = iron.Field()
dependentField.CreateStart(dependentFieldUserNumber,region)
dependentField.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)  
dependentField.MeshDecompositionSet(decomposition)
dependentField.GeometricFieldSet(geometricField) 
dependentField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT) 
# Set the field to have 5 variables: U - dependent; del U/del n - tractions; U1 - strain; U2 - stress; U3 - growth
dependentField.NumberOfVariablesSet(5)
dependentField.VariableTypesSet([iron.FieldVariableTypes.U,iron.FieldVariableTypes.DELUDELN,iron.FieldVariableTypes.U1,iron.FieldVariableTypes.U2,iron.FieldVariableTypes.U3])
dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
dependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,"del U/del n")
dependentField.VariableLabelSet(iron.FieldVariableTypes.U1,"Strain")
dependentField.VariableLabelSet(iron.FieldVariableTypes.U2,"Stress")
dependentField.VariableLabelSet(iron.FieldVariableTypes.U3,"Growth")
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,4)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,4)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U1,6)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U2,6)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U3,3)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.ELEMENT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.ELEMENT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,1,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,2,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,3,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,4,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,5,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,6,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,1,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,2,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,3,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,4,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,5,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,6,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U3,1,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U3,2,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U3,3,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
dependentField.CreateFinish()

# Initialise dependent field from undeformed geometry
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
# Initialise the hydrostatic pressure
iron.Field.ComponentValuesInitialiseDP(
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,-8.0)

# Update the dependent field
dependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
dependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# Create the equations_set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
    iron.EquationsSetTypes.FINITE_ELASTICITY,
    iron.EquationsSetSubtypes.CONSTIT_AND_GROWTH_LAW_IN_CELLML]
equationsSet.CreateStart(equationsSetUserNumber,region,fibreField,
                         equationsSetSpecification,equationsSetFieldUserNumber,equationsSetField)
equationsSet.CreateFinish()

equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
equationsSet.DependentCreateFinish()

# Create the CellML environment for the growth law. Set the rates as known so that we can spatially vary them.
growthCellML = iron.CellML()
growthCellML.CreateStart(growthCellMLUserNumber,region)
growthCellMLIdx = growthCellML.ModelImport("simplegrowth.cellml")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/fibrerate")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/sheetrate")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/normalrate")
growthCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps. Map the lambda's to the U3/growth dependent field variable
growthCellML.FieldMapsCreateStart()
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda1",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U3,1,iron.FieldParameterSetTypes.VALUES)
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda2",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U3,2,iron.FieldParameterSetTypes.VALUES)
growthCellML.CreateCellMLToFieldMap(growthCellMLIdx,"Main/lambda3",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U3,3,iron.FieldParameterSetTypes.VALUES)
growthCellML.FieldMapsCreateFinish()

# Create the CELL models field
growthCellMLModelsField = iron.Field()
growthCellML.ModelsFieldCreateStart(growthCellMLModelsFieldUserNumber,growthCellMLModelsField)
growthCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthModelMap")
growthCellML.ModelsFieldCreateFinish()

# Create the CELL parameters field
growthCellMLParametersField = iron.Field()
growthCellML.ParametersFieldCreateStart(growthCellMLParametersFieldUserNumber,growthCellMLParametersField)
growthCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthParameters")
growthCellML.ParametersFieldCreateFinish()

# Set the parameters so that only the bottom layer of Gauss points has a non-zero growth rate
for yElem in range(1, numberYElements+1):
    for xElem in range(1, numberXElements+1):
        elementNumber = xElem + (yElem-1)*numberXElements
        for yGauss in range(1, numberOfGaussXi+1):
            for xGauss in range(1, numberOfGaussXi+1):
                gaussNumber = xGauss + (yGauss-1)*numberOfGaussXi
                #print 'Setting growth parameter at element ',elementNumber,' and Gauss point number ',gaussNumber
                growthCellMLParametersField.ParameterSetUpdateGaussPoint(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, gaussNumber, elementNumber, 1, xGrowthRate)
                growthCellMLParametersField.ParameterSetUpdateGaussPoint(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, gaussNumber, elementNumber, 2, yGrowthRate)
                growthCellMLParametersField.ParameterSetUpdateGaussPoint(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, gaussNumber, elementNumber, 3, zGrowthRate)
 
# Update the parameters field
growthCellMLParametersField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
growthCellMLParametersField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

# Create the CELL state field
growthCellMLStateField = iron.Field()
growthCellML.StateFieldCreateStart(growthCellMLStateFieldUserNumber,growthCellMLStateField)
growthCellMLStateField.VariableLabelSet(iron.FieldVariableTypes.U,"GrowthState")
growthCellML.StateFieldCreateFinish()

# Create the CellML environment for the consitutative law
constituativeCellML = iron.CellML()
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

# Create CellML <--> OpenCMISS field maps. Map the stress and strain fields.
constituativeCellML.FieldMapsCreateStart()
constituativeCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,1,iron.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/E11",iron.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,2,iron.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/E12",iron.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,3,iron.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/E13",iron.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,4,iron.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/E22",iron.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,5,iron.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/E23",iron.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,6,iron.FieldParameterSetTypes.VALUES,
    constituativeCellMLIdx,"equations/E33",iron.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev11",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,1,iron.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev12",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,2,iron.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev13",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,3,iron.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev22",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,4,iron.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev23",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,5,iron.FieldParameterSetTypes.VALUES)
constituativeCellML.CreateCellMLToFieldMap(constituativeCellMLIdx,"equations/Tdev33",iron.FieldParameterSetTypes.VALUES,
    dependentField,iron.FieldVariableTypes.U2,6,iron.FieldParameterSetTypes.VALUES)
constituativeCellML.FieldMapsCreateFinish()

# Create the CELL models field
constituativeCellMLModelsField = iron.Field()
constituativeCellML.ModelsFieldCreateStart(constituativeCellMLModelsFieldUserNumber,constituativeCellMLModelsField)
constituativeCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstituativeModelMap")
constituativeCellML.ModelsFieldCreateFinish()

# Create the CELL parameters field
constituativeCellMLParametersField = iron.Field()
constituativeCellML.ParametersFieldCreateStart(constituativeCellMLParametersFieldUserNumber,constituativeCellMLParametersField)
constituativeCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstituativeParameters")
constituativeCellML.ParametersFieldCreateFinish()

# Create the CELL intermediate field
constituativeCellMLIntermediateField = iron.Field()
constituativeCellML.IntermediateFieldCreateStart(constituativeCellMLIntermediateFieldUserNumber,constituativeCellMLIntermediateField)
constituativeCellMLIntermediateField.VariableLabelSet(iron.FieldVariableTypes.U,"ConstituativeIntermediate")
constituativeCellML.IntermediateFieldCreateFinish()

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Define the problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.ELASTICITY,
        iron.ProblemTypes.FINITE_ELASTICITY,
        iron.ProblemSubtypes.FINITE_ELASTICITY_WITH_GROWTH_CELLML]
problem.CreateStart(problemUserNumber,problemSpecification)
problem.CreateFinish()

# Create control loops
timeLoop = iron.ControlLoop()
problem.ControlLoopCreateStart()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],timeLoop)
timeLoop.TimesSet(startTime,stopTime,timeIncrement)
problem.ControlLoopCreateFinish()

# Create problem solvers
odeIntegrationSolver = iron.Solver()
nonlinearSolver = iron.Solver()
linearSolver = iron.Solver()
cellMLEvaluationSolver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,odeIntegrationSolver)
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,nonlinearSolver)
nonlinearSolver.outputType = iron.SolverOutputTypes.PROGRESS
nonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
nonlinearSolver.NewtonCellMLSolverGet(cellMLEvaluationSolver)
nonlinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.linearType = iron.LinearSolverTypes.DIRECT
problem.SolversCreateFinish()

# Create nonlinear equations and add equations set to solver equations
nonlinearEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
nonlinearSolver.SolverEquationsGet(nonlinearEquations)
nonlinearEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
nonlinearEquationsSetIndex = nonlinearEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Create CellML equations and add growth and constituative equations to the solvers
growthEquations = iron.CellMLEquations()
constituativeEquations = iron.CellMLEquations()
problem.CellMLEquationsCreateStart()
odeIntegrationSolver.CellMLEquationsGet(growthEquations)
growthEquationsIndex = growthEquations.CellMLAdd(growthCellML)
cellMLEvaluationSolver.CellMLEquationsGet(constituativeEquations)
constituativeEquationsIndex = constituativeEquations.CellMLAdd(constituativeCellML)
problem.CellMLEquationsCreateFinish()

# Prescribe boundary conditions (absolute nodal parameters)
boundaryConditions = iron.BoundaryConditions()
nonlinearEquations.BoundaryConditionsCreateStart(boundaryConditions)

#Set x=0 nodes to built-in
for zNode in range(1, numberZElements+2):
    for yNode in range(1, numberYElements+2):
        nodeNumber = 1+(yNode-1)*(numberXElements+1)+(zNode-1)*(numberXElements+1)*(numberYElements+1)
        #print 'Setting boundary condition at node ',nodeNumber
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)
 
nonlinearEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("CantileverGrowth","FORTRAN")
fields.ElementsExport("CantileverGrowth","FORTRAN")
fields.Finalise()

