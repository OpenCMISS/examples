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

# Define a basis
basis = CMISS.Basis()
basis.CreateStart(basisUserNumber)
basis.type = CMISS.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = 3
basis.interpolationXi = [CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
basis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
basis.CreateFinish()

# Start the creation of a generated mesh in the region
generatedMesh = CMISS.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = CMISS.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [length,width,height]
generatedMesh.numberOfElements = [numberXElements,numberYElements,numberZElements]
# Finish the creation of a generated mesh in the region
mesh = CMISS.Mesh()
generatedMesh.CreateFinish(meshUserNumber,mesh)

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
geometricField.fieldScalingType = CMISS.FieldScalingTypes.ARITHMETIC_MEAN
geometricField.CreateFinish()

# Update the geometric field parameters from generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create a fibre field and attach it to the geometric field
fibreField = CMISS.Field()
fibreField.CreateStart(fibreFieldUserNumber,region)
fibreField.TypeSet(CMISS.FieldTypes.FIBRE)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(CMISS.FieldVariableTypes.U,"Fibre")
fibreField.fieldScalingType = CMISS.FieldScalingTypes.ARITHMETIC_MEAN
fibreField.CreateFinish()

# Create the dependent field
dependentField = CMISS.Field()
dependentField.CreateStart(dependentFieldUserNumber,region)
dependentField.TypeSet(CMISS.FieldTypes.GEOMETRIC_GENERAL)  
dependentField.MeshDecompositionSet(decomposition)
dependentField.GeometricFieldSet(geometricField) 
dependentField.DependentTypeSet(CMISS.FieldDependentTypes.DEPENDENT) 
# Set the field to have 5 variables: U - dependent; del U/del n - tractions; U1 - strain; U2 - stress; U3 - growth
dependentField.NumberOfVariablesSet(5)
dependentField.VariableTypesSet([CMISS.FieldVariableTypes.U,CMISS.FieldVariableTypes.DELUDELN,CMISS.FieldVariableTypes.U1,CMISS.FieldVariableTypes.U2,CMISS.FieldVariableTypes.U3])
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.U,"Dependent")
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.DELUDELN,"del U/del n")
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.U1,"Strain")
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.U2,"Stress")
dependentField.VariableLabelSet(CMISS.FieldVariableTypes.U3,"Growth")
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U,4)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.DELUDELN,4)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U1,6)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U2,6)
dependentField.NumberOfComponentsSet(CMISS.FieldVariableTypes.U3,3)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.U,4,CMISS.FieldInterpolationTypes.ELEMENT_BASED)
dependentField.ComponentInterpolationSet(CMISS.FieldVariableTypes.DELUDELN,4,CMISS.FieldInterpolationTypes.ELEMENT_BASED)
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
# Initialise the hydrostatic pressure
CMISS.Field.ComponentValuesInitialiseDP(
    dependentField,CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,4,-8.0)

# Update the dependent field
dependentField.ParameterSetUpdateStart(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)
dependentField.ParameterSetUpdateFinish(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)

# Create the equations_set
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equationsSetUserNumber,region,fibreField,
    CMISS.EquationsSetClasses.ELASTICITY,
    CMISS.EquationsSetTypes.FINITE_ELASTICITY,
    CMISS.EquationsSetSubtypes.CONSTIT_AND_GROWTH_LAW_IN_CELLML,
    equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
equationsSet.DependentCreateFinish()

# Create the CellML environment for the growth law. Set the rates as known so that we can spatially vary them.
growthCellML = CMISS.CellML()
growthCellML.CreateStart(growthCellMLUserNumber,region)
growthCellMLIdx = growthCellML.ModelImport("simplegrowth.cellml")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/fibrerate")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/sheetrate")
growthCellML.VariableSetAsKnown(growthCellMLIdx,"Main/normalrate")
growthCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps. Map the lambda's to the U3/growth dependent field variable
growthCellML.FieldMapsCreateStart()
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
growthCellMLModelsField.VariableLabelSet(CMISS.FieldVariableTypes.U,"GrowthModelMap")
growthCellML.ModelsFieldCreateFinish()

# Create the CELL parameters field
growthCellMLParametersField = CMISS.Field()
growthCellML.ParametersFieldCreateStart(growthCellMLParametersFieldUserNumber,growthCellMLParametersField)
growthCellMLParametersField.VariableLabelSet(CMISS.FieldVariableTypes.U,"GrowthParameters")
growthCellML.ParametersFieldCreateFinish()

# Set the parameters so that only the bottom layer of Gauss points has a non-zero growth rate
for yElem in range(1, numberYElements+1):
    for xElem in range(1, numberXElements+1):
        elementNumber = xElem + (yElem-1)*numberXElements
        for yGauss in range(1, numberOfGaussXi+1):
            for xGauss in range(1, numberOfGaussXi+1):
                gaussNumber = xGauss + (yGauss-1)*numberOfGaussXi
                #print 'Setting growth parameter at element ',elementNumber,' and Gauss point number ',gaussNumber
                growthCellMLParametersField.ParameterSetUpdateGaussPoint(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, gaussNumber, elementNumber, 1, xGrowthRate)
                growthCellMLParametersField.ParameterSetUpdateGaussPoint(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, gaussNumber, elementNumber, 2, yGrowthRate)
                growthCellMLParametersField.ParameterSetUpdateGaussPoint(CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES, gaussNumber, elementNumber, 3, zGrowthRate)
 
# Update the parameters field
growthCellMLParametersField.ParameterSetUpdateStart(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)
growthCellMLParametersField.ParameterSetUpdateFinish(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES)

# Create the CELL state field
growthCellMLStateField = CMISS.Field()
growthCellML.StateFieldCreateStart(growthCellMLStateFieldUserNumber,growthCellMLStateField)
growthCellMLStateField.VariableLabelSet(CMISS.FieldVariableTypes.U,"GrowthState")
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

# Create CellML <--> OpenCMISS field maps. Map the stress and strain fields.
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
constituativeCellMLModelsField.VariableLabelSet(CMISS.FieldVariableTypes.U,"ConstituativeModelMap")
constituativeCellML.ModelsFieldCreateFinish()

# Create the CELL parameters field
constituativeCellMLParametersField = CMISS.Field()
constituativeCellML.ParametersFieldCreateStart(constituativeCellMLParametersFieldUserNumber,constituativeCellMLParametersField)
constituativeCellMLParametersField.VariableLabelSet(CMISS.FieldVariableTypes.U,"ConstituativeParameters")
constituativeCellML.ParametersFieldCreateFinish()

# Create the CELL intermediate field
constituativeCellMLIntermediateField = CMISS.Field()
constituativeCellML.IntermediateFieldCreateStart(constituativeCellMLIntermediateFieldUserNumber,constituativeCellMLIntermediateField)
constituativeCellMLIntermediateField.VariableLabelSet(CMISS.FieldVariableTypes.U,"ConstituativeIntermediate")
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

#Set x=0 nodes to built-in
for zNode in range(1, numberZElements+2):
    for yNode in range(1, numberYElements+2):
        nodeNumber = 1+(yNode-1)*(numberXElements+1)+(zNode-1)*(numberXElements+1)*(numberYElements+1)
        #print 'Setting boundary condition at node ',nodeNumber
        boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,nodeNumber,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,nodeNumber,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
        boundaryConditions.AddNode(dependentField,CMISS.FieldVariableTypes.U,1,1,nodeNumber,3,CMISS.BoundaryConditionsTypes.FIXED,0.0)
 
nonlinearEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Export results
fields = CMISS.Fields()
fields.CreateRegion(region)
fields.NodesExport("CantileverGrowth","FORTRAN")
fields.ElementsExport("CantileverGrowth","FORTRAN")
fields.Finalise()

