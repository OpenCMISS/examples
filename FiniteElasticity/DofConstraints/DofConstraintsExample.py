#!/usr/bin/env python

#> \file
#> \author Adam Reeve
#> \brief This is an example script showing how to constrain DOFs to be equal in a finite deformation elasticity problem.
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

#> \example FiniteElasticity/DofConstraints/src/DofConstraintsExample.py
## Example script showing how to constrain DOFs to be equal in a finite deformation elasticity problem.
## Gravity loading of a cantilever beam is simulated where the free end of cantilever has DOFs
## constrained so that the end remains oriented in the y-z plane.
#<

from opencmiss import iron

# Problem parameters:
density = 9.0e-4  # in g mm^-3
gravity = [0.0, 0.0, -9.81]  # in m s^-2

numberGlobalElements = [12, 8, 8]
dimensions = [60.0, 40.0, 40.0]
numberOfXi = 3

numberOfLoadIncrements = 3

constitutiveRelation = iron.EquationsSetSubtypes.MOONEY_RIVLIN
c0, c1 = 2.0, 1.0
constitutiveParameters = [c0, c1]
initialHydrostaticPressure = -c0 - 2.0 * c1

# User numbers for identifying OpenCMISS-Iron objects:
coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
generatedMeshUserNumber = 1
meshUserNumber = 1
decompositionUserNumber = 1
equationsSetUserNumber = 1
(geometricFieldUserNumber,
    materialFieldUserNumber,
    dependentFieldUserNumber,
    equationsSetFieldUserNumber,
    sourceFieldUserNumber) = range(1, 6)
problemUserNumber = 1

# Get the number of computational nodes and this computational node number
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Create a 3D rectangular cartesian coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = iron.Region()
region.CreateStart(regionUserNumber, iron.WorldRegion)
region.LabelSet("Region")
region.CoordinateSystemSet(coordinateSystem)
region.CreateFinish()

# Define basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.NumberOfXiSet(numberOfXi)
basis.InterpolationXiSet([
        iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE] * numberOfXi)
basis.QuadratureNumberOfGaussXiSet([2] * numberOfXi)
basis.CreateFinish()

# Start the creation of a generated mesh in the region
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber, region)
generatedMesh.TypeSet(iron.GeneratedMeshTypes.REGULAR)
generatedMesh.BasisSet([basis])
generatedMesh.ExtentSet(dimensions)
generatedMesh.NumberOfElementsSet(numberGlobalElements)
mesh = iron.Mesh()
generatedMesh.CreateFinish(meshUserNumber, mesh)

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U, "Geometry")
geometricField.CreateFinish()

# Update the geometric field parameters from generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create the equations_set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
    iron.EquationsSetTypes.FINITE_ELASTICITY,
    constitutiveRelation]
equationsSet.CreateStart(equationsSetUserNumber, region, geometricField,
                         equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create default materials field
materialField = iron.Field()
equationsSet.MaterialsCreateStart(materialFieldUserNumber, materialField)
equationsSet.MaterialsCreateFinish()

# Create default dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
equationsSet.DependentCreateFinish()

# Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
for component in range(1, 4):
    iron.Field.ParametersToFieldParametersComponentCopy(
        geometricField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, component,
        dependentField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, component)
iron.Field.ComponentValuesInitialiseDP(
    dependentField, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 4, initialHydrostaticPressure)

# Set constitutive parameters
for component, parameter in enumerate(constitutiveParameters, 1):
    iron.Field.ComponentValuesInitialiseDP(
        materialField,iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES,
        component, parameter)

materialField.ComponentValuesInitialise(
    iron.FieldVariableTypes.V, iron.FieldParameterSetTypes.VALUES, 1, density)

#Create the source field with the gravity vector
sourceField = iron.Field()
equationsSet.SourceCreateStart(sourceFieldUserNumber, sourceField)
equationsSet.SourceCreateFinish()

#Set the gravity vector component values
for component in range(1, 4):
    sourceField.ComponentValuesInitialiseDP(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, component, gravity[component - 1])

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
equationsSet.EquationsCreateFinish()

# Define the problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.ELASTICITY,
        iron.ProblemTypes.FINITE_ELASTICITY,
        iron.ProblemSubtypes.NONE]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create the problem control loop
problem.ControlLoopCreateStart()
controlLoop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], controlLoop)
controlLoop.MaximumIterationsSet(numberOfLoadIncrements)
problem.ControlLoopCreateFinish()

# Create problem solver
nonLinearSolver = iron.Solver()
linearSolver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, nonLinearSolver)
nonLinearSolver.outputType = iron.SolverOutputTypes.PROGRESS
nonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS)
nonLinearSolver.NewtonLinearSolverGet(linearSolver)
linearSolver.linearType = iron.LinearSolverTypes.DIRECT
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Prescribe boundary conditions
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

nodes = iron.Nodes()
region.NodesGet(nodes)
eps = 1.0e-10
constrainedNodes = set()
for node in range(1, nodes.NumberOfNodesGet() + 1):
    position = [geometricField.ParameterSetGetNode(
                iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
                1, 1, node, component)
            for component in range(1, 4)]
    # Fix x=0 face
    if abs(position[0]) < eps:
        version = 1
        derivative = 1
        for component in range(1, 4):
            boundaryConditions.AddNode(
                    dependentField, iron.FieldVariableTypes.U,
                    version, derivative, node, component,
                    iron.BoundaryConditionsTypes.FIXED, 0.0)
    # Find nodes to constrain:
    if abs(position[0] - dimensions[0]) < eps:
        constrainedNodes.add(node)

# Constrain nodes at max x end to have the same x component value
version = 1
derivative = 1
component = 1
boundaryConditions.ConstrainNodeDofsEqual(
        dependentField, iron.FieldVariableTypes.U,
        version, derivative, component,
        list(constrainedNodes))

solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Export results
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("Cantilever", "FORTRAN")
fields.ElementsExport("Cantilever", "FORTRAN")
fields.Finalise()
