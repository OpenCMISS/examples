#!/usr/bin/env python
#> \file
#> \author Adam Reeve
#> \brief Test Neumann boundary conditions with Laplace's equation
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
#> Contributor(s): Adam Reeve
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

#> \example ClassicalField/Laplace/NeumannConditions/NeumannConditionsExample.py
## Test Neumann boundary conditions with Laplace's equation on quadratic elements
#<

import sys
import os
from optparse import OptionParser

# Add OpenCMISS python bindings directory to path so we can import it
sys.path.append(os.sep.join((
    os.environ['OPENCMISS_ROOT'], 'cm', 'bindings', 'python')))
from opencmiss import iron

# Read in number of elements
parser = OptionParser()
parser.add_option("-x", "--x-elements", dest="x",type="int", default=3,
        help="Number of x elements")
parser.add_option("-y", "--y-elements", dest="y", type="int", default=3,
        help="Number of y elements")
parser.add_option("-z", "--z-elements", dest="z", type="int", default=0,
        help="Number of z elements, defaults to zero for 2D")
parser.add_option("-i", "--interpolation", dest="interpolation",
        choices=("quadratic-lagrange", "cubic-hermite"),
        default="quadratic-lagrange",
        help="Basis interpolation type")
(opts, args) = parser.parse_args()
if opts.x <= 0:
    raise ValueError("Number of x elements must be > 0")
if opts.y < 0:
    raise ValueError("Number of y elements must be >= 0")
if opts.z < 0:
    raise ValueError("Number of z elements must be >= 0")
numberGlobalElements = (opts.x, opts.y, opts.z)
numberOfXi = sum(ne > 0 for ne in numberGlobalElements)

# Set problem parameters
width = 1.0
if numberOfXi > 1:
    length = 3.0
else:
    length = 0.0
if numberOfXi > 2:
    height = 1.0
else:
    height = 0.0

if opts.interpolation == "quadratic-lagrange":
    interpolationType = iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE
    numGauss = 3
    hasDerivatives = False
elif opts.interpolation == "cubic-hermite":
    interpolationType = iron.BasisInterpolationSpecifications.CUBIC_HERMITE
    numGauss = 4
    hasDerivatives = True

(coordinateSystemUserNumber,
    regionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    meshUserNumber,
    decompositionUserNumber,
    geometricFieldUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    equationsSetUserNumber,
    problemUserNumber) = range(1, 12)

# Output all diagnostics in Neumann integration routine
iron.DiagnosticsSetOn(
        iron.DiagnosticTypes.IN, [1,2,3,4,5], "diagnostics",
        ["BoundaryConditions_NeumannIntegrate",
        "BoundaryConditions_NeumannMatricesInitialise"])

# Get the computational nodes information
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

# Creation a rectangular cartesian coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region in the world region
region = iron.Region()
region.CreateStart(regionUserNumber, iron.WorldRegion)
region.label = "LaplaceRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create a tri-quadratic lagrange basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = numberOfXi
basis.interpolationXi = [interpolationType] * numberOfXi
basis.quadratureNumberOfGaussXi = [numGauss] * numberOfXi
basis.CreateFinish()

# Create a generated regular mesh
generatedMesh = iron.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber, region)
generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [width, length, height]
generatedMesh.numberOfElements = [ne for ne in numberGlobalElements if ne > 0]

mesh = iron.Mesh()
generatedMesh.CreateFinish(meshUserNumber, mesh)

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber, mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.calculateFaces = True
decomposition.calculateLines = True
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber, region)
geometricField.meshDecomposition = decomposition
geometricField.scalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create standard Laplace equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
        iron.EquationsSetTypes.LAPLACE_EQUATION,
        iron.EquationsSetSubtypes.STANDARD_LAPLACE]
equationsSet.CreateStart(equationsSetUserNumber, region, geometricField,
        equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
dependentField.scalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
equationsSet.DependentCreateFinish()

# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
        1, 0.5)

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equationsSet.EquationsCreateFinish()

# Create Laplace problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
        iron.ProblemTypes.LAPLACE_EQUATION,
        iron.ProblemSubtypes.STANDARD_LAPLACE]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.outputType = iron.SolverOutputTypes.PROGRESS
solver.linearType = iron.LinearSolverTypes.ITERATIVE
solver.linearIterativeAbsoluteTolerance = 1.0E-12
solver.linearIterativeRelativeTolerance = 1.0E-12
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Set boundary conditions
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
boundaryConditions.neumannSparsityType = iron.BoundaryConditionSparsityTypes.SPARSE
nodes = iron.Nodes()
region.NodesGet(nodes)
for node in range(1, nodes.numberOfNodes + 1):
    nodeDomain = decomposition.NodeDomainGet(node, 1)
    if nodeDomain != computationalNodeNumber:
        continue
    # get x, y and z positions at node
    position = [
        geometricField.ParameterSetGetNodeDP(
            iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
            1, 1, node, component + 1)
        for component in range(numberOfXi)]
    tol = 1.0e-12
    if abs(position[0] - width) < tol:
        # Fix right nodes at zero
        boundaryConditions.SetNode(dependentField, iron.FieldVariableTypes.U,
                1, 1, node, 1, iron.BoundaryConditionsTypes.FIXED, 0.0)
        if numberOfXi > 1 and hasDerivatives:
            boundaryConditions.SetNode(dependentField,
                    iron.FieldVariableTypes.U, 1, 3, node, 1,
                    iron.BoundaryConditionsTypes.FIXED, 0.0)
        if numberOfXi > 2 and hasDerivatives:
            boundaryConditions.SetNode(dependentField,
                    iron.FieldVariableTypes.U, 1, 5, node, 1,
                    iron.BoundaryConditionsTypes.FIXED, 0.0)
            boundaryConditions.SetNode(dependentField,
                    iron.FieldVariableTypes.U, 1, 7, node, 1,
                    iron.BoundaryConditionsTypes.FIXED, 0.0)
    elif abs(position[0]) < tol:
        # Set Neumann condition of 1 at left side
        boundaryConditions.SetNode(dependentField,
                iron.FieldVariableTypes.DELUDELN, 1, 1, node, 1,
                iron.BoundaryConditionsTypes.NEUMANN_POINT, 1.0)
    elif ((numberOfXi > 1 and
            (abs(position[1]) < tol or abs(position[1] - length) < tol)) or
            (numberOfXi > 2 and
            (abs(position[2]) < tol or abs(position[2] - height) < tol))):
        # Set "integrated only" free conditions at top and bottom, front and back
        # (excluding nodes at left corners) to prevent integrated terms from
        # point values being added here
        boundaryConditions.SetNode(dependentField,
                iron.FieldVariableTypes.DELUDELN, 1, 1, node, 1,
                iron.BoundaryConditionsTypes.NEUMANN_INTEGRATED_ONLY, 0.0)
solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Export results
if numberOfComputationalNodes == 1 and not hasDerivatives:
    # Use FieldML output if we can, but it doesn't support
    # parallel output or cubic Hermite interpolation
    baseName = "laplace"
    dataFormat = "PLAIN_TEXT"
    fml = iron.FieldMLIO()
    fml.OutputCreate(mesh, "", baseName, dataFormat)
    fml.OutputAddFieldNoType(baseName + ".geometric", dataFormat, geometricField,
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    fml.OutputAddFieldNoType(baseName + ".phi", dataFormat, dependentField,
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
    fml.OutputWrite("LaplaceExample.xml")
    fml.Finalise()
else:
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport("Laplace","FORTRAN")
    fields.ElementsExport("Laplace","FORTRAN")
    fields.Finalise()

iron.Finalise()
