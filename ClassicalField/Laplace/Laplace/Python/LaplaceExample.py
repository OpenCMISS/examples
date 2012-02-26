#!/usr/bin/env python

#> \file
#> \author Chris Bradley
#> \brief This is an example script to solve a Laplace problem using openCMISS calls in python.
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

#> \example ClassicalField/Laplace/LaplacePy/LaplaceExample.py
## Example script to solve a Laplace problem using openCMISS calls in python.
## \par Latest Builds:
## \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/Laplace/LaplacePy/build-intel'>Linux Intel Build</a>
## \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/Laplace/LaplacePy/build-gnu'>Linux GNU Build</a>
#<


# Add Python bindings directory to PATH
import sys, os
sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

# Intialise OpenCMISS
from opencmiss import CMISS

# Set problem parameters
height = 1.0
width = 2.0
length = 3.0

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
    problemUserNumber) = range(1,12)

numberGlobalXElements = 5
numberGlobalYElements = 5
numberGlobalZElements = 5

CMISS.DiagnosticsSetOn(CMISS.DiagnosticTypes.IN,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

# Get the computational nodes information
numberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet()
computationalNodeNumber = CMISS.ComputationalNodeNumberGet()

# Creation a RC coordinate system
coordinateSystem = CMISS.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = CMISS.Region()
region.CreateStart(regionUserNumber,CMISS.WorldRegion)
region.label = "LaplaceRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create a tri-linear lagrange basis
basis = CMISS.Basis()
basis.CreateStart(basisUserNumber)
basis.type = CMISS.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = 3
basis.interpolationXi = [CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
basis.quadratureNumberOfGaussXi = [2]*3
basis.CreateFinish()

# Create a generated mesh
generatedMesh = CMISS.GeneratedMesh()
generatedMesh.CreateStart(generatedMeshUserNumber,region)
generatedMesh.type = CMISS.GeneratedMeshTypes.REGULAR
generatedMesh.basis = [basis]
generatedMesh.extent = [width,height,length]
generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements,numberGlobalZElements]

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
geometricField.meshDecomposition = decomposition
geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,3,1)
geometricField.CreateFinish()

# Set geometry from the generated mesh
generatedMesh.GeometricParametersCalculate(geometricField)

# Create standard Laplace equations set
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet()
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        CMISS.EquationsSetClasses.CLASSICAL_FIELD,
        CMISS.EquationsSetTypes.LAPLACE_EQUATION,
        CMISS.EquationsSetSubtypes.STANDARD_LAPLACE,
        equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

# Create dependent field
dependentField = CMISS.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.DOFOrderTypeSet(CMISS.FieldVariableTypes.U,CMISS.FieldDOFOrderTypes.SEPARATED)
dependentField.DOFOrderTypeSet(CMISS.FieldVariableTypes.DELUDELN,CMISS.FieldDOFOrderTypes.SEPARATED)
equationsSet.DependentCreateFinish()

# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,0.5)

# Create equations
equations = CMISS.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = CMISS.EquationsSparsityTypes.SPARSE
equations.outputType = CMISS.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

# Create Laplace problem
problem = CMISS.Problem()
problem.CreateStart(problemUserNumber)
problem.SpecificationSet(CMISS.ProblemClasses.CLASSICAL_FIELD,
        CMISS.ProblemTypes.LAPLACE_EQUATION,
        CMISS.ProblemSubTypes.STANDARD_LAPLACE)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = CMISS.Solver()
problem.SolversCreateStart()
problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],1,solver)
solver.outputType = CMISS.SolverOutputTypes.SOLVER
solver.linearType = CMISS.LinearSolverTypes.ITERATIVE
solver.linearIterativeAbsoluteTolerance = 1.0E-12
solver.linearIterativeRelativeTolerance = 1.0E-12
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

# Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = CMISS.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
firstNodeNumber=1
nodes = CMISS.Nodes()
region.NodesGet(nodes)
lastNodeNumber = nodes.numberOfNodes
firstNodeDomain = decomposition.NodeDomainGet(firstNodeNumber,1)
lastNodeDomain = decomposition.NodeDomainGet(lastNodeNumber,1)
if firstNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField,CMISS.FieldVariableTypes.U,1,1,firstNodeNumber,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
if lastNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField,CMISS.FieldVariableTypes.U,1,1,lastNodeNumber,1,CMISS.BoundaryConditionsTypes.FIXED,1.0)
solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Export results
baseName = "laplace"
dataFormat = "PLAIN_TEXT"
fml = CMISS.FieldMLIO()
fml.OutputCreate(mesh, "", baseName, dataFormat)
fml.OutputAddFieldNoType(baseName+".geometric", dataFormat, geometricField,
    CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES)
fml.OutputAddFieldNoType(baseName+".phi", dataFormat, dependentField,
    CMISS.FieldVariableTypes.U, CMISS.FieldParameterSetTypes.VALUES)
fml.OutputWrite("LaplaceExample.xml")
fml.Finalise()

CMISS.Finalise()
