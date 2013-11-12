#!/usr/bin/env python

#> \file
#> \author David Ladd
#> \brief This is an example script to solve a Navier-Stokes lid driven cavity benchmark problem using openCMISS calls in python.
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

#> \example FluidMechanics/NavierStokes/LidDrivenCavity/LidDriven.py
## Example script to solve a Navier-Stokes lid driven cavity benchmark problem using openCMISS calls in python.
## \par Latest Builds:
#<


# Add Python bindings directory to PATH
import sys, os
sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

# Intialise OpenCMISS
from opencmiss import CMISS
import numpy
import gzip
import pylab
import time
import matplotlib.pyplot as plt


#      Square Lid-Driven Cavity
#
#              v=1
#           >>>>>>>>>>
#          |          |
#          |          |
#     v=0  |          |  v=0
#          |          |
#          |          |
#          ------------
#              v=0


def LidDriven(numberOfElements,cavityDimensions,lidVelocity,viscosity,density,outputFilename):
    """ Sets up the lid driven cavity problem and 
    solves with the provided parameter values"""

    CMISS.DiagnosticsSetOn(CMISS.DiagnosticTypes.IN,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

    # Get the computational nodes information
    numberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet()
    computationalNodeNumber = CMISS.ComputationalNodeNumberGet()

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
     equationsSetUserNumber,
     problemUserNumber) = range(1,14)

    # Creation a RC coordinate system
    coordinateSystem = CMISS.CoordinateSystem()
    coordinateSystem.CreateStart(coordinateSystemUserNumber)
    coordinateSystem.dimension = 2
    coordinateSystem.CreateFinish()

    # Create a region
    region = CMISS.Region()
    region.CreateStart(regionUserNumber,CMISS.WorldRegion)
    region.label = "Cavity"
    region.coordinateSystem = coordinateSystem
    region.CreateFinish()

    # Create a biquadratic lagrange basis
    quadraticBasis = CMISS.Basis()
    quadraticBasis.CreateStart(quadraticBasisUserNumber)
    quadraticBasis.type = CMISS.BasisTypes.LAGRANGE_HERMITE_TP
    quadraticBasis.numberOfXi = 2
    quadraticBasis.interpolationXi = [CMISS.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*2
    quadraticBasis.quadratureNumberOfGaussXi = [2]*2
    quadraticBasis.CreateFinish()

    # Create a bilinear lagrange basis
    linearBasis = CMISS.Basis()
    linearBasis.CreateStart(linearBasisUserNumber)
    linearBasis.type = CMISS.BasisTypes.LAGRANGE_HERMITE_TP
    linearBasis.numberOfXi = 2
    linearBasis.interpolationXi = [CMISS.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*2
    linearBasis.quadratureNumberOfGaussXi = [2]*2
    linearBasis.CreateFinish()

    # Create a generated mesh
    generatedMesh = CMISS.GeneratedMesh()
    generatedMesh.CreateStart(generatedMeshUserNumber,region)
    generatedMesh.type = CMISS.GeneratedMeshTypes.REGULAR
    generatedMesh.basis = [quadraticBasis,linearBasis]
    generatedMesh.extent = cavityDimensions
    generatedMesh.numberOfElements = numberOfElements

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
    geometricField.CreateFinish()

    # Set geometry from the generated mesh
    generatedMesh.GeometricParametersCalculate(geometricField)

    # Create standard Navier-Stokes equations set
    equationsSetField = CMISS.Field()
    equationsSet = CMISS.EquationsSet()
    equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
            CMISS.EquationsSetClasses.FLUID_MECHANICS,
            CMISS.EquationsSetTypes.NAVIER_STOKES_EQUATION,
            CMISS.EquationsSetSubtypes.STATIC_SUPG_NAVIER_STOKES,
            equationsSetFieldUserNumber, equationsSetField)
    equationsSet.CreateFinish()

    # Create dependent field
    dependentField = CMISS.Field()
    equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,1,1)
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,2,1)
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.U,3,2)
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN,1,1)
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN,2,1)
    dependentField.ComponentMeshComponentSet(CMISS.FieldVariableTypes.DELUDELN,3,2)
    dependentField.DOFOrderTypeSet(CMISS.FieldVariableTypes.U,CMISS.FieldDOFOrderTypes.SEPARATED)
    dependentField.DOFOrderTypeSet(CMISS.FieldVariableTypes.DELUDELN,CMISS.FieldDOFOrderTypes.SEPARATED)
    equationsSet.DependentCreateFinish()
    # Initialise dependent field
    dependentField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,0.0)

    # Create materials field
    materialsField = CMISS.Field()
    equationsSet.MaterialsCreateStart(materialsFieldUserNumber,materialsField)
    equationsSet.MaterialsCreateFinish()
    # Initialise materials field parameters
    materialsField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,1,viscosity)
    materialsField.ComponentValuesInitialiseDP(CMISS.FieldVariableTypes.U,CMISS.FieldParameterSetTypes.VALUES,2,density)

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
                             CMISS.ProblemSubTypes.STATIC_NAVIER_STOKES)
    problem.CreateFinish()

    # Create control loops
    problem.ControlLoopCreateStart()
    problem.ControlLoopCreateFinish()

    # Create problem solver
    nonlinearSolver = CMISS.Solver()
    linearSolver = CMISS.Solver()
    problem.SolversCreateStart()
    problem.SolverGet([CMISS.ControlLoopIdentifiers.NODE],1,nonlinearSolver)
    nonlinearSolver.newtonJacobianCalculationType = CMISS.JacobianCalculationTypes.EQUATIONS
    nonlinearSolver.outputType = CMISS.SolverOutputTypes.TIMING
    nonlinearSolver.newtonAbsoluteTolerance = 1.0E-8
    nonlinearSolver.newtonRelativeTolerance = 1.0E-8
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
    # Set values for boundary 
    firstNodeNumber=1
    nodes = CMISS.Nodes()
    region.NodesGet(nodes)
    boundaryTolerance = 1.0e-6
    # Currently issues with getting generated mesh surfaces through python so easier to just loop over all nodes
    for node in range(nodes.numberOfNodes):
        nodeId = node + 1
        nodeNumber = nodes.UserNumberGet(nodeId)
        nodeDomain=decomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            xLocation = geometricField.ParameterSetGetNodeDP(CMISS.FieldVariableTypes.U,
                                                             CMISS.FieldParameterSetTypes.VALUES,
                                                             1,1,nodeNumber,1)
            yLocation = geometricField.ParameterSetGetNodeDP(CMISS.FieldVariableTypes.U,
                                                             CMISS.FieldParameterSetTypes.VALUES,
                                                             1,1,nodeNumber,2)
            # rigid wall (left,right,bottom) conditions: v=0
            if (xLocation < boundaryTolerance or 
                cavityDimensions[0]-xLocation < boundaryTolerance or
                yLocation < boundaryTolerance):
                boundaryConditions.SetNode(dependentField,CMISS.FieldVariableTypes.U,1,1,nodeNumber,1,CMISS.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.SetNode(dependentField,CMISS.FieldVariableTypes.U,1,1,nodeNumber,2,CMISS.BoundaryConditionsTypes.FIXED,0.0)
            # lid (top) conditions: v=v
            elif (cavityDimensions[1]-yLocation < boundaryTolerance):
                boundaryConditions.SetNode(dependentField,CMISS.FieldVariableTypes.U,1,1,nodeNumber,1,CMISS.BoundaryConditionsTypes.FIXED,lidVelocity[0])
                boundaryConditions.SetNode(dependentField,CMISS.FieldVariableTypes.U,1,1,nodeNumber,2,CMISS.BoundaryConditionsTypes.FIXED,lidVelocity[1])
    solverEquations.BoundaryConditionsCreateFinish()

    # Solve the problem
    problem.Solve()

    print("exporting CMGUI data")
    # Export results
    fields = CMISS.Fields()
    fields.CreateRegion(region)
    fields.NodesExport(outputFilename,"FORTRAN")
    fields.ElementsExport(outputFilename,"FORTRAN")
    fields.Finalise()

    CMISS.Finalise()


# Problem defaults
dimensions = [1.0,1.0]
elementDimensions = [10,10]
ReynoldsNumbers = [100] #[100,400,1000,3200,5000]
lidVelocity = [1.0,0.0]
viscosity = 1.0
density = 1.0

for Re in ReynoldsNumbers:
    viscosity = 1.0/Re
    outputDirectory = "./output/Re" + str(Re) + "Dim" +str(elementDimensions[0])+"x" +str(elementDimensions[1]) + "/"
    if not os.path.exists(outputDirectory):
        os.makedirs(outputDirectory)
    outputFile = outputDirectory +"LidDrivenCavity"
    LidDriven(elementDimensions,dimensions,lidVelocity,viscosity,density,outputFile)






