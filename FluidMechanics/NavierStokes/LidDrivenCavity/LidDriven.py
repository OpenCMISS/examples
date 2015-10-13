#!/usr/bin/env python

#> \file
#> \author David Ladd
#> \brief This is an example script to solve a Navier-Stokes lid driven cavity benchmark problem using OpenCMISS calls in python.
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

#> \example FluidMechanics/NavierStokes/LidDrivenCavity/LidDriven.py
## Example script to solve a Navier-Stokes lid driven cavity benchmark problem using OpenCMISS calls in python.
## \par Latest Builds:
#<


# Add Python bindings directory to PATH
import sys, os
sys.path.append(os.sep.join((os.environ['OPENCMISS_ROOT'],'cm','bindings','python')))

import math
import time

# Intialise OpenCMISS
from opencmiss import iron

# Get the computational nodes information
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()

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
 analyticFieldUserNumber,
 problemUserNumber) = range(1,15)

# Creation a RC coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 2
coordinateSystem.CreateFinish()

# Create a region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "Cavity"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create a biquadratic lagrange basis
quadraticBasis = iron.Basis()
quadraticBasis.CreateStart(quadraticBasisUserNumber)
quadraticBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
quadraticBasis.numberOfXi = 2
quadraticBasis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*2
quadraticBasis.quadratureNumberOfGaussXi = [3]*2
quadraticBasis.CreateFinish()

# Create a bilinear lagrange basis
linearBasis = iron.Basis()
linearBasis.CreateStart(linearBasisUserNumber)
linearBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
linearBasis.numberOfXi = 2
linearBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*2
linearBasis.quadratureNumberOfGaussXi = [3]*2
linearBasis.CreateFinish()


def LidDriven(numberOfElements,cavityDimensions,lidVelocity,viscosity,density,
              outputFilename,transient,RBS,fdJacobian,analytic,basisList):
    """ Sets up the lid driven cavity problem and solves with the provided parameter values

          Square Lid-Driven Cavity

                  v=1
               >>>>>>>>>>
             1|          |
              |          |
         v=0  |          |  v=0
              |          |
              |          |
              ------------
             0    v=0    1
    """

    # Create a generated mesh
    generatedMesh = iron.GeneratedMesh()
    generatedMesh.CreateStart(generatedMeshUserNumber,region)
    generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
    generatedMesh.basis = basisList
    generatedMesh.extent = cavityDimensions
    generatedMesh.numberOfElements = numberOfElements

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
    geometricField.meshDecomposition = decomposition
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
    geometricField.CreateFinish()

    # Set geometry from the generated mesh
    generatedMesh.GeometricParametersCalculate(geometricField)

    # Create standard Navier-Stokes equations set
    equationsSetField = iron.Field()
    equationsSet = iron.EquationsSet()
    if RBS:
        equationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                iron.EquationsSetSubtypes.TRANSIENT_RBS_NAVIER_STOKES]
    else:
        equationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                iron.EquationsSetSubtypes.TRANSIENT_NAVIER_STOKES]
    equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
                equationsSetSpecification,equationsSetFieldUserNumber, equationsSetField)
    equationsSet.CreateFinish()

    if RBS:
        # Set max CFL number (default 1.0)
        equationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                      iron.FieldParameterSetTypes.VALUES,2,1.0E20)
        # Set time increment (default 0.0)
        equationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                      iron.FieldParameterSetTypes.VALUES,3,transient[2])
        # Set stabilisation type (default 1.0 = RBS)
        equationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                      iron.FieldParameterSetTypes.VALUES,4,1.0)

    # Create dependent field
    dependentField = iron.Field()
    equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,2)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,1,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,2,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,3,2)
    dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
    dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
    equationsSet.DependentCreateFinish()
    # Initialise dependent field
    dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)

    # Create materials field
    materialsField = iron.Field()
    equationsSet.MaterialsCreateStart(materialsFieldUserNumber,materialsField)
    equationsSet.MaterialsCreateFinish()
    # Initialise materials field parameters
    materialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,viscosity)
    materialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,density)

    # If specified, use a sinusoidal waveform to ramp up lid velocity from 0 to 1
    if analytic:
        # yOffset + amplitude*sin(frequency*time + phaseShift))))
        # Set the time it takes to ramp velocity up to full lid velocity
        rampPeriod = 10.0
        frequency = math.pi/(rampPeriod)
        amplitude = 0.5*lidVelocity[0]
        yOffset = 0.5*lidVelocity[0]
        phaseShift = -math.pi/2.0
        startSine = 0.0
        stopSine = rampPeriod
        analyticField = iron.Field()
        equationsSet.AnalyticCreateStart(iron.NavierStokesAnalyticFunctionTypes.SINUSOID,analyticFieldUserNumber,analyticField)
        equationsSet.AnalyticCreateFinish()
        analyticParameters = [1.0,0.0,0.0,0.0,amplitude,yOffset,frequency,phaseShift,startSine,stopSine]

    # Create equations
    equations = iron.Equations()
    equationsSet.EquationsCreateStart(equations)
    equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
    equations.outputType = iron.EquationsOutputTypes.NONE
    equationsSet.EquationsCreateFinish()

    # Create Navier-Stokes problem
    problem = iron.Problem()
    if RBS:
        problemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
                                iron.ProblemTypes.NAVIER_STOKES_EQUATION,
                                iron.ProblemSubtypes.TRANSIENT_RBS_NAVIER_STOKES]
    else:
        problemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
                                iron.ProblemTypes.NAVIER_STOKES_EQUATION,
                                iron.ProblemSubtypes.TRANSIENT_NAVIER_STOKES]
    problem.CreateStart(problemUserNumber,problemSpecification)
    problem.CreateFinish()

    # Create control loops
    problem.ControlLoopCreateStart()
    controlLoop = iron.ControlLoop()
    problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],controlLoop)
    controlLoop.TimesSet(transient[0],transient[1],transient[2])
    controlLoop.TimeOutputSet(transient[3])
    problem.ControlLoopCreateFinish()

    # Create problem solver
    dynamicSolver = iron.Solver()
    problem.SolversCreateStart()
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,dynamicSolver)
    dynamicSolver.outputType = iron.SolverOutputTypes.NONE
    dynamicSolver.dynamicTheta = [0.5]
    nonlinearSolver = iron.Solver()
    dynamicSolver.DynamicNonlinearSolverGet(nonlinearSolver)
    if fdJacobian:
        nonlinearSolver.newtonJacobianCalculationType = iron.JacobianCalculationTypes.FD
    else:
        nonlinearSolver.newtonJacobianCalculationType = iron.JacobianCalculationTypes.EQUATIONS
    nonlinearSolver.outputType = iron.SolverOutputTypes.NONE
    nonlinearSolver.newtonAbsoluteTolerance = 1.0E-8
    nonlinearSolver.newtonRelativeTolerance = 1.0E-9
    nonlinearSolver.newtonSolutionTolerance = 1.0E-9
    nonlinearSolver.newtonMaximumFunctionEvaluations = 10000
    nonlinearSolver.newtonLineSearchType = iron.NewtonLineSearchTypes.QUADRATIC
    linearSolver = iron.Solver()
    nonlinearSolver.NewtonLinearSolverGet(linearSolver)
    linearSolver.outputType = iron.SolverOutputTypes.NONE
    linearSolver.linearType = iron.LinearSolverTypes.DIRECT
    linearSolver.libraryType = iron.SolverLibraries.MUMPS
    problem.SolversCreateFinish()

    # Create solver equations and add equations set to solver equations
    solver = iron.Solver()
    solverEquations = iron.SolverEquations()
    problem.SolverEquationsCreateStart()
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
    solver.SolverEquationsGet(solverEquations)
    solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
    equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
    problem.SolverEquationsCreateFinish()

    # Create boundary conditions
    boundaryConditions = iron.BoundaryConditions()
    solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
    nodes = iron.Nodes()
    region.NodesGet(nodes)
    print("Total # of nodes: " + str(nodes.numberOfNodes))
    print("Analytic Parameters: " + str(analyticParameters))
    boundaryTolerance = 1.0e-6

    # Currently issues with getting generated mesh surfaces through python so easier to just loop over all nodes
    for node in range(nodes.numberOfNodes):
        nodeId = node + 1
        nodeNumber = nodes.UserNumberGet(nodeId)
        # print('node number: '+ str(nodeNumber))
        # Velocity nodes
        nodeDomain=decomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            xLocation = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES,
                                                             1,1,nodeNumber,1)
            yLocation = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,
                                                             iron.FieldParameterSetTypes.VALUES,
                                                             1,1,nodeNumber,2)
            # rigid wall (left,right,bottom) conditions: v=0
            if (xLocation < boundaryTolerance or 
                cavityDimensions[0]-xLocation < boundaryTolerance or
                yLocation < boundaryTolerance):
                boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
            # lid (top) conditions: v=v
            elif (cavityDimensions[1]-yLocation < boundaryTolerance):
                if not (xLocation < boundaryTolerance or 
                        cavityDimensions[0]-xLocation < boundaryTolerance):
                    if analytic:
                        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                        # Set analytic parameters
                        parameterNumber = 0
                        for parameter in analyticParameters:
                            parameterNumber += 1
                            analyticField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,
                                                                   nodeNumber,parameterNumber,parameter)

                    else:
                        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,lidVelocity[0])
                        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,lidVelocity[1])


    # Pressure node
    nodeNumber = 1 
    nodeDomain=decomposition.NodeDomainGet(nodeNumber,2)
    if (nodeDomain == computationalNodeNumber):
        # bottom left node - reference pressure: p=0
        boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)
        print('pressure node: '+str(nodeNumber))

    solverEquations.BoundaryConditionsCreateFinish()

    # Solve the problem
    print("solving...")
    problem.Solve()

    print("exporting CMGUI data")
    # Export results
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport(outputFilename,"FORTRAN")
    fields.ElementsExport(outputFilename,"FORTRAN")
    fields.Finalise()

    # Clear fields so can run in batch mode on this region
    generatedMesh.Destroy()
    nodes.Destroy()
    mesh.Destroy()
    geometricField.Destroy()
    if RBS:
        equationsSetField.Destroy()
    if analytic:
        analyticField.Destroy()
    dependentField.Destroy()
    materialsField.Destroy()
    equationsSet.Destroy()
    problem.Destroy()


#DOC-START Control Panel
#==========================================================
# P r o b l e m     C o n t r o l
#==========================================================

dimensions = [1.0,1.0]
elementResolutions = [10]
ReynoldsNumbers = [1000]
lidVelocity = [1.0,0.0]
density = 1.0
basisList = [quadraticBasis,linearBasis]
RBSTypes = [True,False]
timeIncrement = 0.1
stopTime = 100.001
fdJacobian = False
analyticLidVelocity = True
# Note: viscosity will be calculated based on specified Reynolds number

#==========================================================
#DOC-END Control Panel

#Check for command line arguments- used for nightly test on GFEM & RBS code
if len(sys.argv) > 1:
    if len(sys.argv) > 4:
        sys.exit('Error: too many arguments- currently only accepting 3 options: "RBS" or "GFEM", time increment, stop time')
    if sys.argv[1] == 'RBS':
        RBSTypes = [True]
    elif sys.argv[1] == 'GFEM':
        RBSTypes = [False]
    else:
        sys.exit('Error: unknown argument for the problem type: choose "RBS" or "GFEM"')
    if len(sys.argv) > 2:
        timeIncrement = float(sys.argv[2])
    if len(sys.argv) > 3:
        stopTime = float(sys.argv[3])

runtimes = []
for elemRes in elementResolutions:
    for Re in ReynoldsNumbers:
        for RBS in RBSTypes:

            elementResolution = [elemRes,elemRes]
            viscosity = density*lidVelocity[0]/Re

            # transient parameters: startTime,stopTime,timeIncrement,outputFrequency
            transient = [0.0,stopTime,timeIncrement,100000]
            if RBS:    
                outputDirectory = "./output/Re" + str(Re) + "Elem" +str(elementResolution[0])+"x" +str(elementResolution[1]) + "_RBS/"
            else:
                outputDirectory = "./output/Re" + str(Re) + "Elem" +str(elementResolution[0])+"x" +str(elementResolution[1]) + "_GFEM/"

            # Create a results directory
            try:
                os.makedirs(outputDirectory)
            except OSError, e:
                if e.errno != 17:
                    raise   

            outputFile = outputDirectory +"LidDrivenCavity"
            # Display solve info
            if computationalNodeNumber == 0:
                print('Starting solve with parameters: ')
                print('---------------------------------------------------------------------')
                print('    Reynolds number (Re): ' + str(Re))
                print('    ElementResolution   : ' + str(elementResolution))
                print('    RBS                 : ' + str(RBS))
                print('    Transient parameters: ' + str(transient))
                print('    FD Jacobian: ' + str(fdJacobian))
            start = time.time()
            LidDriven(elementResolution,dimensions,lidVelocity,viscosity,density,
                      outputFile,transient,RBS,fdJacobian,analyticLidVelocity,basisList)
            end = time.time()
            runtime = end - start
            runtimes.append(runtime)
            if computationalNodeNumber == 0:
                runInfo = ''
                runInfo += 'Successfully finished solve with parameters: ' + '\n'
                runInfo += '---------------------------------------------------------------------' + '\n'
                runInfo += '    Reynolds number (Re): ' + str(Re) + '\n'
                runInfo += '    ElementResolution   : ' + str(elementResolution) + '\n'
                runInfo += '    RBS                 : ' + str(RBS) + '\n'
                runInfo += '    Transient parameters: ' + str(transient) + '\n'
                runInfo += '    FD Jacobian         : ' + str(fdJacobian) + '\n'
                runInfo += '    output results to   : ' + outputFile + '\n'
                runInfo += '    Runtime             : ' + str(runtime) + '\n'
                print(runInfo)
                f = open(outputDirectory+"RunInfo.txt","w")
                f.write(runInfo)
                f.close() 

iron.Finalise()





