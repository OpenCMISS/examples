#!/usr/bin/env python

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

CMISS.DiagnosticsSetOn(CMISS.InDiagType,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

# Get the computational nodes information
numberOfComputationalNodes = CMISS.ComputationalNumberOfNodesGet()
computationalNodeNumber = CMISS.ComputationalNodeNumberGet()

# Creation a RC coordinate system
coordinateSystem = CMISS.CoordinateSystem(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

# Create a region
region = CMISS.Region(regionUserNumber,CMISS.WorldRegion)
region.label = "LaplaceRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Create a tri-linear lagrange basis
basis = CMISS.Basis(basisUserNumber)
basis.type = CMISS.BasisLagrangeHermiteTPType
basis.numberOfXi = 3
basis.interpolationXi = [CMISS.BasisLinearLagrangeInterpolation]*3
basis.quadratureNumberOfGaussXi = [2]*3
basis.CreateFinish()

# Create a generated mesh
generatedMesh = CMISS.GeneratedMesh(generatedMeshUserNumber,region)
generatedMesh.type = CMISS.GeneratedMeshRegularMeshType
generatedMesh.basis = [basis]
generatedMesh.extent = [width,height,length]
generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements,numberGlobalZElements]

mesh = CMISS.Mesh() #Create null mesh type by not passing any arguments
generatedMesh.CreateFinish(meshUserNumber,mesh)

# Create a decomposition for the mesh
decomposition = CMISS.Decomposition(decompositionUserNumber,mesh)
decomposition.type = CMISS.DecompositionCalculatedType
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = CMISS.Field(geometricFieldUserNumber,region)
geometricField.meshDecomposition = decomposition
geometricField.ComponentMeshComponentSet(CMISS.FieldUVariableType,1,1)
geometricField.ComponentMeshComponentSet(CMISS.FieldUVariableType,2,1)
geometricField.ComponentMeshComponentSet(CMISS.FieldUVariableType,3,1)
geometricField.CreateFinish()

# Set geometry from the generated mesh
CMISS.GeneratedMeshGeometricParametersCalculate(geometricField,generatedMesh)

# Create standard Laplace equations set
equationsSetField = CMISS.Field()
equationsSet = CMISS.EquationsSet(equationsSetUserNumber,region,geometricField,CMISS.EquationsSetClassicalFieldClass,
    CMISS.EquationsSetLaplaceEquationType,CMISS.EquationsSetStandardLaplaceSubtype,equationsSetFieldUserNumber,
    equationsSetField)
equationsSet.CreateFinish()

# Create dependent field
dependentField = CMISS.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.DOFOrderTypeSet(CMISS.FieldUVariableType,CMISS.FieldSeparatedComponentDOFOrder)
dependentField.DOFOrderTypeSet(CMISS.FieldDelUDelNVariableType,CMISS.FieldSeparatedComponentDOFOrder)
equationsSet.DependentCreateFinish()

# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(CMISS.FieldUVariableType,CMISS.FieldValuesSetType,1,0.5)

# Create equations
equations = CMISS.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityType = CMISS.EquationsSparseMatrices
equations.OutputType = CMISS.EquationsNoOutput
equationsSet.EquationsCreateFinish()

# Create Laplace problem
problem = CMISS.Problem(problemUserNumber)
problem.SpecificationSet(CMISS.ProblemClassicalFieldClass,CMISS.ProblemLaplaceEquationType,CMISS.ProblemStandardLaplaceSubtype)
problem.CreateFinish()

# Create control loops
problem.ControlLoopCreateStart()
problem.ControlLoopCreateFinish()

# Create problem solver
solver = CMISS.Solver()
problem.SolversCreateStart()
problem.SolverGet([CMISS.ControlLoopNode],1,solver)
solver.OutputType = CMISS.SolverSolverOutput
solver.LinearType = CMISS.SolverLinearIterativeSolveType
solver.LinearIterativeAbsoluteTolerance = 1.0E-12
solver.LinearIterativeRelativeTolerance = 1.0E-12
problem.SolversCreateFinish()

# Create solver equations and add equations set to solver equations
solver = CMISS.Solver()
solverEquations = CMISS.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([CMISS.ControlLoopNode],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.SparsityType = CMISS.SolverEquationsSparseMatrices
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

# Create boundary conditions and set first and last nodes to 0.0 and 1.0
boundaryConditions = CMISS.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)
firstNodeNumber=1
nodes = CMISS.Nodes()
region.NodesGet(nodes)
lastNodeNumber = nodes.NumberOfNodes
firstNodeDomain = decomposition.NodeDomainGet(firstNodeNumber,1)
lastNodeDomain = decomposition.NodeDomainGet(lastNodeNumber,1)
if firstNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField,CMISS.FieldUVariableType,1,1,firstNodeNumber,1,CMISS.BoundaryConditionFixed,0.0)
if lastNodeDomain == computationalNodeNumber:
    boundaryConditions.SetNode(dependentField,CMISS.FieldUVariableType,1,1,lastNodeNumber,1,CMISS.BoundaryConditionFixed,1.0)
solverEquations.BoundaryConditionsCreateFinish()

# Solve the problem
problem.Solve()

# Export results
fields = CMISS.Fields()
CMISS.FieldsTypeCreateRegion(region,fields)
CMISS.FieldIONodesExport(fields,"Laplace","FORTRAN")
CMISS.FieldIOElementsExport(fields,"Laplace","FORTRAN")
CMISS.FieldsTypeFinalise(fields)

CMISS.Finalise()
