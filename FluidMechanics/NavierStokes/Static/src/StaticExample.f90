!> \file
!> $Id: StaticExample.f90 20 2009-10-28 20:22:52Z sebk $
!> \author Sebastian Krittian
!> \brief This is an example program to solve a static Navier-Stokes equation using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example FluidMechanics/NavierStokes/Static/src/StaticExample.f90
!! Example program to solve a static Navier-Stokes equation using OpenCMISS calls.
!! \htmlinclude FluidMechanics/NavierStokes/Static/history.html
!!
!<

!> Main program

PROGRAM NAVIERSTOKESSTATICEXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OPENCMISS
  USE FLUID_MECHANICS_IO_ROUTINES
  USE MPI
  USE FIELDML_INPUT_ROUTINES
  USE FIELDML_UTIL_ROUTINES
  USE FIELDML_OUTPUT_ROUTINES
  USE FIELDML_API

#ifdef WIN32
  USE IFQWINCMISS
#endif

  !
  !================================================================================================================================
  !

  !PROGRAM VARIABLES AND TYPES

  IMPLICIT NONE

  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberNavierStokes=6
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokes=7
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberNavierStokes=8
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberNavierStokes=9
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=10

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverNavierStokesUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesMu=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesRho=2

  INTEGER(CMISSIntg), PARAMETER :: basisNumberTrilinear=1
  INTEGER(CMISSIntg), PARAMETER :: basisNumberTriquadratic=2

  INTEGER(CMISSIntg), PARAMETER :: gaussQuadrature(3) = (/3,3,3/)

  CHARACTER(KIND=C_CHAR), PARAMETER :: NUL = C_NULL_CHAR

  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: inputFilename = "HEX-M2-V2-P1_FE.xml"

  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputDirectory = "output"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputFilename = outputDirectory//"/HEX-M2-V2-P1_FE_out.xml"

CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: basename = "static_navier_stokes"
  !Program types

  !Program variables

  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_PRESSURE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_SPACE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE
!   INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES_NAVIER_STOKES
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES

  INTEGER(CMISSIntg) :: EQUATIONS_NAVIER_STOKES_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: NODE_COUNTER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_NAVIER_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_NAVIER_STOKES

  REAL(CMISSDP) :: INITIAL_FIELD_NAVIER_STOKES(3)
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_NAVIER_STOKES(3)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE
  REAL(CMISSDP) :: MU_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: RHO_PARAM_NAVIER_STOKES

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG
  LOGICAL :: FIXED_WALL_NODES_NAVIER_STOKES_FLAG
  LOGICAL :: INLET_WALL_NODES_NAVIER_STOKES_FLAG

  !CMISS variables

  !Regions
  TYPE(CMISSRegionType) :: Region
  TYPE(CMISSRegionType) :: WorldRegion
  !Coordinate systems
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  !Nodes
  TYPE(CMISSNodesType) :: Nodes
  !Meshes
  TYPE(CMISSMeshType) :: Mesh
  !Decompositions
  TYPE(CMISSDecompositionType) :: Decomposition
  !Fields
  TYPE(CMISSFieldsType) :: Fields
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: DependentFieldNavierStokes
  TYPE(CMISSFieldType) :: MaterialsFieldNavierStokes
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsNavierStokes
  !Equations sets
  TYPE(CMISSEquationsSetType) :: EquationsSetNavierStokes
  !Equations
  TYPE(CMISSEquationsType) :: EquationsNavierStokes
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: NonlinearSolverNavierStokes
  TYPE(CMISSSolverType) :: LinearSolverNavierStokes
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsNavierStokes
  
  !FieldML parsing variables
  TYPE(FieldmlInfoType) :: fieldmlInfo, outputInfo
  
  INTEGER(CMISSIntg) :: spaceHandle, velocityHandle, pressureHandle, meshComponentCount
  
  INTEGER(CMISSIntg) :: domainHandle, nodeCount, elementCount
  INTEGER(CMISSIntg) :: coordinateType, coordinateCount, xiDimensions
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables

  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: err

#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !
  !================================================================================================================================
  !

  !PROBLEM CONTROL PANEL
  
  !Set initial values
  INITIAL_FIELD_NAVIER_STOKES(1)=0.0_CMISSDP
  INITIAL_FIELD_NAVIER_STOKES(2)=0.0_CMISSDP
  INITIAL_FIELD_NAVIER_STOKES(3)=0.0_CMISSDP
  !Set boundary conditions
  FIXED_WALL_NODES_NAVIER_STOKES_FLAG=.TRUE.
  INLET_WALL_NODES_NAVIER_STOKES_FLAG=.TRUE.
  IF(FIXED_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    NUMBER_OF_FIXED_WALL_NODES_NAVIER_STOKES=80
    ALLOCATE(FIXED_WALL_NODES_NAVIER_STOKES(NUMBER_OF_FIXED_WALL_NODES_NAVIER_STOKES))
    FIXED_WALL_NODES_NAVIER_STOKES=(/1,2,3,4,5,7,9,10,11,12,13,14,17,20,24,28,29,30,31,32,33,34,35,37,39, & 
    & 41,44,46,47,48,50,51,52,53,54,57,60,64,65,66,67,68,70,72,74,76,77,78,79,80,83,86, & 
    & 89,90,91,92,93,94,95,97,99,101,102,103,104,105,106,107,108,111,114,115,116,117,118, & 
    & 120,122,123,124,125/)
  ENDIF
  IF(INLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES=9
    ALLOCATE(INLET_WALL_NODES_NAVIER_STOKES(NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES))
    INLET_WALL_NODES_NAVIER_STOKES=(/6,15,16,23,36,42,81,82,96/)
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_NAVIER_STOKES(1)=0.0_CMISSDP
    BOUNDARY_CONDITIONS_NAVIER_STOKES(2)=1.0_CMISSDP
    BOUNDARY_CONDITIONS_NAVIER_STOKES(3)=0.0_CMISSDP
  ENDIF
  !Set material parameters
  MU_PARAM_NAVIER_STOKES=1.0_CMISSDP
  RHO_PARAM_NAVIER_STOKES=1.0_CMISSDP
  !Set interpolation parameters
  BASIS_XI_GAUSS_SPACE=3
  BASIS_XI_GAUSS_VELOCITY=3
  BASIS_XI_GAUSS_PRESSURE=3
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISSSolverNoOutput
  NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISSSolverNoOutput
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_NAVIER_STOKES_OUTPUT=CMISSEquationsNoOutput
  !Set solver parameters
  LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E20 !default: 1.0E5
  MAXIMUM_ITERATIONS=100000 !default: 100000
  RESTART_VALUE=3000 !default: 30
  LINESEARCH_ALPHA=1.0


  !
  !================================================================================================================================
  !

  !INITIALISE OPENCMISS

  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion, err )

  !
  !================================================================================================================================
  !
  
  CALL FieldmlInput_InitializeFromFile( fieldmlInfo, inputFilename, err )

  CALL FieldmlInput_SetDofVariables( fieldmlInfo, "test_mesh.nodal_dofs", "test_mesh.element_dofs", &
    & "test_mesh.constant_dofs", err )
  
  spaceHandle = Fieldml_GetNamedObject( fieldmlInfo%fmlHandle, "test_mesh.coordinates"//NUL )
  velocityHandle = Fieldml_GetNamedObject( fieldmlInfo%fmlHandle, "test_mesh.velocity"//NUL )
  pressureHandle = Fieldml_GetNamedObject( fieldmlInfo%fmlHandle, "test_mesh.pressure"//NUL )

  CALL FieldmlInput_GetMeshInfo( fieldmlInfo, "test_mesh.domain", err )
  
  CALL FieldmlInput_GetCoordinateSystemInfo( fieldmlInfo%fmlHandle, spaceHandle, coordinateType, coordinateCount, err )
  
  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise( CoordinateSystem, err )
  CALL CMISSCoordinateSystemCreateStart( CoordinateSystemUserNumber, CoordinateSystem, err )
  !Set the coordinate system dimension and type
  CALL CMISSCoordinateSystemDimensionSet( CoordinateSystem, coordinateCount, err )
  CALL CMISSCoordinateSystemTypeSet( CoordinateSystem, coordinateType, err )
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish( CoordinateSystem, err )

  !
  !================================================================================================================================
  !

  !REGION

  !Start the creation of a new region
  CALL CMISSRegionTypeInitialise( Region, err )
  CALL CMISSRegionCreateStart( RegionUserNumber, WorldRegion, Region, err )
  !Set the regions coordinate system as defined above
  CALL CMISSRegionCoordinateSystemSet( Region, CoordinateSystem, err )
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish( Region, err )

  !
  !================================================================================================================================
  !

  !BASES
  CALL FieldmlInput_CreateBasis( fieldmlInfo, basisNumberTrilinear, "test_mesh.trilinear_lagrange", gaussQuadrature, err )
  CALL FieldmlInput_CreateBasis( fieldmlInfo, basisNumberTriquadratic, "test_mesh.triquadratic_lagrange", gaussQuadrature, &
    & err )
  
  
  !
  !================================================================================================================================
  !

  !MESH
  
  meshComponentCount = 2
  
  nodeCount = Fieldml_GetEnsembleDomainElementCount( fieldmlInfo%fmlHandle, fieldmlInfo%nodesHandle )
  CALL CMISSNodesTypeInitialise( Nodes, err )
  CALL CMISSNodesCreateStart( Region, nodeCount, nodes, err )
  CALL CMISSNodesCreateFinish( Nodes, err )

  xiDimensions = Fieldml_GetDomainComponentCount( fieldmlInfo%fmlHandle, fieldmlInfo%xiHandle )
  elementCount = Fieldml_GetEnsembleDomainElementCount( fieldmlInfo%fmlHandle, fieldmlInfo%elementsHandle )
  CALL CMISSMeshTypeInitialise( Mesh, err )
  CALL CMISSMeshCreateStart( MeshUserNumber, Region, xiDimensions, Mesh, err )
  CALL CMISSMeshNumberOfElementsSet( Mesh, elementCount, err )
  CALL CMISSMeshNumberOfComponentsSet( Mesh, meshComponentCount, err )
  
  CALL FieldmlInput_CreateMeshComponent( fieldmlInfo, RegionUserNumber, MeshUserNumber, 1, &
    & "test_mesh.template.triquadratic", err )
  CALL FieldmlInput_CreateMeshComponent( fieldmlInfo, RegionUserNumber, MeshUserNumber, 2, &
    & "test_mesh.template.trilinear", err )
  
  MESH_COMPONENT_NUMBER_SPACE = 1
  MESH_COMPONENT_NUMBER_VELOCITY = 1
  MESH_COMPONENT_NUMBER_PRESSURE = 2

  !Finish the creation of the mesh
  CALL CMISSMeshCreateFinish(Mesh, err )

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition, err )
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition, err )
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType, err )
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,DomainUserNumber, err )
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition, err )
  
  CALL FieldmlInput_CreateField( fieldmlInfo, Region, Mesh, Decomposition, GeometricFieldUserNumber, GeometricField, &
    & "test_mesh.coordinates", err )

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for static Navier-Stokes
  CALL CMISSEquationsSetTypeInitialise(EquationsSetNavierStokes, err )
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField,EquationsSetNavierStokes, err )
  !Set the equations set to be a static Navier-Stokes problem
  CALL CMISSEquationsSetSpecificationSet(EquationsSetNavierStokes,CMISSEquationsSetFluidMechanicsClass, &
    & CMISSEquationsSetNavierStokesEquationType,CMISSEquationsSetStaticNavierStokesSubtype, err )
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSetNavierStokes, err )


  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for static Navier-Stokes
  CALL CMISSFieldTypeInitialise(DependentFieldNavierStokes, err )
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetNavierStokes,DependentFieldUserNumberNavierStokes, & 
    & DependentFieldNavierStokes, err )
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,coordinateCount
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY, err )
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY, err )
  ENDDO
  COMPONENT_NUMBER=coordinateCount+1
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE, err )
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE, err )
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetNavierStokes, err )

  !Initialise dependent field
  DO COMPONENT_NUMBER=1,coordinateCount
    CALL CMISSFieldComponentValuesInitialise(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
      & COMPONENT_NUMBER,INITIAL_FIELD_NAVIER_STOKES(COMPONENT_NUMBER), err )
  ENDDO


  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for static Navier-Stokes
  CALL CMISSFieldTypeInitialise(MaterialsFieldNavierStokes, err )
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetNavierStokes,MaterialsFieldUserNumberNavierStokes, & 
    & MaterialsFieldNavierStokes, err )
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetNavierStokes, err )
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & MaterialsFieldUserNumberNavierStokesMu,MU_PARAM_NAVIER_STOKES, err )
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & MaterialsFieldUserNumberNavierStokesRho,RHO_PARAM_NAVIER_STOKES, err )


  !
  !================================================================================================================================
  !

  !EQUATIONS


  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(EquationsNavierStokes, err )
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetNavierStokes,EquationsNavierStokes, err )
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(EquationsNavierStokes,CMISSEquationsSparseMatrices, err )
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(EquationsNavierStokes,EQUATIONS_NAVIER_STOKES_OUTPUT, err )
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetNavierStokes, err )


  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  !Start the creation of the equations set boundary conditions for Navier-Stokes
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsNavierStokes, err )
  CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSetNavierStokes,BoundaryConditionsNavierStokes, err )
  !Set fixed wall nodes
  IF(FIXED_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_FIXED_WALL_NODES_NAVIER_STOKES
      NODE_NUMBER=FIXED_WALL_NODES_NAVIER_STOKES(NODE_COUNTER)
      CONDITION=CMISSBoundaryConditionFixedWall
      DO COMPONENT_NUMBER=1,coordinateCount
        VALUE=0.0_CMISSDP
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsNavierStokes,CMISSFieldUVariableType,CMISSNoGlobalDerivative, & 
          & NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE, err )
      ENDDO
    ENDDO
  ENDIF
  !Set velocity boundary conditions
  IF(INLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES
      NODE_NUMBER=INLET_WALL_NODES_NAVIER_STOKES(NODE_COUNTER)
      CONDITION=CMISSBoundaryConditionInletWall
      DO COMPONENT_NUMBER=1,coordinateCount
        VALUE=BOUNDARY_CONDITIONS_NAVIER_STOKES(COMPONENT_NUMBER)
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsNavierStokes,CMISSFieldUVariableType,CMISSNoGlobalDerivative, & 
          & NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE, err )
      ENDDO
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSetNavierStokes, err )
  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem, err )
  CALL CMISSControlLoopTypeInitialise(ControlLoop, err )
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem, err )
  !Set the problem to be a static Navier-Stokes problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemFluidMechanicsClass,CMISSProblemNavierStokesEquationType, &
    & CMISSProblemStaticNavierStokesSubtype, err )
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem, err )
  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem, err )
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem, err )

  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(NonlinearSolverNavierStokes, err )
  CALL CMISSSolverTypeInitialise(LinearSolverNavierStokes, err )
  CALL CMISSProblemSolversCreateStart(Problem, err )
  !Get the nonlinear static solver
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverNavierStokesUserNumber,NonlinearSolverNavierStokes, err )
  !Set the nonlinear Jacobian type
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(NonlinearSolverNavierStokes,CMISSSolverNewtonJacobianAnalyticCalculated, err )
  !Set the output type
  CALL CMISSSolverOutputTypeSet(NonlinearSolverNavierStokes,NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE, err )
  !Set the solver settings
  CALL CMISSSolverNewtonAbsoluteToleranceSet(NonlinearSolverNavierStokes,ABSOLUTE_TOLERANCE, err )
  CALL CMISSSolverNewtonRelativeToleranceSet(NonlinearSolverNavierStokes,RELATIVE_TOLERANCE, err )
  !Get the nonlinear linear solver
  CALL CMISSSolverNewtonLinearSolverGet(NonlinearSolverNavierStokes,LinearSolverNavierStokes, err )
  !Set the output type
  CALL CMISSSolverOutputTypeSet(LinearSolverNavierStokes,LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE, err )


  !Set the solver settings
  IF(LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG) THEN
    CALL CMISSSolverLinearTypeSet(LinearSolverNavierStokes,CMISSSolverLinearDirectSolveType, err )
    CALL CMISSSolverLibraryTypeSet(LinearSolverNavierStokes,CMISSSolverMUMPSLibrary, err )
  ELSE
    CALL CMISSSolverLinearTypeSet(LinearSolverNavierStokes,CMISSSolverLinearIterativeSolveType, err )
    CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolverNavierStokes,MAXIMUM_ITERATIONS, err )
    CALL CMISSSolverLinearIterativeDivergenceToleranceSet(LinearSolverNavierStokes,DIVERGENCE_TOLERANCE, err )
    CALL CMISSSolverLinearIterativeRelativeToleranceSet(LinearSolverNavierStokes,RELATIVE_TOLERANCE, err )
    CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(LinearSolverNavierStokes,ABSOLUTE_TOLERANCE, err )
    CALL CMISSSolverLinearIterativeGMRESRestartSet(LinearSolverNavierStokes,RESTART_VALUE, err )
  ENDIF
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem, err )

  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL CMISSSolverTypeInitialise(LinearSolverNavierStokes, err )
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsNavierStokes, err )
  CALL CMISSProblemSolverEquationsCreateStart(Problem, err )
  !Get the linear solver equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverNavierStokesUserNumber,LinearSolverNavierStokes, err )
  CALL CMISSSolverSolverEquationsGet(LinearSolverNavierStokes,SolverEquationsNavierStokes, err )
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsNavierStokes,CMISSSolverEquationsSparseMatrices, err )
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsNavierStokes,EquationsSetNavierStokes,EquationsSetIndex, err )
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem, err )

  !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn of PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL CMISSProblemSolve(Problem, err )
  WRITE(*,'(A)') "Problem solved!"

  !
  !================================================================================================================================
  !

  !OUTPUT
  
    CALL FieldmlOutput_InitializeInfo( Region, Mesh, xiDimensions, outputDirectory, basename, outputInfo, err )

    CALL FieldmlOutput_AddField( outputInfo, baseName//".geometric", region, mesh, GeometricField, err )

    domainHandle = Fieldml_GetNamedObject( outputInfo%fmlHandle, "library.velocity.rc.3d"//NUL )
    CALL FieldmlOutput_AddFieldComponents( outputInfo, domainHandle, baseName//".velocity", Mesh, DependentFieldNavierStokes, &
      & (/1,2,3/), err )
    
    domainHandle = Fieldml_GetNamedObject( outputInfo%fmlHandle, "library.pressure"//NUL )
    CALL FieldmlOutput_AddFieldComponents( outputInfo, domainHandle, baseName//".pressure", Mesh, DependentFieldNavierStokes, &
      & (/4/), err )
    
    CALL FieldmlOutput_Write( outputInfo, outputFilename, err )
    
    CALL FieldmlUtil_FinalizeInfo( outputInfo )

  EXPORT_FIELD_IO=.TRUE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL CMISSFieldsTypeInitialise(Fields, err )
    CALL CMISSFieldsTypeCreate(Region,Fields, err )
    CALL CMISSFieldIONodesExport(Fields,"StaticNavierStokes","FORTRAN", err )
    CALL CMISSFieldIOElementsExport(Fields,"StaticNavierStokes","FORTRAN", err )
    CALL CMISSFieldsTypeFinalise(Fields, err )
    WRITE(*,'(A)') "Field exported!"
  ENDIF

  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM NAVIERSTOKESSTATICEXAMPLE
