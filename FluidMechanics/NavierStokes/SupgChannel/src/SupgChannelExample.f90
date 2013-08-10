!> \file
!> $Id: BlockChannelFieldMLExample.f90
!> \author David Ladd
!> \brief This is an example program to demonstrate SUPG stabilization using OpenCMISS calls.
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

!> \example FluidMechanics/NavierStokes/BlockChannelFieldml/src/BlockChannelFieldmlExample.f90
!! Example program to solve a 2D SUPG stabilized problem using OpenCMISS calls.
!! \htmlinclude FluidMechanics/NavierStokes/BlockChannelFieldml/history.html
!!
!<

!> Main program

PROGRAM SupgChannel

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OPENCMISS
  USE FIELDML_API
  USE MPI


#ifdef WIN32
  USE IFQWINCMISS
#endif

  !
  !================================================================================================================================
  !

  !PROGRAM VARIABLES AND TYPES

  IMPLICIT NONE

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField


  !Test program parameters
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberNavierStokes=8
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberNavierStokes=9
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=10

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: SolverNavierStokesUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMu=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberRho=2

  INTEGER(CMISSIntg), PARAMETER :: basisNumberBiquadratic=1
  INTEGER(CMISSIntg), PARAMETER :: basisNumberBilinear=2
  INTEGER(CMISSIntg), PARAMETER :: gaussQuadrature(2) = [2,2]

  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: inputFilename = "input/BlockChannel.xml"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputDirectory = "output"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputFilename = "BlockChannel_out.xml"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: basename = "DynamicBlockChannel"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: dataFormat = "PLAIN_TEXT"

  !Program variables

  INTEGER(CMISSIntg) :: maximumIterations
  INTEGER(CMISSIntg) :: restartValue
  INTEGER(CMISSIntg) :: numberOfFixedWallNodes
  INTEGER(CMISSIntg) :: numberOfInletNodes

  INTEGER(CMISSIntg) :: equationsOutputType
  INTEGER(CMISSIntg) :: componentNumber
  INTEGER(CMISSIntg) :: nodeNumber
  INTEGER(CMISSIntg) :: i
  INTEGER(CMISSIntg) :: condition

  INTEGER, ALLOCATABLE, DIMENSION(:):: fixedWallNodes
  INTEGER, ALLOCATABLE, DIMENSION(:):: inletNodes

  INTEGER(CMISSIntg) :: dynamicSolverOutputFrequency
  INTEGER(CMISSIntg) :: dynamicSolverOutputType
  INTEGER(CMISSIntg) :: nonlinearSolverOutputType
  INTEGER(CMISSIntg) :: linearSolverOutputType

  INTEGER(CMISSIntg) :: EquationsSetSubtype
  INTEGER(CMISSIntg) :: ProblemSubtype

  REAL(CMISSDP) :: initialConditions(2)
  REAL(CMISSDP) :: inletBoundaryConditions(2)
  REAL(CMISSDP) :: divergenceTolerance
  REAL(CMISSDP) :: relativeTolerance
  REAL(CMISSDP) :: absoluteTolerance
  REAL(CMISSDP) :: linesearchAlpha
  REAL(CMISSDP) :: VALUE
  REAL(CMISSDP) :: viscosity
  REAL(CMISSDP) :: density

  REAL(CMISSDP) :: dynamicSolverStartTime
  REAL(CMISSDP) :: dynamicSolverStopTime
  REAL(CMISSDP) :: dynamicSolverTheta
  REAL(CMISSDP) :: dynamicSolverTimeIncrement

  LOGICAL :: directLinearSolverFlag
  LOGICAL :: fixedWallNodesFlag
  LOGICAL :: inletNodesFlag
  LOGICAL :: supgFlag
  LOGICAL :: calculateFacesFlag

  CHARACTER *15 buffer
  CHARACTER *15 arg

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
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: DependentField
  TYPE(CMISSFieldType) :: MaterialsField



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
  TYPE(CMISSSolverType) :: DynamicSolverNavierStokes
  TYPE(CMISSSolverType) :: NonlinearSolverNavierStokes
  TYPE(CMISSSolverType) :: LinearSolverNavierStokes
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsNavierStokes

  !FieldML parsing variables
  TYPE(CMISSFieldMLIOType) :: fieldmlInfo, outputInfo
  INTEGER(CMISSIntg) :: typeHandle
  INTEGER(CMISSIntg) :: coordinateCount

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables

  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,BoundaryNodeDomain
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: Err
  
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

  !Use SUPG stabilisation?
  supgFlag=.TRUE.

  !Set initial values
  initialConditions(1)=0.0_CMISSDP
  initialConditions(2)=0.0_CMISSDP
  !Set default boundary conditions
  inletBoundaryConditions(1)=1.0_CMISSDP
  inletBoundaryConditions(2)=0.0_CMISSDP
  fixedWallNodesFlag=.FALSE.
  inletNodesFlag=.FALSE.
  !Initialize calc faces
  calculateFacesFlag=.FALSE.
  !Set material parameters
  viscosity=0.01_CMISSDP
  density=1.0_CMISSDP

  !Set output types
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  dynamicSolverOutputType=CMISS_SOLVER_PROGRESS_OUTPUT
  linearSolverOutputType=CMISS_SOLVER_NO_OUTPUT
  nonlinearSolverOutputType=CMISS_SOLVER_PROGRESS_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  equationsOutputType=CMISS_EQUATIONS_NO_OUTPUT

  !Set dynamic solver parameters
  dynamicSolverStartTime=0.0_CMISSDP
  dynamicSolverStopTime=10.0001_CMISSDP 
  dynamicSolverTimeIncrement=0.1_CMISSDP
  dynamicSolverTheta=1.0_CMISSDP
  !Set result output parameter (e.g. 1 for every step, 5 for every 5 steps, etc)
  dynamicSolverOutputFrequency=10
  !Set solver parameters
  directLinearSolverFlag=.TRUE.
  relativeTolerance=1.0E-5_CMISSDP !default: 1.0E-05_CMISSDP
  absoluteTolerance=1.0E-6_CMISSDP !default: 1.0E-10_CMISSDP
  divergenceTolerance=1.0E5 !default: 1.0E5
  maximumIterations=100000 !default: 100000
  restartValue=300 !default: 30
  linesearchAlpha=1.0

  !Get command line arguments
  DO i=1,COMMAND_ARGUMENT_COUNT()
    CALL GET_COMMAND_ARGUMENT(i,arg)
    SELECT CASE(arg)
      CASE('-density')
        CALL GET_COMMAND_ARGUMENT(i+1,buffer)
        READ(buffer,*) density
      CASE('-viscosity')
        CALL GET_COMMAND_ARGUMENT(i+1,buffer)
        READ(buffer,*) viscosity
      CASE('-directSolver')
        CALL GET_COMMAND_ARGUMENT(i+1,buffer)
        READ(buffer,*) directLinearSolverFlag
      CASE('-outputFrequency')
        CALL GET_COMMAND_ARGUMENT(i+1,buffer)
        READ(buffer,*) dynamicSolverOutputFrequency
      CASE('-SUPG')
        CALL GET_COMMAND_ARGUMENT(i+1,buffer)
        READ(buffer,*) supgFlag
      CASE('-startTime')
        CALL GET_COMMAND_ARGUMENT(i+1,buffer)
        READ(buffer,*) dynamicSolverStartTime
      CASE('-stopTime')
        CALL GET_COMMAND_ARGUMENT(i+1,buffer)
        READ(buffer,*) dynamicSolverStopTime
      CASE('-timeIncrement')
        CALL GET_COMMAND_ARGUMENT(i+1,buffer)
        READ(buffer,*) dynamicSolverTimeIncrement
      CASE('-inletVelocity')
        CALL GET_COMMAND_ARGUMENT(i+1,buffer)
        READ(buffer,*) inletBoundaryConditions(1)
        CALL GET_COMMAND_ARGUMENT(i+2,buffer)
        READ(buffer,*) inletBoundaryConditions(2)
      CASE DEFAULT
        !do nothing
      END SELECT
  ENDDO 
  WRITE(*,*)' '
  WRITE(*,*)' ************************************* '
  WRITE(*,*)' '
  WRITE(*,*)'-density........', density
  WRITE(*,*)'-viscosity......', viscosity
  WRITE(*,*)'-SUPG.......  ', supgFlag
  WRITE(*,*)'-inletVelocity.......', inletBoundaryConditions
  WRITE(*,*)'-startTime......  ', dynamicSolverStartTime
  WRITE(*,*)'-stopTime.......  ', dynamicSolverStopTime
  WRITE(*,*)'-timeIncrement..  ', dynamicSolverTimeIncrement
  WRITE(*,*)'-outputFrequency..  ', dynamicSolverOutputFrequency
  WRITE(*,*)'-directSolver...  ', directLinearSolverFlag
  WRITE(*,*)' '
  WRITE(*,*)' ************************************* '
  WRITE(*,*)' '
  WRITE(*,*) ' ' 
  !Set boundary conditions
  INQUIRE(FILE="./input/bc/wallNodes.dat", EXIST=fixedWallNodesFlag)
  INQUIRE(FILE="./input/bc/inletNodes.dat", EXIST=inletNodesFlag)
  IF(fixedWallNodesFlag) THEN
    OPEN(UNIT=1, FILE="./input/bc/wallNodes.dat",STATUS='unknown')
    READ(1,*) numberOfFixedWallNodes
    ALLOCATE(fixedWallNodes(numberOfFixedWallNodes))
    READ(1,*) fixedWallNodes(1:numberOfFixedWallNodes)
    CLOSE(1)
  ENDIF
  IF(inletNodesFlag) THEN
     OPEN(UNIT=1, FILE="./input/bc/inletNodes.dat",STATUS='unknown')
    READ(1,*) numberOfInletNodes
    ALLOCATE(inletNodes(numberOfInletNodes))
    READ(1,*) inletNodes(1:numberOfInletNodes)
    CLOSE(1)
  ENDIF

  IF(supgFlag) THEN
    EquationsSetSubtype=CMISS_EQUATIONS_SET_TRANSIENT_SUPG_NAVIER_STOKES_SUBTYPE
    ProblemSubtype=CMISS_PROBLEM_TRANSIENT_SUPG_NAVIER_STOKES_SUBTYPE
  ELSE
    EquationsSetSubtype=CMISS_EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE
    ProblemSubtype=CMISS_PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE
  ENDIF

  !
  !================================================================================================================================
  !

  !INITIALISE OPENCMISS

  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  !
  !================================================================================================================================
  !

  !CHECK COMPUTATIONAL NODE

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !
  !================================================================================================================================
  !

  !INITIALISE FieldML

  CALL CMISSFieldMLIO_Initialise( fieldmlInfo, err ) 
  CALL CMISSFieldML_InputCreateFromFile( inputFilename, fieldmlInfo, err )

  !
  !================================================================================================================================
  !

  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSFieldML_InputCoordinateSystemCreateStart( fieldmlInfo, "BlockChannelMesh.coordinates", CoordinateSystem, &
    & CoordinateSystemUserNumber, err )

  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_DimensionGet( CoordinateSystem, coordinateCount, err )

  !
  !================================================================================================================================
  !

  !REGION

  !Start the creation of a new region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system as defined above
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !NODES
  CALL CMISSFieldML_InputNodesCreateStart( fieldmlInfo, "BlockChannelMesh.nodes.argument", Region, nodes, err )
  CALL CMISSNodes_CreateFinish( Nodes, err )

  !
  !================================================================================================================================
  !

  !BASES
  CALL CMISSFieldML_InputBasisCreateStart( fieldmlInfo, "BlockChannelMesh.biquadratic_lagrange", basisNumberBiquadratic, err )
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet( basisNumberBiquadratic, gaussQuadrature, err )
  CALL CMISSBasis_CreateFinish( basisNumberBiquadratic, err )

  CALL CMISSFieldML_InputBasisCreateStart( fieldmlInfo, "BlockChannelMesh.bilinear_lagrange", basisNumberBilinear, err )
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet( basisNumberBilinear, gaussQuadrature, err )
  CALL CMISSBasis_CreateFinish( basisNumberBilinear, err )

  !
  !================================================================================================================================
  !

  !MESH
  CALL CMISSFieldML_InputMeshCreateStart( fieldmlInfo, "BlockChannelMesh.mesh.argument", Mesh, MeshUserNumber, Region, err )
  CALL CMISSMesh_NumberOfComponentsSet( Mesh, 2, err )

  CALL CMISSFieldML_InputCreateMeshComponent( fieldmlInfo, RegionUserNumber, MeshUserNumber, 1, &
    & "BlockChannelMesh.template.biquadratic", err )
  CALL CMISSFieldML_InputCreateMeshComponent( fieldmlInfo, RegionUserNumber, MeshUserNumber, 2, &
    & "BlockChannelMesh.template.bilinear", err )

  !Finish the creation of the mesh
  CALL CMISSMesh_CreateFinish(Mesh, err )

  !
  !================================================================================================================================
  !

  !Decomposition

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  CALL CMISSDecomposition_CalculateFacesSet(Decomposition,calculateFacesFlag,Err)

  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  CALL CMISSFieldML_InputFieldCreateStart( fieldmlInfo, Region, Decomposition, GeometricFieldUserNumber, GeometricField, &
    & CMISS_FIELD_U_VARIABLE_TYPE, "BlockChannelMesh.coordinates", err )
  CALL CMISSField_CreateFinish( RegionUserNumber, GeometricFieldUserNumber, err )

  CALL CMISSFieldML_InputFieldParametersUpdate( fieldmlInfo, GeometricField, "BlockChannelMesh.node.coordinates", &
    & CMISS_FIELD_U_VARIABLE_TYPE, CMISS_FIELD_VALUES_SET_TYPE, err )
  CALL CMISSField_ParameterSetUpdateStart( GeometricField, CMISS_FIELD_U_VARIABLE_TYPE, CMISS_FIELD_VALUES_SET_TYPE, err )
  CALL CMISSField_ParameterSetUpdateFinish( GeometricField, CMISS_FIELD_U_VARIABLE_TYPE, CMISS_FIELD_VALUES_SET_TYPE, err )

  !
  !================================================================================================================================
  !

  CALL CMISSFieldMLIO_Finalise( fieldmlInfo, err )

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  CALL CMISSEquationsSet_Initialise(EquationsSetNavierStokes,Err)
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField, &
    & CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMISS_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE, &
    & EquationsSetSubtype,EquationsSetFieldUserNumber,EquationsSetField,EquationsSetNavierStokes,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSetNavierStokes,Err)

  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for dynamic Navier-Stokes
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetNavierStokes,DependentFieldUserNumber,DependentField,Err)
  !Set the mesh component to be used by the field components.
  !   Velocity components
  DO componentNumber=1,coordinateCount
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,componentNumber,1,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,componentNumber,1,Err)
  ENDDO
  !   Pressure component
  componentNumber=coordinateCount+1
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,componentNumber,2,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,componentNumber,2,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetNavierStokes,Err)

  !Initialise dependent field
  DO componentNumber=1,coordinateCount
    CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & componentNumber,initialConditions(componentNumber),Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !MATERIALS FIELD

  !Create the equations set materials field variables for static Navier-Stokes
  CALL CMISSField_Initialise(MaterialsField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetNavierStokes,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetNavierStokes,Err)

  ! Materials parameters, viscosity and density
  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberMu,viscosity,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberRho,density,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS

  !Create the equations set equations
  CALL CMISSEquations_Initialise(EquationsNavierStokes,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetNavierStokes,EquationsNavierStokes,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(EquationsNavierStokes,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(EquationsNavierStokes,equationsOutputType,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetNavierStokes,Err)

  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a dynamic Navier-Stokes problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_FLUID_MECHANICS_CLASS,CMISS_PROBLEM_NAVIER_STOKES_EQUATION_TYPE, &
    & ProblemSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoop_TimesSet(ControlLoop,dynamicSolverStartTime,dynamicSolverStopTime, & 
    & dynamicSolverTimeIncrement,Err)
  !Set the output timing
  CALL CMISSControlLoop_TimeOutputSet(ControlLoop,dynamicSolverOutputFrequency,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL CMISSSolver_Initialise(DynamicSolverNavierStokes,Err)
  CALL CMISSSolver_Initialise(NonlinearSolverNavierStokes,Err)
  CALL CMISSSolver_Initialise(LinearSolverNavierStokes,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  !Get the dynamic dymamic solver
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(DynamicSolverNavierStokes,dynamicSolverOutputType,Err)
  !Set theta
  CALL CMISSSolver_DynamicThetaSet(DynamicSolverNavierStokes,dynamicSolverTheta,Err)
  !Get the dynamic nonlinear solver
  CALL CMISSSolver_DynamicNonlinearSolverGet(DynamicSolverNavierStokes,NonlinearSolverNavierStokes,Err)
  !Set the nonlinear Jacobian type
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolverNavierStokes, &
    & CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  !Set the nonlinear solver output type
  CALL CMISSSolver_OutputTypeSet(NonlinearSolverNavierStokes,nonlinearSolverOutputType,Err)
  !Set the solver settings
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(NonlinearSolverNavierStokes,absoluteTolerance,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(NonlinearSolverNavierStokes,relativeTolerance,Err)
  CALL CMISSSolver_NewtonMaximumIterationsSet(NonlinearSolverNavierStokes,maximumIterations,Err)
  !Get the dynamic nonlinear linear solver
  CALL CMISSSolver_NewtonLinearSolverGet(NonlinearSolverNavierStokes,LinearSolverNavierStokes,Err)
  !Set the linear solver output type
  CALL CMISSSolver_OutputTypeSet(LinearSolverNavierStokes,linearSolverOutputType,Err)
  !Set the solver settings
  IF(directLinearSolverFlag) THEN
    CALL CMISSSolver_LinearTypeSet(LinearSolverNavierStokes,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL CMISSSolver_LibraryTypeSet(LinearSolverNavierStokes,CMISS_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL CMISSSolver_LinearTypeSet(LinearSolverNavierStokes,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolverNavierStokes,maximumIterations,Err)
    CALL CMISSSolver_LinearIterativeDivergenceToleranceSet(LinearSolverNavierStokes,divergenceTolerance,Err)
    CALL CMISSSolver_LinearIterativeRelativeToleranceSet(LinearSolverNavierStokes,relativeTolerance,Err)
    CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(LinearSolverNavierStokes,absoluteTolerance,Err)
    CALL CMISSSolver_LinearIterativeGMRESRestartSet(LinearSolverNavierStokes,restartValue,Err)
  ENDIF

  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(DynamicSolverNavierStokes,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsNavierStokes,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  !Get the dynamic solver equations
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  CALL CMISSSolver_SolverEquationsGet(DynamicSolverNavierStokes,SolverEquationsNavierStokes,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsNavierStokes,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !Add in the equations set
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsNavierStokes,EquationsSetNavierStokes,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  !Start the creation of the equations set boundary conditions for Stokes
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsNavierStokes,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsNavierStokes,BoundaryConditionsNavierStokes,Err)
  !Set fixed wall nodes
  IF(fixedWallNodesFlag) THEN
    DO i=1,numberOfFixedWallNodes
      nodeNumber=fixedWallNodes(i)
      condition=CMISS_BOUNDARY_CONDITION_FIXED
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,nodeNumber,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO componentNumber=1,coordinateCount
          VALUE=0.0_CMISSDP
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentField, &
            & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV,nodeNumber,componentNumber,condition,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !Set inlet velocity boundary conditions
  IF(inletNodesFlag) THEN
     DO i=1,numberOfInletNodes
      nodeNumber=inletNodes(i)
      condition=CMISS_BOUNDARY_CONDITION_FIXED
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,nodeNumber,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO componentNumber=1,coordinateCount
          VALUE=inletBoundaryConditions(componentNumber)
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentField, &
            & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV,nodeNumber,componentNumber,condition,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsNavierStokes,Err)


  !
  !================================================================================================================================
  !

  !RUN SOLVER

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL CMISSProblem_Solve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"
! 
  !
  !================================================================================================================================
  !

  !OUTPUT

  ! Export FieldML data
  CALL CMISSFieldMLIO_Initialise( outputInfo, err )

  CALL CMISSFieldML_OutputCreate( Mesh, outputDirectory, basename, dataFormat, outputInfo, err )
  CALL CMISSFieldML_OutputAddImport( outputInfo, "coordinates.rc.2d", typeHandle, err )
  CALL CMISSFieldML_OutputAddField( outputInfo, baseName//".geometric", dataFormat, GeometricField, &
    & CMISS_FIELD_U_VARIABLE_TYPE, CMISS_FIELD_VALUES_SET_TYPE, err )
  CALL CMISSFieldML_OutputAddFieldComponents( outputInfo, typeHandle, baseName//".velocity", dataFormat, &
    & DependentField, [1,2], CMISS_FIELD_U_VARIABLE_TYPE, CMISS_FIELD_VALUES_SET_TYPE, err )
  CALL CMISSFieldML_OutputAddImport( outputInfo, "real.1d", typeHandle, err )
  CALL CMISSFieldML_OutputAddFieldComponents( outputInfo, typeHandle, baseName//".pressure", dataFormat, &
    & DependentField, [3], CMISS_FIELD_U_VARIABLE_TYPE, CMISS_FIELD_VALUES_SET_TYPE, err )
  CALL CMISSFieldML_OutputWrite( outputInfo, outputFilename, err )

  CALL CMISSFieldMLIO_Finalise( outputInfo, err )

  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM SupgChannel
