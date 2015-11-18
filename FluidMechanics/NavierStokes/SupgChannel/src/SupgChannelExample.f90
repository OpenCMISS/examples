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

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE FIELDML_API
#ifndef NOMPIMOD
  USE MPI
#endif


#ifdef WIN32
  USE IFQWINCMISS
#endif

  !
  !================================================================================================================================
  !

  !PROGRAM VARIABLES AND TYPES

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberNavierStokes=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberNavierStokes=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11

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

  REAL(CMISSRP) :: initialConditions(2)
  REAL(CMISSRP) :: inletBoundaryConditions(2)
  REAL(CMISSRP) :: divergenceTolerance
  REAL(CMISSRP) :: relativeTolerance
  REAL(CMISSRP) :: absoluteTolerance
  REAL(CMISSRP) :: linesearchAlpha
  REAL(CMISSRP) :: VALUE
  REAL(CMISSRP) :: viscosity
  REAL(CMISSRP) :: density

  REAL(CMISSRP) :: dynamicSolverStartTime
  REAL(CMISSRP) :: dynamicSolverStopTime
  REAL(CMISSRP) :: dynamicSolverTheta
  REAL(CMISSRP) :: dynamicSolverTimeIncrement

  LOGICAL :: directLinearSolverFlag
  LOGICAL :: fixedWallNodesFlag
  LOGICAL :: inletNodesFlag
  LOGICAL :: supgFlag
  LOGICAL :: calculateFacesFlag

  CHARACTER *15 buffer
  CHARACTER *15 arg

  !CMISS variables

  !Regions
  TYPE(cmfe_RegionType) :: Region
  TYPE(cmfe_RegionType) :: WorldRegion
  !Coordinate systems
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem
  TYPE(cmfe_CoordinateSystemType) :: WorldCoordinateSystem
  !Nodes
  TYPE(cmfe_NodesType) :: Nodes
  !Meshes
  TYPE(cmfe_MeshType) :: Mesh
  !Decompositions
  TYPE(cmfe_DecompositionType) :: Decomposition
  !Field types
  TYPE(cmfe_FieldType) :: GeometricField
  TYPE(cmfe_FieldType) :: EquationsSetField
  TYPE(cmfe_FieldType) :: DependentField
  TYPE(cmfe_FieldType) :: MaterialsField
  !Boundary conditions
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsNavierStokes
  !Equations sets
  TYPE(cmfe_EquationsSetType) :: EquationsSetNavierStokes
  !Equations
  TYPE(cmfe_EquationsType) :: EquationsNavierStokes
  !Problems
  TYPE(cmfe_ProblemType) :: Problem
  !Control loops
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  !Solvers
  TYPE(cmfe_SolverType) :: DynamicSolverNavierStokes
  TYPE(cmfe_SolverType) :: NonlinearSolverNavierStokes
  TYPE(cmfe_SolverType) :: LinearSolverNavierStokes
  !Solver equations
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsNavierStokes

  !FieldML parsing variables
  TYPE(cmfe_FieldMLIOType) :: fieldmlInfo, outputInfo
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
  initialConditions(1)=0.0_CMISSRP
  initialConditions(2)=0.0_CMISSRP
  !Set default boundary conditions
  inletBoundaryConditions(1)=1.0_CMISSRP
  inletBoundaryConditions(2)=0.0_CMISSRP
  fixedWallNodesFlag=.FALSE.
  inletNodesFlag=.FALSE.
  !Initialize calc faces
  calculateFacesFlag=.FALSE.
  !Set material parameters
  viscosity=0.01_CMISSRP
  density=1.0_CMISSRP

  !Set output types
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  dynamicSolverOutputType=CMFE_SOLVER_PROGRESS_OUTPUT
  linearSolverOutputType=CMFE_SOLVER_NO_OUTPUT
  nonlinearSolverOutputType=CMFE_SOLVER_PROGRESS_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  equationsOutputType=CMFE_EQUATIONS_NO_OUTPUT

  !Set dynamic solver parameters
  dynamicSolverStartTime=0.0_CMISSRP
  dynamicSolverStopTime=10.0001_CMISSRP 
  dynamicSolverTimeIncrement=0.1_CMISSRP
  dynamicSolverTheta=1.0_CMISSRP
  !Set result output parameter (e.g. 1 for every step, 5 for every 5 steps, etc)
  dynamicSolverOutputFrequency=10
  !Set solver parameters
  directLinearSolverFlag=.TRUE.
  relativeTolerance=1.0E-5_CMISSRP !default: 1.0E-05_CMISSRP
  absoluteTolerance=1.0E-6_CMISSRP !default: 1.0E-10_CMISSRP
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
      CASE('-outFrequency')
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
  WRITE(*,*)'-outFrequency..  ', dynamicSolverOutputFrequency
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
    EquationsSetSubtype=CMFE_EQUATIONS_SET_TRANSIENT_SUPG_NAVIER_STOKES_SUBTYPE
    ProblemSubtype=CMFE_PROBLEM_TRANSIENT_SUPG_NAVIER_STOKES_SUBTYPE
  ELSE
    EquationsSetSubtype=CMFE_EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE
    ProblemSubtype=CMFE_PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE
  ENDIF

  !
  !================================================================================================================================
  !

  !INITIALISE OPENCMISS

  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  !
  !================================================================================================================================
  !

  !CHECK COMPUTATIONAL NODE

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !
  !================================================================================================================================
  !

  !INITIALISE FieldML

  CALL cmfe_FieldMLIO_Initialise( fieldmlInfo, err ) 
  CALL cmfe_FieldML_InputCreateFromFile( inputFilename, fieldmlInfo, err )

  !
  !================================================================================================================================
  !

  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_FieldML_InputCoordinateSystemCreateStart( fieldmlInfo, "BlockChannelMesh.coordinates", CoordinateSystem, &
    & CoordinateSystemUserNumber, err )

  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_DimensionGet( CoordinateSystem, coordinateCount, err )

  !
  !================================================================================================================================
  !

  !REGION

  !Start the creation of a new region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set region label
  CALL cmfe_Region_LabelSet(Region,"OpenCMISS",Err)
  !Set the regions coordinate system as defined above
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !NODES
  CALL cmfe_FieldML_InputNodesCreateStart( fieldmlInfo, "BlockChannelMesh.nodes.argument", Region, nodes, err )
  CALL cmfe_Nodes_CreateFinish( Nodes, err )

  !
  !================================================================================================================================
  !

  !BASES
  CALL cmfe_FieldML_InputBasisCreateStart( fieldmlInfo, "BlockChannelMesh.biquadratic_lagrange", basisNumberBiquadratic, err )
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet( basisNumberBiquadratic, gaussQuadrature, err )
  CALL cmfe_Basis_CreateFinish( basisNumberBiquadratic, err )

  CALL cmfe_FieldML_InputBasisCreateStart( fieldmlInfo, "BlockChannelMesh.bilinear_lagrange", basisNumberBilinear, err )
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet( basisNumberBilinear, gaussQuadrature, err )
  CALL cmfe_Basis_CreateFinish( basisNumberBilinear, err )

  !
  !================================================================================================================================
  !

  !MESH
  CALL cmfe_FieldML_InputMeshCreateStart( fieldmlInfo, "BlockChannelMesh.mesh.argument", Mesh, MeshUserNumber, Region, err )
  CALL cmfe_Mesh_NumberOfComponentsSet( Mesh, 2, err )

  CALL cmfe_FieldML_InputCreateMeshComponent( fieldmlInfo, RegionUserNumber, MeshUserNumber, 1, &
    & "BlockChannelMesh.template.biquadratic", err )
  CALL cmfe_FieldML_InputCreateMeshComponent( fieldmlInfo, RegionUserNumber, MeshUserNumber, 2, &
    & "BlockChannelMesh.template.bilinear", err )

  !Finish the creation of the mesh
  CALL cmfe_Mesh_CreateFinish(Mesh, err )

  !
  !================================================================================================================================
  !

  !Decomposition

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  CALL cmfe_Decomposition_CalculateFacesSet(Decomposition,calculateFacesFlag,Err)

  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  CALL cmfe_FieldML_InputFieldCreateStart( fieldmlInfo, Region, Decomposition, GeometricFieldUserNumber, GeometricField, &
    & CMFE_FIELD_U_VARIABLE_TYPE, "BlockChannelMesh.coordinates", err )
  CALL cmfe_Field_CreateFinish( RegionUserNumber, GeometricFieldUserNumber, err )

  CALL cmfe_FieldML_InputFieldParametersUpdate( fieldmlInfo, GeometricField, "BlockChannelMesh.node.coordinates", &
    & CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err )
  CALL cmfe_Field_ParameterSetUpdateStart( GeometricField, CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err )
  CALL cmfe_Field_ParameterSetUpdateFinish( GeometricField, CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err )

  !
  !================================================================================================================================
  !

  CALL cmfe_FieldMLIO_Finalise( fieldmlInfo, err )

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  CALL cmfe_EquationsSet_Initialise(EquationsSetNavierStokes,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField, &
    & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMFE_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE,EquationsSetSubtype], &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSetNavierStokes,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetNavierStokes,Err)

  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for dynamic Navier-Stokes
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetNavierStokes,DependentFieldUserNumber,DependentField,Err)
  !Set the mesh component to be used by the field components.
  !   Velocity components
  DO componentNumber=1,coordinateCount
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,componentNumber,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,componentNumber,1,Err)
  ENDDO
  !   Pressure component
  componentNumber=coordinateCount+1
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,componentNumber,2,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,componentNumber,2,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetNavierStokes,Err)

  !Initialise dependent field
  DO componentNumber=1,coordinateCount
    CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
      & componentNumber,initialConditions(componentNumber),Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !MATERIALS FIELD

  !Create the equations set materials field variables for static Navier-Stokes
  CALL cmfe_Field_Initialise(MaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetNavierStokes,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetNavierStokes,Err)

  ! Materials parameters, viscosity and density
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberMu,viscosity,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberRho,density,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(EquationsNavierStokes,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetNavierStokes,EquationsNavierStokes,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(EquationsNavierStokes,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(EquationsNavierStokes,equationsOutputType,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetNavierStokes,Err)

  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_FLUID_MECHANICS_CLASS,CMFE_PROBLEM_NAVIER_STOKES_EQUATION_TYPE, &
    & ProblemSubtype],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,dynamicSolverStartTime,dynamicSolverStopTime, & 
    & dynamicSolverTimeIncrement,Err)
  !Set the output timing
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,dynamicSolverOutputFrequency,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(DynamicSolverNavierStokes,Err)
  CALL cmfe_Solver_Initialise(NonlinearSolverNavierStokes,Err)
  CALL cmfe_Solver_Initialise(LinearSolverNavierStokes,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  !Get the dynamic dymamic solver
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(DynamicSolverNavierStokes,dynamicSolverOutputType,Err)
  !Set theta
  CALL cmfe_Solver_DynamicThetaSet(DynamicSolverNavierStokes,dynamicSolverTheta,Err)
  !Get the dynamic nonlinear solver
  CALL cmfe_Solver_DynamicNonlinearSolverGet(DynamicSolverNavierStokes,NonlinearSolverNavierStokes,Err)
  !Set the nonlinear Jacobian type
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(NonlinearSolverNavierStokes, &
    & CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  !Set the nonlinear solver output type
  CALL cmfe_Solver_OutputTypeSet(NonlinearSolverNavierStokes,nonlinearSolverOutputType,Err)
  !Set the solver settings
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(NonlinearSolverNavierStokes,absoluteTolerance,Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(NonlinearSolverNavierStokes,relativeTolerance,Err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(NonlinearSolverNavierStokes,maximumIterations,Err)
  !Get the dynamic nonlinear linear solver
  CALL cmfe_Solver_NewtonLinearSolverGet(NonlinearSolverNavierStokes,LinearSolverNavierStokes,Err)
  !Set the linear solver output type
  CALL cmfe_Solver_OutputTypeSet(LinearSolverNavierStokes,linearSolverOutputType,Err)
  !Set the solver settings
  IF(directLinearSolverFlag) THEN
    CALL cmfe_Solver_LinearTypeSet(LinearSolverNavierStokes,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LibraryTypeSet(LinearSolverNavierStokes,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(LinearSolverNavierStokes,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolverNavierStokes,maximumIterations,Err)
    CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(LinearSolverNavierStokes,divergenceTolerance,Err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolverNavierStokes,relativeTolerance,Err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolverNavierStokes,absoluteTolerance,Err)
    CALL cmfe_Solver_LinearIterativeGMRESRestartSet(LinearSolverNavierStokes,restartValue,Err)
  ENDIF

  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(DynamicSolverNavierStokes,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsNavierStokes,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the dynamic solver equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  CALL cmfe_Solver_SolverEquationsGet(DynamicSolverNavierStokes,SolverEquationsNavierStokes,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsNavierStokes,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsNavierStokes,EquationsSetNavierStokes,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  !Start the creation of the equations set boundary conditions for Stokes
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsNavierStokes,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsNavierStokes,BoundaryConditionsNavierStokes,Err)
  !Set fixed wall nodes
  IF(fixedWallNodesFlag) THEN
    DO i=1,numberOfFixedWallNodes
      nodeNumber=fixedWallNodes(i)
      condition=CMFE_BOUNDARY_CONDITION_FIXED
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeNumber,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO componentNumber=1,coordinateCount
          VALUE=0.0_CMISSRP
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentField, &
            & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV,nodeNumber,componentNumber,condition,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !Set inlet velocity boundary conditions
  IF(inletNodesFlag) THEN
     DO i=1,numberOfInletNodes
      nodeNumber=inletNodes(i)
      condition=CMFE_BOUNDARY_CONDITION_FIXED
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,nodeNumber,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO componentNumber=1,coordinateCount
          VALUE=inletBoundaryConditions(componentNumber)
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentField, &
            & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV,nodeNumber,componentNumber,condition,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsNavierStokes,Err)


  !
  !================================================================================================================================
  !

  !RUN SOLVER

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL cmfe_Problem_Solve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"
! 
  !
  !================================================================================================================================
  !

  !OUTPUT

  ! Export FieldML data
  CALL cmfe_FieldMLIO_Initialise( outputInfo, err )

  CALL cmfe_FieldML_OutputCreate( Mesh, outputDirectory, basename, dataFormat, outputInfo, err )
  CALL cmfe_FieldML_OutputAddImport( outputInfo, "coordinates.rc.2d", typeHandle, err )
  CALL cmfe_FieldML_OutputAddField( outputInfo, baseName//".geometric", dataFormat, GeometricField, &
    & CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err )
  CALL cmfe_FieldML_OutputAddFieldComponents( outputInfo, typeHandle, baseName//".velocity", dataFormat, &
    & DependentField, [1,2], CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err )
  CALL cmfe_FieldML_OutputAddImport( outputInfo, "real.1d", typeHandle, err )
  CALL cmfe_FieldML_OutputAddFieldComponents( outputInfo, typeHandle, baseName//".pressure", dataFormat, &
    & DependentField, [3], CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE, err )
  CALL cmfe_FieldML_OutputWrite( outputInfo, outputFilename, err )

  CALL cmfe_FieldMLIO_Finalise( outputInfo, err )

  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM SupgChannel
