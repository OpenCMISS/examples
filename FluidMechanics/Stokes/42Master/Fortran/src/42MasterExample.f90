!> \file
!> \author Sebastian Krittian
!> \brief This is a master example program to solve static, dynamic and ALE Stokes equations using OpenCMISS calls.
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

!> \example FluidMechanics/Stokes/RoutineCheck/Static/src/StaticExample.f90
!! Example program to solve a static Stokes equation using OpenCMISS calls.
!! \htmlinclude FluidMechanics/Stokes/RoutineCheck/Static/history.html
!!
!<

!> Main program

PROGRAM STOKESMASTEREXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OPENCMISS
  USE FLUID_MECHANICS_IO_ROUTINES
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

  INTEGER(CMFEIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
!   TYPE(cmfe_FieldType) :: EquationsSetField


  !Test program parameters

  INTEGER(CMFEIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMFEIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMFEIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMFEIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMFEIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMFEIntg), PARAMETER :: DependentFieldUserNumberStokes=6
  INTEGER(CMFEIntg), PARAMETER :: DependentFieldUserNumberMovingMesh=42
  INTEGER(CMFEIntg), PARAMETER :: MaterialsFieldUserNumberStokes=8
  INTEGER(CMFEIntg), PARAMETER :: MaterialsFieldUserNumberMovingMesh=9
  INTEGER(CMFEIntg), PARAMETER :: IndependentFieldUserNumberStokes=10
  INTEGER(CMFEIntg), PARAMETER :: IndependentFieldUserNumberMovingMesh=11
  INTEGER(CMFEIntg), PARAMETER :: EquationsSetUserNumberStokes=12
  INTEGER(CMFEIntg), PARAMETER :: EquationsSetUserNumberMovingMesh=13
  INTEGER(CMFEIntg), PARAMETER :: EquationsSetFieldUserNumberStokes=22
  INTEGER(CMFEIntg), PARAMETER :: EquationsSetFieldUserNumberMovingMesh=23
  INTEGER(CMFEIntg), PARAMETER :: ProblemUserNumber=14
  INTEGER(CMFEIntg), PARAMETER :: AnalyticFieldUserNumberStokes=42

  INTEGER(CMFEIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMFEIntg), PARAMETER :: SolverStokesUserNumber=1
  !should be 1 and the Stokes solver 2
  INTEGER(CMFEIntg), PARAMETER :: SolverMovingMeshUserNumber=1
  INTEGER(CMFEIntg), PARAMETER :: MaterialsFieldUserNumberStokesMu=1
  INTEGER(CMFEIntg), PARAMETER :: MaterialsFieldUserNumberStokesRho=2
  INTEGER(CMFEIntg), PARAMETER :: MaterialsFieldUserNumberMovingMeshK=1

  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMFEIntg) :: NUMBER_OF_DIMENSIONS
  
  INTEGER(CMFEIntg) :: BASIS_TYPE
  INTEGER(CMFEIntg) :: BASIS_NUMBER_SPACE
  INTEGER(CMFEIntg) :: BASIS_NUMBER_VELOCITY
  INTEGER(CMFEIntg) :: BASIS_NUMBER_PRESSURE
  INTEGER(CMFEIntg) :: BASIS_GAUSS_SPACE
  INTEGER(CMFEIntg) :: BASIS_GAUSS_VELOCITY
  INTEGER(CMFEIntg) :: BASIS_GAUSS_PRESSURE
  INTEGER(CMFEIntg) :: BASIS_XI_INTERPOLATION_SPACE
  INTEGER(CMFEIntg) :: BASIS_XI_INTERPOLATION_VELOCITY
  INTEGER(CMFEIntg) :: BASIS_XI_INTERPOLATION_PRESSURE
  INTEGER(CMFEIntg) :: MESH_NUMBER_OF_COMPONENTS
  INTEGER(CMFEIntg) :: MESH_COMPONENT_NUMBER_SPACE
  INTEGER(CMFEIntg) :: MESH_COMPONENT_NUMBER_VELOCITY
  INTEGER(CMFEIntg) :: MESH_COMPONENT_NUMBER_PRESSURE
  INTEGER(CMFEIntg) :: NUMBER_OF_NODES_SPACE
  INTEGER(CMFEIntg) :: NUMBER_OF_NODES_VELOCITY
  INTEGER(CMFEIntg) :: NUMBER_OF_NODES_PRESSURE
  INTEGER(CMFEIntg) :: NUMBER_OF_ELEMENT_NODES_SPACE
  INTEGER(CMFEIntg) :: NUMBER_OF_ELEMENT_NODES_VELOCITY
  INTEGER(CMFEIntg) :: NUMBER_OF_ELEMENT_NODES_PRESSURE
  INTEGER(CMFEIntg) :: TOTAL_NUMBER_OF_NODES
  INTEGER(CMFEIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMFEIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMFEIntg) :: RESTART_VALUE
!   INTEGER(CMFEIntg) :: MPI_IERROR
  INTEGER(CMFEIntg) :: NUMBER_OF_FIXED_WALL_NODES_STOKES
  INTEGER(CMFEIntg) :: NUMBER_OF_MOVED_WALL_NODES_STOKES
  INTEGER(CMFEIntg) :: NUMBER_OF_INLET_WALL_NODES_STOKES
  INTEGER(CMFEIntg) :: NUMBER_OF_FIXED_WALL_NODES_MOVING_MESH
  INTEGER(CMFEIntg) :: NUMBER_OF_MOVED_WALL_NODES_MOVING_MESH

  INTEGER(CMFEIntg) :: EQUATIONS_STOKES_OUTPUT
  INTEGER(CMFEIntg) :: EQUATIONS_MOVING_MESH_OUTPUT
  INTEGER(CMFEIntg) :: COMPONENT_NUMBER
  INTEGER(CMFEIntg) :: NODE_NUMBER
  INTEGER(CMFEIntg) :: ELEMENT_NUMBER
  INTEGER(CMFEIntg) :: NODE_COUNTER
  INTEGER(CMFEIntg) :: CONDITION
  INTEGER(CMFEIntg) :: I
  INTEGER(CMFEIntg) :: ANALYTIC_TYPE

  INTEGER(CMFEIntg) :: DYNAMIC_SOLVER_STOKES_INPUT_OPTION=2
  INTEGER(CMFEIntg) :: DYNAMIC_SOLVER_STOKES_OUTPUT_FREQUENCY
  INTEGER(CMFEIntg) :: DYNAMIC_SOLVER_STOKES_OUTPUT_TYPE
  INTEGER(CMFEIntg) :: LINEAR_SOLVER_STOKES_OUTPUT_TYPE
  INTEGER(CMFEIntg) :: LINEAR_SOLVER_MOVING_MESH_OUTPUT_TYPE

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: MOVED_WALL_NODES_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_MOVING_MESH
  INTEGER, ALLOCATABLE, DIMENSION(:):: MOVED_WALL_NODES_MOVING_MESH

  REAL(CMFEDP) :: INITIAL_FIELD_STOKES(3)
  REAL(CMFEDP) :: INITIAL_FIELD_MOVING_MESH(3)
  REAL(CMFEDP) :: BOUNDARY_CONDITIONS_STOKES(3)
  REAL(CMFEDP) :: BOUNDARY_CONDITIONS_MOVING_MESH(3)
  REAL(CMFEDP) :: DIVERGENCE_TOLERANCE
  REAL(CMFEDP) :: RELATIVE_TOLERANCE
  REAL(CMFEDP) :: ABSOLUTE_TOLERANCE
  REAL(CMFEDP) :: LINESEARCH_ALPHA
  REAL(CMFEDP) :: VALUE
  REAL(CMFEDP) :: K_PARAM_MOVING_MESH=1.0_CMFEDP
  REAL(CMFEDP) :: MU_PARAM_STOKES
  REAL(CMFEDP) :: RHO_PARAM_STOKES

  REAL(CMFEDP) :: DYNAMIC_SOLVER_STOKES_START_TIME
  REAL(CMFEDP) :: DYNAMIC_SOLVER_STOKES_STOP_TIME
  REAL(CMFEDP) :: DYNAMIC_SOLVER_STOKES_THETA
  REAL(CMFEDP) :: DYNAMIC_SOLVER_STOKES_TIME_INCREMENT

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_STOKES_DIRECT_FLAG
  LOGICAL :: LINEAR_SOLVER_MOVING_MESH_DIRECT_FLAG
  LOGICAL :: FIXED_WALL_NODES_STOKES_FLAG
  LOGICAL :: MOVED_WALL_NODES_STOKES_FLAG
  LOGICAL :: INLET_WALL_NODES_STOKES_FLAG
  LOGICAL :: ALE_SOLVER_STOKES_FLAG
  LOGICAL :: ANALYTIC_FLAG
  LOGICAL :: DYNAMIC_SOLVER_STOKES_FLAG
  LOGICAL :: FIXED_WALL_NODES_MOVING_MESH_FLAG
  LOGICAL :: MOVED_WALL_NODES_MOVING_MESH_FLAG

  CHARACTER *15 BUFFER
  CHARACTER *15 ARG
  CHARACTER *15 OUTPUT_STRING

  !CMISS variables

  !Regions
  TYPE(cmfe_RegionType) :: Region
  TYPE(cmfe_RegionType) :: WorldRegion
  !Coordinate systems
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem
  TYPE(cmfe_CoordinateSystemType) :: WorldCoordinateSystem
  !Basis
  TYPE(cmfe_BasisType) :: BasisSpace
  TYPE(cmfe_BasisType) :: BasisVelocity
  TYPE(cmfe_BasisType) :: BasisPressure
  !Nodes
  TYPE(cmfe_NodesType) :: Nodes
  !Elements
  TYPE(cmfe_MeshElementsType) :: MeshElementsSpace
  TYPE(cmfe_MeshElementsType) :: MeshElementsVelocity
  TYPE(cmfe_MeshElementsType) :: MeshElementsPressure
  !Meshes
  TYPE(cmfe_MeshType) :: Mesh
  !Decompositions
  TYPE(cmfe_DecompositionType) :: Decomposition
  !Fields
  TYPE(cmfe_FieldsType) :: Fields
  !Field types
  TYPE(cmfe_FieldType) :: GeometricField
  TYPE(cmfe_FieldType) :: AnalyticFieldStokes
  TYPE(cmfe_FieldType) :: DependentFieldStokes
  TYPE(cmfe_FieldType) :: DependentFieldMovingMesh
  TYPE(cmfe_FieldType) :: MaterialsFieldStokes
  TYPE(cmfe_FieldType) :: MaterialsFieldMovingMesh
  TYPE(cmfe_FieldType) :: IndependentFieldStokes
  TYPE(cmfe_FieldType) :: IndependentFieldMovingMesh
  TYPE(cmfe_FieldType) :: EquationsSetFieldStokes
  TYPE(cmfe_FieldType) :: EquationsSetFieldMovingMesh
  !Boundary conditions
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsStokes
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsMovingMesh
  !Equations sets
  TYPE(cmfe_EquationsSetType) :: EquationsSetStokes
  TYPE(cmfe_EquationsSetType) :: EquationsSetMovingMesh
  !Equations
  TYPE(cmfe_EquationsType) :: EquationsStokes
  TYPE(cmfe_EquationsType) :: EquationsMovingMesh
  !Problems
  TYPE(cmfe_ProblemType) :: Problem
  !Control loops
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  !Solvers
  TYPE(cmfe_SolverType) :: DynamicSolverStokes
  TYPE(cmfe_SolverType) :: LinearSolverStokes
  TYPE(cmfe_SolverType) :: LinearSolverMovingMesh
  !Solver equations
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsStokes
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsMovingMesh

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables

  INTEGER(CMFEIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,BoundaryNodeDomain
  INTEGER(CMFEIntg) :: EquationsSetIndex
  INTEGER(CMFEIntg) :: Err
  
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

  !Import cmHeart mesh information
  CALL FLUID_MECHANICS_IO_READ_CMHEART(CM,Err)  
  BASIS_NUMBER_SPACE=CM%ID_M
  BASIS_NUMBER_VELOCITY=CM%ID_V
  BASIS_NUMBER_PRESSURE=CM%ID_P
  NUMBER_OF_DIMENSIONS=CM%D
  BASIS_TYPE=CM%IT_T
  BASIS_XI_INTERPOLATION_SPACE=CM%IT_M
  BASIS_XI_INTERPOLATION_VELOCITY=CM%IT_V
  BASIS_XI_INTERPOLATION_PRESSURE=CM%IT_P
  NUMBER_OF_NODES_SPACE=CM%N_M
  NUMBER_OF_NODES_VELOCITY=CM%N_V
  NUMBER_OF_NODES_PRESSURE=CM%N_P
  TOTAL_NUMBER_OF_NODES=CM%N_T
  TOTAL_NUMBER_OF_ELEMENTS=CM%E_T
  NUMBER_OF_ELEMENT_NODES_SPACE=CM%EN_M
  NUMBER_OF_ELEMENT_NODES_VELOCITY=CM%EN_V
  NUMBER_OF_ELEMENT_NODES_PRESSURE=CM%EN_P
  !Set initial values
  INITIAL_FIELD_STOKES(1)=0.0_CMFEDP
  INITIAL_FIELD_STOKES(2)=0.0_CMFEDP
  INITIAL_FIELD_STOKES(3)=0.0_CMFEDP
  INITIAL_FIELD_MOVING_MESH(1)=0.0_CMFEDP
  INITIAL_FIELD_MOVING_MESH(2)=0.0_CMFEDP
  INITIAL_FIELD_MOVING_MESH(3)=0.0_CMFEDP
  
  !Set defaults
  !Set material parameters
  MU_PARAM_STOKES=1.0_CMFEDP
  RHO_PARAM_STOKES=1.0_CMFEDP
  LINEAR_SOLVER_STOKES_DIRECT_FLAG=.FALSE.
  LINEAR_SOLVER_MOVING_MESH_DIRECT_FLAG=.FALSE.
  DYNAMIC_SOLVER_STOKES_FLAG=.FALSE.
  ALE_SOLVER_STOKES_FLAG=.FALSE.
  ANALYTIC_FLAG=.FALSE.
  IF(NUMBER_OF_DIMENSIONS==2) THEN
    ANALYTIC_TYPE=CMFE_EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1
  ELSE
    ANALYTIC_TYPE=CMFE_EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1
  ENDIF
  OUTPUT_STRING='DEFAULT'
  DYNAMIC_SOLVER_STOKES_START_TIME=0.0_CMFEDP
  DYNAMIC_SOLVER_STOKES_STOP_TIME=3.0_CMFEDP 
  DYNAMIC_SOLVER_STOKES_TIME_INCREMENT=1.0_CMFEDP
  !Set initial boundary conditions
  BOUNDARY_CONDITIONS_STOKES=0.0_CMFEDP

  !Get command line arguments
  DO I=1,COMMAND_ARGUMENT_COUNT()
    CALL GET_COMMAND_ARGUMENT(I,ARG)
    SELECT CASE(ARG)
      CASE('-density')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ(BUFFER,*) RHO_PARAM_STOKES
      CASE('-viscosity')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ(BUFFER,*) MU_PARAM_STOKES
      CASE('-directsolver')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ(BUFFER,*) LINEAR_SOLVER_STOKES_DIRECT_FLAG
      CASE('-dynamic')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ(BUFFER,*) DYNAMIC_SOLVER_STOKES_FLAG
      CASE('-ALE')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ(BUFFER,*) ALE_SOLVER_STOKES_FLAG
      CASE('-analytic')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ(BUFFER,*) ANALYTIC_FLAG
      CASE('-analytictype')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ(BUFFER,*) ANALYTIC_TYPE
      CASE('-analyticoutput')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ(BUFFER,*) OUTPUT_STRING
      CASE('-starttime')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ(BUFFER,*) DYNAMIC_SOLVER_STOKES_START_TIME
      CASE('-stoptime')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ(BUFFER,*) DYNAMIC_SOLVER_STOKES_STOP_TIME
      CASE('-timeincrement')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ(BUFFER,*) DYNAMIC_SOLVER_STOKES_TIME_INCREMENT
      CASE('-velocity')
        CALL GET_COMMAND_ARGUMENT(I+1,BUFFER)
        READ(BUFFER,*) BOUNDARY_CONDITIONS_STOKES(1)
        CALL GET_COMMAND_ARGUMENT(I+2,BUFFER)
        READ(BUFFER,*) BOUNDARY_CONDITIONS_STOKES(2)
        CALL GET_COMMAND_ARGUMENT(I+3,BUFFER)
        READ(BUFFER,*) BOUNDARY_CONDITIONS_STOKES(3)
      CASE DEFAULT
        !do nothing
      END SELECT
  ENDDO 
  IF(NUMBER_OF_DIMENSIONS==2) THEN
    IF(DYNAMIC_SOLVER_STOKES_FLAG) THEN
      ANALYTIC_TYPE=CMFE_EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5
    ELSE
      ANALYTIC_TYPE=CMFE_EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1
    ENDIF
  ELSE
    IF(DYNAMIC_SOLVER_STOKES_FLAG) THEN
      ANALYTIC_TYPE=CMFE_EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5
    ELSE
      ANALYTIC_TYPE=CMFE_EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1
    ENDIF
  ENDIF

  WRITE(*,*)' '
  WRITE(*,*)' ************************************* '
  WRITE(*,*)' '
  WRITE(*,*)'-density........', RHO_PARAM_STOKES
  WRITE(*,*)'-viscosity......', MU_PARAM_STOKES
  WRITE(*,*)'-analytic.......  ', ANALYTIC_FLAG
  IF(.NOT.ANALYTIC_FLAG) THEN
    WRITE(*,*)'-velocity.......', BOUNDARY_CONDITIONS_STOKES
  ELSE
    IF(ANALYTIC_TYPE==CMFE_EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1.AND.NUMBER_OF_DIMENSIONS==2) THEN
      WRITE(*,*)'  -analytictype...', 'CMFE_EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1'
   ELSE IF(ANALYTIC_TYPE==CMFE_EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2.AND.NUMBER_OF_DIMENSIONS==2) THEN
      WRITE(*,*)'  -analytictype...', 'CMFE_EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2'
   ELSE IF(ANALYTIC_TYPE==CMFE_EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3.AND.NUMBER_OF_DIMENSIONS==2) THEN
      WRITE(*,*)'  -analytictype...', 'CMFE_EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3'
   ELSE IF(ANALYTIC_TYPE==CMFE_EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.AND.NUMBER_OF_DIMENSIONS==2) THEN
      WRITE(*,*)'  -analytictype...', 'CMFE_EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4'
    ELSE IF(ANALYTIC_TYPE==CMFE_EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5.AND.NUMBER_OF_DIMENSIONS==2) THEN
      WRITE(*,*)'  -analytictype...', 'CMFE_EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5'
    ELSE IF(ANALYTIC_TYPE==CMFE_EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1.AND.NUMBER_OF_DIMENSIONS==3) THEN
      WRITE(*,*)'  -analytictype...', 'CMFE_EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1'
    ELSE IF(ANALYTIC_TYPE==CMFE_EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2.AND.NUMBER_OF_DIMENSIONS==3) THEN
      WRITE(*,*)'  -analytictype...', 'CMFE_EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2'
    ELSE IF(ANALYTIC_TYPE==CMFE_EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3.AND.NUMBER_OF_DIMENSIONS==3) THEN
      WRITE(*,*)'  -analytictype...', 'CMFE_EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3'
    ELSE IF(ANALYTIC_TYPE==CMFE_EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.AND.NUMBER_OF_DIMENSIONS==3) THEN
      WRITE(*,*)'  -analytictype...', 'CMFE_EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4'
    ELSE IF(ANALYTIC_TYPE==CMFE_EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5.AND.NUMBER_OF_DIMENSIONS==3) THEN
      WRITE(*,*)'  -analytictype...', 'CMFE_EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5'
    ENDIF
   ENDIF 
  WRITE(*,*) ' '
  WRITE(*,*)'-dynamic........  ', DYNAMIC_SOLVER_STOKES_FLAG
  IF(DYNAMIC_SOLVER_STOKES_FLAG) THEN
    WRITE(*,*)'  -starttime......  ', DYNAMIC_SOLVER_STOKES_START_TIME
    WRITE(*,*)'  -stoptime.......  ', DYNAMIC_SOLVER_STOKES_STOP_TIME
    WRITE(*,*)'  -timeincrement..  ', DYNAMIC_SOLVER_STOKES_TIME_INCREMENT
    WRITE(*,*)'  -ALE............  ', ALE_SOLVER_STOKES_FLAG
    WRITE(*,*) ' ' 
  ENDIF
  WRITE(*,*)'-directsolver...  ', LINEAR_SOLVER_STOKES_DIRECT_FLAG
  WRITE(*,*)' '
  WRITE(*,*)' ************************************* '
  WRITE(*,*)' '
  WRITE(*,*) ' ' 
  !Set boundary conditions
  INQUIRE(FILE="./input/bc/FIXED_WALL", EXIST=FIXED_WALL_NODES_STOKES_FLAG)
  INQUIRE(FILE="./input/bc/FREE_INLET", EXIST=INLET_WALL_NODES_STOKES_FLAG)
  INQUIRE(FILE="./input/bc/MOVED_WALL", EXIST=MOVED_WALL_NODES_STOKES_FLAG)
  INQUIRE(FILE="./input/bc/FIXED_MESH", EXIST=FIXED_WALL_NODES_MOVING_MESH_FLAG)
  INQUIRE(FILE="./input/bc/MOVED_MESH", EXIST=MOVED_WALL_NODES_MOVING_MESH_FLAG)
  INITIAL_FIELD_STOKES(1)=BOUNDARY_CONDITIONS_STOKES(1)
  INITIAL_FIELD_STOKES(2)=BOUNDARY_CONDITIONS_STOKES(2)
  INITIAL_FIELD_STOKES(3)=BOUNDARY_CONDITIONS_STOKES(3)
  IF(FIXED_WALL_NODES_STOKES_FLAG) THEN
    OPEN(UNIT=1, FILE="./input/bc/FIXED_WALL",STATUS='unknown')
    READ(1,*) NUMBER_OF_FIXED_WALL_NODES_STOKES
    ALLOCATE(FIXED_WALL_NODES_STOKES(NUMBER_OF_FIXED_WALL_NODES_STOKES))
    READ(1,*) FIXED_WALL_NODES_STOKES(1:NUMBER_OF_FIXED_WALL_NODES_STOKES)
    CLOSE(1)
  ENDIF
  IF(MOVED_WALL_NODES_STOKES_FLAG) THEN
    OPEN(UNIT=1, FILE="./input/bc/MOVED_WALL",STATUS='unknown')
    READ(1,*) NUMBER_OF_MOVED_WALL_NODES_STOKES
    ALLOCATE(MOVED_WALL_NODES_STOKES(NUMBER_OF_MOVED_WALL_NODES_STOKES))
    READ(1,*) MOVED_WALL_NODES_STOKES(1:NUMBER_OF_MOVED_WALL_NODES_STOKES)
    CLOSE(1)
  ENDIF
  IF(INLET_WALL_NODES_STOKES_FLAG) THEN
    OPEN(UNIT=1, FILE="./input/bc/FREE_INLET",STATUS='unknown')
    READ(1,*) NUMBER_OF_INLET_WALL_NODES_STOKES
    ALLOCATE(INLET_WALL_NODES_STOKES(NUMBER_OF_INLET_WALL_NODES_STOKES))
    READ(1,*) INLET_WALL_NODES_STOKES(1:NUMBER_OF_INLET_WALL_NODES_STOKES)
    CLOSE(1)
  ENDIF
  IF(FIXED_WALL_NODES_MOVING_MESH_FLAG) THEN
    OPEN(UNIT=1, FILE="./input/bc/FIXED_MESH",STATUS='unknown')
    READ(1,*) NUMBER_OF_FIXED_WALL_NODES_MOVING_MESH
    ALLOCATE(FIXED_WALL_NODES_MOVING_MESH(NUMBER_OF_FIXED_WALL_NODES_MOVING_MESH))
    READ(1,*) FIXED_WALL_NODES_MOVING_MESH(1:NUMBER_OF_FIXED_WALL_NODES_MOVING_MESH)
    CLOSE(1)
  ENDIF
  IF(MOVED_WALL_NODES_MOVING_MESH_FLAG) THEN
    OPEN(UNIT=1, FILE="./input/bc/MOVED_MESH",STATUS='unknown')
    READ(1,*) NUMBER_OF_MOVED_WALL_NODES_MOVING_MESH
    ALLOCATE(MOVED_WALL_NODES_MOVING_MESH(NUMBER_OF_MOVED_WALL_NODES_MOVING_MESH))
    READ(1,*) MOVED_WALL_NODES_MOVING_MESH(1:NUMBER_OF_MOVED_WALL_NODES_MOVING_MESH)
    CLOSE(1)
    BOUNDARY_CONDITIONS_MOVING_MESH(1)=0.0_CMFEDP
    BOUNDARY_CONDITIONS_MOVING_MESH(2)=0.0_CMFEDP
    BOUNDARY_CONDITIONS_MOVING_MESH(3)=0.0_CMFEDP
  ENDIF
  !Set interpolation parameters
  BASIS_GAUSS_SPACE=4
  BASIS_GAUSS_VELOCITY=4
  BASIS_GAUSS_PRESSURE=4
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  DYNAMIC_SOLVER_STOKES_OUTPUT_TYPE=CMFE_SOLVER_NO_OUTPUT
  LINEAR_SOLVER_MOVING_MESH_OUTPUT_TYPE=CMFE_SOLVER_NO_OUTPUT
  LINEAR_SOLVER_STOKES_OUTPUT_TYPE=CMFE_SOLVER_NO_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_STOKES_OUTPUT=CMFE_EQUATIONS_NO_OUTPUT
  !Set time parameter
  DYNAMIC_SOLVER_STOKES_THETA=1.0_CMFEDP
  !Set result output parameter
  DYNAMIC_SOLVER_STOKES_OUTPUT_FREQUENCY=1
  !Set solver parameters
  RELATIVE_TOLERANCE=1.0E-10_CMFEDP !default: 1.0E-05_CMFEDP
  ABSOLUTE_TOLERANCE=1.0E-10_CMFEDP !default: 1.0E-10_CMFEDP
  DIVERGENCE_TOLERANCE=1.0E20 !default: 1.0E5
  MAXIMUM_ITERATIONS=100000 !default: 100000
  RESTART_VALUE=3000 !default: 30
  LINESEARCH_ALPHA=1.0

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

  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !REGION

  !Start the creation of a new region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system as defined above
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !BASES

  !Start the creation of new bases
  MESH_NUMBER_OF_COMPONENTS=1
  CALL cmfe_Basis_Initialise(BasisSpace,Err)
  CALL cmfe_Basis_CreateStart(BASIS_NUMBER_SPACE,BasisSpace,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL cmfe_Basis_TypeSet(BasisSpace,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL cmfe_Basis_NumberOfXiSet(BasisSpace,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  IF(NUMBER_OF_DIMENSIONS==2) THEN
    CALL cmfe_Basis_InterpolationXiSet(BasisSpace,(/BASIS_XI_INTERPOLATION_SPACE,BASIS_XI_INTERPOLATION_SPACE/),Err)
    IF(BASIS_TYPE/=CMFE_BASIS_SIMPLEX_TYPE) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisSpace,(/BASIS_GAUSS_SPACE,BASIS_GAUSS_SPACE/),Err)
    ELSE
      CALL cmfe_Basis_QuadratureOrderSet(BasisSpace,BASIS_GAUSS_SPACE+1,Err)
    ENDIF
  ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
    CALL cmfe_Basis_InterpolationXiSet(BasisSpace,(/BASIS_XI_INTERPOLATION_SPACE,BASIS_XI_INTERPOLATION_SPACE, & 
      & BASIS_XI_INTERPOLATION_SPACE/),Err)                         
    IF(BASIS_TYPE/=CMFE_BASIS_SIMPLEX_TYPE) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisSpace,(/BASIS_GAUSS_SPACE,BASIS_GAUSS_SPACE,BASIS_GAUSS_SPACE/), & 
        & Err)
    ELSE
      CALL cmfe_Basis_QuadratureOrderSet(BasisSpace,BASIS_GAUSS_SPACE+1,Err)
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(BasisSpace,Err)
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisVelocity=BasisSpace
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new velocity basis
    CALL cmfe_Basis_Initialise(BasisVelocity,Err)
    !Start the creation of a basis
    CALL cmfe_Basis_CreateStart(BASIS_NUMBER_VELOCITY,BasisVelocity,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL cmfe_Basis_TypeSet(BasisVelocity,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL cmfe_Basis_NumberOfXiSet(BasisVelocity,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisVelocity,(/BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY/),Err)
      IF(BASIS_TYPE/=CMFE_BASIS_SIMPLEX_TYPE) THEN 
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisVelocity,(/BASIS_GAUSS_VELOCITY,BASIS_GAUSS_VELOCITY/),Err)
      ELSE
        CALL cmfe_Basis_QuadratureOrderSet(BasisVelocity,BASIS_GAUSS_VELOCITY+1,Err)
      ENDIF
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisVelocity,(/BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY, & 
        & BASIS_XI_INTERPOLATION_VELOCITY/),Err)                         
      IF(BASIS_TYPE/=CMFE_BASIS_SIMPLEX_TYPE) THEN
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisVelocity,(/BASIS_GAUSS_VELOCITY,BASIS_GAUSS_VELOCITY, & 
          & BASIS_GAUSS_VELOCITY/),Err)
      ELSE
        CALL cmfe_Basis_QuadratureOrderSet(BasisVelocity,BASIS_GAUSS_VELOCITY+1,Err)
      ENDIF
    ENDIF
    !Finish the creation of the basis
    CALL cmfe_Basis_CreateFinish(BasisVelocity,Err)
  ENDIF
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisPressure=BasisSpace
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_VELOCITY) THEN
    BasisPressure=BasisVelocity
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new pressure basis
    CALL cmfe_Basis_Initialise(BasisPressure,Err)
    !Start the creation of a basis
    CALL cmfe_Basis_CreateStart(BASIS_NUMBER_PRESSURE,BasisPressure,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL cmfe_Basis_TypeSet(BasisPressure,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL cmfe_Basis_NumberOfXiSet(BasisPressure,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisPressure,(/BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE/),Err)
      IF(BASIS_TYPE/=CMFE_BASIS_SIMPLEX_TYPE) THEN
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisPressure,(/BASIS_GAUSS_PRESSURE,BASIS_GAUSS_PRESSURE/),Err)
      ELSE
        CALL cmfe_Basis_QuadratureOrderSet(BasisPressure,BASIS_GAUSS_PRESSURE+1,Err)
      ENDIF
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisPressure,(/BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE, & 
        & BASIS_XI_INTERPOLATION_PRESSURE/),Err)                         
      IF(BASIS_TYPE/=CMFE_BASIS_SIMPLEX_TYPE) THEN
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisPressure,(/BASIS_GAUSS_PRESSURE,BASIS_GAUSS_PRESSURE, & 
          & BASIS_GAUSS_PRESSURE/),Err)
      ELSE
        CALL cmfe_Basis_QuadratureOrderSet(BasisPressure,BASIS_GAUSS_PRESSURE+1,Err)
      ENDIF
    ENDIF
    !Finish the creation of the basis
    CALL cmfe_Basis_CreateFinish(BasisPressure,Err)
  ENDIF

  !
  !================================================================================================================================
  !

  !MESH

  !Start the creation of mesh nodes
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Nodes_CreateStart(Region,TOTAL_NUMBER_OF_NODES,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,MESH_NUMBER_OF_COMPONENTS,Err)
  !Specify spatial mesh component
  CALL cmfe_MeshElements_Initialise(MeshElementsSpace,Err)
  CALL cmfe_MeshElements_Initialise(MeshElementsVelocity,Err)
  CALL cmfe_MeshElements_Initialise(MeshElementsPressure,Err)
  MESH_COMPONENT_NUMBER_SPACE=1
  MESH_COMPONENT_NUMBER_VELOCITY=1
  MESH_COMPONENT_NUMBER_PRESSURE=1
  CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_SPACE,BasisSpace,MeshElementsSpace,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL cmfe_MeshElements_NodesSet(MeshElementsSpace,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_SPACE),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(MeshElementsSpace,Err)
  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsVelocity=MeshElementsSpace
  ELSE
    MESH_COMPONENT_NUMBER_VELOCITY=MESH_COMPONENT_NUMBER_SPACE+1
    CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_VELOCITY,BasisVelocity,MeshElementsVelocity,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL cmfe_MeshElements_NodesSet(MeshElementsVelocity,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_VELOCITY),Err)
    ENDDO
    CALL cmfe_MeshElements_CreateFinish(MeshElementsVelocity,Err)
  ENDIF
  !Specify pressure mesh component
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsPressure=MeshElementsSpace
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_SPACE
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_VELOCITY) THEN
    MeshElementsPressure=MeshElementsVelocity
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_VELOCITY
  ELSE
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_VELOCITY+1
    CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_PRESSURE,BasisPressure,MeshElementsPressure,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL cmfe_MeshElements_NodesSet(MeshElementsPressure,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_PRESSURE),Err)
    ENDDO
    CALL cmfe_MeshElements_CreateFinish(MeshElementsPressure,Err)
  ENDIF
  !Finish the creation of the mesh
  CALL cmfe_Mesh_CreateFinish(Mesh,Err)

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the field type
  CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL cmfe_Field_ScalingTypeSet(GeometricField,CMFE_FIELD_NO_SCALING,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDDO
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_SPACE
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
          & 1,CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
      ENDIF
    ENDDO
  ENDDO
  CALL cmfe_Field_ParameterSetUpdateStart(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for static/dynamic Stokes
  CALL cmfe_EquationsSet_Initialise(EquationsSetStokes,Err)
  CALL cmfe_Field_Initialise(EquationsSetFieldStokes,Err)
  !Set the equations set to be a static/dynamic Stokes problem
  IF(DYNAMIC_SOLVER_STOKES_FLAG) THEN
    IF(ALE_SOLVER_STOKES_FLAG) THEN
      CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberStokes,Region,GeometricField, &
        & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
        & CMFE_EQUATIONS_SET_STOKES_EQUATION_TYPE,CMFE_EQUATIONS_SET_ALE_STOKES_SUBTYPE],EquationsSetFieldUserNumber, &
        & EquationsSetFieldStokes,EquationsSetStokes,Err)
    ELSE
      CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberStokes,Region,GeometricField, &
        & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
        & CMFE_EQUATIONS_SET_STOKES_EQUATION_TYPE,CMFE_EQUATIONS_SET_TRANSIENT_STOKES_SUBTYPE],EquationsSetFieldUserNumber, &
        & EquationsSetFieldStokes,EquationsSetStokes,Err)
    ENDIF
  ELSE
    CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberStokes,Region,GeometricField, &
      & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
      & CMFE_EQUATIONS_SET_STOKES_EQUATION_TYPE,CMFE_EQUATIONS_SET_STATIC_STOKES_SUBTYPE],EquationsSetFieldUserNumber, &
      & EquationsSetFieldStokes,EquationsSetStokes,Err)
  ENDIF
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetStokes,Err)

  IF(ALE_SOLVER_STOKES_FLAG) THEN
    !Create the equations set for moving mesh
    CALL cmfe_Field_Initialise(EquationsSetFieldMovingMesh,Err)
    CALL cmfe_EquationsSet_Initialise(EquationsSetMovingMesh,Err)
    !Set the equations set to be a moving mesh problem
    CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberMovingMesh,Region,GeometricField, &
      & [CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
      & CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMFE_EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE],&
      & EquationsSetFieldUserNumberMovingMesh,EquationsSetFieldMovingMesh,EquationsSetMovingMesh,Err)
    !Finish creating the equations set
    CALL cmfe_EquationsSet_CreateFinish(EquationsSetMovingMesh,Err)
  ENDIF


  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for static/dynamic Stokes
  CALL cmfe_Field_Initialise(DependentFieldStokes,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetStokes,DependentFieldUserNumberStokes, & 
    & DependentFieldStokes,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldStokes,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldStokes,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
  ENDDO
  COMPONENT_NUMBER=NUMBER_OF_DIMENSIONS+1
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldStokes,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldStokes,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetStokes,Err)

  !Initialise dependent field
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
      & COMPONENT_NUMBER,INITIAL_FIELD_STOKES(COMPONENT_NUMBER),Err)
  ENDDO

  IF(ALE_SOLVER_STOKES_FLAG) THEN
    !Create the equations set dependent field variables for moving mesh
    CALL cmfe_Field_Initialise(DependentFieldMovingMesh,Err)
    CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetMovingMesh,DependentFieldUserNumberMovingMesh, & 
      & DependentFieldMovingMesh,Err)
    !Set the mesh component to be used by the field components.
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
        & MESH_COMPONENT_NUMBER_SPACE,Err)
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldMovingMesh,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
        & MESH_COMPONENT_NUMBER_SPACE,Err)
    ENDDO
    !Finish the equations set dependent field variables
    CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetMovingMesh,Err)
 
    !Initialise dependent field
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      CALL cmfe_Field_ComponentValuesInitialise(DependentFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & COMPONENT_NUMBER,INITIAL_FIELD_MOVING_MESH(COMPONENT_NUMBER),Err)
    ENDDO
  ENDIF

  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for static/dynamic Stokes
  CALL cmfe_Field_Initialise(MaterialsFieldStokes,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetStokes,MaterialsFieldUserNumberStokes, & 
    & MaterialsFieldStokes,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetStokes,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberStokesMu,MU_PARAM_STOKES,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberStokesRho,RHO_PARAM_STOKES,Err)

  IF(ALE_SOLVER_STOKES_FLAG) THEN
    !Create the equations set materials field variables for moving mesh
    CALL cmfe_Field_Initialise(MaterialsFieldMovingMesh,Err)
    CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetMovingMesh,MaterialsFieldUserNumberMovingMesh, & 
      & MaterialsFieldMovingMesh,Err)
    !Finish the equations set materials field variables
    CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetMovingMesh,Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
      & MaterialsFieldUserNumberMovingMeshK,K_PARAM_MOVING_MESH,Err)
  ENDIF

  !
  !================================================================================================================================
  !

  !INDEPENDENT FIELDS

  IF(ALE_SOLVER_STOKES_FLAG) THEN
    !Create the equations set independent field variables for ALE Stokes
    CALL cmfe_Field_Initialise(IndependentFieldStokes,Err)
    CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetStokes,IndependentFieldUserNumberStokes, & 
      & IndependentFieldStokes,Err)
    !Set the mesh component to be used by the field components.
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      CALL cmfe_Field_ComponentMeshComponentSet(InDependentFieldStokes,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
        & MESH_COMPONENT_NUMBER_SPACE,Err)
    ENDDO
    !Finish the equations set independent field variables
    CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetStokes,Err)
    !Create the equations set independent field variables for moving mesh
    CALL cmfe_Field_Initialise(IndependentFieldMovingMesh,Err)
    CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetMovingMesh,IndependentFieldUserNumberMovingMesh, & 
      & IndependentFieldMovingMesh,Err)
    !Set the mesh component to be used by the field components.
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      CALL cmfe_Field_ComponentMeshComponentSet(InDependentFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
        & MESH_COMPONENT_NUMBER_SPACE,Err)
    ENDDO
    !Finish the equations set independent field variables
    CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetMovingMesh,Err)
  ENDIF

  !
  !================================================================================================================================
  !

  IF(ANALYTIC_FLAG) THEN
  !ANALYTIC FIELDS
  !Create the equations set analytic field variables for static/dynamic Stokes
    CALL cmfe_Field_Initialise(AnalyticFieldStokes,Err)
    IF(NUMBER_OF_DIMENSIONS==2) THEN  
      CALL cmfe_EquationsSet_AnalyticCreateStart(EquationsSetStokes,ANALYTIC_TYPE,AnalyticFieldUserNumberStokes, &
        & AnalyticFieldStokes,Err)
    ELSE
      CALL cmfe_EquationsSet_AnalyticCreateStart(EquationsSetStokes,ANALYTIC_TYPE,AnalyticFieldUserNumberStokes, &
        & AnalyticFieldStokes,Err)
    ENDIF
    !Finish the equations set analytic field variables
    CALL cmfe_EquationsSet_AnalyticCreateFinish(EquationsSetStokes,Err)
  ENDIF 

  !
  !================================================================================================================================
  !

  !EQUATIONS

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(EquationsStokes,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetStokes,EquationsStokes,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(EquationsStokes,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  IF(DYNAMIC_SOLVER_STOKES_FLAG) THEN
    !Set the equations lumping type
    CALL cmfe_Equations_LumpingTypeSet(EquationsStokes,CMFE_EQUATIONS_UNLUMPED_MATRICES,Err)
  ENDIF
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(EquationsStokes,EQUATIONS_STOKES_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetStokes,Err)

  IF(ALE_SOLVER_STOKES_FLAG) THEN
    !Create the equations set equations
    CALL cmfe_Equations_Initialise(EquationsMovingMesh,Err)
    CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetMovingMesh,EquationsMovingMesh,Err)
    !Set the equations matrices sparsity type
    CALL cmfe_Equations_SparsityTypeSet(EquationsMovingMesh,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
    !Set the equations set output
    CALL cmfe_Equations_OutputTypeSet(EquationsMovingMesh,EQUATIONS_MOVING_MESH_OUTPUT,Err)
    !Finish the equations set equations
    CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetMovingMesh,Err)
  ENDIF

  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  !Set the problem to be a static/dynamic Stokes problem
  IF(DYNAMIC_SOLVER_STOKES_FLAG) THEN
    IF(ALE_SOLVER_STOKES_FLAG) THEN
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_FLUID_MECHANICS_CLASS,CMFE_PROBLEM_STOKES_EQUATION_TYPE, &
        & CMFE_PROBLEM_ALE_STOKES_SUBTYPE],Problem,Err)
    ELSE  
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_FLUID_MECHANICS_CLASS,CMFE_PROBLEM_STOKES_EQUATION_TYPE, &
        & CMFE_PROBLEM_TRANSIENT_STOKES_SUBTYPE],Problem,Err)
    ENDIF
  ELSE
    CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_FLUID_MECHANICS_CLASS,CMFE_PROBLEM_STOKES_EQUATION_TYPE, &
      & CMFE_PROBLEM_STATIC_STOKES_SUBTYPE],Problem,Err)
  ENDIF
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  IF(DYNAMIC_SOLVER_STOKES_FLAG) THEN
    !Get the control loop
    CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
   !Set the times
    CALL cmfe_ControlLoop_TimesSet(ControlLoop,DYNAMIC_SOLVER_STOKES_START_TIME,DYNAMIC_SOLVER_STOKES_STOP_TIME, & 
      & DYNAMIC_SOLVER_STOKES_TIME_INCREMENT,Err)
    !Set the output timing
    CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,DYNAMIC_SOLVER_STOKES_OUTPUT_FREQUENCY,Err)
    IF(ALE_SOLVER_STOKES_FLAG) THEN
      CALL cmfe_ControlLoop_TimeInputSet(ControlLoop,DYNAMIC_SOLVER_STOKES_INPUT_OPTION,Err)
    ENDIF
  ENDIF
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)


  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(LinearSolverMovingMesh,Err)
  CALL cmfe_Solver_Initialise(DynamicSolverStokes,Err)
  CALL cmfe_Solver_Initialise(LinearSolverStokes,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  IF(ALE_SOLVER_STOKES_FLAG) THEN
    !Get the moving mesh solver
    CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverMovingMeshUserNumber,LinearSolverMovingMesh,Err)
    !Set the output type
    !CALL cmfe_Solver_OutputTypeSet(LinearSolverMovingMesh,LINEAR_SOLVER_MOVING_MESH_OUTPUT_TYPE,Err)
    CALL cmfe_Solver_OutputTypeSet(LinearSolverMovingMesh,4,Err)
    !Set the solver settings
    IF(LINEAR_SOLVER_MOVING_MESH_DIRECT_FLAG) THEN
      CALL cmfe_Solver_LinearTypeSet(LinearSolverMovingMesh,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
      CALL cmfe_Solver_LibraryTypeSet(LinearSolverMovingMesh,CMFE_SOLVER_MUMPS_LIBRARY,Err)
    ELSE
      CALL cmfe_Solver_LinearTypeSet(LinearSolverMovingMesh,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
      CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolverMovingMesh,MAXIMUM_ITERATIONS,Err)
      CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(LinearSolverMovingMesh,DIVERGENCE_TOLERANCE,Err)
      CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolverMovingMesh,RELATIVE_TOLERANCE,Err)
      CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolverMovingMesh,ABSOLUTE_TOLERANCE,Err)
      CALL cmfe_Solver_LinearIterativeGMRESRestartSet(LinearSolverMovingMesh,RESTART_VALUE,Err)
    ENDIF
  ENDIF

  IF(DYNAMIC_SOLVER_STOKES_FLAG) THEN
    !Get the dynamic dymamic solver
    CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverStokesUserNumber,DynamicSolverStokes,Err)
    !Set the output type
    CALL cmfe_Solver_OutputTypeSet(DynamicSolverStokes,DYNAMIC_SOLVER_STOKES_OUTPUT_TYPE,Err)
    !Set theta
    CALL cmfe_Solver_DynamicThetaSet(DynamicSolverStokes,DYNAMIC_SOLVER_STOKES_THETA,Err)
  !   CALL cmfe_SolverDynamicDynamicSet(DynamicSolverStokes,.TRUE.,Err)
    !Get the dynamic linear solver
    CALL cmfe_Solver_DynamicLinearSolverGet(DynamicSolverStokes,LinearSolverStokes,Err)
    !Set the output type
    CALL cmfe_Solver_OutputTypeSet(LinearSolverStokes,LINEAR_SOLVER_STOKES_OUTPUT_TYPE,Err)
    !Set the solver settings
  ELSE
    !Get the linear static solver
    CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverStokesUserNumber,LinearSolverStokes,Err)
    !Set the output type
    CALL cmfe_Solver_OutputTypeSet(LinearSolverStokes,LINEAR_SOLVER_STOKES_OUTPUT_TYPE,Err)
  ENDIF
  !Set the solver settings
  IF(LINEAR_SOLVER_STOKES_DIRECT_FLAG) THEN
    CALL cmfe_Solver_LinearTypeSet(LinearSolverStokes,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LibraryTypeSet(LinearSolverStokes,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(LinearSolverStokes,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolverStokes,MAXIMUM_ITERATIONS,Err)
    CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(LinearSolverStokes,DIVERGENCE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolverStokes,RELATIVE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolverStokes,ABSOLUTE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeGMRESRestartSet(LinearSolverStokes,RESTART_VALUE,Err)
  ENDIF
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS
  CALL cmfe_Solver_Initialise(DynamicSolverStokes,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsStokes,Err)
  CALL cmfe_Solver_Initialise(LinearSolverMovingMesh,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsMovingMesh,Err)
  CALL cmfe_Solver_Initialise(LinearSolverStokes,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsStokes,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)

  IF(DYNAMIC_SOLVER_STOKES_FLAG) THEN
    !Start the creation of the problem solver equations
    IF(ALE_SOLVER_STOKES_FLAG) THEN
      !Get the linear solver equations
      CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverMovingMeshUserNumber,LinearSolverMovingMesh,Err)
      CALL cmfe_Solver_SolverEquationsGet(LinearSolverMovingMesh,SolverEquationsMovingMesh,Err)
      !Set the solver equations sparsity
      CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsMovingMesh,CMFE_SOLVER_SPARSE_MATRICES,Err)
      !Add in the equations set
      CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsMovingMesh,EquationsSetMovingMesh,EquationsSetIndex,Err)
      !Finish the creation of the problem solver equations
      !Get the dynamic solver equations
      CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverStokesUserNumber,DynamicSolverStokes,Err)
      CALL cmfe_Solver_SolverEquationsGet(DynamicSolverStokes,SolverEquationsStokes,Err)
    ELSE
      !Get the dynamic solver equations
      CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverStokesUserNumber,DynamicSolverStokes,Err)
      CALL cmfe_Solver_SolverEquationsGet(DynamicSolverStokes,SolverEquationsStokes,Err)
    ENDIF
  ELSE
    !Start the creation of the problem solver equations
    !Get the linear solver equations
    CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverStokesUserNumber,LinearSolverStokes,Err)
    CALL cmfe_Solver_SolverEquationsGet(LinearSolverStokes,SolverEquationsStokes,Err)
  ENDIF
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsStokes,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsStokes,EquationsSetStokes,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  !Start the creation of the equations set boundary conditions for Stokes
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsStokes,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsStokes,BoundaryConditionsStokes,Err)
  IF(.NOT.ANALYTIC_FLAG) THEN
    !Set fixed wall nodes
    IF(FIXED_WALL_NODES_STOKES_FLAG) THEN
      DO NODE_COUNTER=1,NUMBER_OF_FIXED_WALL_NODES_STOKES
        NODE_NUMBER=FIXED_WALL_NODES_STOKES(NODE_COUNTER)
        CONDITION=CMFE_BOUNDARY_CONDITION_FIXED_WALL
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
        IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
          DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
            VALUE=0.0_CMFEDP
            CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsStokes,DependentFieldStokes,CMFE_FIELD_U_VARIABLE_TYPE,1, &
              & CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    !Set moved wall nodes
    IF(MOVED_WALL_NODES_STOKES_FLAG) THEN
      DO NODE_COUNTER=1,NUMBER_OF_MOVED_WALL_NODES_STOKES
        NODE_NUMBER=MOVED_WALL_NODES_STOKES(NODE_COUNTER)
        CONDITION=CMFE_BOUNDARY_CONDITION_MOVED_WALL
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
        IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
          DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
            VALUE=0.0_CMFEDP
            CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsStokes,DependentFieldStokes,CMFE_FIELD_U_VARIABLE_TYPE,1, &
              & CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    !Set velocity boundary conditions
    IF(INLET_WALL_NODES_STOKES_FLAG) THEN
      DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_STOKES
        NODE_NUMBER=INLET_WALL_NODES_STOKES(NODE_COUNTER)
        CONDITION=CMFE_BOUNDARY_CONDITION_FIXED_INLET
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
        IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
          DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
            VALUE=BOUNDARY_CONDITIONS_STOKES(COMPONENT_NUMBER)
            CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsStokes,DependentFieldStokes,CMFE_FIELD_U_VARIABLE_TYPE,1, &
              & CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
          ENDDO
        ENDIF
      ENDDO
    ENDIF

    IF(ALE_SOLVER_STOKES_FLAG) THEN
      !Start the creation of the equations set boundary conditions for moving mesh
      CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsMovingMesh,Err)
      CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsMovingMesh,BoundaryConditionsMovingMesh,Err)
      !Set fixed wall nodes
      IF(FIXED_WALL_NODES_MOVING_MESH_FLAG) THEN
        DO NODE_COUNTER=1,NUMBER_OF_FIXED_WALL_NODES_MOVING_MESH
          NODE_NUMBER=FIXED_WALL_NODES_MOVING_MESH(NODE_COUNTER)
          CONDITION=CMFE_BOUNDARY_CONDITION_FIXED_WALL
          CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
          IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
            DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
              VALUE=0.0_CMFEDP
              CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
                & CMFE_FIELD_U_VARIABLE_TYPE,1, &
                & CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      !Set moved wall nodes
      IF(MOVED_WALL_NODES_MOVING_MESH_FLAG) THEN
        DO NODE_COUNTER=1,NUMBER_OF_MOVED_WALL_NODES_MOVING_MESH
          NODE_NUMBER=MOVED_WALL_NODES_MOVING_MESH(NODE_COUNTER)
          CONDITION=CMFE_BOUNDARY_CONDITION_MOVED_WALL
          CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
          IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
            DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
              VALUE=BOUNDARY_CONDITIONS_MOVING_MESH(COMPONENT_NUMBER)
              CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
                & CMFE_FIELD_U_VARIABLE_TYPE,1, &
                & CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      !Finish the creation of the equations set boundary conditions
      CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsMovingMesh,Err)
    ENDIF
  ELSE
    !Set up the boundary conditions as per the analytic solution
    CALL cmfe_SolverEquations_BoundaryConditionsAnalytic(SolverEquationsStokes,Err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsStokes,Err)

  !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn of PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL cmfe_Problem_Solve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"

  !
  !================================================================================================================================
  !

  !OUTPUT

  IF(ANALYTIC_FLAG) THEN
    !Output Analytic analysis
    CALL cmfe_AnalyticAnalysisOutput(DependentFieldStokes,OUTPUT_STRING,Err)
  ENDIF

  EXPORT_FIELD_IO=.FALSE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"MasterStokes","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"MasterStokes","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF
  
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM STOKESMASTEREXAMPLE
