!> \file
!> \author Sebastian Krittian
!> \brief This is an example program to solve a ALE Stokes equation using OpenCMISS calls.
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

!> \example FluidMechanics/Stokes/RoutineCheck/ALE/src/ALEExample.f90
!! Example program to solve a ALE Stokes equation using OpenCMISS calls.
!! \htmlinclude FluidMechanics/Stokes/RoutineCheck/ALE/history.html
!!
!<

!> Main program

PROGRAM STOKESALEEXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OpenCMISS
  USE OpenCMISS_Iron
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


  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberStokes=6
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberMovingMesh=42
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberStokes=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMovingMesh=9
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberStokes=10
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberMovingMesh=11
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberStokes=12
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberMovingMesh=13
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberStokes=22
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberMovingMesh=23
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=14

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverMovingMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverStokesUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberStokesMu=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberStokesRho=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMovingMeshK=1

  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS
  
  INTEGER(CMISSIntg) :: BASIS_TYPE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_SPACE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_PRESSURE
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_SPACE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_PRESSURE
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES_STOKES
  INTEGER(CMISSIntg) :: NUMBER_OF_MOVED_WALL_NODES_STOKES
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_WALL_NODES_STOKES
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES_MOVING_MESH
  INTEGER(CMISSIntg) :: NUMBER_OF_MOVED_WALL_NODES_MOVING_MESH

  INTEGER(CMISSIntg) :: EQUATIONS_STOKES_OUTPUT
  INTEGER(CMISSIntg) :: EQUATIONS_MOVING_MESH_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_COUNTER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_STOKES_INPUT_OPTION
  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_STOKES_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_MOVING_MESH_OUTPUT_TYPE

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: MOVED_WALL_NODES_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_MOVING_MESH
  INTEGER, ALLOCATABLE, DIMENSION(:):: MOVED_WALL_NODES_MOVING_MESH

  REAL(CMISSRP) :: INITIAL_FIELD_STOKES(3)
  REAL(CMISSRP) :: INITIAL_FIELD_MOVING_MESH(3)
  REAL(CMISSRP) :: BOUNDARY_CONDITIONS_STOKES(3)
  REAL(CMISSRP) :: BOUNDARY_CONDITIONS_MOVING_MESH(3)
  REAL(CMISSRP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSRP) :: RELATIVE_TOLERANCE
  REAL(CMISSRP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSRP) :: LINESEARCH_ALPHA
  REAL(CMISSRP) :: VALUE
  REAL(CMISSRP) :: K_PARAM_MOVING_MESH
  REAL(CMISSRP) :: MU_PARAM_STOKES
  REAL(CMISSRP) :: RHO_PARAM_STOKES

  REAL(CMISSRP) :: DYNAMIC_SOLVER_STOKES_START_TIME
  REAL(CMISSRP) :: DYNAMIC_SOLVER_STOKES_STOP_TIME
  REAL(CMISSRP) :: DYNAMIC_SOLVER_STOKES_THETA
  REAL(CMISSRP) :: DYNAMIC_SOLVER_STOKES_TIME_INCREMENT

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_STOKES_DIRECT_FLAG
  LOGICAL :: LINEAR_SOLVER_MOVING_MESH_DIRECT_FLAG
  LOGICAL :: FIXED_WALL_NODES_STOKES_FLAG
  LOGICAL :: MOVED_WALL_NODES_STOKES_FLAG
  LOGICAL :: INLET_WALL_NODES_STOKES_FLAG
  LOGICAL :: FIXED_WALL_NODES_MOVING_MESH_FLAG
  LOGICAL :: MOVED_WALL_NODES_MOVING_MESH_FLAG


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
  INITIAL_FIELD_STOKES(1)=0.0_CMISSRP
  INITIAL_FIELD_STOKES(2)=0.0_CMISSRP
  INITIAL_FIELD_STOKES(3)=0.0_CMISSRP
  INITIAL_FIELD_MOVING_MESH(1)=0.0_CMISSRP
  INITIAL_FIELD_MOVING_MESH(2)=0.0_CMISSRP
  INITIAL_FIELD_MOVING_MESH(3)=0.0_CMISSRP
  !Set boundary conditions
  FIXED_WALL_NODES_STOKES_FLAG=.TRUE.
  MOVED_WALL_NODES_STOKES_FLAG=.TRUE.
  INLET_WALL_NODES_STOKES_FLAG=.FALSE.
  FIXED_WALL_NODES_MOVING_MESH_FLAG=.TRUE.
  MOVED_WALL_NODES_MOVING_MESH_FLAG=.TRUE.
  IF(FIXED_WALL_NODES_STOKES_FLAG) THEN
    NUMBER_OF_FIXED_WALL_NODES_STOKES=16
    ALLOCATE(FIXED_WALL_NODES_STOKES(NUMBER_OF_FIXED_WALL_NODES_STOKES))
    FIXED_WALL_NODES_STOKES=[46,47,48,53,57,64,65,68,72,106,107,111,117,118,122,125]
  ENDIF
  IF(MOVED_WALL_NODES_STOKES_FLAG) THEN
    NUMBER_OF_MOVED_WALL_NODES_STOKES=73
    ALLOCATE(MOVED_WALL_NODES_STOKES(NUMBER_OF_MOVED_WALL_NODES_STOKES))
    MOVED_WALL_NODES_STOKES=[3,4,7,10,11,12,13,17,20,24,29,31,33,34,35,39,41,44,50,51,52,54,60,66,67,70,74,78,79,83, &
      & 86,90,91,92,93,95,99,101,103,104,105,108,114,115,116,120,123,124,1,2,5,6,9,14,15,16,23,28,30,32,36, & 
      & 37,42,76,77,80,81,82,89,94,96,97,102]
  ENDIF
  IF(INLET_WALL_NODES_STOKES_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_STOKES=25
    ALLOCATE(INLET_WALL_NODES_STOKES(NUMBER_OF_INLET_WALL_NODES_STOKES))
    INLET_WALL_NODES_STOKES=[46,47,48,49,53,57,58,59,63,64,65,68,71,72,75,106,107,111,112,113,117,118,121,122,125]
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_STOKES(1)=0.0_CMISSRP
    BOUNDARY_CONDITIONS_STOKES(2)=0.0_CMISSRP
    BOUNDARY_CONDITIONS_STOKES(3)=0.0_CMISSRP
  ENDIF
  IF(FIXED_WALL_NODES_MOVING_MESH_FLAG) THEN
    NUMBER_OF_FIXED_WALL_NODES_MOVING_MESH=25
    ALLOCATE(FIXED_WALL_NODES_MOVING_MESH(NUMBER_OF_FIXED_WALL_NODES_MOVING_MESH))
    FIXED_WALL_NODES_MOVING_MESH=[46,47,48,49,53,57,58,59,63,64,65,68,71,72,75,106,107,111,112,113,117,118,121,122,125]
  ENDIF
  IF(MOVED_WALL_NODES_MOVING_MESH_FLAG) THEN
    NUMBER_OF_MOVED_WALL_NODES_MOVING_MESH=25
    ALLOCATE(MOVED_WALL_NODES_MOVING_MESH(NUMBER_OF_MOVED_WALL_NODES_MOVING_MESH))
    MOVED_WALL_NODES_MOVING_MESH=[1,2,5,6,9,14,15,16,23,28,30,32,36,37,42,76,77,80,81,82,89,94,96,97,102]
    BOUNDARY_CONDITIONS_MOVING_MESH(1)=0.0_CMISSRP
    BOUNDARY_CONDITIONS_MOVING_MESH(2)=0.0_CMISSRP
    BOUNDARY_CONDITIONS_MOVING_MESH(3)=0.0_CMISSRP
  ENDIF
  !Set material parameters
  MU_PARAM_STOKES=1.0_CMISSRP
  RHO_PARAM_STOKES=1.0_CMISSRP
  K_PARAM_MOVING_MESH=1.0_CMISSRP
  !Set interpolation parameters
  BASIS_XI_GAUSS_SPACE=3
  BASIS_XI_GAUSS_VELOCITY=3
  BASIS_XI_GAUSS_PRESSURE=3
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_MOVING_MESH_OUTPUT_TYPE=CMFE_SOLVER_NO_OUTPUT
  DYNAMIC_SOLVER_STOKES_OUTPUT_TYPE=CMFE_SOLVER_NO_OUTPUT
  LINEAR_SOLVER_STOKES_OUTPUT_TYPE=CMFE_SOLVER_NO_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_STOKES_OUTPUT=CMFE_EQUATIONS_NO_OUTPUT
  EQUATIONS_MOVING_MESH_OUTPUT=CMFE_EQUATIONS_NO_OUTPUT
  !Set time parameter
  DYNAMIC_SOLVER_STOKES_START_TIME=0.0_CMISSRP
  DYNAMIC_SOLVER_STOKES_STOP_TIME=5.0_CMISSRP 
  DYNAMIC_SOLVER_STOKES_TIME_INCREMENT=1.0_CMISSRP
  DYNAMIC_SOLVER_STOKES_THETA=1.0_CMISSRP
  !Set result output parameter
  DYNAMIC_SOLVER_STOKES_OUTPUT_FREQUENCY=1
  !Set bc input option
  DYNAMIC_SOLVER_STOKES_INPUT_OPTION=2
  !Set solver parameters
  LINEAR_SOLVER_MOVING_MESH_DIRECT_FLAG=.TRUE.
  LINEAR_SOLVER_STOKES_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-05_CMISSRP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-10_CMISSRP
  DIVERGENCE_TOLERANCE=1.0E20 !default: 1.0E5
  MAXIMUM_ITERATIONS=100000 !default: 100000
  RESTART_VALUE=3000 !default: 30
  LINESEARCH_ALPHA=1.0

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
    CALL cmfe_Basis_InterpolationXiSet(BasisSpace,[BASIS_XI_INTERPOLATION_SPACE,BASIS_XI_INTERPOLATION_SPACE],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisSpace,[BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE],Err)
  ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
    CALL cmfe_Basis_InterpolationXiSet(BasisSpace,[BASIS_XI_INTERPOLATION_SPACE,BASIS_XI_INTERPOLATION_SPACE, & 
      & BASIS_XI_INTERPOLATION_SPACE],Err)                         
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisSpace,[BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE],Err)
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
      CALL cmfe_Basis_InterpolationXiSet(BasisVelocity,[BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY],Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisVelocity,[BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY],Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisVelocity,[BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY, & 
        & BASIS_XI_INTERPOLATION_VELOCITY],Err)                         
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisVelocity,[BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY, & 
        & BASIS_XI_GAUSS_VELOCITY],Err)
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
      CALL cmfe_Basis_InterpolationXiSet(BasisPressure,[BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE],Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisPressure,[BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE],Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisPressure,[BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE, & 
        & BASIS_XI_INTERPOLATION_PRESSURE],Err)                         
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisPressure,[BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE, & 
        & BASIS_XI_GAUSS_PRESSURE],Err)
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

  !Create the equations set for ALE Stokes
  CALL cmfe_Field_Initialise(EquationsSetFieldStokes,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetStokes,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberStokes,Region,GeometricField, &
    & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMFE_EQUATIONS_SET_STOKES_EQUATION_TYPE,CMFE_EQUATIONS_SET_ALE_STOKES_SUBTYPE], &
    & EquationsSetFieldUserNumberStokes,EquationsSetFieldStokes,EquationsSetStokes,Err)
  !Set the equations set to be a ALE Stokes problem
!   CALL cmfe_EquationsSet_SpecificationSet(EquationsSetStokes,CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
!     & CMFE_EQUATIONS_SET_STOKES_EQUATION_TYPE,CMFE_EQUATIONS_SET_ALE_STOKES_SUBTYPE,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetStokes,Err)

  !Create the equations set for moving mesh
  CALL cmfe_Field_Initialise(EquationsSetFieldMovingMesh,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetMovingMesh,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberMovingMesh,Region,GeometricField, &
    & [CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS,CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE],EquationsSetFieldUserNumberMovingMesh,EquationsSetFieldMovingMesh, &
    & EquationsSetMovingMesh,Err)
  !Set the equations set to be a moving mesh problem
!   CALL cmfe_EquationsSet_SpecificationSet(EquationsSetMovingMesh,CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
!     & CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMFE_EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetMovingMesh,Err)


  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for ALE Stokes
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


  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for ALE Stokes
  CALL cmfe_Field_Initialise(MaterialsFieldStokes,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetStokes,MaterialsFieldUserNumberStokes, & 
    & MaterialsFieldStokes,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetStokes,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberStokesMu,MU_PARAM_STOKES,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldStokes,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberStokesRho,RHO_PARAM_STOKES,Err)
  !Create the equations set materials field variables for moving mesh
  CALL cmfe_Field_Initialise(MaterialsFieldMovingMesh,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetMovingMesh,MaterialsFieldUserNumberMovingMesh, & 
    & MaterialsFieldMovingMesh,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetMovingMesh,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberMovingMeshK,K_PARAM_MOVING_MESH,Err)

  !
  !================================================================================================================================
  !

  !INDEPENDENT FIELDS

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

  !
  !================================================================================================================================
  !

  !EQUATIONS


  !Create the equations set equations
  CALL cmfe_Equations_Initialise(EquationsStokes,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetStokes,EquationsStokes,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(EquationsStokes,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations lumping type
  CALL cmfe_Equations_LumpingTypeSet(EquationsStokes,CMFE_EQUATIONS_UNLUMPED_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(EquationsStokes,EQUATIONS_STOKES_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetStokes,Err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(EquationsMovingMesh,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetMovingMesh,EquationsMovingMesh,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(EquationsMovingMesh,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(EquationsMovingMesh,EQUATIONS_MOVING_MESH_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetMovingMesh,Err)

  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_FLUID_MECHANICS_CLASS,CMFE_PROBLEM_STOKES_EQUATION_TYPE, &
    & CMFE_PROBLEM_ALE_STOKES_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,DYNAMIC_SOLVER_STOKES_START_TIME,DYNAMIC_SOLVER_STOKES_STOP_TIME, & 
    & DYNAMIC_SOLVER_STOKES_TIME_INCREMENT,Err)
  !Set the output timing
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,DYNAMIC_SOLVER_STOKES_OUTPUT_FREQUENCY,Err)
  CALL cmfe_ControlLoop_TimeInputSet(ControlLoop,DYNAMIC_SOLVER_STOKES_INPUT_OPTION,Err)
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
  !Get the dynamic ALE solver
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverStokesUserNumber,DynamicSolverStokes,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(DynamicSolverStokes,DYNAMIC_SOLVER_STOKES_OUTPUT_TYPE,Err)
  !Set theta
  CALL cmfe_Solver_DynamicThetaSet(DynamicSolverStokes,DYNAMIC_SOLVER_STOKES_THETA,Err)
!   CALL cmfe_SolverDynamicALESet(DynamicSolverStokes,.TRUE.,Err)
  !Get the dynamic linear solver
  CALL cmfe_Solver_DynamicLinearSolverGet(DynamicSolverStokes,LinearSolverStokes,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(LinearSolverStokes,LINEAR_SOLVER_STOKES_OUTPUT_TYPE,Err)
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

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(LinearSolverMovingMesh,Err)
  CALL cmfe_Solver_Initialise(DynamicSolverStokes,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsMovingMesh,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsStokes,Err)

  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
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
  !Set fixed wall nodes
  IF(FIXED_WALL_NODES_STOKES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_FIXED_WALL_NODES_STOKES
      NODE_NUMBER=FIXED_WALL_NODES_STOKES(NODE_COUNTER)
      CONDITION=CMFE_BOUNDARY_CONDITION_FIXED_WALL
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
          VALUE=0.0_CMISSRP
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
          VALUE=0.0_CMISSRP
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
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsStokes,Err)
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
          VALUE=0.0_CMISSRP
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE, &
            & 1, &
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
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE, &
            & 1, &
            & CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsMovingMesh,Err)

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

  EXPORT_FIELD_IO=.TRUE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"ALEStokes","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"ALEStokes","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF

  !Finialise CMISS
!   CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

 
  STOP

END PROGRAM STOKESALEEXAMPLE
