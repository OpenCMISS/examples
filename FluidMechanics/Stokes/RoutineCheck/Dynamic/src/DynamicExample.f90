!> \file
!> \author Sebastian Krittian
!> \brief This is an example program to solve a dynamic Stokes equation using OpenCMISS calls.
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

!> \example FluidMechanics/Stokes/RoutineCheck/Dynamic/src/DynamicExample.f90
!! Example program to solve a dynamic Stokes equation using OpenCMISS calls.
!! \htmlinclude FluidMechanics/Stokes/RoutineCheck/Dynamic/history.html
!!
!<

!> Main program

PROGRAM STOKESDYNAMICEXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OPENCMISS
  USE FLUID_MECHANICS_IO_ROUTINES
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
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberStokes=6
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberStokes=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberStokes=8
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=9

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: SolverStokesUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberStokesMu=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberStokesRho=2

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
!   INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES_STOKES
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_WALL_NODES_STOKES

  INTEGER(CMISSIntg) :: EQUATIONS_STOKES_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_COUNTER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_STOKES_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_STOKES_OUTPUT_TYPE

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_STOKES


  REAL(CMISSDP) :: INITIAL_FIELD_STOKES(3)
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_STOKES(3)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE
  REAL(CMISSDP) :: MU_PARAM_STOKES
  REAL(CMISSDP) :: RHO_PARAM_STOKES

  REAL(CMISSDP) :: DYNAMIC_SOLVER_STOKES_START_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_STOKES_STOP_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_STOKES_THETA
  REAL(CMISSDP) :: DYNAMIC_SOLVER_STOKES_TIME_INCREMENT

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_STOKES_DIRECT_FLAG
  LOGICAL :: FIXED_WALL_NODES_STOKES_FLAG
  LOGICAL :: INLET_WALL_NODES_STOKES_FLAG

  !CMISS variables

  !Regions
  TYPE(CMISSRegionType) :: Region
  TYPE(CMISSRegionType) :: WorldRegion
  !Coordinate systems
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  !Basis
  TYPE(CMISSBasisType) :: BasisSpace
  TYPE(CMISSBasisType) :: BasisVelocity
  TYPE(CMISSBasisType) :: BasisPressure
  !Nodes
  TYPE(CMISSNodesType) :: Nodes
  !Elements
  TYPE(CMISSMeshElementsType) :: MeshElementsSpace
  TYPE(CMISSMeshElementsType) :: MeshElementsVelocity
  TYPE(CMISSMeshElementsType) :: MeshElementsPressure
  !Meshes
  TYPE(CMISSMeshType) :: Mesh
  !Decompositions
  TYPE(CMISSDecompositionType) :: Decomposition
  !Fields
  TYPE(CMISSFieldsType) :: Fields
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: DependentFieldStokes
  TYPE(CMISSFieldType) :: MaterialsFieldStokes
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsStokes
  !Equations sets
  TYPE(CMISSEquationsSetType) :: EquationsSetStokes
  !Equations
  TYPE(CMISSEquationsType) :: EquationsStokes
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: DynamicSolverStokes
  TYPE(CMISSSolverType) :: LinearSolverStokes
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsStokes

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

  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

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
  INITIAL_FIELD_STOKES(1)=0.0_CMISSDP
  INITIAL_FIELD_STOKES(2)=0.0_CMISSDP
  INITIAL_FIELD_STOKES(3)=0.0_CMISSDP
  !Set boundary conditions
  FIXED_WALL_NODES_STOKES_FLAG=.TRUE.
  INLET_WALL_NODES_STOKES_FLAG=.TRUE.
  IF(FIXED_WALL_NODES_STOKES_FLAG) THEN
    NUMBER_OF_FIXED_WALL_NODES_STOKES=80
    ALLOCATE(FIXED_WALL_NODES_STOKES(NUMBER_OF_FIXED_WALL_NODES_STOKES))
    FIXED_WALL_NODES_STOKES=(/1,2,3,4,5,7,9,10,11,12,13,14,17,20,24,28,29,30,31,32,33,34,35,37,39, & 
    & 41,44,46,47,48,50,51,52,53,54,57,60,64,65,66,67,68,70,72,74,76,77,78,79,80,83,86, & 
    & 89,90,91,92,93,94,95,97,99,101,102,103,104,105,106,107,108,111,114,115,116,117,118, & 
    & 120,122,123,124,125/)
  ENDIF
  IF(INLET_WALL_NODES_STOKES_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_STOKES=9
    ALLOCATE(INLET_WALL_NODES_STOKES(NUMBER_OF_INLET_WALL_NODES_STOKES))
    INLET_WALL_NODES_STOKES=(/6,15,16,23,36,42,81,82,96/)
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_STOKES(1)=0.0_CMISSDP
    BOUNDARY_CONDITIONS_STOKES(2)=1.0_CMISSDP
    BOUNDARY_CONDITIONS_STOKES(3)=0.0_CMISSDP
  ENDIF
  !Set material parameters
  MU_PARAM_STOKES=1.0_CMISSDP
  RHO_PARAM_STOKES=1.0_CMISSDP
  !Set interpolation parameters
  BASIS_XI_GAUSS_SPACE=3
  BASIS_XI_GAUSS_VELOCITY=3
  BASIS_XI_GAUSS_PRESSURE=3
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  DYNAMIC_SOLVER_STOKES_OUTPUT_TYPE=CMISSSolverNoOutput
  LINEAR_SOLVER_STOKES_OUTPUT_TYPE=CMISSSolverNoOutput
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_STOKES_OUTPUT=CMISSEquationsNoOutput
  !Set time parameter
  DYNAMIC_SOLVER_STOKES_START_TIME=0.0_CMISSDP
  DYNAMIC_SOLVER_STOKES_STOP_TIME=5.0_CMISSDP 
  DYNAMIC_SOLVER_STOKES_TIME_INCREMENT=1.0_CMISSDP
  DYNAMIC_SOLVER_STOKES_THETA=1.0_CMISSDP
  !Set result output parameter
  DYNAMIC_SOLVER_STOKES_OUTPUT_FREQUENCY=1
  !Set solver parameters
  LINEAR_SOLVER_STOKES_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E20 !default: 1.0E5
  MAXIMUM_ITERATIONS=100000 !default: 100000
  RESTART_VALUE=3000 !default: 30
  LINESEARCH_ALPHA=1.0

  !
  !================================================================================================================================
  !

  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !REGION

  !Start the creation of a new region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system as defined above
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !BASES

  !Start the creation of new bases
  MESH_NUMBER_OF_COMPONENTS=1
  CALL CMISSBasisTypeInitialise(BasisSpace,Err)
  CALL CMISSBasisCreateStart(BASIS_NUMBER_SPACE,BasisSpace,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasisTypeSet(BasisSpace,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL CMISSBasisNumberOfXiSet(BasisSpace,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  IF(NUMBER_OF_DIMENSIONS==2) THEN
    CALL CMISSBasisInterpolationXiSet(BasisSpace,(/BASIS_XI_INTERPOLATION_SPACE,BASIS_XI_INTERPOLATION_SPACE/),Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisSpace,(/BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE/),Err)
  ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
    CALL CMISSBasisInterpolationXiSet(BasisSpace,(/BASIS_XI_INTERPOLATION_SPACE,BASIS_XI_INTERPOLATION_SPACE, & 
      & BASIS_XI_INTERPOLATION_SPACE/),Err)                         
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisSpace,(/BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE/),Err)
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(BasisSpace,Err)
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisVelocity=BasisSpace
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new velocity basis
    CALL CMISSBasisTypeInitialise(BasisVelocity,Err)
    !Start the creation of a basis
    CALL CMISSBasisCreateStart(BASIS_NUMBER_VELOCITY,BasisVelocity,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasisTypeSet(BasisVelocity,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL CMISSBasisNumberOfXiSet(BasisVelocity,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL CMISSBasisInterpolationXiSet(BasisVelocity,(/BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY/),Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisVelocity,(/BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY/),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL CMISSBasisInterpolationXiSet(BasisVelocity,(/BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY, & 
        & BASIS_XI_INTERPOLATION_VELOCITY/),Err)                         
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisVelocity,(/BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY, & 
        & BASIS_XI_GAUSS_VELOCITY/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(BasisVelocity,Err)
  ENDIF
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisPressure=BasisSpace
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_VELOCITY) THEN
    BasisPressure=BasisVelocity
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new pressure basis
    CALL CMISSBasisTypeInitialise(BasisPressure,Err)
    !Start the creation of a basis
    CALL CMISSBasisCreateStart(BASIS_NUMBER_PRESSURE,BasisPressure,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasisTypeSet(BasisPressure,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL CMISSBasisNumberOfXiSet(BasisPressure,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL CMISSBasisInterpolationXiSet(BasisPressure,(/BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE/),Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisPressure,(/BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE/),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL CMISSBasisInterpolationXiSet(BasisPressure,(/BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE, & 
        & BASIS_XI_INTERPOLATION_PRESSURE/),Err)                         
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisPressure,(/BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE, & 
        & BASIS_XI_GAUSS_PRESSURE/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(BasisPressure,Err)
  ENDIF

  !
  !================================================================================================================================
  !

  !MESH

  !Start the creation of mesh nodes
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSNodesCreateStart(Region,TOTAL_NUMBER_OF_NODES,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL CMISSMeshNumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL CMISSMeshNumberOfComponentsSet(Mesh,MESH_NUMBER_OF_COMPONENTS,Err)
  !Specify spatial mesh component
  CALL CMISSMeshElementsTypeInitialise(MeshElementsSpace,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsVelocity,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsPressure,Err)
  MESH_COMPONENT_NUMBER_SPACE=1
  MESH_COMPONENT_NUMBER_VELOCITY=1
  MESH_COMPONENT_NUMBER_PRESSURE=1
  CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_SPACE,BasisSpace,MeshElementsSpace,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL CMISSMeshElementsNodesSet(MeshElementsSpace,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_SPACE),Err)
  ENDDO
  CALL CMISSMeshElementsCreateFinish(MeshElementsSpace,Err)
  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsVelocity=MeshElementsSpace
  ELSE
    MESH_COMPONENT_NUMBER_VELOCITY=MESH_COMPONENT_NUMBER_SPACE+1
    CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_VELOCITY,BasisVelocity,MeshElementsVelocity,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL CMISSMeshElementsNodesSet(MeshElementsVelocity,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_VELOCITY),Err)
    ENDDO
    CALL CMISSMeshElementsCreateFinish(MeshElementsVelocity,Err)
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
    CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_PRESSURE,BasisPressure,MeshElementsPressure,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL CMISSMeshElementsNodesSet(MeshElementsPressure,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_PRESSURE),Err)
    ENDDO
    CALL CMISSMeshElementsCreateFinish(MeshElementsPressure,Err)
  ENDIF
  !Finish the creation of the mesh
  CALL CMISSMeshCreateFinish(Mesh,Err)

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the field type
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL CMISSFieldScalingTypeSet(GeometricField,CMISSFieldNoScaling,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDDO
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_SPACE
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
          & 1,CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
      ENDIF
    ENDDO
  ENDDO
  CALL CMISSFieldParameterSetUpdateStart(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)



  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for dynamic Stokes
  CALL CMISSEquationsSetTypeInitialise(EquationsSetStokes,Err)
    CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberStokes,Region,GeometricField,CMISSEquationsSetFluidMechanicsClass, &
    & CMISSEquationsSetStokesEquationType,CMISSEquationsSetTransientStokesSubtype,EquationsSetFieldUserNumber,EquationsSetField, &
    & EquationsSetStokes,Err)
  !Set the equations set to be a dynamic Stokes problem
  
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSetStokes,Err)


  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for dynamic Stokes
  CALL CMISSFieldTypeInitialise(DependentFieldStokes,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetStokes,DependentFieldUserNumberStokes, & 
    & DependentFieldStokes,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldStokes,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldStokes,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
  ENDDO
  COMPONENT_NUMBER=NUMBER_OF_DIMENSIONS+1
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldStokes,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldStokes,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_PRESSURE,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetStokes,Err)

  !Initialise dependent field
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSFieldComponentValuesInitialise(DependentFieldStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
      & COMPONENT_NUMBER,INITIAL_FIELD_STOKES(COMPONENT_NUMBER),Err)
  ENDDO


  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for dynamic Stokes
  CALL CMISSFieldTypeInitialise(MaterialsFieldStokes,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetStokes,MaterialsFieldUserNumberStokes, & 
    & MaterialsFieldStokes,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetStokes,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & MaterialsFieldUserNumberStokesMu,MU_PARAM_STOKES,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & MaterialsFieldUserNumberStokesRho,RHO_PARAM_STOKES,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS


  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(EquationsStokes,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetStokes,EquationsStokes,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(EquationsStokes,CMISSEquationsSparseMatrices,Err)
  !Set the equations lumping type
  CALL CMISSEquationsLumpingTypeSet(EquationsStokes,CMISSEquationsUnlumpedMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(EquationsStokes,EQUATIONS_STOKES_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetStokes,Err)


  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a dynamic Stokes problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemFluidMechanicsClass,CMISSProblemStokesEquationType, &
    & CMISSProblemTransientStokesSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,DYNAMIC_SOLVER_STOKES_START_TIME,DYNAMIC_SOLVER_STOKES_STOP_TIME, & 
    & DYNAMIC_SOLVER_STOKES_TIME_INCREMENT,Err)
  !Set the output timing
  CALL CMISSControlLoopTimeOutputSet(ControlLoop,DYNAMIC_SOLVER_STOKES_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(DynamicSolverStokes,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverStokes,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  !Get the dynamic dymamic solver
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverStokesUserNumber,DynamicSolverStokes,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(DynamicSolverStokes,DYNAMIC_SOLVER_STOKES_OUTPUT_TYPE,Err)
  !Set theta
  CALL CMISSSolverDynamicThetaSet(DynamicSolverStokes,DYNAMIC_SOLVER_STOKES_THETA,Err)
!   CALL CMISSSolverDynamicDynamicSet(DynamicSolverStokes,.TRUE.,Err)
  !Get the dynamic linear solver
  CALL CMISSSolverDynamicLinearSolverGet(DynamicSolverStokes,LinearSolverStokes,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(LinearSolverStokes,LINEAR_SOLVER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_STOKES_DIRECT_FLAG) THEN
    CALL CMISSSolverLinearTypeSet(LinearSolverStokes,CMISSSolverLinearDirectSolveType,Err)
    CALL CMISSSolverLibraryTypeSet(LinearSolverStokes,CMISSSolverMUMPSLibrary,Err)
  ELSE
    CALL CMISSSolverLinearTypeSet(LinearSolverStokes,CMISSSolverLinearIterativeSolveType,Err)
    CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolverStokes,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolverLinearIterativeDivergenceToleranceSet(LinearSolverStokes,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeRelativeToleranceSet(LinearSolverStokes,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(LinearSolverStokes,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeGMRESRestartSet(LinearSolverStokes,RESTART_VALUE,Err)
  ENDIF
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL CMISSSolverTypeInitialise(DynamicSolverStokes,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsStokes,Err)

  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the dynamic solver equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverStokesUserNumber,DynamicSolverStokes,Err)
  CALL CMISSSolverSolverEquationsGet(DynamicSolverStokes,SolverEquationsStokes,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsStokes,CMISSSolverEquationsSparseMatrices,Err)
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsStokes,EquationsSetStokes,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  !Start the creation of the equations set boundary conditions for Stokes
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsStokes,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsStokes,BoundaryConditionsStokes,Err)
  !Set fixed wall nodes
  IF(FIXED_WALL_NODES_STOKES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_FIXED_WALL_NODES_STOKES
      NODE_NUMBER=FIXED_WALL_NODES_STOKES(NODE_COUNTER)
      CONDITION=CMISSBoundaryConditionFixedWall
      CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
          VALUE=0.0_CMISSDP
          CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsStokes,DependentFieldStokes,CMISSFieldUVariableType,1, &
            & CMISSNoGlobalDerivative, &
            & NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !Set velocity boundary conditions
  IF(INLET_WALL_NODES_STOKES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_STOKES
      NODE_NUMBER=INLET_WALL_NODES_STOKES(NODE_COUNTER)
      CONDITION=CMISSBoundaryConditionInletWall
      CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
          VALUE=BOUNDARY_CONDITIONS_STOKES(COMPONENT_NUMBER)
          CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsStokes,DependentFieldStokes,CMISSFieldUVariableType,1, &
            & CMISSNoGlobalDerivative, &
            & NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsStokes,Err)

  !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn of PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL CMISSProblemSolve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"

  !
  !================================================================================================================================
  !

  !OUTPUT

  EXPORT_FIELD_IO=.TRUE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"DynamicStokes","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"DynamicStokes","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF
  
  !Finialise CMISS
!   CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM STOKESDYNAMICEXAMPLE
