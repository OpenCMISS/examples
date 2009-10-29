!> \file
!> $Id: ALEExample.f90 20 2009-10-28 20:22:52Z sebk $
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

!> \example FluidMechanics/Stokes/ALE/src/ALEExample.f90
!! Example program to solve a ALE Stokes equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FluidMechanics/Stokes/ALE/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FluidMechanics/Stokes/ALE/build-intel'>Linux GNU Build</a>
!!
!<

!> Main program

PROGRAM STOKESALEEXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OPENCMISS
  USE FLUID_MECHANICS_IO_ROUTINES
  USE MPI

#ifdef WIN32
  USE IFQWIN
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
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberStokes=6
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberMesh=7
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberStokes=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMesh=9
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberStokes=10
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberMesh=11
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberStokes=12
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberMesh=13
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=15

  INTEGER(CMISSIntg), PARAMETER :: SolverMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverStokesUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberStokesMu=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberStokesRho=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMeshK=1

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
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES
  INTEGER(CMISSIntg) :: NUMBER_OF_MOVED_WALL_NODES_STOKES
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES_MESH
  INTEGER(CMISSIntg) :: NUMBER_OF_MOVED_WALL_NODES_MESH

  INTEGER(CMISSIntg) :: EQUATIONS_STOKES_OUTPUT
  INTEGER(CMISSIntg) :: EQUATIONS_MESH_OUTPUT

  INTEGER(CMISSIntg),PARAMETER :: NO_DERIVATIVE_NUMBER=1

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES
  INTEGER, ALLOCATABLE, DIMENSION(:):: DOF_FIXED_INDICES_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: DOF_FIXED_CONDITION_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: BC_MOVED_WALL_NODES_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: DOF_MOVED_INDICES_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: DOF_MOVED_CONDITION_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_MESH
  INTEGER, ALLOCATABLE, DIMENSION(:):: DOF_FIXED_INDICES_MESH
  INTEGER, ALLOCATABLE, DIMENSION(:):: DOF_FIXED_CONDITION_MESH
  INTEGER, ALLOCATABLE, DIMENSION(:):: BC_MOVED_WALL_NODES_MESH
  INTEGER, ALLOCATABLE, DIMENSION(:):: DOF_MOVED_INDICES_MESH
  INTEGER, ALLOCATABLE, DIMENSION(:):: DOF_MOVED_CONDITION_MESH


  REAL(CMISSDP) :: INITIAL_FIELD_STOKES(3)
  REAL(CMISSDP) :: INITIAL_FIELD_MESH(3)
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_STOKES(3)
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_MESH(3)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA

  REAL(CMISSDP),ALLOCATABLE, DIMENSION(:):: DOF_FIXED_VALUES_STOKES
  REAL(CMISSDP),ALLOCATABLE, DIMENSION(:):: DOF_MOVED_VALUES_STOKES
  REAL(CMISSDP),ALLOCATABLE, DIMENSION(:):: DOF_FIXED_VALUES_MESH
  REAL(CMISSDP),ALLOCATABLE, DIMENSION(:):: DOF_MOVED_VALUES_MESH


  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: EXPORT_FIELD_FLUID_IO

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
  TYPE(CMISSFieldType) :: DependentFieldMesh
  TYPE(CMISSFieldType) :: MaterialsFieldStokes
  TYPE(CMISSFieldType) :: MaterialsFieldMesh
  TYPE(CMISSFieldType) :: IndependentFieldStokes
  TYPE(CMISSFieldType) :: IndependentFieldMesh
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsStokes
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsMesh
  !Equations sets
  TYPE(CMISSEquationsSetType) :: EquationsSetStokes
  TYPE(CMISSEquationsSetType) :: EquationsSetMesh
  !Equations
  TYPE(CMISSEquationsType) :: EquationsStokes
  TYPE(CMISSEquationsType) :: EquationsMesh
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: DynamicSolverStokes
  TYPE(CMISSSolverType) :: LinearSolverStokes
  TYPE(CMISSSolverType) :: LinearSolverMesh
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsStokes
  TYPE(CMISSSolverEquationsType) :: SolverEquationsMesh

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
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

  !Import cmHeart mesh information
  CALL FLUID_MECHANICS_IO_READ_CMHEART(CM,Err)  
  BASIS_NUMBER_SPACE=CM%ID_M
  BASIS_NUMBER_VELOCITY=CM%ID_V
  BASIS_NUMBER_PRESSURE=CM%ID_P
  MESH_COMPONENT_NUMBER_SPACE=CM%ID_M
  MESH_COMPONENT_NUMBER_VELOCITY=CM%ID_V
  MESH_COMPONENT_NUMBER_PRESSURE=CM%ID_P
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
  INITIAL_FIELD_MESH(1)=0.0_CMISSDP
  INITIAL_FIELD_MESH(2)=0.0_CMISSDP
  INITIAL_FIELD_MESH(3)=0.0_CMISSDP
  !Set boundary conditions
  FIXED_WALL_NODES_STOKES_FLAG=.TRUE.
  MOVED_WALL_NODES_STOKES_FLAG=.TRUE.
  INLET_WALL_NODES_STOKES_FLAG=.FALSE.
  FIXED_WALL_NODES_MESH_FLAG=.TRUE.
  MOVED_WALL_NODES_MESH_FLAG=.TRUE.
  IF(FIXED_WALL_NODES_STOKES_FLAG) THEN
    NUMBER_OF_FIXED_WALL_NODES_STOKES=16
    ALLOCATE(FIXED_WALL_NODES_STOKES(NUMBER_OF_FIXED_WALL_NODES_STOKES))
    FIXED_WALL_NODES_STOKES=(/46,47,48,53,57,64,65,68,72,106,107,111,117,118,122,125/)
  ENDIF
  IF(MOVED_WALL_NODES_STOKES_FLAG) THEN
    NUMBER_OF_MOVED_WALL_NODES_STOKES=73
    ALLOCATE(MOVED_WALL_NODES_STOKES(NUMBER_OF_MOVED_WALL_NODES_STOKES))
    MOVED_WALL_NODES_STOKES=(/3,4,7,10,11,12,13,17,20,24,29,31,33,34,35,39,41,44,50,51,52,54,60,66,67,70,74,78,79,83, &
      & 86,90,91,92,93,95,99,101,103,104,105,108,114,115,116,120,123,124,1,2,5,6,9,14,15,16,23,28,30,32,36, & 
      & 37,42,76,77,80,81,82,89,94,96,97,102/)
  ENDIF
  IF(INLET_WALL_NODES_STOKES_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_STOKES=25
    ALLOCATE(INLET_WALL_NODES_STOKES(NUMBER_OF_INLET_NODES_STOKES))
    INLET_WALL_NODES_STOKES=(/46,47,48,49,53,57,58,59,63,64,65,68,71,72,75,106,107,111,112,113,117,118,121,122,125/)
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_STOKES(1)=0.0_CMISSDP
    BOUNDARY_CONDITIONS_STOKES(2)=0.0_CMISSDP
    BOUNDARY_CONDITIONS_STOKES(3)=0.0_CMISSDP
  ENDIF
  IF(FIXED_WALL_NODES_MESH_FLAG) THEN
    NUMBER_OF_FIXED_WALL_NODES_MESH=25
    ALLOCATE(FIXED_WALL_NODES_MESH(NUMBER_OF_FIXED_NODES_MESH))
    FIXED_WALL_NODES_MESH=(/46,47,48,49,53,57,58,59,63,64,65,68,71,72,75,106,107,111,112,113,117,118,121,122,125/)
  ENDIF
  IF(MOVED_WALL_NODES_MESH_FLAG) THEN
    NUMBER_OF_MOVED_WALL_NODES_MESH=25
    ALLOCATE(MOVED_WALL_NODES_MESH(NUMBER_OF_MOVED_WALL_NODES_MESH))
    MOVED_WALL_NODES_MESH=(/1,2,5,6,9,14,15,16,23,28,30,32,36,37,42,76,77,80,81,82,89,94,96,97,102/)
    BOUNDARY_CONDITIONS_MESH(1)=0.0_CMISSDP
    BOUNDARY_CONDITIONS_MESH(2)=0.0_CMISSDP
    BOUNDARY_CONDITIONS_MESH(3)=0.0_CMISSDP
  ENDIF
  !Set material parameters
  MU_PARAM_STOKES=1.0_CMISSDP
  RHO_PARAM_STOKES=1.0_CMISSDP
  K_PARAM_MESH=1.0_CMISSDP
  !Set interpolation parameters
  BASIS_XI_GAUSS_SPACE=3
  BASIS_XI_GAUSS_VELOCITY=3
  BASIS_XI_GAUSS_PRESSURE=3
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_MESH_OUTPUT_TYPE=CMISSSolverNoOutput
  DYNAMIC_SOLVER_STOKES_OUTPUT_TYPE=CMISSSolverNoOutput
  LINEAR_SOLVER_STOKES_OUTPUT_TYPE=CMISSSolverNoOutput
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_STOKES_OUTPUT=CMISSEquationsNoOutput
  EQUATIONS_MESH_OUTPUT=CMISSEquationsNoOutput
  !Set time parameter
  DYNAMIC_SOLVER_STOKES_START_TIME=0.0_CMISSDP
  DYNAMIC_SOLVER_STOKES_STOP_TIME=4.0_CMISSDP 
  DYNAMIC_SOLVER_STOKES_TIME_INCREMENT=1.0_CMISSDP
  DYNAMIC_SOLVER_STOKES_THETA=2.0_CMISSDP/3.0_CMISSDP
  !Set solver parameters
  LINEAR_SOLVER_MESH_DIRECT_FLAG=.FALSE.
  LINEAR_SOLVER_STOKES_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-14_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-14_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E20 !default: 1.0E5
  MAXIMUM_ITERATIONS=100000 !default: 100000
  RESTART_VALUE=300 !default: 30
  LINESEARCH_ALPHA=1.0


  !
  !================================================================================================================================
  !

  !INITIALISE OPENCMISS

  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !
  !================================================================================================================================
  !

  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL CMISSNumberOfDimensionsSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
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
    CALL CMISSBasisQuadratureNumberOfGaussXiGet(BasisSpace,(BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE),Err)
  ELSE IF(NUMBER_OF_DIMENSIONS==3)
    CALL CMISSBasisInterpolationXiSet(BasisSpace,(/BASIS_XI_INTERPOLATION_SPACE,BASIS_XI_INTERPOLATION_SPACE, & 
      & BASIS_XI_INTERPOLATION_SPACE/),Err)                         
    CALL CMISSBasisQuadratureNumberOfGaussXiGet(BasisSpace,(BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE),Err)
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
      CALL CMISSBasisQuadratureNumberOfGaussXiGet(BasisVelocity,(BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3)
      CALL CMISSBasisInterpolationXiSet(BasisVelocity,(/BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY, & 
        & BASIS_XI_INTERPOLATION_VELOCITY/),Err)                         
      CALL CMISSBasisQuadratureNumberOfGaussXiGet(BasisVelocity,(BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY, & 
        & BASIS_XI_GAUSS_VELOCITY),Err)
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
      CALL CMISSBasisQuadratureNumberOfGaussXiGet(BasisPressure,(BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3)
      CALL CMISSBasisInterpolationXiSet(BasisPressure,(/BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE, & 
        & BASIS_XI_INTERPOLATION_PRESSURE/),Err)                         
      CALL CMISSBasisQuadratureNumberOfGaussXiGet(BasisPressure,(BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE, & 
        & BASIS_XI_GAUSS_PRESSURE),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(BasisPressure,Err)
  ENDIF

  !
  !================================================================================================================================
  !

  !MESH

  !Start the creation of mesh nodes
  CALL CMISSNodesCreateStart(Region,TOTAL_NUMBER_OF_NODES,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL CMISSMeshNumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL CMISSMeshNumberOfComponentsSet(Mesh,MESH_NUMBER_OF_COMPONENTS,Err)
  !Specify spatial mesh component
  CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_SPACE,BasisSpace,MeshElementsSpace,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL CMISSMeshElementsNodesSet(MeshElementsSpace,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_SPACE),Err)
  ENDDO
  CALL CMISSMeshElementsCreateFinish(MeshElementsSpace,Err)
  !Specify velocity mesh component
  IF(MESH_COMPONENT_NUMBER_VELOCITY==MESH_COMPONENT_NUMBER_SPACE) THEN
    MeshElementsVelocity=MeshElementsSpace
  ELSE
    CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_VELOCITY,BasisVelocity,MeshElementsVelocity,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL CMISSMeshElementsNodesSet(MeshElementsVelocity,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_VELOCITY),Err)
    ENDDO
    CALL CMISSMeshElementsCreateFinish(MeshElementsVelocity,Err)
  ENDIF
  !Specify pressure mesh component
  IF(MESH_COMPONENT_NUMBER_PRESSURE==MESH_COMPONENT_NUMBER_SPACE) THEN
    MeshElementsPressure=MeshElementsSpace
  ELSE IF(MESH_COMPONENT_NUMBER_PRESSURE==MESH_COMPONENT_NUMBER_VELOCITY) THEN
    MeshElementsPressure=MeshElementsVelocity
  ELSE
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

  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the field type
  CALL CMISSFieldTypeSetObj(GeometricField,CMISSFieldGeometricType,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL CMISSFieldScalingTypeSetObj(GeometricField,CMISSFieldNoScaling,Err)
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
      CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,NO_DERIVATIVE_NUMBER, &
        NODE_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSFieldParameterSetUpdateStart(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,DomainUserNumber,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for ALE Stokes
  CALL CMISSEquationsSetTypeInitialise(EquationsSetStokes,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberStokes,Region,GeometricField,EquationsSetStokes,Err)
  !Set the equations set to be a ALE Stokes problem
  CALL CMISSEquationsSetSpecificationSet(EquationsSetStokes,CMISSEquationsSetFluidMechanicsClass, &
    & CMISSEquationsSetStokesEquationType,CMISSEquationsSetALEStokesSubtype,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSetStokes,Err)
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
  !Create the equations set for moving mesh
  CALL CMISSEquationsSetTypeInitialise(EquationsSetMesh,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberMesh,Region,GeometricField,EquationsSetMesh,Err)
  !Set the equations set to be a moving mesh problem
  CALL CMISSEquationsSetSpecificationSet(EquationsSetMesh,CMISSEquationsSetClassicalFieldClass, &
    & CMISSEquationsSetLaplaceEquationType,CMISSEquationsSetMeshLaplaceSubtype,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSetMesh,Err)
  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(EquationsMesh,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetMesh,EquationsMesh,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(EquationsMesh,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(EquationsMesh,EQUATIONS_MESH_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetMesh,Err)


  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for ALE Stokes
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
  !Create the equations set dependent field variables for moving mesh
  CALL CMISSFieldTypeInitialise(DependentFieldMesh,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetMesh,DependentFieldUserNumberMesh, & 
    & DependentFieldMesh,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldMesh,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldMesh,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
  ENDDO
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetMesh,Err)
  !Initialise dependent field
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSFieldComponentValuesInitialise(DependentFieldMesh,CMISSFieldUVariableType,CMISSFieldValuesSetType,COMPONENT_NUMBER, & 
      & INITIAL_FIELD_MESH(COMPONENT_NUMBER),Err)
  ENDDO


  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for ALE Stokes
  CALL CMISSFieldTypeInitialise(MaterialsFieldStokes,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetStokes,MaterialsFieldUserNumberStokes, & 
    & MaterialsFieldStokes,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetStokes,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & MaterialsFieldUserNumberStokesMu,MU_PARAM_STOKES,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & MaterialsFieldUserNumberStokesRho,RHO_PARAM_STOKES,Err)
  !Create the equations set materials field variables for moving mesh
  CALL CMISSFieldTypeInitialise(MaterialsFieldMesh,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetMesh,MaterialsFieldUserNumberMesh, & 
    & MaterialsFieldMesh,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetMesh,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & MaterialsFieldUserNumberMeshK,K_PARAM_MESH,Err)

  !
  !================================================================================================================================
  !

  !INDEPENDENT FIELDS

  !Create the equations set independent field variables for ALE Stokes
  CALL CMISSFieldTypeInitialise(IndependentFieldStokes,Err)
  CALL CMISSEquationsSetIndependentCreateStart(EquationsSetStokes,IndependentFieldUserNumberStokes, & 
    & IndependentFieldStokes,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSFieldComponentMeshComponentSet(InDependentFieldStokes,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDDO
  !Finish the equations set independent field variables
  CALL CMISSEquationsSetIndependentCreateFinish(EquationsSetStokes,Err)
  !Create the equations set independent field variables for moving mesh
  CALL CMISSFieldTypeInitialise(IndependentFieldMesh,Err)
  CALL CMISSEquationsSetIndependentCreateStart(EquationsSetMesh,IndependentFieldUserNumberMesh, & 
    & IndependentFieldMesh,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSFieldComponentMeshComponentSet(InDependentFieldMesh,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDDO
  !Finish the equations set independent field variables
  CALL CMISSEquationsSetIndependentCreateFinish(EquationsSetMesh,Err)


  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  !Start the creation of the equations set boundary conditions for Stokes
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsStokes,Err)
  CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSetStokes,BoundaryConditionsStokes,Err)
  !Set fixed wall nodes
  IF(FIXED_WALL_NODES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_FIXED_WALL_NODES
      NODE_NUMBER=FIXED_WALL_NODES(NODE_COUNTER)
      CONDITION=BOUNDARY_CONDITION_FIXED_WALL
      DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=0.0_DPCMISS
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsStokes,CMISSFieldUVariableType,NO_DERIVATIVE_NUMBER,NODE_NUMBER, & 
          & COMPONENT_NUMBER,CONDITION,VALUE,Err)
      ENDDO
    ENDDO
  ENDIF
  !Set moved wall nodes
  IF(MOVED_WALL_NODES_STOKES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_MOVED_WALL_NODES_STOKES
      NODE_NUMBER=MOVED_WALL_NODES_STOKES(NODE_COUNTER)
      CONDITION=BOUNDARY_CONDITION_MOVED_WALL
      DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=0.0_DPCMISS
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsStokes,CMISSFieldUVariableType,NO_DERIVATIVE_NUMBER,NODE_NUMBER, & 
          & COMPONENT_NUMBER,CONDITION,VALUE,Err)
      ENDDO
    ENDDO
  ENDIF
  !Set velocity boundary conditions
  IF(INLET_WALL_NODES_STOKES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_NODES_STOKES
      NODE_NUMBER=INLET_NODES_STOKES(NODE_COUNTER)
      CONDITION=BOUNDARY_CONDITION_INLET_WALL
      DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=BOUNDARY_CONDITIONS_STOKES(COMPONENT_NUMBER)
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsStokes,CMISSFieldUVariableType,NO_DERIVATIVE_NUMBER,NODE_NUMBER, & 
          & COMPONENT_NUMBER,CONDITION,VALUE,Err)
      ENDDO
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSetStokes,Err)
  !Start the creation of the equations set boundary conditions for moving mesh
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsMesh,Err)
  CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSetMesh,BoundaryConditionsMesh,Err)
  !Set fixed wall nodes
  IF(FIXED_WALL_NODES_MESH_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_FIXED_WALL_NODES_MESH
      NODE_NUMBER=FIXED_WALL_NODES_MESH(NODE_COUNTER)
      CONDITION=BOUNDARY_CONDITION_FIXED_WALL
      DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=0.0_DPCMISS
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsMesh,CMISSFieldUVariableType,NO_DERIVATIVE_NUMBER,NODE_NUMBER, & 
          & COMPONENT_NUMBER,CONDITION,VALUE,Err)
      ENDDO
    ENDDO
  ENDIF
  !Set moved wall nodes
  IF(MOVED_WALL_NODES_MESH_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_MOVED_WALL_NODES_MESH
      NODE_NUMBER=MOVED_WALL_NODES_MESH(NODE_COUNTER)
      CONDITION=BOUNDARY_CONDITION_MOVED_WALL
      DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=BOUNDARY_CONDITIONS_MESH(COMPONENT_NUMBER)
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsMesh,CMISSFieldUVariableType,NO_DERIVATIVE_NUMBER,NODE_NUMBER, & 
          & COMPONENT_NUMBER,CONDITION,VALUE,Err)
      ENDDO
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSetMesh,Err)

  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a ALE Stokes problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemFluidMechanicsClass,CMISSProblemStokesEquationType, &
    & CMISSProblemALEStokesSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,CMISSControlLoop,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(CMISSControlLoop,DYNAMIC_SOLVER_STOKES_START_TIME,DYNAMIC_SOLVER_STOKES_STOP_TIME, & 
    & DYNAMIC_SOLVER_STOKES_TIME_INCREMENT,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  !Get the moving mesh solver
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverMeshUserNumber,LinearSolverMesh,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(LinearSolverMesh,LINEAR_SOLVER_MESH_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_MESH_DIRECT_FLAG) THEN
    CALL CMISSSolverLinearTypeSet(LinearSolverMesh,CMISSSolverLinearDirectSolveType,Err)
    CALL CMISSSolverLinearDirectTypeSet(LinearSolverMesh,CMISSSolverDirectMumps,Err)
  ELSE
    CALL CMISSSolverLinearTypeSet(LinearSolverMesh,CMISSSolverLinearIterativeSolveType,Err)
    CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolverMesh,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolverLinearIterativeDivergenceToleranceSet(LinearSolverMesh,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeRelativeToleranceSet(LinearSolverMesh,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(LinearSolverMesh,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeGMRESRestartSet(LinearSolverMesh,RESTART_VALUE,Err)
  ENDIF
  !Get the dynamic ALE solver
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverStokesUserNumber,DynamicSolverStokes,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(DynamicSolverStokes,DYNAMIC_SOLVER_STOKES_OUTPUT_TYPE,Err)
  !Set theta
  CALL CMISSSolverDynamicThetaSet(DynamicSolverStokes,DYNAMIC_SOLVER_STOKES_THETA,Err)
  CALL CMISSSolverDynamicALESet(DynamicSolverStokes,.TRUE.,Err)
  !Get the dynamic linear solver
  CALL CMISSSolverDynamicLinearSolverGet(DynamicSolverStokes,LinearSolverStokes,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(LinearSolverStokes,LINEAR_SOLVER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_MESH_DIRECT_FLAG) THEN
    CALL CMISSSolverLinearTypeSet(LinearSolverStokes,CMISSSolverLinearDirectSolveType,Err)
    CALL CMISSSolverLinearDirectTypeSet(LinearSolverStokes,CMISSSolverDirectMumps,Err)
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
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the linear solver equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverMeshUserNumber,LinearSolverMesh,Err)
  CALL CMISSSolverSolverEquationsGet(LinearSolverMesh,SolverEquationsMesh,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsMesh,CMISSSolverEquationsSparseMatrices,Err)
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsMesh,EquationsSetMesh,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
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

  !RUN SOLVERS

  !Turn of PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  CALL CMISSProblemSolve(Problem,Err)
  WRITE(*,'(A)') "Solving problem..."
  CALL PROBLEM_SOLVE(PROBLEM,ERR,ERROR,*999) 
  WRITE(*,'(A)') "Problem solved!"

  !
  !================================================================================================================================
  !

  !OUTPUT

  EXPORT_FIELD_IO=.FALSE.
  IF(EXPORT_FIELD_IO) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"ALEStokes","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"ALEStokes","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
  ENDIF
  

  FILE="cmgui"

  EXPORT_FIELD_FLUID_IO=.TRUE.
  IF(EXPORT_FIELD_FLUID_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL FLUID_MECHANICS_IO_WRITE_CMGUI(Region,"cmgui",Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF


  !Finialise CMISS
  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM STOKESALEEXAMPLE
