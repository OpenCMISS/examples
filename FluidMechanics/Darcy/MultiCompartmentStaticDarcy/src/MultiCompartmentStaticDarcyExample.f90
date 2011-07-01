!> \file
!> \author Christian Michler
!> \brief This is an example program to solve a static Darcy equation using openCMISS calls.
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

!> \example FluidMechanics/Darcy/Static/src/StaticExample.f90
!! Example program to solve a static Darcy equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FluidMechanics/Darcy/Static/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FluidMechanics/Darcy/Static/build-intel'>Linux GNU Build</a>
!!
!<

! ! 
! !  This example considers a static Darcy problem.
! ! 

!> Main program

PROGRAM DARCYSTATICEXAMPLE

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

  !Test program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberDarcy=6
  INTEGER(CMISSIntg) :: MaterialsFieldUserNumberDarcy
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=9

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverDarcyUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPorosity=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPermOverVis=2


  REAL(CMISSDP), PARAMETER :: PI=3.141592653589793238462643383279502884197_CMISSDP


  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS
  
  INTEGER(CMISSIntg) :: BASIS_TYPE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_PRESSURE
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_GEOMETRY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_GEOMETRY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_GEOMETRY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_PRESSURE
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE
!   INTEGER(CMISSIntg) :: MPI_IERROR

  INTEGER(CMISSIntg) :: EQUATIONS_DARCY_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: LINEAR_SOLVER_DARCY_OUTPUT_TYPE

  REAL(CMISSDP) :: COORD_X, COORD_Y, COORD_Z
  REAL(CMISSDP) :: DOMAIN_X1, DOMAIN_X2, DOMAIN_Y1, DOMAIN_Y2, DOMAIN_Z1, DOMAIN_Z2
  REAL(CMISSDP) :: GEOMETRY_TOLERANCE

  REAL(CMISSDP) :: INITIAL_FIELD_DARCY(3)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: VALUE
  REAL(CMISSDP) :: POROSITY_PARAM_DARCY, PERM_OVER_VIS_PARAM_DARCY

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_DARCY_DIRECT_FLAG

  !CMISS variables

  !Regions
  TYPE(CMISSRegionType) :: Region
  TYPE(CMISSRegionType) :: WorldRegion
  !Coordinate systems
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  !Basis
  TYPE(CMISSBasisType) :: BasisGeometry
  TYPE(CMISSBasisType) :: BasisVelocity
  TYPE(CMISSBasisType) :: BasisPressure
  !Nodes
  TYPE(CMISSNodesType) :: Nodes
  !Elements
  TYPE(CMISSMeshElementsType) :: MeshElementsGeometry
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
  TYPE(CMISSFieldType) :: DependentFieldDarcy
  TYPE(CMISSFieldType), ALLOCATABLE :: MaterialsFieldDarcy(:)
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsDarcy
  !Equations sets
!   TYPE(CMISSEquationsSetType) :: EquationsSetDarcy
  TYPE(CMISSEquationsSetType), ALLOCATABLE :: EquationsSetDarcy(:)
  TYPE(CMISSEquationsSetType) :: EquationsSetMatProperties
  !Equations
  TYPE(CMISSEquationsType), ALLOCATABLE :: EquationsDarcy(:)
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: LinearSolverDarcy
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsDarcy

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg), ALLOCATABLE :: EquationsSetIndices(:)
  INTEGER(CMISSIntg) :: Err
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:) :: VariableTypes


  INTEGER(CMISSIntg) :: DIAG_LEVEL_LIST(5)
  CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(1) !,TIMING_ROUTINE_LIST(1)


  TYPE(CMISSFieldType), ALLOCATABLE :: EquationsSetFieldDarcy(:)
  INTEGER(CMISSIntg) :: EquationsSetFieldDarcyUserNumber
  INTEGER(CMISSIntg) :: EquationsSetUserNumberDarcy
  INTEGER(CMISSIntg) :: icomp,Ncompartments,num_var,icompartment

  
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

  !PROBLEM CONTROL PANEL

  !Import cmHeart mesh information
  CALL FLUID_MECHANICS_IO_READ_CMHEART(CM,Err)  
  BASIS_NUMBER_GEOMETRY=CM%ID_M
  BASIS_NUMBER_VELOCITY=CM%ID_V
  BASIS_NUMBER_PRESSURE=CM%ID_P
  NUMBER_OF_DIMENSIONS=CM%D
  BASIS_TYPE=CM%IT_T
  BASIS_XI_INTERPOLATION_GEOMETRY=CM%IT_M
  BASIS_XI_INTERPOLATION_VELOCITY=CM%IT_V
  BASIS_XI_INTERPOLATION_PRESSURE=CM%IT_P
  NUMBER_OF_NODES_GEOMETRY=CM%N_M
  NUMBER_OF_NODES_VELOCITY=CM%N_V
  NUMBER_OF_NODES_PRESSURE=CM%N_P
  TOTAL_NUMBER_OF_NODES=CM%N_T
  TOTAL_NUMBER_OF_ELEMENTS=CM%E_T
  NUMBER_OF_ELEMENT_NODES_GEOMETRY=CM%EN_M
  NUMBER_OF_ELEMENT_NODES_VELOCITY=CM%EN_V
  NUMBER_OF_ELEMENT_NODES_PRESSURE=CM%EN_P
  !Set domain dimensions
  DOMAIN_X1 = -5.0_CMISSDP
  DOMAIN_X2 =  5.0_CMISSDP
  DOMAIN_Y1 = -5.0_CMISSDP
  DOMAIN_Y2 =  5.0_CMISSDP
  DOMAIN_Z1 = -5.0_CMISSDP
  DOMAIN_Z2 =  5.0_CMISSDP
  !Set geometric tolerance
  GEOMETRY_TOLERANCE = 1.0E-12_CMISSDP
  !Set initial values
  INITIAL_FIELD_DARCY(1)=0.0_CMISSDP
  INITIAL_FIELD_DARCY(2)=0.0_CMISSDP
  INITIAL_FIELD_DARCY(3)=0.0_CMISSDP
  !Set material parameters
  POROSITY_PARAM_DARCY=0.3_CMISSDP
  PERM_OVER_VIS_PARAM_DARCY=0.8_CMISSDP 
  !Set number of Gauss points (Mind that also material field may be interpolated)
  BASIS_XI_GAUSS_GEOMETRY=3 !4
  BASIS_XI_GAUSS_VELOCITY=3 !4
  BASIS_XI_GAUSS_PRESSURE=3 !4
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_DARCY_OUTPUT_TYPE=CMISSSolverSolverOutput
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DARCY_OUTPUT=CMISSEquationsMatrixOutput
  !Set solver parameters
  LINEAR_SOLVER_DARCY_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E5_CMISSDP !default: 1.0E5
  MAXIMUM_ITERATIONS=10000_CMISSIntg !default: 100000
  RESTART_VALUE=3000_CMISSIntg !default: 30

  Ncompartments = 4_CMISSIntg

  !
  !================================================================================================================================
  !

  !Set diagnostics

  DIAG_LEVEL_LIST(1)=1
  DIAG_LEVEL_LIST(2)=2
  DIAG_LEVEL_LIST(3)=3
  DIAG_LEVEL_LIST(4)=4
  DIAG_LEVEL_LIST(5)=5

  DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_FINITE_ELEMENT_CALCULATE"

  !CMISSAllDiagType/CMISSInDiagType/CMISSFromDiagType
!   CALL CMISSDiagnosticsSetOn(CMISSInDiagType,DIAG_LEVEL_LIST,"DarcyDiagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISSAllTimingType/CMISSInTimingType/CMISSFromTimingType
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

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

  !Start the creation of new bases: Geometry
  MESH_NUMBER_OF_COMPONENTS=1
  CALL CMISSBasisTypeInitialise(BasisGeometry,Err)
  CALL CMISSBasisCreateStart(BASIS_NUMBER_GEOMETRY,BasisGeometry,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasisTypeSet(BasisGeometry,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL CMISSBasisNumberOfXiSet(BasisGeometry,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  IF(NUMBER_OF_DIMENSIONS==2) THEN
    CALL CMISSBasisInterpolationXiSet(BasisGeometry,(/BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY/),Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisGeometry,(/BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY/),Err)
  ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
    CALL CMISSBasisInterpolationXiSet(BasisGeometry,(/BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY, & 
      & BASIS_XI_INTERPOLATION_GEOMETRY/),Err)                         
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisGeometry,(/BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY, &
      & BASIS_XI_GAUSS_GEOMETRY/),Err)
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(BasisGeometry,Err)
  !
  !Start the creation of another basis: Velocity
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisVelocity=BasisGeometry
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
  !
  !Start the creation of another basis: Pressure
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisPressure=BasisGeometry
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
  CALL CMISSNodesCreateStart(Region,TOTAL_NUMBER_OF_NODES,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL CMISSMeshNumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL CMISSMeshNumberOfComponentsSet(Mesh,MESH_NUMBER_OF_COMPONENTS,Err)
  !Specify spatial mesh component
  CALL CMISSMeshElementsTypeInitialise(MeshElementsGeometry,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsVelocity,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsPressure,Err)
  MESH_COMPONENT_NUMBER_GEOMETRY=1
  MESH_COMPONENT_NUMBER_VELOCITY=1
  MESH_COMPONENT_NUMBER_PRESSURE=1
  CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_GEOMETRY,BasisGeometry,MeshElementsGeometry,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL CMISSMeshElementsNodesSet(MeshElementsGeometry,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_GEOMETRY),Err)
  ENDDO
  CALL CMISSMeshElementsCreateFinish(MeshElementsGeometry,Err)
  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    MeshElementsVelocity=MeshElementsGeometry
  ELSE
    MESH_COMPONENT_NUMBER_VELOCITY=MESH_COMPONENT_NUMBER_GEOMETRY+1
    CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_VELOCITY,BasisVelocity,MeshElementsVelocity,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL CMISSMeshElementsNodesSet(MeshElementsVelocity,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_VELOCITY),Err)
    ENDDO
    CALL CMISSMeshElementsCreateFinish(MeshElementsVelocity,Err)
  ENDIF
  !Specify pressure mesh component
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    MeshElementsPressure=MeshElementsGeometry
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_GEOMETRY
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
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,DomainUserNumber,Err)
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
      & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
  ENDDO

  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1, & 
        & CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSFieldParameterSetUpdateStart(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  IF(.NOT.ALLOCATED(EquationsSetFieldDarcy)) ALLOCATE(EquationsSetFieldDarcy(Ncompartments))
  IF(.NOT.ALLOCATED(EquationsSetDarcy)) ALLOCATE(EquationsSetDarcy(Ncompartments))
  IF(.NOT.ALLOCATED(MaterialsFieldDarcy)) ALLOCATE(MaterialsFieldDarcy(Ncompartments))
  IF(.NOT.ALLOCATED(EquationsSetIndices)) ALLOCATE(EquationsSetIndices(Ncompartments))
  IF(.NOT.ALLOCATED(EquationsDarcy)) ALLOCATE(EquationsDarcy(Ncompartments))


  DO icomp = 1,Ncompartments
    EquationsSetFieldDarcyUserNumber = 100_CMISSIntg + icomp
    EquationsSetUserNumberDarcy = 200_CMISSIntg + icomp

    !Create the equations set field for Darcy - *** MULTIPLE COMPARTMENTS ***
    CALL CMISSFieldTypeInitialise(EquationsSetFieldDarcy(icomp),Err)

!   CALL CMISSFieldCreateStart(FieldFibreSolidUserNumber,Region,EquationsSetFieldDarcy,Err)
!   CALL CMISSFieldTypeSet(EquationsSetFieldDarcy,CMISSFieldFibreType,Err)
!   CALL CMISSFieldMeshDecompositionSet(EquationsSetFieldDarcy,Decomposition,Err)        
!   CALL CMISSFieldGeometricFieldSet(EquationsSetFieldDarcy,GeometricFieldSolid,Err)
!   CALL CMISSFieldNumberOfVariablesSet(EquationsSetFieldDarcy,FieldFibreSolidNumberOfVariables,Err)
!   CALL CMISSFieldNumberOfComponentsSet(EquationsSetFieldDarcy,CMISSFieldUVariableType,FieldFibreSolidNumberOfComponents,Err)  
!   CALL CMISSFieldComponentMeshComponentSet(EquationsSetFieldDarcy,CMISSFieldUVariableType,1,SolidMeshComponenetNumber,Err)
!   CALL CMISSFieldComponentMeshComponentSet(EquationsSetFieldDarcy,CMISSFieldUVariableType,2,SolidMeshComponenetNumber,Err)
!   CALL CMISSFieldComponentMeshComponentSet(EquationsSetFieldDarcy,CMISSFieldUVariableType,3,SolidMeshComponenetNumber,Err)
!   CALL CMISSFieldCreateFinish(EquationsSetFieldDarcy,Err)

    !These two lines have to go into each example file:
    CALL CMISSEquationsSetTypeInitialise(EquationsSetDarcy(icomp),Err)
    CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberDarcy,Region,GeometricField, &
      & CMISSEquationsSetFluidMechanicsClass,CMISSEquationsSetDarcyEquationType,CMISSEquationsSetMultiCompartmentDarcySubtype, &
      & EquationsSetFieldDarcyUserNumber,EquationsSetFieldDarcy(icomp),EquationsSetDarcy(icomp), &
      & Err)

!Create materials_field for coupling coefficients acccording to create_equations_set_field
!  for first material field (auto-)create it, subsequently just pass it in as we do for the
!  shared dependent field between daryc and elaasticy


    CALL CMISSEquationsSetCreateFinish(EquationsSetDarcy(icomp),Err)

    CALL CMISSFieldParameterSetUpdateConstant(EquationsSetFieldDarcy(icomp),CMISSFieldUVariableType, &
      & CMISSFieldValuesSetType,1,icomp,Err)
    CALL CMISSFieldParameterSetUpdateConstant(EquationsSetFieldDarcy(icomp),CMISSFieldUVariableType, &
      & CMISSFieldValuesSetType,2,Ncompartments,Err)
  ENDDO
  write (*,*) "Equations set fields now created."
  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS
! 
!     !Initialise dependent field (velocity components)
!     DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!       CALL CMISSFieldComponentValuesInitialise(DependentFieldDarcy(icomp),CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
!         & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
!     ENDDO
!   ENDDO

   CALL CMISSFieldTypeInitialise(DependentFieldDarcy,Err)
 

    CALL CMISSFieldCreateStart(DependentFieldUserNumberDarcy,Region,DependentFieldDarcy,Err)

    CALL CMISSFieldTypeSet(DependentFieldDarcy,CMISSFieldGeneralType,Err)  
    CALL CMISSFieldMeshDecompositionSet(DependentFieldDarcy,Decomposition,Err)
    CALL CMISSFieldGeometricFieldSet(DependentFieldDarcy,GeometricField,Err) 
    CALL CMISSFieldDependentTypeSet(DependentFieldDarcy,CMISSFieldDependentType,Err) 

    CALL CMISSFieldNumberOfVariablesSet(DependentFieldDarcy,2*Ncompartments,Err) 
    !create two variables for each compartment
    ALLOCATE(VariableTypes(2*Ncompartments))
    DO num_var=1,Ncompartments
       VariableTypes(2*num_var-1)=CMISSFieldUVariableType+(CMISSFieldNumberOfVariableSubtypes*(num_var-1))
       VariableTypes(2*num_var)=CMISSFieldDelUDelNVariableType+(CMISSFieldNumberOfVariableSubtypes*(num_var-1))
    ENDDO
    CALL CMISSFieldVariableTypesSet(DependentFieldDarcy,VariableTypes,Err) 

    DO icompartment=1,2*Ncompartments
      !set dimension type
      CALL CMISSFieldDimensionSet(DependentFieldDarcy,VariableTypes(icompartment), &
         & CMISSFieldVectorDimensionType,Err)
      CALL CMISSFieldNumberOfComponentsSet(DependentFieldDarcy,VariableTypes(icompartment),NUMBER_OF_DIMENSIONS+1,Err)
    ENDDO
    DO icompartment=1,2*Ncompartments
     DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,VariableTypes(icompartment),COMPONENT_NUMBER, & 
         & MESH_COMPONENT_NUMBER_VELOCITY,Err)
     ENDDO
    ENDDO
    DO icompartment=1,2*Ncompartments
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,VariableTypes(icompartment),NUMBER_OF_DIMENSIONS+1, & 
         & MESH_COMPONENT_NUMBER_PRESSURE,Err)
    ENDDO
    CALL CMISSFieldCreateFinish(DependentFieldDarcy,Err)
    
    DO icomp=1,Ncompartments
     
      CALL CMISSEquationsSetDependentCreateStart(EquationsSetDarcy(icomp),DependentFieldUserNumberDarcy,DependentFieldDarcy,Err)
      CALL CMISSEquationsSetDependentCreateFinish(EquationsSetDarcy(icomp),Err)

    ENDDO



  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  DO icomp = 1,Ncompartments
    !Create the equations set materials field variables for Static Darcy
    MaterialsFieldUserNumberDarcy = 400+icomp
    CALL CMISSFieldTypeInitialise(MaterialsFieldDarcy(icomp),Err)
    CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetDarcy(icomp),MaterialsFieldUserNumberDarcy, & 
      & MaterialsFieldDarcy(icomp),Err)
    !Finish the equations set materials field variables
    CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetDarcy(icomp),Err)
    CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy(icomp),CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
      & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
    CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy(icomp),CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
      & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !EQUATIONS

  DO icomp = 1,Ncompartments
    !Create the equations set equations
    CALL CMISSEquationsTypeInitialise(EquationsDarcy(icomp),Err)
    CALL CMISSEquationsSetEquationsCreateStart(EquationsSetDarcy(icomp),EquationsDarcy(icomp),Err)
    !Set the equations matrices sparsity type
    CALL CMISSEquationsSparsityTypeSet(EquationsDarcy(icomp),CMISSEquationsSparseMatrices,Err)
    !Set the equations set output
    CALL CMISSEquationsOutputTypeSet(EquationsDarcy(icomp),EQUATIONS_DARCY_OUTPUT,Err)
    !Finish the equations set equations
    CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetDarcy(icomp),Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a Static Darcy problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemFluidMechanicsClass,CMISSProblemDarcyEquationType, &
    & CMISSProblemStandardDarcySubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(LinearSolverDarcy,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  !Get the Darcy solver
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverDarcyUserNumber,LinearSolverDarcy,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(LinearSolverDarcy,LINEAR_SOLVER_DARCY_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_DARCY_DIRECT_FLAG) THEN
    CALL CMISSSolverLinearTypeSet(LinearSolverDarcy,CMISSSolverLinearDirectSolveType,Err)
    CALL CMISSSolverLibraryTypeSet(LinearSolverDarcy,CMISSSolverMUMPSLibrary,Err)
  ELSE
    CALL CMISSSolverLinearTypeSet(LinearSolverDarcy,CMISSSolverLinearIterativeSolveType,Err)
    CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolverDarcy,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolverLinearIterativeDivergenceToleranceSet(LinearSolverDarcy,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeRelativeToleranceSet(LinearSolverDarcy,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(LinearSolverDarcy,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeGMRESRestartSet(LinearSolverDarcy,RESTART_VALUE,Err)
  ENDIF
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL CMISSSolverTypeInitialise(LinearSolverDarcy,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsDarcy,Err)

  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the Darcy solver equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverDarcyUserNumber,LinearSolverDarcy,Err)
  CALL CMISSSolverSolverEquationsGet(LinearSolverDarcy,SolverEquationsDarcy,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsDarcy,CMISSSolverEquationsSparseMatrices,Err)

  DO icomp=1,Ncompartments
    CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsDarcy,EquationsSetDarcy(icomp), &
      & EquationsSetIndices(icomp),Err) !extend EquationsSetDarcy(icomp),EquationsSetIndices(icomp) to be arrays
  ENDDO

  !
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsDarcy,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)
  DO icomp = 1,Ncompartments

!     !--- BCs on normal velocity only
!     CONDITION = CMISSBoundaryConditionMovedWall
! 
!     IF( CM%D==2_CMISSIntg ) THEN
!       DO NODE_NUMBER=1_CMISSIntg,NUMBER_OF_NODES_GEOMETRY
!         COORD_X = CM%N(NODE_NUMBER,1_CMISSIntg)
!         COORD_Y = CM%N(NODE_NUMBER,2_CMISSIntg)
! 
!         IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) ) THEN
!           !x-velocity
!           VALUE = 1.0_CMISSDP
!           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMISSFieldUVariableType,CMISSNoGlobalDerivative, & 
!             & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
!           !x-velocity
!           VALUE = 1.0_CMISSDP
!           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMISSFieldUVariableType,CMISSNoGlobalDerivative, & 
!             & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
!           !y-velocity
!           VALUE = 2.0_CMISSDP
!           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMISSFieldUVariableType,CMISSNoGlobalDerivative, & 
!             & NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
!           !y-velocity
!           VALUE = 2.0_CMISSDP
!           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMISSFieldUVariableType,CMISSNoGlobalDerivative, & 
!             & NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!       END DO
!     ELSE IF( CM%D==3_CMISSIntg ) THEN
!       DO NODE_NUMBER=1_CMISSIntg,NUMBER_OF_NODES_GEOMETRY
!         COORD_X = CM%N(NODE_NUMBER,1_CMISSIntg)
!         COORD_Y = CM%N(NODE_NUMBER,2_CMISSIntg)
!         COORD_Z = CM%N(NODE_NUMBER,3_CMISSIntg)
! 
!         IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) ) THEN
!           !x-velocity
!           VALUE = 2.0_CMISSDP
!           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMISSFieldUVariableType,CMISSNoGlobalDerivative, & 
!             & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
!           !x-velocity
!           VALUE = 1.0_CMISSDP
!           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMISSFieldUVariableType,CMISSNoGlobalDerivative, & 
!             & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
!           !y-velocity
!           VALUE = 1.0_CMISSDP
!           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMISSFieldUVariableType,CMISSNoGlobalDerivative, & 
!             & NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
!           !y-velocity
!           VALUE = 2.0_CMISSDP
!           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMISSFieldUVariableType,CMISSNoGlobalDerivative, & 
!             & NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_Z-DOMAIN_Z1) < GEOMETRY_TOLERANCE) ) THEN
!           !z-velocity
!           VALUE = 1.0_CMISSDP
!           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMISSFieldUVariableType,CMISSNoGlobalDerivative, & 
!             & NODE_NUMBER,3_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!         !
!         IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) ) THEN
!           !z-velocity
!           VALUE = 1.0_CMISSDP
!           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldDarcy,CMISSFieldUVariableType,CMISSNoGlobalDerivative, & 
!             & NODE_NUMBER,3_CMISSIntg,CONDITION,VALUE,Err)
!         END IF
!       END DO
!     END IF

  ENDDO
  !Finish the creation of the equations set boundary conditions for Darcy
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsDarcy,Err)

  !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn off PETSc error handling
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
    CALL CMISSFieldIONodesExport(Fields,"StaticDarcy","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"StaticDarcy","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF

  !Finialise CMISS
!   CALL CMISSFinalise(Err)

  IF (ALLOCATED(EquationsSetFieldDarcy)) DEALLOCATE(EquationsSetFieldDarcy)
  IF (ALLOCATED(EquationsSetDarcy)) DEALLOCATE(EquationsSetDarcy)
  IF (ALLOCATED(MaterialsFieldDarcy)) DEALLOCATE(MaterialsFieldDarcy)
  IF (ALLOCATED(EquationsSetIndices)) DEALLOCATE(EquationsSetIndices)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM DARCYSTATICEXAMPLE
