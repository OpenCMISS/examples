!> \file
!> \author Christian Michler, Adam Reeve
!> \brief This is an example program to solve a coupled Finite Elastiticity Darcy equation using openCMISS calls.
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

!> \example MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDrivenDarcy/src/IncompressibleElasticityDrivenDarcyExample.f90
!! Example program to solve coupled FiniteElasticityDarcy equations using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDrivenDarcy/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDrivenDarcy/build-intel'>Linux GNU Build</a>
!!
!<

! !
! !  This example considers a coupled Finite Elasticity Darcy problem
! !

!> Main program

PROGRAM FINITEELASTICITYDARCYIOEXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OPENCMISS
!   USE FLUID_MECHANICS_IO_ROUTINES
  USE IOSTUFF
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

  REAL(CMISSDP), PARAMETER :: Y_DIM=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: X_DIM=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: Z_DIM=1.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: CubicBasisUserNumber=3

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcy=8
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberDarcy=12
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberDarcy=22

  INTEGER(CMISSIntg), PARAMETER :: ControlLoopSolidNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopFluidNumber=2
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopSubiterationNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverSolidIndex=1
  INTEGER(CMISSIntg), PARAMETER :: SolverDarcyIndex=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPorosity=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPermOverVis=2

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=3

  !Program types

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS

!   INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS

  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE

  INTEGER(CMISSIntg) :: EQUATIONS_DARCY_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER, NODE_NUMBER

  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_DARCY_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_DARCY_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_MAT_PROPERTIES_OUTPUT_TYPE

  REAL(CMISSDP) :: GEOMETRY_TOLERANCE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SOLID
  REAL(CMISSDP) :: INITIAL_FIELD_DARCY(4)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE
  REAL(CMISSDP) :: POROSITY_PARAM_DARCY, PERM_OVER_VIS_PARAM_DARCY
  REAL(CMISSDP) :: INNER_PRESSURE,OUTER_PRESSURE

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
!   TYPE(CMISSBasisType) :: CubicBasis, QuadraticBasis, LinearBasis, Bases(2)
  TYPE(CMISSBasisType),allocatable :: Bases(:)
  !Meshes
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh

  !Decompositions
  TYPE(CMISSDecompositionType) :: Decomposition
  !Fields
  TYPE(CMISSFieldsType) :: Fields
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: MaterialsFieldDarcy
  TYPE(CMISSFieldType) :: EquationsSetFieldDarcy
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsDarcy
  !Equations sets
  TYPE(CMISSEquationsSetType) :: EquationsSetDarcy
  !Equations
  TYPE(CMISSEquationsType) :: EquationsDarcy
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: DynamicSolverDarcy
  TYPE(CMISSSolverType) :: LinearSolverDarcy
  TYPE(CMISSSolverType) :: LinearSolverSolid
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsDarcy
  ! nodes and elements
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSMeshElementsType),allocatable :: Elements(:)

  !Other variables
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face1Nodes(:),Face2Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face3Nodes(:),Face4Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face5Nodes(:),Face6Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face7Nodes(:),Face8Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face9Nodes(:),Face10Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face11Nodes(:),Face12Nodes(:)
  INTEGER(CMISSIntg) :: FaceXi(6)
  INTEGER(CMISSIntg) :: NN,NODE,NodeDomain
  REAL(CMISSDP) :: XCoord,YCoord,ZCoord
  LOGICAL :: X_FIXED,Y_FIXED,X_OKAY,Y_OKAY

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables

  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: Err


  INTEGER(CMISSIntg) :: DIAG_LEVEL_LIST(5)
!   CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(8) !,TIMING_ROUTINE_LIST(1)
  CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(1) !,TIMING_ROUTINE_LIST(1)

  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !

  !Program variables and types (finite elasticity part)

  !Test program parameters

  INTEGER(CMISSIntg) :: SolidMeshComponenetNumber

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometrySolidUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometrySolidNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometrySolidNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreSolidUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreSolidNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreSolidNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialSolidUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialSolidNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialSolidNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentSolidUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentSolidNumberOfVariables=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentSolidNumberOfComponents=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentFluidNumberOfComponents=4  !(u,v,w,m)

  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldDarcyUserNumber=15

  INTEGER(CMISSIntg), PARAMETER :: EquationSetSolidUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldSolidUserNumber=25

  INTEGER(CMISSIntg), PARAMETER :: SolidDisplMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolidLagrMultMeshComponentNumber=2
  INTEGER(CMISSIntg), PARAMETER :: SolidGeometryMeshComponentNumber=SolidDisplMeshComponentNumber

  INTEGER(CMISSIntg), PARAMETER :: DarcyVelMeshComponentNumber=SolidLagrMultMeshComponentNumber
  INTEGER(CMISSIntg), PARAMETER :: DarcyMassIncreaseMeshComponentNumber=SolidLagrMultMeshComponentNumber
!   INTEGER(CMISSIntg), PARAMETER :: DarcyGeometryMeshComponentNumber=SolidDisplMeshComponentNumber

  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=32
  !Program types
  !Program variables

  REAL(CMISSDP) :: DYNAMIC_SOLVER_DARCY_START_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_DARCY_STOP_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_DARCY_THETA
  REAL(CMISSDP) :: DYNAMIC_SOLVER_DARCY_TIME_INCREMENT

  !CMISS variables

  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsSolid
  TYPE(CMISSEquationsType) :: EquationsSolid
  TYPE(CMISSEquationsSetType) :: EquationsSetSolid
  TYPE(CMISSFieldType) :: GeometricFieldSolid,FibreFieldSolid,MaterialFieldSolid
  TYPE(CMISSFieldType) :: DependentFieldSolid,EquationsSetFieldSolid
  TYPE(CMISSFieldType) :: IndependentFieldDarcy
  TYPE(CMISSSolverType) :: SolverSolid
  TYPE(CMISSSolverEquationsType) :: SolverEquationsSolid

  INTEGER(CMISSIntg),allocatable :: surface_lin_inner(:),surface_lin_outer(:),surface_lin_base(:)
  INTEGER(CMISSIntg),allocatable :: surface_quad_inner(:),surface_quad_outer(:),surface_quad_base(:)

  !End - Program variables and types (finite elasticity part)

  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !


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
  NUMBER_GLOBAL_X_ELEMENTS=1
  NUMBER_GLOBAL_Y_ELEMENTS=1
  NUMBER_GLOBAL_Z_ELEMENTS=3

  IF(NUMBER_GLOBAL_Z_ELEMENTS==0)THEN
    NUMBER_OF_DIMENSIONS=2
  ELSE
    NUMBER_OF_DIMENSIONS=3
  ENDIF
  !PROBLEM CONTROL PANEL

!   BASIS_XI_INTERPOLATION_SOLID=CMISSBasisLinearLagrangeInterpolation
  BASIS_XI_INTERPOLATION_SOLID=CMISSBasisQuadraticLagrangeInterpolation
  !Set geometric tolerance
  GEOMETRY_TOLERANCE = 1.0E-12_CMISSDP
  !Set initial values
  INITIAL_FIELD_DARCY(1)=0.0_CMISSDP
  INITIAL_FIELD_DARCY(2)=0.0_CMISSDP
  INITIAL_FIELD_DARCY(3)=0.0_CMISSDP
  INITIAL_FIELD_DARCY(4)=0.0_CMISSDP
  !Set material parameters
  POROSITY_PARAM_DARCY=0.1_CMISSDP
  PERM_OVER_VIS_PARAM_DARCY=1.0_CMISSDP
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE=CMISSSolverProgressOutput
  LINEAR_SOLVER_DARCY_OUTPUT_TYPE=CMISSSolverSolverOutput
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DARCY_OUTPUT=CMISSEquationsNoOutput

  !Set time parameter
  DYNAMIC_SOLVER_DARCY_START_TIME=0.0_CMISSDP
  DYNAMIC_SOLVER_DARCY_TIME_INCREMENT=1.0e-4_CMISSDP
  DYNAMIC_SOLVER_DARCY_STOP_TIME=100_CMISSIntg * DYNAMIC_SOLVER_DARCY_TIME_INCREMENT
  DYNAMIC_SOLVER_DARCY_THETA=1.0_CMISSDP !2.0_CMISSDP/3.0_CMISSDP
  !Set result output parameter
  DYNAMIC_SOLVER_DARCY_OUTPUT_FREQUENCY=1
  !Set solver parameters
  LINEAR_SOLVER_DARCY_DIRECT_FLAG=.TRUE.
  RELATIVE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E5_CMISSDP !default: 1.0E5
  MAXIMUM_ITERATIONS=10000_CMISSIntg !default: 100000
  RESTART_VALUE=30_CMISSIntg !default: 30
  LINESEARCH_ALPHA=1.0_CMISSDP


  !
  !================================================================================================================================
  !

  !INITIALISE OPENCMISS

  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  !
  !================================================================================================================================
  !

  !Set diagnostics

  DIAG_LEVEL_LIST(1)=1
  DIAG_LEVEL_LIST(2)=2
  DIAG_LEVEL_LIST(3)=3
  DIAG_LEVEL_LIST(4)=4
  DIAG_LEVEL_LIST(5)=5

  DIAG_ROUTINE_LIST(1)="WRITE_IP_INFO"
!   DIAG_ROUTINE_LIST(2)="FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR"

  !CMISSAllDiagType/CMISSInDiagType/CMISSFromDiagType
  CALL CMISSDiagnosticsSetOn(CMISSInDiagType,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISSAllTimingType/CMISSInTimingType/CMISSFromTimingType
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

  !
  !================================================================================================================================
  !

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains = NumberOfComputationalNodes
  write(*,*) "NumberOfDomains = ",NumberOfDomains

  !
  !================================================================================================================================
  !

  !COORDINATE SYSTEM

  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !REGION
  !For a volume-coupled problem, solid and fluid are based in the same region

  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

call READ_MESH('input/pigheart-mesh',MeshUserNumber,Region, Mesh,Bases,Nodes,Elements)
!   !BASES
!   !Define basis functions
!   CALL CMISSBasisTypeInitialise(LinearBasis,Err)
!   CALL CMISSBasisCreateStart(LinearBasisUserNumber,LinearBasis,Err)
!   CALL CMISSBasisQuadratureNumberOfGaussXiSet(LinearBasis, &
!     & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
!   !CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err)
!   CALL CMISSBasisCreateFinish(LinearBasis,Err)
! 
!   CALL CMISSBasisTypeInitialise(QuadraticBasis,Err)
!   CALL CMISSBasisCreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
!   CALL CMISSBasisInterpolationXiSet(QuadraticBasis,(/CMISSBasisQuadraticLagrangeInterpolation, &
!     & CMISSBasisQuadraticLagrangeInterpolation,CMISSBasisQuadraticLagrangeInterpolation/),Err)
!   CALL CMISSBasisQuadratureNumberOfGaussXiSet(QuadraticBasis, &
!     & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
!   !CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err)
!   CALL CMISSBasisCreateFinish(QuadraticBasis,Err)
! 
!   CALL CMISSBasisTypeInitialise(CubicBasis,Err)
!   CALL CMISSBasisCreateStart(CubicBasisUserNumber,CubicBasis,Err)
!   CALL CMISSBasisInterpolationXiSet(CubicBasis,(/CMISSBasisCubicLagrangeInterpolation, &
!     & CMISSBasisCubicLagrangeInterpolation,CMISSBasisCubicLagrangeInterpolation/),Err)
!   CALL CMISSBasisQuadratureNumberOfGaussXiSet(CubicBasis, &
!     & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
!   !CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(CubicBasis,.true.,Err) !Enable 3D interpolation on faces
!   CALL CMISSBasisCreateFinish(CubicBasis,Err)
! 
!   !LinearBasis/QuadraticBasis/CubicBasis
!   Bases(1)=QuadraticBasis
!   Bases(2)=LinearBasis
! 
!   !Start the creation of a generated mesh in the region
!   CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
!   CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
!   CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
!   CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Bases,Err)
!   IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
!     CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/X_DIM,Y_DIM/),Err)
!     CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/),Err)
!   ELSE
!     CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/X_DIM,Y_DIM,Z_DIM/),Err)
!     CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
!       & NUMBER_GLOBAL_Z_ELEMENTS/),Err)
!   ENDIF
!   CALL CMISSMeshTypeInitialise(Mesh,Err)
!   CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
! -------------------------------------------------------------------------------
  

  !GEOMETRIC FIELD

  !Create a decomposition:
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecompositionCalculateFacesSet(Decomposition,.true.,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)
  CALL CMISSFieldNumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,3,Err)  
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)
CALL READ_NODES('input/pigheart-nodes',GeometricField)
!   CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create a decomposition

  !Create a field to put the geometry (defualt is geometry)

  SolidMeshComponenetNumber = SolidGeometryMeshComponentNumber

  CALL CMISSFieldTypeInitialise(GeometricFieldSolid,Err)
  CALL CMISSFieldCreateStart(FieldGeometrySolidUserNumber,Region,GeometricFieldSolid,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricFieldSolid,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricFieldSolid,CMISSFieldGeometricType,Err)
  CALL CMISSFieldNumberOfVariablesSet(GeometricFieldSolid,FieldGeometrySolidNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricFieldSolid,CMISSFieldUVariableType,FieldGeometrySolidNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,1,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,2,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,3,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldCreateFinish(GeometricFieldSolid,Err)
  !Set the mesh component to be used by the field components.
CALL READ_NODES('input/pigheart-nodes',GeometricFieldSolid)
!  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricFieldSolid,GeneratedMesh,Err)

  !Create a fibre field and attach it to the geometric field
  CALL CMISSFieldTypeInitialise(FibreFieldSolid,Err)
  CALL CMISSFieldCreateStart(FieldFibreSolidUserNumber,Region,FibreFieldSolid,Err)
  CALL CMISSFieldTypeSet(FibreFieldSolid,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreFieldSolid,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(FibreFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSFieldNumberOfVariablesSet(FibreFieldSolid,FieldFibreSolidNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreFieldSolid,CMISSFieldUVariableType,FieldFibreSolidNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,1,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,2,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,3,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(FibreFieldSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for ALE Darcy
  CALL CMISSFieldTypeInitialise(EquationsSetFieldDarcy,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSetDarcy,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberDarcy,Region,GeometricField,CMISSEquationsSetFluidMechanicsClass, &
    & CMISSEquationsSetDarcyEquationType,CMISSEquationsSetIncompressibleElasticityDrivenDarcySubtype,&
    & EquationsSetFieldUserNumberDarcy,EquationsSetFieldDarcy,EquationsSetDarcy,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSetDarcy,Err)

  !Create the equations set for the solid
  CALL CMISSFieldTypeInitialise(EquationsSetFieldSolid,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSetSolid,Err)
  CALL CMISSEquationsSetCreateStart(EquationSetSolidUserNumber,Region,FibreFieldSolid,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetIncompressibleElasticityDrivenDarcySubtype,&
    & EquationsSetFieldSolidUserNumber,EquationsSetFieldSolid,EquationsSetSolid,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSetSolid,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid Materials Field

  !Create a material field and attach it to the geometric field
  CALL CMISSFieldTypeInitialise(MaterialFieldSolid,Err)
  !
  CALL CMISSFieldCreateStart(FieldMaterialSolidUserNumber,Region,MaterialFieldSolid,Err)
  !
  CALL CMISSFieldTypeSet(MaterialFieldSolid,CMISSFieldMaterialType,Err)
  CALL CMISSFieldMeshDecompositionSet(MaterialFieldSolid,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(MaterialFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialFieldSolid,FieldMaterialSolidNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialFieldSolid,CMISSFieldUVariableType,FieldMaterialSolidNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldSolid,CMISSFieldUVariableType,1,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldSolid,CMISSFieldUVariableType,2,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldSolid,CMISSFieldUVariableType,3,SolidGeometryMeshComponentNumber,Err)
  !
  CALL CMISSFieldCreateFinish(MaterialFieldSolid,Err)

  !Set material parameters
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2.0_CMISSDP,Err)
!   CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2.0e3_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,6.0_CMISSDP,Err)
!   CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,33.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,10.0_CMISSDP,Err)


  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetSolid,FieldMaterialSolidUserNumber,MaterialFieldSolid,Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------


  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create a dependent field with four variables (U, DelUDelN = solid, V, DelVDelN = Darcy) and four components
  CALL CMISSFieldTypeInitialise(DependentFieldSolid,Err)
  !
  CALL CMISSFieldCreateStart(FieldDependentSolidUserNumber,Region,DependentFieldSolid,Err)
  !
  CALL CMISSFieldTypeSet(DependentFieldSolid,CMISSFieldGeneralType,Err)
  CALL CMISSFieldMeshDecompositionSet(DependentFieldSolid,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(DependentFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSFieldDependentTypeSet(DependentFieldSolid,CMISSFieldDependentType,Err)
  CALL CMISSFieldNumberOfVariablesSet(DependentFieldSolid,FieldDependentSolidNumberOfVariables,Err)
  CALL CMISSFieldVariableTypesSet(DependentFieldSolid,(/CMISSFieldUVariableType, &
    & CMISSFieldDelUDelNVariableType,CMISSFieldVVariableType,CMISSFieldDelVDelNVariableType/),Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldUVariableType,FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldVVariableType,FieldDependentFluidNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,FieldDependentFluidNumberOfComponents,Err)
  !
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,1,SolidDisplMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,2,SolidDisplMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,3,SolidDisplMeshComponentNumber,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,4,CMISSFieldNodeBasedInterpolation,Err)
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,4,CMISSFieldElementBasedInterpolation,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,4,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,4,SolidLagrMultMeshComponentNumber,Err)
  !
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,1,SolidDisplMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,2,SolidDisplMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,3,SolidDisplMeshComponentNumber,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,4, &
    & CMISSFieldNodeBasedInterpolation,Err)
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,4, &
!     & CMISSFieldElementBasedInterpolation,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,4,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,4,SolidLagrMultMeshComponentNumber, &
    & Err)

  !For this equation type, MESH_COMPONENT_NUMBER_PRESSURE is actually the mass increase component as the pressure is taken from the solid equations
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,1,DarcyVelMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,2,DarcyVelMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,3,DarcyVelMeshComponentNumber,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,4,MESH_COMPONENT_NUMBER_PRESSURE,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,4,DarcyMassIncreaseMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,DarcyVelMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,2,DarcyVelMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,3,DarcyVelMeshComponentNumber,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,4,MESH_COMPONENT_NUMBER_PRESSURE,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,4, &
    & DarcyMassIncreaseMeshComponentNumber,Err)

  CALL CMISSFieldCreateFinish(DependentFieldSolid,Err)
  !
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetSolid,FieldDependentSolidUserNumber,DependentFieldSolid,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------


  !Create the equations set dependent field variables for ALE Darcy
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetDarcy,FieldDependentSolidUserNumber, & ! ??? UserNumber ???
    & DependentFieldSolid,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetDarcy,Err)

  !Initialise dependent field (velocity components,mass increase)
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS+1
    CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldVVariableType,CMISSFieldValuesSetType, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !INDEPENDENT FIELD Darcy for storing BC flags

  CALL CMISSFieldTypeInitialise(IndependentFieldDarcy,Err)
  CALL CMISSEquationsSetIndependentCreateStart(EquationsSetDarcy,IndependentFieldDarcyUserNumber, &
    & IndependentFieldDarcy,Err)

  CALL CMISSFieldComponentMeshComponentSet(IndependentFieldDarcy,CMISSFieldUVariableType,1,DarcyVelMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(IndependentFieldDarcy,CMISSFieldUVariableType,2,DarcyVelMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(IndependentFieldDarcy,CMISSFieldUVariableType,3,DarcyVelMeshComponentNumber,Err)
!   CALL CMISSFieldComponentMeshComponentSet(IndependentFieldDarcy,CMISSFieldUVariableType,4,DarcyMassIncreaseMeshComponentNumber,Err)

  CALL CMISSEquationsSetIndependentCreateFinish(EquationsSetDarcy,Err)

  CALL CMISSFieldComponentValuesInitialise(IndependentFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(IndependentFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(IndependentFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,0.0_CMISSDP,Err)

  !
  !================================================================================================================================


  !
  !================================================================================================================================
  !

call READ_FIELD('input/pigheart-material',MaterialsFieldUserNumberDarcy,Region,GeometricField,MaterialsFieldDarcy)
CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetDarcy,MaterialsFieldUserNumberDarcy,MaterialsFieldDarcy,Err)
CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetDarcy,Err)
!   !MATERIALS FIELDS
! 
!   !Create the equations set materials field variables for ALE Darcy
!   CALL CMISSFieldTypeInitialise(MaterialsFieldDarcy,Err)
!   CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetDarcy,MaterialsFieldUserNumberDarcy, &
!     & MaterialsFieldDarcy,Err)
!   !Finish the equations set materials field variables
!   CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetDarcy,Err)
!   CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
!     & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
!   CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
!     & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SET EQUATIONS

  !Darcy
  CALL CMISSEquationsTypeInitialise(EquationsDarcy,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetDarcy,EquationsDarcy,Err)
  CALL CMISSEquationsSparsityTypeSet(EquationsDarcy,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(EquationsDarcy,EQUATIONS_DARCY_OUTPUT,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetDarcy,Err)

  !Solid
  CALL CMISSEquationsTypeInitialise(EquationsSolid,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetSolid,EquationsSolid,Err)
  CALL CMISSEquationsSparsityTypeSet(EquationsSolid,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(EquationsSolid,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 1,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,0.0_CMISSDP, &
    & Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  !
  !================================================================================================================================
  !

  !PROBLEMS

  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemMultiPhysicsClass,CMISSProblemFiniteElasticityDarcyType, &
    & CMISSProblemQuasistaticElasticityTransientDarcySubtype,Err)
  CALL CMISSProblemCreateFinish(Problem,Err)

  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
!   CALL CMISSControlLoopMaximumIterationsSet(ControlLoop,2,Err)
  CALL CMISSControlLoopTimesSet(ControlLoop,DYNAMIC_SOLVER_DARCY_START_TIME,DYNAMIC_SOLVER_DARCY_STOP_TIME, &
    & DYNAMIC_SOLVER_DARCY_TIME_INCREMENT,Err)
  CALL CMISSControlLoopTimeOutputSet(ControlLoop,DYNAMIC_SOLVER_DARCY_OUTPUT_FREQUENCY,Err)
!   CALL CMISSControlLoopOutputTypeSet(ControlLoop,CMISSControlLoopProgressOutput,Err)
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVERS

  CALL CMISSSolverTypeInitialise(SolverSolid,Err)
  CALL CMISSSolverTypeInitialise(DynamicSolverDarcy,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverDarcy,Err)

  CALL CMISSProblemSolversCreateStart(Problem,Err)

  ! Solid
  CALL CMISSProblemSolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopSolidNumber,CMISSControlLoopNode/), &
    & SolverSolidIndex,SolverSolid,Err)
  CALL CMISSSolverOutputTypeSet(SolverSolid,CMISSSolverProgressOutput,Err)
!   CALL CMISSSolverNewtonJacobianCalculationTypeSet(SolverSolid,CMISSSolverNewtonJacobianFDCalculated,Err)
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(SolverSolid,CMISSSolverNewtonJacobianAnalyticCalculated,Err)

  CALL CMISSSolverNewtonAbsoluteToleranceSet(SolverSolid,ABSOLUTE_TOLERANCE,Err)
  CALL CMISSSolverNewtonRelativeToleranceSet(SolverSolid,RELATIVE_TOLERANCE,Err)
  CALL CMISSSolverNewtonMaximumIterationsSet(SolverSolid,MAXIMUM_ITERATIONS,Err)

!   CALL CMISSSolverNonLinearTypeSet(SolverSolid,CMISSSolverNonlinearNewton,Err)
!   CALL CMISSSolverLibraryTypeSet(SolverSolid,CMISSSolverPETScLibrary,Err)

  CALL CMISSSolverNewtonLinearSolverGet(SolverSolid,LinearSolverSolid,Err)
  CALL CMISSSolverLinearTypeSet(LinearSolverSolid,CMISSSolverLinearDirectSolveType,Err)


  !Darcy
  CALL CMISSProblemSolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMISSControlLoopNode/), &
    & SolverDarcyIndex,DynamicSolverDarcy,Err)
  CALL CMISSSolverOutputTypeSet(DynamicSolverDarcy,DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE,Err)
  CALL CMISSSolverDynamicThetaSet(DynamicSolverDarcy,DYNAMIC_SOLVER_DARCY_THETA,Err)
!   CALL CMISSSolverDynamicDynamicSet(DynamicSolverDarcy,.TRUE.,Err)
  CALL CMISSSolverDynamicLinearSolverGet(DynamicSolverDarcy,LinearSolverDarcy,Err)
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

  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  CALL CMISSSolverTypeInitialise(SolverSolid,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverDarcy,Err)

  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsSolid,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsDarcy,Err)

  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !
  !Get the finite elasticity solver equations
  CALL CMISSProblemSolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopSolidNumber,CMISSControlLoopNode/), &
    & SolverSolidIndex,SolverSolid,Err)
  CALL CMISSSolverSolverEquationsGet(SolverSolid,SolverEquationsSolid,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsSolid,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsSolid,EquationsSetSolid,EquationsSetIndex,Err)
  !
  !Get the Darcy solver equations
  CALL CMISSProblemSolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMISSControlLoopNode/), &
    & SolverDarcyIndex,LinearSolverDarcy,Err)
  CALL CMISSSolverSolverEquationsGet(LinearSolverDarcy,SolverEquationsDarcy,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsDarcy,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsDarcy,EquationsSetDarcy,EquationsSetIndex,Err)
  !
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !------------------------------------
  ! Grab surfaces to assign boundary conditions

  CALL READ_SURFACE('input/surface-lin-inner',surface_lin_inner)
  CALL READ_SURFACE('input/surface-lin-outer',surface_lin_outer)
  CALL READ_SURFACE('input/surface-lin-base',surface_lin_base)
  CALL READ_SURFACE('input/surface-quad-inner',surface_quad_inner)
  CALL READ_SURFACE('input/surface-quad-outer',surface_quad_outer)
  CALL READ_SURFACE('input/surface-quad-base',surface_quad_base)

  !------------------------------------
  ! ASSIGN BOUNDARY CONDITIONS - SOLID (absolute nodal parameters)
  !Solid is computed in absolute position, rather than displacement. Thus BCs for absolute position
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsSolid,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsSolid,BoundaryConditionsSolid,Err)

  write(*,*)'surface_quad_base = ',surface_quad_base

  INNER_PRESSURE = 0.01_CMISSDP
  OUTER_PRESSURE = 0.0_CMISSDP

  ! ASSIGN BOUNDARY CONDITIONS
  !Fix base of the ellipsoid in z direction
  DO NN=1,SIZE(surface_quad_base,1)
    NODE=surface_quad_base(NN)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetGetNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NODE,3,ZCoord,Err)
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldUVariableType,1,1,NODE,3, &
        & CMISSBoundaryConditionFixed,ZCoord,Err)
    ENDIF
  ENDDO

!   !Apply inner surface pressure
!   !NOTE: Surface pressure goes into pressure_values_set_type of the DELUDELN type
!   DO NN=1,SIZE(surface_quad_inner,1)
!     NODE=surface_quad_inner(NN)
!     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
!     IF(NodeDomain==ComputationalNodeNumber) THEN
!       COMPONENT_NUMBER = 3  ! Does it matter which number ??? It used to be linked to the normal ... Check this !!!
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldDelUDelNVariableType,1,NODE,COMPONENT_NUMBER, &
!         & CMISSBoundaryConditionPressure,INNER_PRESSURE,Err)
!     ENDIF
!   ENDDO
!
!   !Apply outer surface pressure
!   DO NN=1,SIZE(surface_quad_outer,1)
!     NODE=surface_quad_outer(NN)
!     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
!     IF(NodeDomain==ComputationalNodeNumber) THEN
!       COMPONENT_NUMBER = 3  ! Does it matter which number ??? It used to be linked to the normal ... Check this !!!
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldDelUDelNVariableType,1,NODE,COMPONENT_NUMBER, &
!         & CMISSBoundaryConditionPressure,OUTER_PRESSURE,Err)
!     ENDIF
!   ENDDO

  !Fix more nodes at the base to stop free body motion: 600 in x, 584 in y
  X_FIXED=.FALSE.
  Y_FIXED=.FALSE.
  DO NN=1,SIZE(surface_quad_base,1)
    NODE=surface_quad_base(NN)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetGetNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NODE,1,XCoord,Err)
      CALL CMISSFieldParameterSetGetNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NODE,2,YCoord,Err)
!       IF(ABS(XCoord)<1.0E-6_CMISSDP) THEN
      IF(NODE==600) THEN
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldUVariableType,1,1,NODE,1, &
          & CMISSBoundaryConditionFixed,XCoord,Err)
        WRITE(*,*) "FIXING NODE",NODE,"IN X DIRECTION"
        X_FIXED=.TRUE.
      ENDIF
!       IF(ABS(YCoord)<1.0E-6_CMISSDP) THEN
      IF(NODE==584) THEN
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldUVariableType,1,1,NODE,2, &
          & CMISSBoundaryConditionFixed,YCoord,Err)
        WRITE(*,*) "FIXING NODE",NODE,"IN Y DIRECTION"
        Y_FIXED=.TRUE.
    ENDIF
    ENDIF
  ENDDO
!   CALL MPI_REDUCE(X_FIXED,X_OKAY,1,MPI_LOGICAL,MPI_LOR,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_REDUCE(Y_FIXED,Y_OKAY,1,MPI_LOGICAL,MPI_LOR,0,MPI_COMM_WORLD,MPI_IERROR)
  IF(ComputationalNodeNumber==0) THEN
!     IF(.NOT.(X_OKAY.AND.Y_OKAY)) THEN
    IF(.NOT.(X_FIXED.AND.Y_FIXED)) THEN
      WRITE(*,*) "Free body motion could not be prevented!"
      CALL CMISSFinalise(Err)
    STOP
    ENDIF
  ENDIF

  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsSolid,Err)

  !
  !================================================================================================================================
  !

  !BCs Darcy
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsDarcy,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)

    !In 'generated_mesh_routines.f90/GENERATED_MESH_ELLIPSOID_SURFACE_GET' there is a bug:
    !  BASIS=>ELLIPSOID_MESH%BASES(MESH_COMPONENT)%PTR does not account for the fact that:
    !  in 'generated_mesh_routines.f90/GENERATED_MESH_ELLIPSOID_CREATE_FINISH' the following is done:
    !  CALL MESH_NUMBER_OF_COMPONENTS_SET(GENERATED_MESH%MESH,SIZE(ELLIPSOID_MESH%BASES)/2,ERR,ERROR,*999)
    !Temporary work around, until bug fix:


    !  I N N E R   S U R F A C E

    write(*,*)'surface_lin_inner = ',surface_lin_inner

    !Set all inner surface nodes impermeable
    !MIND: CMISSFieldDelVDelNVariableType -> RHS invoked in DARCY_EQUATION_FINITE_ELEMENT_CALCULATE
    !      CMISSBoundaryConditionImpermeableWall
    DO NN=1,SIZE(surface_lin_inner,1)
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 1
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,surface_lin_inner(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", surface_lin_inner(NN)
!
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 2
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,surface_lin_inner(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", surface_lin_inner(NN)
!
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 3
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,surface_lin_inner(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", surface_lin_inner(NN)

!       VALUE = 1.0_CMISSDP
!       COMPONENT_NUMBER = ABS(InnerNormalXi)
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,surface_lin_inner(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", surface_lin_inner(NN)

      NODE_NUMBER = surface_lin_inner(NN)
      COMPONENT_NUMBER = 1
      VALUE = 1.0_CMISSDP
      CALL CMISSFieldParameterSetUpdateNode(IndependentFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType,1, &
        & CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO


    !  O U T E R   S U R F A C E

    write(*,*)'surface_lin_outer = ',surface_lin_outer

    !Set all outer surface nodes impermeable
    DO NN=1,SIZE(surface_lin_outer,1)
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 1
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,surface_lin_outer(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", surface_lin_outer(NN)
!
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 2
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,surface_lin_outer(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", surface_lin_outer(NN)
!
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 3
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,surface_lin_outer(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", surface_lin_outer(NN)

!       VALUE = 1.0_CMISSDP
!       COMPONENT_NUMBER = ABS(OuterNormalXi)
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,surface_lin_outer(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", surface_lin_outer(NN)

      NODE_NUMBER = surface_lin_outer(NN)
      COMPONENT_NUMBER = 1
      VALUE = 1.0_CMISSDP
      CALL CMISSFieldParameterSetUpdateNode(IndependentFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType,1, &
        & CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO


    !  T O P   S U R F A C E

    write(*,*)'surface_lin_base = ',surface_lin_base

    !Set all top surface nodes to Darcy inflow BC
    DO NN=1,SIZE(surface_lin_base,1)
!       VALUE = +1.0_CMISSDP  ! Mind the sign !
      VALUE = 0.0_CMISSDP  ! Mind the sign !
      COMPONENT_NUMBER = 3
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,1, &
        & surface_lin_base(NN), &
        & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING TOP DARCY BC TO NODE", surface_lin_base(NN)
    ENDDO

    DO NN=1,SIZE(surface_lin_inner,1)
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 1
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,InnerSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)
!
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 2
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,InnerSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)
!
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 3
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,InnerSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)

!       VALUE = 1.0_CMISSDP
!       COMPONENT_NUMBER = ABS(InnerNormalXi)
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,InnerSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)

      NODE_NUMBER = surface_lin_inner(NN)
      !COMPONENT_NUMBER = ABS(InnerNormalXi)
      COMPONENT_NUMBER=1_CMISSIntg
      VALUE = 1.0_CMISSDP
      CALL CMISSFieldParameterSetUpdateNode(IndependentFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType,1, &
        & CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)

    ENDDO
    DO NN=1,SIZE(surface_lin_outer,1)
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 1
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)
!
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 2
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)
!
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 3
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)

!       VALUE = 1.0_CMISSDP
!       COMPONENT_NUMBER = ABS(OuterNormalXi)
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)

      NODE_NUMBER = surface_lin_outer(NN)
      !COMPONENT_NUMBER = ABS(OuterNormalXi)
      COMPONENT_NUMBER=1_CMISSIntg
      VALUE = 1.0_CMISSDP
      CALL CMISSFieldParameterSetUpdateNode(IndependentFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType,1, &
        & CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO



  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsDarcy,Err)

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

  EXPORT_FIELD_IO=.FALSE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"FiniteElasticityDarcy","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"FiniteElasticityDarcy","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF


  !Finialise CMISS
!   CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM FINITEELASTICITYDARCYIOEXAMPLE
