!> \file
!> \author Christian Michler, Adam Reeve, Andrew Cookson
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
! !  This example considers a coupled Finite Elasticity Multi-Compartment Darcy problem
! !

!> Main program

PROGRAM FINITEELASTICITYMULTICOMPDARCYIOEXAMPLE

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
  REAL(CMISSDP), PARAMETER :: Z_DIM=3.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: CubicBasisUserNumber=3

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberDarcy=6
  INTEGER(CMISSIntg) :: MaterialsFieldUserNumberDarcy
  INTEGER(CMISSIntg) :: EquationsSetUserNumberDarcy
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=14
  INTEGER(CMISSIntg) :: EquationsSetFieldUserNumberDarcy
  INTEGER(CMISSIntg) :: icompartment,Ncompartments,num_var,componentnum,Nparams
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

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS

!   INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  !TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE

  INTEGER(CMISSIntg) :: EQUATIONS_DARCY_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_DARCY_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_DARCY_OUTPUT_TYPE

  REAL(CMISSDP) :: COORD_X, COORD_Y, COORD_Z
  REAL(CMISSDP) :: DOMAIN_X1, DOMAIN_X2, DOMAIN_Y1, DOMAIN_Y2, DOMAIN_Z1, DOMAIN_Z2
  REAL(CMISSDP) :: GEOMETRY_TOLERANCE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SOLID
  REAL(CMISSDP) :: INITIAL_FIELD_DARCY(4)
  REAL(CMISSDP) :: INITIAL_FIELD_SOLID(4)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
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
  !TYPE(CMISSBasisType) :: CubicBasis, QuadraticBasis, LinearBasis, Bases(2)
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
  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: MaterialsFieldDarcy
  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: EquationsSetFieldDarcy
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsDarcy
  !Equations sets
  TYPE(CMISSEquationsSetType), ALLOCATABLE, DIMENSION(:) :: EquationsSetDarcy
  !Equations
  TYPE(CMISSEquationsType), ALLOCATABLE, DIMENSION(:) :: EquationsDarcy
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: DynamicSolverDarcy
  TYPE(CMISSSolverType) :: LinearSolverDarcy
!   TYPE(CMISSSolverType) :: LinearSolverSolid
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
  LOGICAL :: X_FIXED,Y_FIXED !,X_OKAY,Y_OKAY

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables

  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: Err
  !Array containing the field variable types that will be used (for ease of incorporating inside a loop)
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:) :: VariableTypes
  REAL(CMISSDP), ALLOCATABLE, DIMENSION(:,:) :: CouplingCoeffs,ConstitutiveParams

  INTEGER(CMISSIntg) :: DIAG_LEVEL_LIST(5)
!   CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(8) !,TIMING_ROUTINE_LIST(1)
  CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(1) !,TIMING_ROUTINE_LIST(1)

  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !

  !Program variables and types (finite elasticity part)

  !Test program parameters

  INTEGER(CMISSIntg) :: SolidMeshComponenetNumber

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometrySolidUserNumber=51
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometrySolidNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometrySolidNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreSolidUserNumber=52
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreSolidNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreSolidNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialSolidUserNumber=53
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialSolidNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialSolidNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentSolidUserNumber=54
  INTEGER(CMISSIntg) :: FieldDependentSolidNumberOfVariables
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentSolidNumberOfComponents=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentFluidNumberOfComponents=4  !(u,v,w,m)

  INTEGER(CMISSIntg), PARAMETER :: EquationSetSolidUserNumber=55
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

!   TYPE(CMISSBasisType) :: BasisSolid
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsSolid
  TYPE(CMISSEquationsType) :: EquationsSolid
  TYPE(CMISSEquationsSetType) :: EquationsSetSolid
  TYPE(CMISSFieldType) :: GeometricFieldSolid,FibreFieldSolid,MaterialFieldSolid,DependentFieldSolid,EquationsSetFieldSolid
  TYPE(CMISSSolverType) :: SolverSolid
  TYPE(CMISSSolverEquationsType) :: SolverEquationsSolid
  TYPE(CMISSMeshElementsType) :: MeshElementsSolid

  !End - Program variables and types (finite elasticity part)

  INTEGER(CMISSIntg) :: dummy_counter

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

  !INITIALISE OPENCMISS

  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  !
  !================================================================================================================================
  !

  !PROBLEM CONTROL PANEL
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
  DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE=CMISSSolverSolverMatrixOutput
  LINEAR_SOLVER_DARCY_OUTPUT_TYPE=CMISSSolverNoOutput
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DARCY_OUTPUT=CMISSEquationsNoOutput

  !Set time parameter
  DYNAMIC_SOLVER_DARCY_START_TIME=0.0_CMISSDP
  DYNAMIC_SOLVER_DARCY_TIME_INCREMENT=1.0e-3_CMISSDP
  DYNAMIC_SOLVER_DARCY_STOP_TIME=2_CMISSIntg * DYNAMIC_SOLVER_DARCY_TIME_INCREMENT
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

  icompartment =1_CMISSIntg
  Ncompartments=2_CMISSIntg
  !
  !================================================================================================================================
  !

! !   !Set diagnostics
! ! 
! !   DIAG_LEVEL_LIST(1)=1
! !   DIAG_LEVEL_LIST(2)=2
! !   DIAG_LEVEL_LIST(3)=3
! !   DIAG_LEVEL_LIST(4)=4
! !   DIAG_LEVEL_LIST(5)=5
! ! 
! ! !   DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_FINITE_ELEMENT_CALCULATE"
! ! !   DIAG_ROUTINE_LIST(2)="DARCY_EQUATION_PRE_SOLVE_STORE_REFERENCE_DATA"
! ! !   DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_PRE_SOLVE_GET_SOLID_DISPLACEMENT"
! ! !   DIAG_ROUTINE_LIST(2)="DARCY_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH"
! ! !   DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_PRE_SOLVE_GET_SOLID_DISPLACEMENT"
! ! !   DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS"
! !   DIAG_ROUTINE_LIST(1)="FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE"
! ! !   DIAG_ROUTINE_LIST(2)="FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR"
! ! !   DIAG_ROUTINE_LIST(3)="EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION"
! ! !   DIAG_ROUTINE_LIST(5)="DARCY_EQUATION_PRE_SOLVE_MAT_PROPERTIES"
! ! !   DIAG_ROUTINE_LIST(6)="DATA_FITTING_FINITE_ELEMENT_CALCULATE"
! ! !   DIAG_ROUTINE_LIST(7)="FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE"
! ! !   DIAG_ROUTINE_LIST(8)="FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE"
! ! !   DIAG_ROUTINE_LIST(1)="PROBLEM_SOLVER_EQUATIONS_SOLVE"
! ! !   DIAG_ROUTINE_LIST(1)="SOLVER_NEWTON_SOLVE"
! ! !   DIAG_ROUTINE_LIST(2)="SOLVER_NEWTON_LINESEARCH_SOLVE"
! ! !   DIAG_ROUTINE_LIST(1)="SOLVER_SOLUTION_UPDATE"
! ! !   DIAG_ROUTINE_LIST(1)="FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE"


!Set diagnostics

  DIAG_LEVEL_LIST(1)=1
  DIAG_LEVEL_LIST(2)=2
  DIAG_LEVEL_LIST(3)=3
  DIAG_LEVEL_LIST(4)=4
  DIAG_LEVEL_LIST(5)=5

  !DIAG_ROUTINE_LIST(1)="WRITE_IP_INFO"
!   DIAG_ROUTINE_LIST(2)="FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR"
  DIAG_ROUTINE_LIST(1)="FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE"

  !CMISSAllDiagType/CMISSInDiagType/CMISSFromDiagType
  CALL CMISSDiagnosticsSetOn(CMISSInDiagType,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISSAllTimingType/CMISSInTimingType/CMISSFromTimingType
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

!   !CMISSAllDiagType/CMISSInDiagType/CMISSFromDiagType
!    CALL CMISSDiagnosticsSetOn(CMISSInDiagType,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISSAllTimingType/CMISSInTimingType/CMISSFromTimingType
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)
  ALLOCATE (EquationsSetDarcy(Ncompartments))
  ALLOCATE (EquationsSetFieldDarcy(Ncompartments))
  ALLOCATE (MaterialsFieldDarcy(Ncompartments))
  ALLOCATE (EquationsDarcy(Ncompartments))
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
call READ_MESH('input/example-mesh-3els',MeshUserNumber,Region, Mesh,Bases,Nodes,Elements)
  !BASES
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
! !   Bases(1)=CubicBasis
! !   Bases(2)=QuadraticBasis

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

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  !Create a decomposition:
  !All mesh components (associated with G.Projection / Darcy / solid) share the same decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  !Finish the decomposition
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
CALL READ_NODES('input/example-nodes-3els',GeometricField)
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
CALL READ_NODES('input/example-nodes-3els',GeometricFieldSolid)
!   CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricFieldSolid,GeneratedMesh,Err)
!   !Set the scaling to use
!   CALL CMISSFieldScalingTypeSet(GeometricFieldSolid,CMISSFieldNoScaling,Err)

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
  DO icompartment = 1,Ncompartments
    EquationsSetFieldUserNumberDarcy = 100_CMISSIntg+icompartment
    EquationsSetUserNumberDarcy = 200_CMISSIntg+icompartment
  !Create the equations set for ALE Darcy
    CALL CMISSFieldTypeInitialise(EquationsSetFieldDarcy(icompartment),Err)
    CALL CMISSEquationsSetTypeInitialise(EquationsSetDarcy(icompartment),Err)
    CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberDarcy,Region,GeometricField,CMISSEquationsSetFluidMechanicsClass, & 
      & CMISSEquationsSetDarcyEquationType,CMISSEquationsSetIncompressibleElastMultiCompDarcySubtype,&
      & EquationsSetFieldUserNumberDarcy,EquationsSetFieldDarcy(icompartment),EquationsSetDarcy(icompartment),Err)
    !Finish creating the equations set
    CALL CMISSEquationsSetCreateFinish(EquationsSetDarcy(icompartment),Err)
    !Set the values for the equations set field to be the current compartment number (1 - N), and the total number of compartments (N)
    CALL CMISSFieldParameterSetUpdateConstant(EquationsSetFieldDarcy(icompartment),CMISSFieldUVariableType, &
      & CMISSFieldValuesSetType,1,icompartment,Err)
    CALL CMISSFieldParameterSetUpdateConstant(EquationsSetFieldDarcy(icompartment),CMISSFieldUVariableType, &
      & CMISSFieldValuesSetType,2,Ncompartments,Err)
  ENDDO

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create the equations_set
  CALL CMISSFieldTypeInitialise(EquationsSetFieldSolid,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSetSolid,Err)
  CALL CMISSEquationsSetCreateStart(EquationSetSolidUserNumber,Region,FibreFieldSolid,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetIncompressibleElastMultiCompDarcySubtype,&
    & EquationsSetFieldSolidUserNumber,EquationsSetFieldSolid,EquationsSetSolid,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSetSolid,Err)
  !Set the values for the equations set field to be the current compartment number (O for the finite elasticity equations_set), and the total number of compartments (N)
  !Need to store number of compartments, as finite elasticity uses this to calculate the total mass increase for the constiutive law
  CALL CMISSFieldParameterSetUpdateConstant(EquationsSetFieldSolid,CMISSFieldUVariableType, &
     & CMISSFieldValuesSetType,1,0_CMISSIntg,Err)
  CALL CMISSFieldParameterSetUpdateConstant(EquationsSetFieldSolid,CMISSFieldUVariableType, &
     & CMISSFieldValuesSetType,2,Ncompartments,Err)


  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------
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


  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create a dependent field with two variables and four components
  CALL CMISSFieldTypeInitialise(DependentFieldSolid,Err)
  !
  CALL CMISSFieldCreateStart(FieldDependentSolidUserNumber,Region,DependentFieldSolid,Err)
  !
  CALL CMISSFieldTypeSet(DependentFieldSolid,CMISSFieldGeneralType,Err)
  CALL CMISSFieldMeshDecompositionSet(DependentFieldSolid,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(DependentFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSFieldDependentTypeSet(DependentFieldSolid,CMISSFieldDependentType,Err)
  !Create 2N+2 number of variables - 2 for solid, 2N for N Darcy compartments
  FieldDependentSolidNumberOfVariables=2*Ncompartments+2
  CALL CMISSFieldNumberOfVariablesSet(DependentFieldSolid,FieldDependentSolidNumberOfVariables,Err)
  !create two variables for each compartment
  ALLOCATE(VariableTypes(2*Ncompartments+2))
  DO num_var=1,Ncompartments+1
     VariableTypes(2*num_var-1)=CMISSFieldUVariableType+(CMISSFieldNumberOfVariableSubtypes*(num_var-1))
     VariableTypes(2*num_var)=CMISSFieldDelUDelNVariableType+(CMISSFieldNumberOfVariableSubtypes*(num_var-1))
  ENDDO
  CALL CMISSFieldVariableTypesSet(DependentFieldSolid,VariableTypes,Err) 
!   CALL CMISSFieldVariableTypesSet(DependentFieldSolid,(/CMISSFieldUVariableType, &
!     & CMISSFieldDelUDelNVariableType,CMISSFieldVVariableType,CMISSFieldDelVDelNVariableType/),Err)
    CALL CMISSFieldDimensionSet(DependentFieldSolid,CMISSFieldUVariableType, &
       & CMISSFieldVectorDimensionType,Err)
    CALL CMISSFieldDimensionSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType, &
       & CMISSFieldVectorDimensionType,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldUVariableType,FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,FieldDependentSolidNumberOfComponents,Err)
  DO icompartment=3,2*Ncompartments+2
    CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,VariableTypes(icompartment),FieldDependentFluidNumberOfComponents,Err)
  ENDDO
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,1,CMISSFieldNodeBasedInterpolation,Err)
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,2,CMISSFieldNodeBasedInterpolation,Err)
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,3,CMISSFieldNodeBasedInterpolation,Err)

  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,1,SolidDisplMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,2,SolidDisplMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,3,SolidDisplMeshComponentNumber,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,4,CMISSFieldNodeBasedInterpolation,Err)
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,4,CMISSFieldElementBasedInterpolation,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,4,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,4,SolidLagrMultMeshComponentNumber,Err)
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,1, &
!     & CMISSFieldNodeBasedInterpolation,Err)
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,2, &
!     & CMISSFieldNodeBasedInterpolation,Err)
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,3, &
!     & CMISSFieldNodeBasedInterpolation,Err)
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
  !loop over the number of compartments
  DO icompartment=3,2*Ncompartments+2
!     CALL CMISSFieldDimensionSet(DependentFieldSolid,VariableTypes(icompartment), &
!        & CMISSFieldVectorDimensionType,Err)
    !CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,VariableTypes(icompartment),FieldDependentFluidNumberOfComponents,Err)
    DO componentnum=1,FieldDependentFluidNumberOfComponents-1
    !set dimension type
!     CALL CMISSFieldDimensionSet(DependentField,VariableTypes(icompartment), &
!        & CMISSFieldScalarDimensionType,Err)
      CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,VariableTypes(icompartment),componentnum, &
       & CMISSFieldNodeBasedInterpolation,Err)
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,VariableTypes(icompartment),componentnum, & 
         & DarcyVelMeshComponentNumber,Err)
    ENDDO
      CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,VariableTypes(icompartment), &
       & FieldDependentFluidNumberOfComponents, &
       & CMISSFieldNodeBasedInterpolation,Err)
!     CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,VariableTypes(icompartment), &
!       & FieldDependentFluidNumberOfComponents,MESH_COMPONENT_NUMBER_PRESSURE,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,VariableTypes(icompartment), &
      & FieldDependentFluidNumberOfComponents,DarcyMassIncreaseMeshComponentNumber,Err)
    
  ENDDO

!   CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldVVariableType,FieldDependentFluidNumberOfComponents,Err)
!   CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,FieldDependentFluidNumberOfComponents,Err)
!   !For this equation type, MESH_COMPONENT_NUMBER_PRESSURE is actually the mass increase component as the pressure is taken from the solid equations
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,1,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,2,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,3,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,4,MESH_COMPONENT_NUMBER_PRESSURE,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,2,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,3,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,4,MESH_COMPONENT_NUMBER_PRESSURE,Err)

  !
  CALL CMISSFieldCreateFinish(DependentFieldSolid,Err)
  !
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetSolid,FieldDependentSolidUserNumber,DependentFieldSolid,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetSolid,Err)

!   !Initialise dependent field (solid displacement and pressure)
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS !+1
!     CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
!       & COMPONENT_NUMBER,INITIAL_FIELD_SOLID(COMPONENT_NUMBER),Err)
!   ENDDO

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  DO icompartment = 1,Ncompartments
    CALL CMISSEquationsSetDependentCreateStart(EquationsSetDarcy(icompartment),FieldDependentSolidUserNumber,&
      & DependentFieldSolid,Err)
    CALL CMISSEquationsSetDependentCreateFinish(EquationsSetDarcy(icompartment),Err)
  ENDDO


  !Create the equations set dependent field variables for ALE Darcy
!   CALL CMISSFieldTypeInitialise(DependentFieldDarcy,Err)
!   CALL CMISSEquationsSetDependentCreateStart(EquationsSetDarcy,DependentFieldUserNumberDarcy, & ! ??? UserNumber ???
!   CALL CMISSEquationsSetDependentCreateStart(EquationsSetDarcy,FieldDependentSolidUserNumber, & ! ??? UserNumber ???
!     & DependentFieldSolid,Err)
! !   !Set the mesh component to be used by the field components.
! !   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,CMISSFieldUVariableType,COMPONENT_NUMBER, &
! !       & MESH_COMPONENT_NUMBER_VELOCITY,Err)
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, &
! !       & MESH_COMPONENT_NUMBER_VELOCITY,Err)
! !   ENDDO
! !   COMPONENT_NUMBER=NUMBER_OF_DIMENSIONS+1
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,CMISSFieldUVariableType,COMPONENT_NUMBER, &
! !       & MESH_COMPONENT_NUMBER_PRESSURE,Err)
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, &
! !       & MESH_COMPONENT_NUMBER_PRESSURE,Err)
!   !Finish the equations set dependent field variables
!   CALL CMISSEquationsSetDependentCreateFinish(EquationsSetDarcy,Err)

  !Initialise dependent field (velocity components,pressure,mass increase)
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS+1
    CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldVVariableType,CMISSFieldValuesSetType, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
    CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldU1VariableType,CMISSFieldValuesSetType, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
  ENDDO
  !
  !================================================================================================================================
  !
  ALLOCATE(CouplingCoeffs(Ncompartments,Ncompartments))
  IF(Ncompartments==2)THEN
    CouplingCoeffs(1,1)=0.0E-01_CMISSDP
!     CouplingCoeffs(1,2)=-1.0E-04_CMISSDP
!     CouplingCoeffs(2,1)=-1.0E-04_CMISSDP
    CouplingCoeffs(1,2)=0.0E-01_CMISSDP
    CouplingCoeffs(2,1)=0.0E-01_CMISSDP
    CouplingCoeffs(2,2)=0.0E-01_CMISSDP
  ELSE IF(Ncompartments==3)THEN
    CouplingCoeffs(1,1)=1.0E-02_CMISSDP
    CouplingCoeffs(1,2)=1.0E-02_CMISSDP
    CouplingCoeffs(1,3)=0.0E-02_CMISSDP
    CouplingCoeffs(2,1)=1.0E-02_CMISSDP
    CouplingCoeffs(2,2)=2.0E-02_CMISSDP
    CouplingCoeffs(2,3)=1.0E-02_CMISSDP
    CouplingCoeffs(3,1)=0.0E-02_CMISSDP
    CouplingCoeffs(3,2)=1.0E-02_CMISSDP
    CouplingCoeffs(3,3)=1.0E-02_CMISSDP
  ELSE IF(Ncompartments==4)THEN
    CouplingCoeffs(1,1)=0.0E-02_CMISSDP
    CouplingCoeffs(1,2)=0.0E-02_CMISSDP
    CouplingCoeffs(1,3)=0.0E-02_CMISSDP
    CouplingCoeffs(1,4)=0.0E-02_CMISSDP
    CouplingCoeffs(2,1)=0.0E-02_CMISSDP
    CouplingCoeffs(2,2)=0.0E-02_CMISSDP
    CouplingCoeffs(2,3)=0.0E-02_CMISSDP
    CouplingCoeffs(2,4)=0.0E-02_CMISSDP
    CouplingCoeffs(3,1)=0.0E-02_CMISSDP
    CouplingCoeffs(3,2)=0.0E-02_CMISSDP
    CouplingCoeffs(3,3)=0.0E-02_CMISSDP
    CouplingCoeffs(3,4)=0.0E-02_CMISSDP
    CouplingCoeffs(4,1)=0.0E-02_CMISSDP
    CouplingCoeffs(4,2)=0.0E-02_CMISSDP
    CouplingCoeffs(4,3)=0.0E-02_CMISSDP
    CouplingCoeffs(4,4)=0.0E-02_CMISSDP
  ELSE IF(Ncompartments==5)THEN
    CouplingCoeffs(1,1)=0.0E-02_CMISSDP
    CouplingCoeffs(1,2)=0.0E-02_CMISSDP
    CouplingCoeffs(1,3)=0.0E-02_CMISSDP
    CouplingCoeffs(1,4)=0.0E-02_CMISSDP
    CouplingCoeffs(1,5)=0.0E-02_CMISSDP
    CouplingCoeffs(2,1)=0.0E-02_CMISSDP
    CouplingCoeffs(2,2)=0.0E-02_CMISSDP
    CouplingCoeffs(2,3)=0.0E-02_CMISSDP
    CouplingCoeffs(2,4)=0.0E-02_CMISSDP
    CouplingCoeffs(2,5)=0.0E-02_CMISSDP
    CouplingCoeffs(3,1)=0.0E-02_CMISSDP
    CouplingCoeffs(3,2)=0.0E-02_CMISSDP
    CouplingCoeffs(3,3)=0.0E-02_CMISSDP
    CouplingCoeffs(3,4)=0.0E-02_CMISSDP
    CouplingCoeffs(3,5)=0.0E-02_CMISSDP
    CouplingCoeffs(4,1)=0.0E-02_CMISSDP
    CouplingCoeffs(4,2)=0.0E-02_CMISSDP
    CouplingCoeffs(4,3)=0.0E-02_CMISSDP
    CouplingCoeffs(4,4)=0.0E-02_CMISSDP
    CouplingCoeffs(4,5)=0.0E-02_CMISSDP
    CouplingCoeffs(5,1)=0.0E-02_CMISSDP
    CouplingCoeffs(5,2)=0.0E-02_CMISSDP
    CouplingCoeffs(5,3)=0.0E-02_CMISSDP
    CouplingCoeffs(5,4)=0.0E-02_CMISSDP
    CouplingCoeffs(5,5)=0.0E-02_CMISSDP
  ELSE
    write(*,*) "Can't initialise coupling coefficients array."
  ENDIF
  !Define the material parameters for each compartments' constitutive law (for determining pressure)
  Nparams=3
  ALLOCATE(ConstitutiveParams(Ncompartments,Nparams))
  IF(Ncompartments==2)THEN
    ConstitutiveParams(1,1)=10.0E-01_CMISSDP
    ConstitutiveParams(1,2)=10.0E-01_CMISSDP
    ConstitutiveParams(1,3)=10.0E-01_CMISSDP
    ConstitutiveParams(2,1)=10.0E-01_CMISSDP
    ConstitutiveParams(2,2)=10.0E-01_CMISSDP
    ConstitutiveParams(2,3)=10.0E-01_CMISSDP
  ELSE IF(Ncompartments==3)THEN
    ConstitutiveParams(1,1)=1.0E-02_CMISSDP
    ConstitutiveParams(1,2)=1.0E-02_CMISSDP
    ConstitutiveParams(1,3)=0.0E-02_CMISSDP
    ConstitutiveParams(2,1)=1.0E-02_CMISSDP
    ConstitutiveParams(2,2)=2.0E-02_CMISSDP
    ConstitutiveParams(2,3)=1.0E-02_CMISSDP
    ConstitutiveParams(3,1)=0.0E-02_CMISSDP
    ConstitutiveParams(3,2)=1.0E-02_CMISSDP
    ConstitutiveParams(3,3)=1.0E-02_CMISSDP
  ELSE IF(Ncompartments==4)THEN
    ConstitutiveParams(1,1)=0.0E-02_CMISSDP
    ConstitutiveParams(1,2)=0.0E-02_CMISSDP
    ConstitutiveParams(1,3)=0.0E-02_CMISSDP
    ConstitutiveParams(2,1)=0.0E-02_CMISSDP
    ConstitutiveParams(2,2)=0.0E-02_CMISSDP
    ConstitutiveParams(2,3)=0.0E-02_CMISSDP
    ConstitutiveParams(3,1)=0.0E-02_CMISSDP
    ConstitutiveParams(3,2)=0.0E-02_CMISSDP
    ConstitutiveParams(3,3)=0.0E-02_CMISSDP
    ConstitutiveParams(4,1)=0.0E-02_CMISSDP
    ConstitutiveParams(4,2)=0.0E-02_CMISSDP
    ConstitutiveParams(4,3)=0.0E-02_CMISSDP
  ELSE IF(Ncompartments==5)THEN
    ConstitutiveParams(1,1)=0.0E-02_CMISSDP
    ConstitutiveParams(1,2)=0.0E-02_CMISSDP
    ConstitutiveParams(1,3)=0.0E-02_CMISSDP
    ConstitutiveParams(2,1)=0.0E-02_CMISSDP
    ConstitutiveParams(2,2)=0.0E-02_CMISSDP
    ConstitutiveParams(2,3)=0.0E-02_CMISSDP
    ConstitutiveParams(3,1)=0.0E-02_CMISSDP
    ConstitutiveParams(3,2)=0.0E-02_CMISSDP
    ConstitutiveParams(3,3)=0.0E-02_CMISSDP
    ConstitutiveParams(4,1)=0.0E-02_CMISSDP
    ConstitutiveParams(4,2)=0.0E-02_CMISSDP
    ConstitutiveParams(4,3)=0.0E-02_CMISSDP
    ConstitutiveParams(5,1)=0.0E-02_CMISSDP
    ConstitutiveParams(5,2)=0.0E-02_CMISSDP
    ConstitutiveParams(5,3)=0.0E-02_CMISSDP
  ELSE
    write(*,*) "Can't initialise constitutive parameters array."
  ENDIF
  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS
  !Auto-created field contains a U variable type to store the diffusion coefficient(s)
  !It also contains a V variable type to store the coupling coefficients 
! call READ_FIELD('input/example-field-darcy',MaterialsFieldUserNumberDarcy,Region,GeometricField,MaterialsFieldDarcy)
!   DO icompartment = 1,Ncompartments
!     MaterialsFieldUserNumberDarcy = 400+icompartment
!     CALL CMISSFieldTypeInitialise(MaterialsFieldDarcy(icompartment),Err)
!     CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetDarcy(icompartment),MaterialsFieldUserNumberDarcy,&
!          & MaterialsFieldDarcy(icompartment),Err)
!     CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetDarcy(icompartment),Err)
!   END DO
!   DO icompartment = 1,Ncompartments
!     CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMISSFieldUVariableType, &
!       & CMISSFieldValuesSetType, &
!       & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
!     CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMISSFieldUVariableType, &
!       & CMISSFieldValuesSetType, &
!       & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)
!   END DO
!   DO icompartment = 1, Ncompartments
!     DO COMPONENT_NUMBER=1, Ncompartments
!       CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMISSFieldVVariableType, &
!          & CMISSFieldValuesSetType,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
! !         CALL CMISSFieldParameterSetUpdateConstant(MaterialsFieldDarcy(icompartment),CMISSFieldVVariableType, &
! !           & CMISSFieldValuesSetType,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
!     END DO
!   END DO
!   DO icompartment = 1, Ncompartments
!     DO COMPONENT_NUMBER=1,Nparams
!       CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMISSFieldU1VariableType, &
!          & CMISSFieldValuesSetType,COMPONENT_NUMBER,ConstitutiveParams(icompartment,COMPONENT_NUMBER),Err)
! !         CALL CMISSFieldParameterSetUpdateConstant(MaterialsFieldDarcy(icompartment),CMISSFieldVVariableType, &
! !           & CMISSFieldValuesSetType,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
!     END DO
!   END DO
  DO icompartment = 1,Ncompartments
     MaterialsFieldUserNumberDarcy = 400+icompartment
call READ_FIELD('input/new-format-field-darcy',MaterialsFieldUserNumberDarcy,Region,GeometricField, &
   & MaterialsFieldDarcy(icompartment))
    CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetDarcy(icompartment),MaterialsFieldUserNumberDarcy,&
         & MaterialsFieldDarcy(icompartment),Err)
    CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetDarcy(icompartment),Err)
  ENDDO

  !Create the equations set materials field variables for ALE Darcy
!   CALL CMISSFieldTypeInitialise(MaterialsFieldDarcy,Err)
!   CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetDarcy,MaterialsFieldUserNumberDarcy, &
!     & MaterialsFieldDarcy,Err)
!   !Finish the equations set materials field variables
!   CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetDarcy,Err)
!   CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
!     & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
!   CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
!     & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)
  !Create the equations set materials field variables for deformation-dependent material properties

  !
  !================================================================================================================================
  !

  !EQUATIONS

 
  DO icompartment=1,Ncompartments
    !Create the equations set equations
    CALL CMISSEquationsTypeInitialise(EquationsDarcy(icompartment),Err)
    CALL CMISSEquationsSetEquationsCreateStart(EquationsSetDarcy(icompartment),EquationsDarcy(icompartment),Err)
    !Set the equations matrices sparsity type
    CALL CMISSEquationsSparsityTypeSet(EquationsDarcy(icompartment),CMISSEquationsSparseMatrices,Err)
  !   !Set the equations lumping type
  !   CALL CMISSEquationsLumpingTypeSet(EquationsDarcy,CMISSEquationsUnlumpedMatrices,Err)
    !Set the equations set output
    CALL CMISSEquationsOutputTypeSet(EquationsDarcy(icompartment),EQUATIONS_DARCY_OUTPUT,Err)
  !Finish the equations set equations
    CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetDarcy(icompartment),Err)
  ENDDO
  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(EquationsSolid,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetSolid,EquationsSolid,Err)
  CALL CMISSEquationsSparsityTypeSet(EquationsSolid,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(EquationsSolid,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

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

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(SolverSolid,Err)
  CALL CMISSSolverTypeInitialise(DynamicSolverDarcy,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverDarcy,Err)

  CALL CMISSProblemSolversCreateStart(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Get the finite elasticity solver
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

!   CALL CMISSSolverNewtonLinearSolverGet(SolverSolid,LinearSolverSolid,Err)
!   CALL CMISSSolverLinearTypeSet(LinearSolverSolid,CMISSSolverLinearDirectSolveType,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  !Get the Darcy solver
  CALL CMISSProblemSolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMISSControlLoopNode/), &
    & SolverDarcyIndex,DynamicSolverDarcy,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(DynamicSolverDarcy,DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE,Err)
  !Set theta
  CALL CMISSSolverDynamicThetaSet(DynamicSolverDarcy,DYNAMIC_SOLVER_DARCY_THETA,Err)
  !CALL CMISSSolverDynamicDynamicSet(DynamicSolverDarcy,.TRUE.,Err)
  !Get the dynamic linear solver
  CALL CMISSSolverDynamicLinearSolverGet(DynamicSolverDarcy,LinearSolverDarcy,Err)
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
  DO icompartment=1,Ncompartments
    CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsDarcy,EquationsSetDarcy(icompartment),EquationsSetIndex,Err)
  ENDDO
  !
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)
!   ENDDO




  !Prescribe boundary conditions (absolute nodal parameters)
  !Solid is computed in absolute position, rather than displacement. Thus BCs for absolute position
   ! ASSIGN BOUNDARY CONDITIONS - SOLID (absolute nodal parameters)
  !Solid is computed in absolute position, rather than displacement. Thus BCs for absolute position
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsSolid,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsSolid,BoundaryConditionsSolid,Err)

  !Get surfaces
!   CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISSGeneratedMeshRegularFrontSurface, &
!     & Face1Nodes,FaceXi(1),Err)
!   CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISSGeneratedMeshRegularBackSurface, &
!     & Face2Nodes,FaceXi(2),Err)
!   CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISSGeneratedMeshRegularRightSurface, &
!     & Face3Nodes,FaceXi(3),Err)
!   CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISSGeneratedMeshRegularLeftSurface, &
!     & Face4Nodes,FaceXi(4),Err)
!   CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISSGeneratedMeshRegularTopSurface, &
!     & Face5Nodes,FaceXi(5),Err)
!   CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISSGeneratedMeshRegularBottomSurface, &
!     & Face6Nodes,FaceXi(6),Err)

! ! write(*,'("Face1Nodes=[",100I3,"]")') Face1Nodes
! ! write(*,'("Face2Nodes=[",100I3,"]")') Face2Nodes
! ! write(*,'("Face3Nodes=[",100I3,"]")') Face3Nodes
! ! write(*,'("Face4Nodes=[",100I3,"]")') Face4Nodes
! ! write(*,'("Face5Nodes=[",100I3,"]")') Face5Nodes
! ! write(*,'("Face6Nodes=[",100I3,"]")') Face6Nodes
! ! write(*,'("FaceXi=[",100I3,"]")') FaceXi
!
allocate(Face1Nodes(21),Face2Nodes(21),Face3Nodes(21),Face4Nodes(21),Face5Nodes(9),Face6Nodes(9))
Face1Nodes=[  1, 17,  2, 22, 23, 24,  5, 31,  6, 36, 37, 38,  9, 45, 10, 50, 51, 52, 13, 59, 14]
Face2Nodes=[  3, 21,  4, 28, 29, 30,  7, 35,  8, 42, 43, 44, 11, 49, 12, 56, 57, 58, 15, 63, 16]
Face3Nodes=[  2, 20,  4, 24, 27, 30,  6, 34,  8, 38, 41, 44, 10, 48, 12, 52, 55, 58, 14, 62, 16]
Face4Nodes=[  1, 18,  3, 22, 25, 28,  5, 32,  7, 36, 39, 42,  9, 46, 11, 50, 53, 56, 13, 60, 15]
Face5Nodes=[ 13, 59, 14, 60, 61, 62, 15, 63, 16]
Face6Nodes=[  1, 17,  2, 18, 19, 20,  3, 21,  4]
FaceXi=[ -2,  2,  1, -1,  3, -3]

  ! Fix the bottom in z direction
  DO NN=1,SIZE(Face6Nodes,1)
    NODE=Face6Nodes(NN)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetGetNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,NODE,3,ZCoord,Err)
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldUVariableType,1,NODE,3, &
        & CMISSBoundaryConditionFixed,ZCoord,Err)
      WRITE(*,*) "FIXING NODE",NODE,"AT BOTTOM IN Z DIRECTION"
    ENDIF
  ENDDO

  ! Fix the top in z direction
  DO NN=1,SIZE(Face5Nodes,1)
    NODE=Face5Nodes(NN)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetGetNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,NODE,3,ZCoord,Err)
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldUVariableType,1,NODE,3, &
        & CMISSBoundaryConditionFixed,ZCoord,Err)
      WRITE(*,*) "FIXING NODE",NODE,"AT TOP IN Z DIRECTION"
    ENDIF
  ENDDO

  !Fix more nodes at the bottom to stop free body motion
  X_FIXED=.FALSE.
  Y_FIXED=.FALSE.
  DO NN=1,SIZE(Face6Nodes,1)
    NODE=Face6Nodes(NN)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetGetNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,NODE,1,XCoord,Err)
      CALL CMISSFieldParameterSetGetNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,NODE,2,YCoord,Err)

      !Fix Origin displacement in x and y (z already fixed)
      IF(ABS(XCoord)<1.0E-6_CMISSDP.AND.ABS(YCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldUVariableType,1,NODE,1, &
          & CMISSBoundaryConditionFixed,XCoord,Err)
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldUVariableType,1,NODE,2, &
          & CMISSBoundaryConditionFixed,YCoord,Err)
        WRITE(*,*) "FIXING ORIGIN NODE",NODE,"IN X AND Y DIRECTION"
        X_FIXED=.TRUE.
        Y_FIXED=.TRUE.
      ENDIF

      !Fix nodal displacements at (X_DIM,0) in y
      IF(ABS(XCoord - X_DIM)<1.0E-6_CMISSDP .AND. ABS(YCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldUVariableType,1,NODE,2, &
          & CMISSBoundaryConditionFixed,YCoord,Err)
        WRITE(*,*) "FIXING NODES",NODE,"AT (X_DIM,0) IN Y DIRECTION"
        Y_FIXED=.TRUE.
      ENDIF

      !Fix nodal displacements at (0,Y_DIM) in x
      IF(ABS(XCoord)<1.0E-6_CMISSDP .AND. ABS(YCoord - Y_DIM)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISSFieldUVariableType,1,NODE,1, &
          & CMISSBoundaryConditionFixed,XCoord,Err)
        WRITE(*,*) "FIXING NODES",NODE,"AT (0,Y_DIM) IN X DIRECTION"
        X_FIXED=.TRUE.
      ENDIF

    ENDIF
  ENDDO
!   CALL MPI_REDUCE(X_FIXED,X_OKAY,1,MPI_LOGICAL,MPI_LOR,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_REDUCE(Y_FIXED,Y_OKAY,1,MPI_LOGICAL,MPI_LOR,0,MPI_COMM_WORLD,MPI_IERROR)
!   IF(ComputationalNodeNumber==0) THEN
!     IF(.NOT.(X_OKAY.AND.Y_OKAY)) THEN
!       WRITE(*,*) "Free body motion could not be prevented!"
!       CALL CMISSFinalise(Err)
!       STOP
!     ENDIF
!   ENDIF

  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsSolid,Err)


  !
  !------------------------------------
  ! ASSIGN BOUNDARY CONDITIONS - FLUID
    !Get surfaces
!     CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISSGeneratedMeshRegularFrontSurface, &
!       & Face7Nodes,FaceXi(1),Err)
!     CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISSGeneratedMeshRegularBackSurface, &
!       & Face8Nodes,FaceXi(2),Err)
!     CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISSGeneratedMeshRegularRightSurface, &
!       & Face9Nodes,FaceXi(3),Err)
!     CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISSGeneratedMeshRegularLeftSurface, &
!       & Face10Nodes,FaceXi(4),Err)
!     CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISSGeneratedMeshRegularTopSurface, &
!       & Face11Nodes,FaceXi(5),Err)
!     CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISSGeneratedMeshRegularBottomSurface, &
!       & Face12Nodes,FaceXi(6),Err)


allocate(Face7Nodes(8),Face8Nodes(8),Face9Nodes(8),Face10Nodes(8),Face11Nodes(4),Face12Nodes(4))
Face7Nodes=[  1,  2,  5,  6,  9, 10, 13, 14]
Face8Nodes=[  3,  4,  7,  8, 11, 12, 15, 16]
Face9Nodes=[  2,  4,  6,  8, 10, 12, 14, 16]
Face10Nodes=[  1,  3,  5,  7,  9, 11, 13, 15]
Face11Nodes=[ 13, 14, 15, 16]
Face12Nodes=[  1,  2,  3,  4]
FaceXi=[ -2,  2,  1, -1,  3, -3]

  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsDarcy,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)
  DO icompartment=1,Ncompartments

    IF(icompartment==1) THEN
    ! At the top impose Darcy velocity in z direction
    DO NN=1,SIZE(Face11Nodes,1)
      NODE=Face11Nodes(NN)
! !     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = -1.0_CMISSDP
        COMPONENT_NUMBER = 3
        write(*,*)'Marker 0'
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,NODE, &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        WRITE(*,*) "SPECIFIED INFLOW AT NODE",NODE,"IN Z DIRECTION"

! !       CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,NODE,1,XCoord,Err)
! !       CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,NODE,2,YCoord,Err)
! !       CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,NODE,3,ZCoord,Err)
! !       WRITE(*,*) "XCoord, YCoord, ZCoord = ",XCoord, YCoord, ZCoord
! !     ENDIF
    ENDDO

    !All other faces are impermeable
    DO NN=1,SIZE(Face7Nodes,1)
      NODE=Face7Nodes(NN)
! !     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 1
        write(*,*)'Marker 1'
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,NODE, &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN X DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face8Nodes,1)
      NODE=Face8Nodes(NN)
! !     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 1
        write(*,*)'Marker 2'
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,NODE, &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN X DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face9Nodes,1)
      NODE=Face9Nodes(NN)
! !     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 2
        write(*,*)'Marker 3'
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,NODE, &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Y DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face10Nodes,1)
      NODE=Face10Nodes(NN)
! !     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 2
        write(*,*)'Marker 4'
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,NODE, &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Y DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face12Nodes,1)
      NODE=Face12Nodes(NN)
! !     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 3
        write(*,*)'Marker 5'
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,NODE, &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Z DIRECTION"
! !     ENDIF
    ENDDO
    ENDIF
    !Compartment TWO boundary conditions!
    ! At the top impose Darcy velocity in z direction
    IF(icompartment==2) THEN
    DO NN=1,SIZE(Face11Nodes,1)
      NODE=Face11Nodes(NN)
! !     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = -1.0_CMISSDP
        COMPONENT_NUMBER = 3
        write(*,*)'Marker 0'
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldU1VariableType,1,NODE, &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        WRITE(*,*) "SPECIFIED INFLOW AT NODE",NODE,"IN Z DIRECTION"

! !       CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,NODE,1,XCoord,Err)
! !       CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,NODE,2,YCoord,Err)
! !       CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,NODE,3,ZCoord,Err)
! !       WRITE(*,*) "XCoord, YCoord, ZCoord = ",XCoord, YCoord, ZCoord
! !     ENDIF
    ENDDO
!     DO NN=1,SIZE(Face11Nodes,1)
!       NODE=Face11Nodes(NN)
! ! !     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! ! !     IF(NodeDomain==ComputationalNodeNumber) THEN
!         VALUE = 0.0_CMISSDP
!         COMPONENT_NUMBER = 3
!         write(*,*)'Marker 0'
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldU1VariableType,1,NODE,COMPONENT_NUMBER,&
!           & CMISSBoundaryConditionFixed,VALUE,Err)
!         WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Z DIRECTION"
!
! ! !       CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,NODE,1,XCoord,Err)
! ! !       CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,NODE,2,YCoord,Err)
! ! !       CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,NODE,3,ZCoord,Err)
! ! !       WRITE(*,*) "XCoord, YCoord, ZCoord = ",XCoord, YCoord, ZCoord
! ! !     ENDIF
!     ENDDO

    !All other faces are impermeable
    DO NN=1,SIZE(Face7Nodes,1)
      NODE=Face7Nodes(NN)
! !     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 1
        write(*,*)'Marker 1'
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldU1VariableType,1,NODE, &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN X DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face8Nodes,1)
      NODE=Face8Nodes(NN)
! !     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 1
        write(*,*)'Marker 2'
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldU1VariableType,1,NODE, &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN X DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face9Nodes,1)
      NODE=Face9Nodes(NN)
! !     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 2
        write(*,*)'Marker 3'
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldU1VariableType,1,NODE, &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Y DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face10Nodes,1)
      NODE=Face10Nodes(NN)
! !     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 2
        write(*,*)'Marker 4'
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldU1VariableType,1,NODE, &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Y DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face12Nodes,1)
      NODE=Face12Nodes(NN)
! !     CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 3
        write(*,*)'Marker 5'
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldU1VariableType,1,NODE, &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Z DIRECTION"
! !     ENDIF
    ENDDO
    ENDIF

  !BOUNDARY CONDITIONS
  !Start the creation of the equations set boundary conditions for Darcy
!   DO icompartment=1,Ncompartments  
!     CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsDarcy,DependentFieldSolid,Err)
!     CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSetDarcy(icompartment),BoundaryConditionsDarcy,DependentFieldSolid,Err)
!     CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSetDarcy(icompartment),Err)
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

  WRITE(*,*)'dummy_counter = ',dummy_counter

  EXPORT_FIELD_IO=.FALSE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"FiniteElasticityMultiCompDarcy","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"FiniteElasticityMultiCompDarcy","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF


  !Finialise CMISS
!   CALL CMISSFinalise(Err)

  IF (ALLOCATED(EquationsSetFieldDarcy)) DEALLOCATE(EquationsSetFieldDarcy)
  IF (ALLOCATED(EquationsSetDarcy)) DEALLOCATE(EquationsSetDarcy)
  IF (ALLOCATED(MaterialsFieldDarcy)) DEALLOCATE(MaterialsFieldDarcy)
  IF (ALLOCATED(BoundaryConditionsDarcy)) DEALLOCATE(BoundaryConditionsDarcy)
  IF (ALLOCATED(EquationsDarcy)) DEALLOCATE(EquationsDarcy)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM FINITEELASTICITYMULTICOMPDARCYIOEXAMPLE
