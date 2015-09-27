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

  !CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

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

!   BASIS_XI_INTERPOLATION_SOLID=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  BASIS_XI_INTERPOLATION_SOLID=CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
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
  DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE=CMISS_SOLVER_MATRIX_OUTPUT
  LINEAR_SOLVER_DARCY_OUTPUT_TYPE=CMISS_SOLVER_NO_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DARCY_OUTPUT=CMISS_EQUATIONS_NO_OUTPUT

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

  !CMISS_ALL_DIAG_TYPE/CMISS_IN_DIAG_TYPE/CMISS_FROM_DIAG_TYPE
  CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISS_ALL_TIMING_TYPE/CMISS_IN_TIMING_TYPE/CMISS_FROM_TIMING_TYPE
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

!   !CMISS_ALL_DIAG_TYPE/CMISS_IN_DIAG_TYPE/CMISS_FROM_DIAG_TYPE
!    CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISS_ALL_TIMING_TYPE/CMISS_IN_TIMING_TYPE/CMISS_FROM_TIMING_TYPE
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

  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !REGION
  !For a volume-coupled problem, solid and fluid are based in the same region

  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !
  !================================================================================================================================
  !
call READ_MESH('input/example-mesh-3els',MeshUserNumber,Region, Mesh,Bases,Nodes,Elements)
  !BASES
!   !Define basis functions
!   CALL CMISSBasis_Initialise(LinearBasis,Err)
!   CALL CMISSBasis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
!   CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearBasis, &
!     & (/CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME/),Err)
!   !CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err)
!   CALL CMISSBasis_CreateFinish(LinearBasis,Err)
! 
!   CALL CMISSBasis_Initialise(QuadraticBasis,Err)
!   CALL CMISSBasis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
!   CALL CMISSBasis_InterpolationXiSet(QuadraticBasis,(/CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
!     & CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION/),Err)
!   CALL CMISSBasis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
!     & (/CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME/),Err)
!   !CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err)
!   CALL CMISSBasis_CreateFinish(QuadraticBasis,Err)
! 
!   CALL CMISSBasis_Initialise(CubicBasis,Err)
!   CALL CMISSBasis_CreateStart(CubicBasisUserNumber,CubicBasis,Err)
!   CALL CMISSBasis_InterpolationXiSet(CubicBasis,(/CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
!     & CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION/),Err)
!   CALL CMISSBasis_QuadratureNumberOfGaussXiSet(CubicBasis, &
!     & (/CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME/),Err)
!   !CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(CubicBasis,.true.,Err) !Enable 3D interpolation on faces
!   CALL CMISSBasis_CreateFinish(CubicBasis,Err)
! 
!   !LinearBasis/QuadraticBasis/CubicBasis
!   Bases(1)=QuadraticBasis
!   Bases(2)=LinearBasis
! !   Bases(1)=CubicBasis
! !   Bases(2)=QuadraticBasis

!   !Start the creation of a generated mesh in the region
!   CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
!   CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
!   CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
!   CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,Bases,Err)
!   IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
!     CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/X_DIM,Y_DIM/),Err)
!     CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/),Err)
!   ELSE
!     CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/X_DIM,Y_DIM,Z_DIM/),Err)
!     CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
!       & NUMBER_GLOBAL_Z_ELEMENTS/),Err)
!   ENDIF
!   CALL CMISSMesh_Initialise(Mesh,Err)
!   CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  !Create a decomposition:
  !All mesh components (associated with G.Projection / Darcy / solid) share the same decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSField_TypeSet(GeometricField,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,Err)  
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,1,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,2,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSField_CreateFinish(GeometricField,Err)
CALL READ_NODES('input/example-nodes-3els',GeometricField)
!   CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create a decomposition

  !Create a field to put the geometry (defualt is geometry)

  SolidMeshComponenetNumber = SolidGeometryMeshComponentNumber

  CALL CMISSField_Initialise(GeometricFieldSolid,Err)
  CALL CMISSField_CreateStart(FieldGeometrySolidUserNumber,Region,GeometricFieldSolid,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricFieldSolid,Decomposition,Err)
  CALL CMISSField_TypeSet(GeometricFieldSolid,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(GeometricFieldSolid,FieldGeometrySolidNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometrySolidNumberOfComponents,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,2,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,3,SolidMeshComponenetNumber,Err)
  CALL CMISSField_CreateFinish(GeometricFieldSolid,Err)
  !Set the mesh component to be used by the field components.
CALL READ_NODES('input/example-nodes-3els',GeometricFieldSolid)
!   CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricFieldSolid,Err)
!   !Set the scaling to use
!   CALL CMISSField_ScalingTypeSet(GeometricFieldSolid,CMISS_FIELD_NO_SCALING,Err)

  !Create a fibre field and attach it to the geometric field
  CALL CMISSField_Initialise(FibreFieldSolid,Err)
  CALL CMISSField_CreateStart(FieldFibreSolidUserNumber,Region,FibreFieldSolid,Err)
  CALL CMISSField_TypeSet(FibreFieldSolid,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreFieldSolid,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(FibreFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSField_NumberOfVariablesSet(FibreFieldSolid,FieldFibreSolidNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,FieldFibreSolidNumberOfComponents,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,2,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,3,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSField_CreateFinish(FibreFieldSolid,Err)

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
    CALL CMISSField_Initialise(EquationsSetFieldDarcy(icompartment),Err)
    CALL CMISSEquationsSet_Initialise(EquationsSetDarcy(icompartment),Err)
    CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumberDarcy,Region,GeometricField, &
      & CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
      & CMISS_EQUATIONS_SET_DARCY_EQUATION_TYPE,CMISS_EQUATIONS_SET_INCOMPRESS_ELAST_MULTI_COMP_DARCY_SUBTYPE,&
      & EquationsSetFieldUserNumberDarcy,EquationsSetFieldDarcy(icompartment),EquationsSetDarcy(icompartment),Err)
    !Finish creating the equations set
    CALL CMISSEquationsSet_CreateFinish(EquationsSetDarcy(icompartment),Err)
    !Set the values for the equations set field to be the current compartment number (1 - N), and the total number of compartments (N)
    CALL CMISSField_ParameterSetUpdateConstant(EquationsSetFieldDarcy(icompartment),CMISS_FIELD_U_VARIABLE_TYPE, &
      & CMISS_FIELD_VALUES_SET_TYPE,1,icompartment,Err)
    CALL CMISSField_ParameterSetUpdateConstant(EquationsSetFieldDarcy(icompartment),CMISS_FIELD_U_VARIABLE_TYPE, &
      & CMISS_FIELD_VALUES_SET_TYPE,2,Ncompartments,Err)
  ENDDO

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create the equations_set
  CALL CMISSField_Initialise(EquationsSetFieldSolid,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSetSolid,Err)
  CALL CMISSEquationsSet_CreateStart(EquationSetSolidUserNumber,Region,FibreFieldSolid,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_INCOMPRESS_ELAST_MULTI_COMP_DARCY_SUBTYPE,&
    & EquationsSetFieldSolidUserNumber,EquationsSetFieldSolid,EquationsSetSolid,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSetSolid,Err)
  !Set the values for the equations set field to be the current compartment number (O for the finite elasticity equations_set), and the total number of compartments (N)
  !Need to store number of compartments, as finite elasticity uses this to calculate the total mass increase for the constiutive law
  CALL CMISSField_ParameterSetUpdateConstant(EquationsSetFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE, &
     & CMISS_FIELD_VALUES_SET_TYPE,1,0_CMISSIntg,Err)
  CALL CMISSField_ParameterSetUpdateConstant(EquationsSetFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE, &
     & CMISS_FIELD_VALUES_SET_TYPE,2,Ncompartments,Err)


  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid Materials Field

  !Create a material field and attach it to the geometric field
  CALL CMISSField_Initialise(MaterialFieldSolid,Err)
  !
  CALL CMISSField_CreateStart(FieldMaterialSolidUserNumber,Region,MaterialFieldSolid,Err)
  !
  CALL CMISSField_TypeSet(MaterialFieldSolid,CMISS_FIELD_MATERIAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(MaterialFieldSolid,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(MaterialFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialFieldSolid,FieldMaterialSolidNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,FieldMaterialSolidNumberOfComponents,Err)
  CALL CMISSField_ComponentMeshComponentSet(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,2,SolidGeometryMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,3,SolidGeometryMeshComponentNumber,Err)
  !
  CALL CMISSField_CreateFinish(MaterialFieldSolid,Err)

  !Set material parameters
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & 2.0_CMISSDP,Err)
!   CALL CMISSField_ComponentValuesInitialise(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2.0e3_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2, &
    & 6.0_CMISSDP,Err)
!   CALL CMISSField_ComponentValuesInitialise(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,33.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3, &
    & 10.0_CMISSDP,Err)

  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetSolid,FieldMaterialSolidUserNumber,MaterialFieldSolid,Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create a dependent field with two variables and four components
  CALL CMISSField_Initialise(DependentFieldSolid,Err)
  !
  CALL CMISSField_CreateStart(FieldDependentSolidUserNumber,Region,DependentFieldSolid,Err)
  !
  CALL CMISSField_TypeSet(DependentFieldSolid,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(DependentFieldSolid,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(DependentFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSField_DependentTypeSet(DependentFieldSolid,CMISS_FIELD_DEPENDENT_TYPE,Err)
  !Create 2N+2 number of variables - 2 for solid, 2N for N Darcy compartments
  FieldDependentSolidNumberOfVariables=2*Ncompartments+2
  CALL CMISSField_NumberOfVariablesSet(DependentFieldSolid,FieldDependentSolidNumberOfVariables,Err)
  !create two variables for each compartment
  ALLOCATE(VariableTypes(2*Ncompartments+2))
  DO num_var=1,Ncompartments+1
     VariableTypes(2*num_var-1)=CMISS_FIELD_U_VARIABLE_TYPE+(CMISS_FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
     VariableTypes(2*num_var)=CMISS_FIELD_DELUDELN_VARIABLE_TYPE+(CMISS_FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
  ENDDO
  CALL CMISSField_VariableTypesSet(DependentFieldSolid,VariableTypes,Err) 
!   CALL CMISSField_VariableTypesSet(DependentFieldSolid,(/CMISS_FIELD_U_VARIABLE_TYPE, &
!     & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_DELVDELN_VARIABLE_TYPE/),Err)
    CALL CMISSField_DimensionSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE, &
       & CMISS_FIELD_VECTOR_DIMENSION_TYPE,Err)
    CALL CMISSField_DimensionSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE, &
       & CMISS_FIELD_VECTOR_DIMENSION_TYPE,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE, &
    & FieldDependentSolidNumberOfComponents,Err)
  DO icompartment=3,2*Ncompartments+2
    CALL CMISSField_NumberOfComponentsSet(DependentFieldSolid,VariableTypes(icompartment),FieldDependentFluidNumberOfComponents,Err)
  ENDDO
!   CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
!   CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
!   CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,3,CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)

  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,SolidDisplMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,2,SolidDisplMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,3,SolidDisplMeshComponentNumber,Err)
  CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,4, &
    & CMISS_FIELD_NODE_BASED_INTERPOLATION, &
    & Err)
!   CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
!   CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,4,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,4,SolidLagrMultMeshComponentNumber,Err)
!   CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1, &
!     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
!   CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2, &
!     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
!   CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3, &
!     & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1, &
    & SolidDisplMeshComponentNumber, &
    & Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2, &
    & SolidDisplMeshComponentNumber, &
    & Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3, &
    & SolidDisplMeshComponentNumber, &
    & Err)
  CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
    & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
!   CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
!     & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
!   CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
    & SolidLagrMultMeshComponentNumber, &
    & Err)
  !loop over the number of compartments
  DO icompartment=3,2*Ncompartments+2
!     CALL CMISSField_DimensionSet(DependentFieldSolid,VariableTypes(icompartment), &
!        & CMISS_FIELD_VECTOR_DIMENSION_TYPE,Err)
    !CALL CMISSField_NumberOfComponentsSet(DependentFieldSolid,VariableTypes(icompartment),FieldDependentFluidNumberOfComponents,Err)
    DO componentnum=1,FieldDependentFluidNumberOfComponents-1
    !set dimension type
!     CALL CMISSField_DimensionSet(DependentField,VariableTypes(icompartment), &
!        & CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
      CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,VariableTypes(icompartment),componentnum, &
       & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,VariableTypes(icompartment),componentnum, & 
         & DarcyVelMeshComponentNumber,Err)
    ENDDO
      CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,VariableTypes(icompartment), &
       & FieldDependentFluidNumberOfComponents, &
       & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
!     CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,VariableTypes(icompartment), &
!       & FieldDependentFluidNumberOfComponents,MESH_COMPONENT_NUMBER_PRESSURE,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,VariableTypes(icompartment), &
      & FieldDependentFluidNumberOfComponents,DarcyMassIncreaseMeshComponentNumber,Err)
    
  ENDDO

!   CALL CMISSField_NumberOfComponentsSet(DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,FieldDependentFluidNumberOfComponents,Err)
!   CALL CMISSField_NumberOfComponentsSet(DependentFieldSolid,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,FieldDependentFluidNumberOfComponents,Err)
!   !For this equation type, MESH_COMPONENT_NUMBER_PRESSURE is actually the mass increase component as the pressure is taken from the solid equations
!   CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,3,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,4,MESH_COMPONENT_NUMBER_PRESSURE,Err)
!   CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,3,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,4,MESH_COMPONENT_NUMBER_PRESSURE,Err)

  !
  CALL CMISSField_CreateFinish(DependentFieldSolid,Err)
  !
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetSolid,FieldDependentSolidUserNumber,DependentFieldSolid,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetSolid,Err)

!   !Initialise dependent field (solid displacement and pressure)
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS !+1
!     CALL CMISSField_ComponentValuesInitialise(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
!       & COMPONENT_NUMBER,INITIAL_FIELD_SOLID(COMPONENT_NUMBER),Err)
!   ENDDO

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  DO icompartment = 1,Ncompartments
    CALL CMISSEquationsSet_DependentCreateStart(EquationsSetDarcy(icompartment),FieldDependentSolidUserNumber,&
      & DependentFieldSolid,Err)
    CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetDarcy(icompartment),Err)
  ENDDO


  !Create the equations set dependent field variables for ALE Darcy
!   CALL CMISSField_Initialise(DependentFieldDarcy,Err)
!   CALL CMISSEquationsSet_DependentCreateStart(EquationsSetDarcy,DependentFieldUserNumberDarcy, & ! ??? UserNumber ???
!   CALL CMISSEquationsSet_DependentCreateStart(EquationsSetDarcy,FieldDependentSolidUserNumber, & ! ??? UserNumber ???
!     & DependentFieldSolid,Err)
! !   !Set the mesh component to be used by the field components.
! !   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
! !     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, &
! !       & MESH_COMPONENT_NUMBER_VELOCITY,Err)
! !     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDarcy,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, &
! !       & MESH_COMPONENT_NUMBER_VELOCITY,Err)
! !   ENDDO
! !   COMPONENT_NUMBER=NUMBER_OF_DIMENSIONS+1
! !     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, &
! !       & MESH_COMPONENT_NUMBER_PRESSURE,Err)
! !     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDarcy,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, &
! !       & MESH_COMPONENT_NUMBER_PRESSURE,Err)
!   !Finish the equations set dependent field variables
!   CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetDarcy,Err)

  !Initialise dependent field (velocity components,pressure,mass increase)
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS+1
    CALL CMISSField_ComponentValuesInitialise(DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
    CALL CMISSField_ComponentValuesInitialise(DependentFieldSolid,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
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
!     CALL CMISSField_Initialise(MaterialsFieldDarcy(icompartment),Err)
!     CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetDarcy(icompartment),MaterialsFieldUserNumberDarcy,&
!          & MaterialsFieldDarcy(icompartment),Err)
!     CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetDarcy(icompartment),Err)
!   END DO
!   DO icompartment = 1,Ncompartments
!     CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMISS_FIELD_U_VARIABLE_TYPE, &
!       & CMISS_FIELD_VALUES_SET_TYPE, &
!       & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
!     CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMISS_FIELD_U_VARIABLE_TYPE, &
!       & CMISS_FIELD_VALUES_SET_TYPE, &
!       & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)
!   END DO
!   DO icompartment = 1, Ncompartments
!     DO COMPONENT_NUMBER=1, Ncompartments
!       CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMISS_FIELD_V_VARIABLE_TYPE, &
!          & CMISS_FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
! !         CALL CMISSField_ParameterSetUpdateConstant(MaterialsFieldDarcy(icompartment),CMISS_FIELD_V_VARIABLE_TYPE, &
! !           & CMISS_FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
!     END DO
!   END DO
!   DO icompartment = 1, Ncompartments
!     DO COMPONENT_NUMBER=1,Nparams
!       CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMISS_FIELD_U1_VARIABLE_TYPE, &
!          & CMISS_FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,ConstitutiveParams(icompartment,COMPONENT_NUMBER),Err)
! !         CALL CMISSField_ParameterSetUpdateConstant(MaterialsFieldDarcy(icompartment),CMISS_FIELD_V_VARIABLE_TYPE, &
! !           & CMISS_FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
!     END DO
!   END DO
  DO icompartment = 1,Ncompartments
     MaterialsFieldUserNumberDarcy = 400+icompartment
call READ_FIELD('input/new-format-field-darcy',MaterialsFieldUserNumberDarcy,Region,GeometricField, &
   & MaterialsFieldDarcy(icompartment))
    CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetDarcy(icompartment),MaterialsFieldUserNumberDarcy,&
         & MaterialsFieldDarcy(icompartment),Err)
    CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetDarcy(icompartment),Err)
  ENDDO

  !Create the equations set materials field variables for ALE Darcy
!   CALL CMISSField_Initialise(MaterialsFieldDarcy,Err)
!   CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetDarcy,MaterialsFieldUserNumberDarcy, &
!     & MaterialsFieldDarcy,Err)
!   !Finish the equations set materials field variables
!   CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetDarcy,Err)
!   CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
!     & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
!   CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
!     & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)
  !Create the equations set materials field variables for deformation-dependent material properties

  !
  !================================================================================================================================
  !

  !EQUATIONS

 
  DO icompartment=1,Ncompartments
    !Create the equations set equations
    CALL CMISSEquations_Initialise(EquationsDarcy(icompartment),Err)
    CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetDarcy(icompartment),EquationsDarcy(icompartment),Err)
    !Set the equations matrices sparsity type
    CALL CMISSEquations_SparsityTypeSet(EquationsDarcy(icompartment),CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !   !Set the equations lumping type
  !   CALL CMISSEquations_LumpingTypeSet(EquationsDarcy,CMISS_EQUATIONS_UNLUMPED_MATRICES,Err)
    !Set the equations set output
    CALL CMISSEquations_OutputTypeSet(EquationsDarcy(icompartment),EQUATIONS_DARCY_OUTPUT,Err)
  !Finish the equations set equations
    CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetDarcy(icompartment),Err)
  ENDDO
  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create the equations set equations
  CALL CMISSEquations_Initialise(EquationsSolid,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetSolid,EquationsSolid,Err)
  CALL CMISSEquations_SparsityTypeSet(EquationsSolid,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(EquationsSolid,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  !
  !================================================================================================================================
  !

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  CALL CMISSField_ComponentValuesInitialise(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4, &
    & 0.0_CMISSDP, &
    & Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  !
  !================================================================================================================================
  !

  !PROBLEMS

  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_MULTI_PHYSICS_CLASS,CMISS_PROBLEM_FINITE_ELASTICITY_DARCY_TYPE, &
    & CMISS_PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,Err)
  CALL CMISSProblem_CreateFinish(Problem,Err)

  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
!   CALL CMISSControlLoop_MaximumIterationsSet(ControlLoop,2,Err)
  CALL CMISSControlLoop_TimesSet(ControlLoop,DYNAMIC_SOLVER_DARCY_START_TIME,DYNAMIC_SOLVER_DARCY_STOP_TIME, &
    & DYNAMIC_SOLVER_DARCY_TIME_INCREMENT,Err)
  CALL CMISSControlLoop_TimeOutputSet(ControlLoop,DYNAMIC_SOLVER_DARCY_OUTPUT_FREQUENCY,Err)
!   CALL CMISSControlLoop_OutputTypeSet(ControlLoop,CMISS_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)



  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL CMISSSolver_Initialise(SolverSolid,Err)
  CALL CMISSSolver_Initialise(DynamicSolverDarcy,Err)
  CALL CMISSSolver_Initialise(LinearSolverDarcy,Err)

  CALL CMISSProblem_SolversCreateStart(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Get the finite elasticity solver
  CALL CMISSProblem_SolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopSolidNumber,CMISS_CONTROL_LOOP_NODE/), &
    & SolverSolidIndex,SolverSolid,Err)
  CALL CMISSSolver_OutputTypeSet(SolverSolid,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
!   CALL CMISSSolver_NewtonJacobianCalculationTypeSet(SolverSolid,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(SolverSolid,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)

  CALL CMISSSolver_NewtonAbsoluteToleranceSet(SolverSolid,ABSOLUTE_TOLERANCE,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(SolverSolid,RELATIVE_TOLERANCE,Err)
  CALL CMISSSolver_NewtonMaximumIterationsSet(SolverSolid,MAXIMUM_ITERATIONS,Err)

!   CALL CMISSSolverNonLinearTypeSet(SolverSolid,CMISS_SOLVER_NONLINEAR_NEWTON,Err)
!   CALL CMISSSolver_LibraryTypeSet(SolverSolid,CMISS_SOLVER_PETSC_LIBRARY,Err)

!   CALL CMISSSolver_NewtonLinearSolverGet(SolverSolid,LinearSolverSolid,Err)
!   CALL CMISSSolver_LinearTypeSet(LinearSolverSolid,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  !Get the Darcy solver
  CALL CMISSProblem_SolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMISS_CONTROL_LOOP_NODE/), &
    & SolverDarcyIndex,DynamicSolverDarcy,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(DynamicSolverDarcy,DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE,Err)
  !Set theta
  CALL CMISSSolver_DynamicThetaSet(DynamicSolverDarcy,DYNAMIC_SOLVER_DARCY_THETA,Err)
  !CALL CMISSSolverDynamicDynamicSet(DynamicSolverDarcy,.TRUE.,Err)
  !Get the dynamic linear solver
  CALL CMISSSolver_DynamicLinearSolverGet(DynamicSolverDarcy,LinearSolverDarcy,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_DARCY_DIRECT_FLAG) THEN
    CALL CMISSSolver_LinearTypeSet(LinearSolverDarcy,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL CMISSSolver_LibraryTypeSet(LinearSolverDarcy,CMISS_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL CMISSSolver_LinearTypeSet(LinearSolverDarcy,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolverDarcy,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolver_LinearIterativeDivergenceToleranceSet(LinearSolverDarcy,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeRelativeToleranceSet(LinearSolverDarcy,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(LinearSolverDarcy,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeGMRESRestartSet(LinearSolverDarcy,RESTART_VALUE,Err)
  ENDIF

  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(SolverSolid,Err)
  CALL CMISSSolver_Initialise(LinearSolverDarcy,Err)

  CALL CMISSSolverEquations_Initialise(SolverEquationsSolid,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsDarcy,Err)

  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  !
  !Get the finite elasticity solver equations
  CALL CMISSProblem_SolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopSolidNumber,CMISS_CONTROL_LOOP_NODE/), &
    & SolverSolidIndex,SolverSolid,Err)
  CALL CMISSSolver_SolverEquationsGet(SolverSolid,SolverEquationsSolid,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsSolid,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsSolid,EquationsSetSolid,EquationsSetIndex,Err)
  !
  !Get the Darcy solver equations
  CALL CMISSProblem_SolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMISS_CONTROL_LOOP_NODE/), &
    & SolverDarcyIndex,LinearSolverDarcy,Err)
  CALL CMISSSolver_SolverEquationsGet(LinearSolverDarcy,SolverEquationsDarcy,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsDarcy,CMISS_SOLVER_SPARSE_MATRICES,Err)
  DO icompartment=1,Ncompartments
    CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsDarcy,EquationsSetDarcy(icompartment),EquationsSetIndex,Err)
  ENDDO
  !
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)
!   ENDDO




  !Prescribe boundary conditions (absolute nodal parameters)
  !Solid is computed in absolute position, rather than displacement. Thus BCs for absolute position
   ! ASSIGN BOUNDARY CONDITIONS - SOLID (absolute nodal parameters)
  !Solid is computed in absolute position, rather than displacement. Thus BCs for absolute position
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsSolid,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsSolid,BoundaryConditionsSolid,Err)

  !Get surfaces
!   CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_FRONT_SURFACE, &
!     & Face1Nodes,FaceXi(1),Err)
!   CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_BACK_SURFACE, &
!     & Face2Nodes,FaceXi(2),Err)
!   CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_RIGHT_SURFACE, &
!     & Face3Nodes,FaceXi(3),Err)
!   CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE, &
!     & Face4Nodes,FaceXi(4),Err)
!   CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_TOP_SURFACE, &
!     & Face5Nodes,FaceXi(5),Err)
!   CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_BOTTOM_SURFACE, &
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
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetGetNode(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,3, &
        & ZCoord,Err)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,NODE,3, &
        & CMISS_BOUNDARY_CONDITION_FIXED,ZCoord,Err)
      WRITE(*,*) "FIXING NODE",NODE,"AT BOTTOM IN Z DIRECTION"
    ENDIF
  ENDDO

  ! Fix the top in z direction
  DO NN=1,SIZE(Face5Nodes,1)
    NODE=Face5Nodes(NN)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetGetNode(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,3, &
        & ZCoord,Err)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,NODE,3, &
        & CMISS_BOUNDARY_CONDITION_FIXED,ZCoord,Err)
      WRITE(*,*) "FIXING NODE",NODE,"AT TOP IN Z DIRECTION"
    ENDIF
  ENDDO

  !Fix more nodes at the bottom to stop free body motion
  X_FIXED=.FALSE.
  Y_FIXED=.FALSE.
  DO NN=1,SIZE(Face6Nodes,1)
    NODE=Face6Nodes(NN)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetGetNode(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,1, &
        & XCoord,Err)
      CALL CMISSField_ParameterSetGetNode(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,2, &
        & YCoord,Err)

      !Fix Origin displacement in x and y (z already fixed)
      IF(ABS(XCoord)<1.0E-6_CMISSDP.AND.ABS(YCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,NODE,1, &
          & CMISS_BOUNDARY_CONDITION_FIXED,XCoord,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,NODE,2, &
          & CMISS_BOUNDARY_CONDITION_FIXED,YCoord,Err)
        WRITE(*,*) "FIXING ORIGIN NODE",NODE,"IN X AND Y DIRECTION"
        X_FIXED=.TRUE.
        Y_FIXED=.TRUE.
      ENDIF

      !Fix nodal displacements at (X_DIM,0) in y
      IF(ABS(XCoord - X_DIM)<1.0E-6_CMISSDP .AND. ABS(YCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,NODE,2, &
          & CMISS_BOUNDARY_CONDITION_FIXED,YCoord,Err)
        WRITE(*,*) "FIXING NODES",NODE,"AT (X_DIM,0) IN Y DIRECTION"
        Y_FIXED=.TRUE.
      ENDIF

      !Fix nodal displacements at (0,Y_DIM) in x
      IF(ABS(XCoord)<1.0E-6_CMISSDP .AND. ABS(YCoord - Y_DIM)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,NODE,1, &
          & CMISS_BOUNDARY_CONDITION_FIXED,XCoord,Err)
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

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsSolid,Err)


  !
  !------------------------------------
  ! ASSIGN BOUNDARY CONDITIONS - FLUID
    !Get surfaces
!     CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_FRONT_SURFACE, &
!       & Face7Nodes,FaceXi(1),Err)
!     CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_BACK_SURFACE, &
!       & Face8Nodes,FaceXi(2),Err)
!     CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_RIGHT_SURFACE, &
!       & Face9Nodes,FaceXi(3),Err)
!     CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE, &
!       & Face10Nodes,FaceXi(4),Err)
!     CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_TOP_SURFACE, &
!       & Face11Nodes,FaceXi(5),Err)
!     CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_BOTTOM_SURFACE, &
!       & Face12Nodes,FaceXi(6),Err)


allocate(Face7Nodes(8),Face8Nodes(8),Face9Nodes(8),Face10Nodes(8),Face11Nodes(4),Face12Nodes(4))
Face7Nodes=[  1,  2,  5,  6,  9, 10, 13, 14]
Face8Nodes=[  3,  4,  7,  8, 11, 12, 15, 16]
Face9Nodes=[  2,  4,  6,  8, 10, 12, 14, 16]
Face10Nodes=[  1,  3,  5,  7,  9, 11, 13, 15]
Face11Nodes=[ 13, 14, 15, 16]
Face12Nodes=[  1,  2,  3,  4]
FaceXi=[ -2,  2,  1, -1,  3, -3]

  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsDarcy,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)
  DO icompartment=1,Ncompartments

    IF(icompartment==1) THEN
    ! At the top impose Darcy velocity in z direction
    DO NN=1,SIZE(Face11Nodes,1)
      NODE=Face11Nodes(NN)
! !     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = -1.0_CMISSDP
        COMPONENT_NUMBER = 3
        write(*,*)'Marker 0'
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,NODE, &
          & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        WRITE(*,*) "SPECIFIED INFLOW AT NODE",NODE,"IN Z DIRECTION"

! !       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,1,XCoord,Err)
! !       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,2,YCoord,Err)
! !       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,3,ZCoord,Err)
! !       WRITE(*,*) "XCoord, YCoord, ZCoord = ",XCoord, YCoord, ZCoord
! !     ENDIF
    ENDDO

    !All other faces are impermeable
    DO NN=1,SIZE(Face7Nodes,1)
      NODE=Face7Nodes(NN)
! !     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 1
        write(*,*)'Marker 1'
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,NODE, &
          & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN X DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face8Nodes,1)
      NODE=Face8Nodes(NN)
! !     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 1
        write(*,*)'Marker 2'
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,NODE, &
          & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN X DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face9Nodes,1)
      NODE=Face9Nodes(NN)
! !     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 2
        write(*,*)'Marker 3'
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,NODE, &
          & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Y DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face10Nodes,1)
      NODE=Face10Nodes(NN)
! !     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 2
        write(*,*)'Marker 4'
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,NODE, &
          & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Y DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face12Nodes,1)
      NODE=Face12Nodes(NN)
! !     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 3
        write(*,*)'Marker 5'
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,NODE, &
          & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Z DIRECTION"
! !     ENDIF
    ENDDO
    ENDIF
    !Compartment TWO boundary conditions!
    ! At the top impose Darcy velocity in z direction
    IF(icompartment==2) THEN
    DO NN=1,SIZE(Face11Nodes,1)
      NODE=Face11Nodes(NN)
! !     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = -1.0_CMISSDP
        COMPONENT_NUMBER = 3
        write(*,*)'Marker 0'
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_U1_VARIABLE_TYPE,1,NODE, &
          & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        WRITE(*,*) "SPECIFIED INFLOW AT NODE",NODE,"IN Z DIRECTION"

! !       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,1,XCoord,Err)
! !       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,2,YCoord,Err)
! !       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,3,ZCoord,Err)
! !       WRITE(*,*) "XCoord, YCoord, ZCoord = ",XCoord, YCoord, ZCoord
! !     ENDIF
    ENDDO
!     DO NN=1,SIZE(Face11Nodes,1)
!       NODE=Face11Nodes(NN)
! ! !     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! ! !     IF(NodeDomain==ComputationalNodeNumber) THEN
!         VALUE = 0.0_CMISSDP
!         COMPONENT_NUMBER = 3
!         write(*,*)'Marker 0'
!         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_U1_VARIABLE_TYPE,1,NODE,COMPONENT_NUMBER,&
!           & CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!         WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Z DIRECTION"
!
! ! !       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,1,XCoord,Err)
! ! !       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,2,YCoord,Err)
! ! !       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,3,ZCoord,Err)
! ! !       WRITE(*,*) "XCoord, YCoord, ZCoord = ",XCoord, YCoord, ZCoord
! ! !     ENDIF
!     ENDDO

    !All other faces are impermeable
    DO NN=1,SIZE(Face7Nodes,1)
      NODE=Face7Nodes(NN)
! !     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 1
        write(*,*)'Marker 1'
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_U1_VARIABLE_TYPE,1,NODE, &
          & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN X DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face8Nodes,1)
      NODE=Face8Nodes(NN)
! !     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 1
        write(*,*)'Marker 2'
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_U1_VARIABLE_TYPE,1,NODE, &
          & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN X DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face9Nodes,1)
      NODE=Face9Nodes(NN)
! !     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 2
        write(*,*)'Marker 3'
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_U1_VARIABLE_TYPE,1,NODE, &
          & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Y DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face10Nodes,1)
      NODE=Face10Nodes(NN)
! !     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 2
        write(*,*)'Marker 4'
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_U1_VARIABLE_TYPE,1,NODE, &
          & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Y DIRECTION"
! !     ENDIF
    ENDDO

    DO NN=1,SIZE(Face12Nodes,1)
      NODE=Face12Nodes(NN)
! !     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
! !     IF(NodeDomain==ComputationalNodeNumber) THEN
        VALUE = 0.0_CMISSDP
        COMPONENT_NUMBER = 3
        write(*,*)'Marker 5'
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_U1_VARIABLE_TYPE,1,NODE, &
          & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Z DIRECTION"
! !     ENDIF
    ENDDO
    ENDIF

  !BOUNDARY CONDITIONS
  !Start the creation of the equations set boundary conditions for Darcy
!   DO icompartment=1,Ncompartments  
!     CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsDarcy,DependentFieldSolid,Err)
!     CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSetDarcy(icompartment),BoundaryConditionsDarcy,DependentFieldSolid,Err)
!     CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSetDarcy(icompartment),Err)
  ENDDO
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsDarcy,Err)

  !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn of PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL CMISSProblem_Solve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"


  !
  !================================================================================================================================
  !

  !OUTPUT

  WRITE(*,*)'dummy_counter = ',dummy_counter

  EXPORT_FIELD_IO=.FALSE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"FiniteElasticityMultiCompDarcy","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"FiniteElasticityMultiCompDarcy","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
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
