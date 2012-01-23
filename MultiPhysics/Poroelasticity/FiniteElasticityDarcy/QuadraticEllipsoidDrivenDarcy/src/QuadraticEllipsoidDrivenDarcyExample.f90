!> \file
!> \author Christian Michler
!> \brief This is an example program to solve a coupled Finite Elastiticity Darcy equation on a cylindrical geometry using openCMISS calls.
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
!> Contributor(s): Christian Michler, Jack Lee
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

!> \example MultiPhysics/Poroelasticity/FiniteElasticityDarcy/QuadraticEllipsoidDrivenDarcy/src/QuadraticEllipsoidDrivenDarcyExample.f90
!! Example program to solve coupled FiniteElasticityDarcy equations using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/Poroelasticity/FiniteElasticityDarcy/QuadraticEllipsoidDrivenDarcy/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/Poroelasticity/FiniteElasticityDarcy/QuadraticEllipsoidDrivenDarcy/build-intel'>Linux GNU Build</a>
!!
!<

! !
! !  This example considers a coupled Finite Elasticity Darcy problem on a cylindrical geometry
! !

!> Main program

PROGRAM QUADRATICELLIPSOIDDRIVENDARCYEXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OPENCMISS
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

  REAL(CMISSDP), PARAMETER :: PI=3.14159_CMISSDP
  REAL(CMISSDP), PARAMETER :: LONG_AXIS=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: SHORT_AXIS=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WALL_THICKNESS=0.5_CMISSDP
  REAL(CMISSDP), PARAMETER :: CUTOFF_ANGLE=1.5708_CMISSDP
  REAL(CMISSDP), PARAMETER :: FIBRE_SLOPE_INTERSECTION=1.73205_CMISSDP !Slope of fibres in base endocardium = 60 degrees
  REAL(CMISSDP), PARAMETER :: FIBRE_SLOPE_CHANGE=-3.4641_CMISSDP !Slope change of fibres from 60 to -60 degrees in transmural direction 
  REAL(CMISSDP), PARAMETER :: SHEET_SLOPE_BASE_ENDO=1.0_CMISSDP !Slope of sheet at base endocardium 

  REAL(CMISSDP), PARAMETER :: INNER_PRESSURE=2.0_CMISSDP  !Positive is compressive
  REAL(CMISSDP), PARAMETER :: OUTER_PRESSURE=0.0_CMISSDP !Positive is compressive
  REAL(CMISSDP), PARAMETER :: C1= 2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: C2= 6.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: C3=10.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: NumberGlobalXElements=4  ! X ==NUMBER_GLOBAL_CIRCUMFERENTIAL_ELEMENTS
  INTEGER(CMISSIntg), PARAMETER :: NumberGlobalYElements=4  ! Y ==NUMBER_GLOBAL_LONGITUDINAL_ELEMENTS
  INTEGER(CMISSIntg), PARAMETER :: NumberGlobalZElements=2  ! Z ==NUMBER_GLOBAL_TRANSMURAL_ELEMENTS
  INTEGER(CMISSIntg) :: NumberOfDomains

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticCollapsedBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: LinearCollapsedBasisUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DerivativeUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensions=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2
  INTEGER(CMISSIntg), PARAMETER :: TotalNumberOfElements=1

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumberSolid=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumberDarcy=11
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialNumberOfComponents=3 !2

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentNumberOfVariables=4 !2
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentSolidNumberOfComponents=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentFluidNumberOfComponents=4

  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldDarcyUserNumber=15

  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumberSolid=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberSolid=13
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1


  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcy=8
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberDarcy=12
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberDarcy=22
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldDarcyUserNumber=42

  INTEGER(CMISSIntg), PARAMETER :: ControlLoopSolidNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopFluidNumber=2
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopSubiterationNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverSolidIndex=1
  INTEGER(CMISSIntg), PARAMETER :: SolverDarcyIndex=1

  !Program types


  !Program variables

  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  REAL(CMISSDP) :: FibreFieldAngle(3) 
  REAL(CMISSDP) :: nu,theta,omega,XI3,XI3delta,XI2delta, zero
  INTEGER(CMISSIntg) ::i,j,k,component_idx,node_idx,TOTAL_NUMBER_NODES_XI(3)
  !For grabbing surfaces
  INTEGER(CMISSIntg) :: InnerNormalXi,OuterNormalXi,TopNormalXi
  INTEGER(CMISSIntg), ALLOCATABLE :: InnerSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: OuterSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: TopSurfaceNodes(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: TopSurfaceNodesDarcyVel(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: InnerSurfaceNodesDarcyVel(:)
  INTEGER(CMISSIntg), ALLOCATABLE :: OuterSurfaceNodesDarcyVel(:)
  INTEGER(CMISSIntg) :: NN,NODE,NodeDomain
  REAL(CMISSDP) :: XCoord,YCoord,ZCoord
  LOGICAL :: X_FIXED,Y_FIXED,X_OKAY,Y_OKAY

  INTEGER(CMISSIntg) :: GeometricFieldDarcyMeshComponentNumber, DarcyVelMeshComponentNumber, DarcyMassIncreaseMeshComponentNumber
  INTEGER(CMISSIntg) :: MeshComponentNumber_dummy

  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS

  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE

  INTEGER(CMISSIntg) :: EQUATIONS_DARCY_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER, NODE_NUMBER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_DARCY_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_DARCY_OUTPUT_TYPE

  REAL(CMISSDP) :: COORD_X, COORD_Y, COORD_Z
  REAL(CMISSDP) :: DOMAIN_X1, DOMAIN_X2, DOMAIN_Y1, DOMAIN_Y2, DOMAIN_Z1, DOMAIN_Z2
  REAL(CMISSDP) :: GEOMETRY_TOLERANCE
  INTEGER(CMISSIntg) :: EDGE_COUNT
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SOLID
  REAL(CMISSDP) :: INITIAL_FIELD_DARCY(4)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE
  REAL(CMISSDP) :: POROSITY_PARAM_DARCY, PERM_OVER_VIS_PARAM_DARCY

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_DARCY_DIRECT_FLAG

  !CMISS variables

  TYPE(CMISSBasisType) :: QuadraticBasis,QuadraticCollapsedBasis,LinearBasis,LinearCollapsedBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSetSolid
  TYPE(CMISSFieldType) :: GeometricFieldSolid,FibreFieldSolid,MaterialFieldSolid
  TYPE(CMISSFieldType) :: DependentField,EquationsSetFieldSolid

  TYPE(CMISSFieldType) :: IndependentFieldDarcy

  TYPE(CMISSFieldType) :: GeometricFieldDarcy

  TYPE(CMISSFieldType) :: MaterialsFieldDarcy
  TYPE(CMISSFieldType) :: EquationsSetFieldDarcy
  TYPE(CMISSFieldType) :: SourceFieldDarcy
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsDarcy
  !Equations sets
  TYPE(CMISSEquationsSetType) :: EquationsSetDarcy
  !Equations
  TYPE(CMISSEquationsType) :: EquationsDarcy

  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: SolverSolid,LinearSolverSolid
  TYPE(CMISSSolverEquationsType) :: SolverEquationsSolid
  TYPE(CMISSControlLoopType) :: ControlLoop

  !Solvers
  TYPE(CMISSSolverType) :: DynamicSolverDarcy
  TYPE(CMISSSolverType) :: LinearSolverDarcy
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsDarcy

  REAL(CMISSDP) :: DYNAMIC_SOLVER_DARCY_START_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_DARCY_STOP_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_DARCY_THETA
  REAL(CMISSDP) :: DYNAMIC_SOLVER_DARCY_TIME_INCREMENT

  INTEGER(CMISSIntg),allocatable :: ElementUserNodes(:)
  INTEGER(CMISSIntg) :: NUMBER_USER_ELEMENT_NODES, ELEMENT_NUMBER

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables
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
  DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE=CMISS_SOLVER_PROGRESS_OUTPUT
  LINEAR_SOLVER_DARCY_OUTPUT_TYPE=CMISS_SOLVER_SOLVER_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DARCY_OUTPUT=CMISS_EQUATIONS_NO_OUTPUT

  !Set time parameter
  DYNAMIC_SOLVER_DARCY_START_TIME=1.0E-3_CMISSDP
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


  !LinearMeshComponentNumber/QuadraticMeshComponentNumber
  DarcyVelMeshComponentNumber = LinearMeshComponentNumber
  DarcyMassIncreaseMeshComponentNumber = LinearMeshComponentNumber
!   GeometricFieldDarcyMeshComponentNumber = DarcyVelMeshComponentNumber
  GeometricFieldDarcyMeshComponentNumber = QuadraticMeshComponentNumber


  !
  !================================================================================================================================
  !

  !Intialise cmiss
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
  CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,(/1,2,3,4,5/),"Diagnostics",(/"PROBLEM_FINITE_ELEMENT_CALCULATE"/),Err)

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
!   CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  NumberOfDomains=NumberOfComputationalNodes

  !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_TypeSet(CoordinateSystem,CMISS_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL CMISSCoordinateSystem_OriginSet(CoordinateSystem,(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/),Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !Create a region and assign the CS to the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !Define basis functions - tri-linear Lagrange and tri-Quadratic Lagrange, each with collapsed variant
    !Quadratic Basis
  CALL CMISSBasis_Initialise(QuadraticBasis,Err)
  CALL CMISSBasis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasis_InterpolationXiSet(QuadraticBasis,(/CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & (/CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME/),Err)
!     & (/CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME/),Err)
  CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err) !Have to do this
  CALL CMISSBasis_CreateFinish(QuadraticBasis,Err)
  
    !Collapsed Quadratic Basis
  CALL CMISSBasis_Initialise(QuadraticCollapsedBasis,Err)
  CALL CMISSBasis_CreateStart(QuadraticCollapsedBasisUserNumber,QuadraticCollapsedBasis,Err)
  CALL CMISSBasis_TypeSet(QuadraticCollapsedBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(QuadraticCollapsedBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasis_InterpolationXiSet(QuadraticCollapsedBasis,(/CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
       & CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION/),Err)
  CALL CMISSBasis_CollapsedXiSet(QuadraticCollapsedBasis,(/CMISS_BASIS_XI_COLLAPSED, &
       & CMISS_BASIS_COLLAPSED_AT_XI0,CMISS_BASIS_NOT_COLLAPSED/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(QuadraticCollapsedBasis, &
       & (/CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME/),Err)  
!     & (/CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME/),Err)
  CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(QuadraticCollapsedBasis,.true.,Err) !Have to do this
  CALL CMISSBasis_CreateFinish(QuadraticCollapsedBasis,Err)

    !Linear Basis
  CALL CMISSBasis_Initialise(LinearBasis,Err)
  CALL CMISSBasis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearBasis, &
    & (/CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME/),Err)
!     & (/CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME/),Err)
  CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err) !Have to do this (unused) due to field_interp setup
  CALL CMISSBasis_CreateFinish(LinearBasis,Err)

    !Collapsed Linear Basis
  CALL CMISSBasis_Initialise(LinearCollapsedBasis,Err)
  CALL CMISSBasis_CreateStart(LinearCollapsedBasisUserNumber,LinearCollapsedBasis,Err)
  CALL CMISSBasis_TypeSet(LinearCollapsedBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL CMISSBasis_NumberOfXiSet(LinearCollapsedBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasis_InterpolationXiSet(LinearCollapsedBasis,(/CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
       & CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION/),Err)
  CALL CMISSBasis_CollapsedXiSet(LinearCollapsedBasis,(/CMISS_BASIS_XI_COLLAPSED,CMISS_BASIS_COLLAPSED_AT_XI0, &
    & CMISS_BASIS_NOT_COLLAPSED/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearCollapsedBasis, &
       & (/CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME,CMISS_BASIS_MID_QUADRATURE_SCHEME/),Err)
!     & (/CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME/),Err)
  CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(LinearCollapsedBasis,.true.,Err) !Have to do this (unused) due to field_interp setup
  CALL CMISSBasis_CreateFinish(LinearCollapsedBasis,Err)

  !
  !================================================================================================================================
  !

  !Start the creation of a generated ellipsoid mesh
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up an ellipsoid mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_ELLIPSOID_MESH_TYPE,Err)
  !Set the quadratic and linear bases
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,[QuadraticBasis,QuadraticCollapsedBasis,LinearBasis,LinearCollapsedBasis],Err)
  !Define the mesh on the region
  CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/LONG_AXIS,SHORT_AXIS,WALL_THICKNESS,CUTOFF_ANGLE/),Err)
  CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,(/NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements/),Err)
  
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !
  !================================================================================================================================
  !

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecomposition_CalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !
  !================================================================================================================================
  !

  ! --- GeometricFieldSolid ---
  !Create a field to put the geometry (default is geometry)
  CALL CMISSField_Initialise(GeometricFieldSolid,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumberSolid,Region,GeometricFieldSolid,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricFieldSolid,Decomposition,Err)
  CALL CMISSField_TypeSet(GeometricFieldSolid,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(GeometricFieldSolid,FieldGeometryNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_CreateFinish(GeometricFieldSolid,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricFieldSolid,Err)

  !
  !================================================================================================================================
  !

  ! --- GeometricFieldDarcy ---
  !Create a field to put the geometry (default is geometry)
  CALL CMISSField_Initialise(GeometricFieldDarcy,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumberDarcy,Region,GeometricFieldDarcy,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricFieldDarcy,Decomposition,Err)
  CALL CMISSField_TypeSet(GeometricFieldDarcy,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(GeometricFieldDarcy,FieldGeometryNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricFieldDarcyMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,2, &
    & GeometricFieldDarcyMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,3, &
    & GeometricFieldDarcyMeshComponentNumber,Err)
  CALL CMISSField_CreateFinish(GeometricFieldDarcy,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricFieldDarcy,Err)

  !
  !================================================================================================================================
  !

  !Create a fibre field and attach it to the geometric field
  CALL CMISSField_Initialise(FibreFieldSolid,Err)
  CALL CMISSField_CreateStart(FieldFibreUserNumber,Region,FibreFieldSolid,Err)
  CALL CMISSField_TypeSet(FibreFieldSolid,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreFieldSolid,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(FibreFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSField_NumberOfVariablesSet(FibreFieldSolid,FieldFibreNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)  
  CALL CMISSField_ComponentMeshComponentSet(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_CreateFinish(FibreFieldSolid,Err)

  !Set Fibre directions (this block is parallel-untested)
  node_idx=0  
  !This is valid only for quadratic basis functions
  TOTAL_NUMBER_NODES_XI(1)=NumberGlobalXElements*2
  TOTAL_NUMBER_NODES_XI(2)=NumberGlobalYElements*2+1
  TOTAL_NUMBER_NODES_XI(3)=NumberGlobalZElements*2+1

  XI2delta=(PI-CUTOFF_ANGLE)/(TOTAL_NUMBER_NODES_XI(2)-1)
  XI3=0
  XI3delta=(1.0)/(TOTAL_NUMBER_NODES_XI(3)-1)
  zero=0
  DO k=1, TOTAL_NUMBER_NODES_XI(3)
    !Apex nodes
    j=1
    i=1
    node_idx=node_idx+1
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      FibreFieldAngle=(/zero,zero,zero/) 
      DO component_idx=1,FieldFibreNumberOfComponents
        CALL CMISSField_ParameterSetUpdateNode(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
          & DerivativeUserNumber,node_idx,component_idx,FibreFieldAngle(component_idx),Err)
      ENDDO
    ENDIF
    theta=atan(FIBRE_SLOPE_CHANGE*XI3+FIBRE_SLOPE_INTERSECTION)
    DO j=2, TOTAL_NUMBER_NODES_XI(2) 
      nu=PI-XI2delta*(j-1)
      omega=PI/2+cos(2*nu)*atan(SHEET_SLOPE_BASE_ENDO)*(-2*XI3+1)
      DO i=1, TOTAL_NUMBER_NODES_XI(1)
        node_idx=node_idx+1
        CALL CMISSDecomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          FibreFieldAngle=(/theta,zero,omega/)
          DO component_idx=1,FieldFibreNumberOfComponents
            CALL CMISSField_ParameterSetUpdateNode(FibreFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
              & DerivativeUserNumber, node_idx,component_idx,FibreFieldAngle(component_idx),Err)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    XI3=XI3+XI3delta
  ENDDO

  !Create a material field and attach it to the geometric field
  CALL CMISSField_Initialise(MaterialFieldSolid,Err)
  CALL CMISSField_CreateStart(FieldMaterialUserNumber,Region,MaterialFieldSolid,Err)
  CALL CMISSField_TypeSet(MaterialFieldSolid,CMISS_FIELD_MATERIAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(MaterialFieldSolid,Decomposition,Err)        
  CALL CMISSField_GeometricFieldSet(MaterialFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSField_NumberOfVariablesSet(MaterialFieldSolid,FieldMaterialNumberOfVariables,Err)
  CALL CMISSField_NumberOfComponentsSet(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponents,Err)  
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,2,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_ComponentInterpolationSet(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,3,CMISS_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL CMISSField_CreateFinish(MaterialFieldSolid,Err)

  !Set Mooney-Rivlin constants 
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,C1,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,C2,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,C3,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for ALE Darcy
  CALL CMISSField_Initialise(EquationsSetFieldDarcy,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSetDarcy,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumberDarcy,Region,GeometricFieldDarcy, &
    & CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & CMISS_EQUATIONS_SET_DARCY_EQUATION_TYPE,CMISS_EQUATIONS_SET_INCOMPRESS_ELASTICITY_DRIVEN_DARCY_SUBTYPE,&
    & EquationsSetFieldUserNumberDarcy,EquationsSetFieldDarcy,EquationsSetDarcy,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSetDarcy,Err)

  !Create the equations set for the solid
  CALL CMISSField_Initialise(EquationsSetFieldSolid,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSetSolid,Err)
  CALL CMISSEquationsSet_CreateStart(EquationSetUserNumberSolid,Region,FibreFieldSolid,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_INCOMPRESS_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
    & EquationsSetFieldUserNumberSolid,EquationsSetFieldSolid,EquationsSetSolid,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create a dependent field with four variables (U, DelUDelN = solid, V, DelVDelN = Darcy) and four components
  !Solid: The U, DelUDelN variables have 4 components (3 displacement, 1 pressure)
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSField_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL CMISSField_TypeSet(DependentField,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(DependentField,GeometricFieldSolid,Err)
  CALL CMISSField_DependentTypeSet(DependentField,CMISS_FIELD_DEPENDENT_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)

  CALL CMISSField_VariableTypesSet(DependentField,(/CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_DELVDELN_VARIABLE_TYPE/),Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_V_VARIABLE_TYPE,FieldDependentFluidNumberOfComponents,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,FieldDependentFluidNumberOfComponents,Err)

  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)

  !Darcy: The V, DelVDelN variables have 4 components (3 velocities, 1 mass increase)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_V_VARIABLE_TYPE,1,DarcyVelMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_V_VARIABLE_TYPE,2,DarcyVelMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_V_VARIABLE_TYPE,3,DarcyVelMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_V_VARIABLE_TYPE,4,DarcyMassIncreaseMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1,DarcyVelMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,2,DarcyVelMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,3,DarcyVelMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,4, &
    & DarcyMassIncreaseMeshComponentNumber,Err)

  CALL CMISSField_ScalingTypeSet(DependentField,CMISS_FIELD_UNIT_SCALING,Err)

  CALL CMISSField_CreateFinish(DependentField,Err)

  !
  !================================================================================================================================
  !

  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetSolid,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetSolid,Err)

  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetSolid,FieldMaterialUserNumber,MaterialFieldSolid,Err)  
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetDarcy,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetDarcy,Err)

  DO COMPONENT_NUMBER=1,FieldDependentFluidNumberOfComponents
    CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !INDEPENDENT FIELD Darcy for storing BC flags

  CALL CMISSField_Initialise(IndependentFieldDarcy,Err)
  CALL CMISSEquationsSet_IndependentCreateStart(EquationsSetDarcy,IndependentFieldDarcyUserNumber, &
    & IndependentFieldDarcy,Err)

  CALL CMISSField_ComponentMeshComponentSet(IndependentFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,1,DarcyVelMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(IndependentFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,2,DarcyVelMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(IndependentFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,3,DarcyVelMeshComponentNumber,Err)
!   CALL CMISSField_ComponentMeshComponentSet(IndependentFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,4,DarcyMassIncreaseMeshComponentNumber,Err)

  CALL CMISSEquationsSet_IndependentCreateFinish(EquationsSetDarcy,Err)

  CALL CMISSField_ComponentValuesInitialise(IndependentFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2, &
    & 0.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3, &
    & 0.0_CMISSDP,Err)

  !
  !================================================================================================================================
  !

  !Create the equations set materials field variables for ALE Darcy
  CALL CMISSField_Initialise(MaterialsFieldDarcy,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetDarcy,MaterialsFieldUserNumberDarcy, &
    & MaterialsFieldDarcy,Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetDarcy,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,POROSITY_PARAM_DARCY,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,PERM_OVER_VIS_PARAM_DARCY,Err)

  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,0.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 4,0.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 5,PERM_OVER_VIS_PARAM_DARCY,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 6,0.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 7,PERM_OVER_VIS_PARAM_DARCY,Err)

  !
  !================================================================================================================================
  !

  !Source field
  CALL CMISSField_Initialise(SourceFieldDarcy,Err)
  CALL CMISSEquationsSet_SourceCreateStart(EquationsSetDarcy,SourceFieldDarcyUserNumber,SourceFieldDarcy,Err)
  CALL CMISSEquationsSet_SourceCreateFinish(EquationsSetDarcy,Err)

!   ELEMENT_NUMBER = 5
!   COMPONENT_NUMBER = 4
!   VALUE = 4.2_CMISSDP
! !   CALL CMISSField_ParameterSetUpdateElement(RegionUserNumber,SourceFieldDarcyUserNumber,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
! !     & ELEMENT_NUMBER,COMPONENT_NUMBER,VALUE,Err)
! 
!   NUMBER_USER_ELEMENT_NODES = 27  !hardcoding is bad - but how to access the number of element nodes ?
!                                   !- there is no CMISS library call available for this
!                                   !- traversing the structure of 'CMISSMeshElementsType' does not work either,
!                                   !  since certain members are private
! 
!   allocate( ElementUserNodes(NUMBER_USER_ELEMENT_NODES) )
! 
!   CALL CMISSMeshElements_NodesGet(RegionUserNumber,MeshUserNumber,QuadraticMeshComponentNumber,ELEMENT_NUMBER, &
!     & ElementUserNodes,Err)
! 
!   DO NN=1,NUMBER_USER_ELEMENT_NODES
!     NODE_NUMBER = ElementUserNodes(NN)
!     CALL CMISSField_ParameterSetUpdateNode(SourceFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
!       & 1,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
!   ENDDO

  !
  !================================================================================================================================
  !

  !EQUATIONS SET EQUATIONS

  !Darcy
  CALL CMISSEquations_Initialise(EquationsDarcy,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetDarcy,EquationsDarcy,Err)
  CALL CMISSEquations_SparsityTypeSet(EquationsDarcy,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(EquationsDarcy,EQUATIONS_DARCY_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetDarcy,Err)

  !Solid
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetSolid,Equations,Err)
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
!   CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,-14.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,0.0_CMISSDP, &
    & Err)

  !
  !================================================================================================================================
  !

  !Define the problem
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_MULTI_PHYSICS_CLASS,CMISS_PROBLEM_FINITE_ELASTICITY_DARCY_TYPE, &
    & CMISS_PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,Err)
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !Create the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
!   CALL CMISSControlLoop_MaximumIterationsSet(ControlLoop,1,Err)  ! this one sets the increment loop counter
  CALL CMISSControlLoop_TimesSet(ControlLoop,DYNAMIC_SOLVER_DARCY_START_TIME,DYNAMIC_SOLVER_DARCY_STOP_TIME, &
    & DYNAMIC_SOLVER_DARCY_TIME_INCREMENT,Err)
  CALL CMISSControlLoop_TimeOutputSet(ControlLoop,DYNAMIC_SOLVER_DARCY_OUTPUT_FREQUENCY,Err)
!   CALL CMISSControlLoop_OutputTypeSet(ControlLoop,CMISS_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !Create the problem solvers
  CALL CMISSSolver_Initialise(SolverSolid,Err)
  CALL CMISSSolver_Initialise(LinearSolverSolid,Err)
  CALL CMISSSolver_Initialise(DynamicSolverDarcy,Err)
  CALL CMISSSolver_Initialise(LinearSolverDarcy,Err)

  CALL CMISSProblem_SolversCreateStart(Problem,Err)

  ! Solid
  CALL CMISSProblem_SolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopSolidNumber,CMISS_CONTROL_LOOP_NODE/), &
    & SolverSolidIndex,SolverSolid,Err)
  CALL CMISSSolver_OutputTypeSet(SolverSolid,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
!   CALL CMISSSolver_NewtonJacobianCalculationTypeSet(SolverSolid,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(SolverSolid,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)

!   CALL CMISSSolverNonLinearTypeSet(SolverSolid,CMISS_SOLVER_NONLINEAR_NEWTON,Err)
!   CALL CMISSSolver_LibraryTypeSet(SolverSolid,CMISS_SOLVER_PETSC_LIBRARY,Err)

  CALL CMISSSolver_NewtonAbsoluteToleranceSet(SolverSolid,ABSOLUTE_TOLERANCE,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(SolverSolid,RELATIVE_TOLERANCE,Err)
  CALL CMISSSolver_NewtonMaximumIterationsSet(SolverSolid,MAXIMUM_ITERATIONS,Err)

  CALL CMISSSolver_NewtonLinearSolverGet(SolverSolid,LinearSolverSolid,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolverSolid,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)

  !Darcy
  CALL CMISSProblem_SolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMISS_CONTROL_LOOP_NODE/), &
    & SolverDarcyIndex,DynamicSolverDarcy,Err)
  CALL CMISSSolver_OutputTypeSet(DynamicSolverDarcy,DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE,Err)
  CALL CMISSSolver_DynamicThetaSet(DynamicSolverDarcy,DYNAMIC_SOLVER_DARCY_THETA,Err)
!   CALL CMISSSolverDynamicDynamicSet(DynamicSolverDarcy,.TRUE.,Err)
  CALL CMISSSolver_DynamicLinearSolverGet(DynamicSolverDarcy,LinearSolverDarcy,Err)
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

CALL CMISSProblem_SolversCreateFinish(Problem,Err)


  !
  !================================================================================================================================
  !

  !Create the problem solver equations
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
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsDarcy,EquationsSetDarcy,EquationsSetIndex,Err)
  !
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsSolid,BoundaryConditions,Err)

  !Grab the list of nodes on inner, outer and top surfaces
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_ELLIPSOID_TOP_SURFACE,TopSurfaceNodes,TopNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_ELLIPSOID_INNER_SURFACE,InnerSurfaceNodes,InnerNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_ELLIPSOID_OUTER_SURFACE,OuterSurfaceNodes,OuterNormalXi,Err)

  write(*,*)'TopSurfaceNodes = ',TopSurfaceNodes

  ! ASSIGN BOUNDARY CONDITIONS
  !Fix base of the ellipsoid in z direction
  DO NN=1,SIZE(TopSurfaceNodes,1)
    NODE=TopSurfaceNodes(NN)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetGetNode(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE,3, &
        & ZCoord,Err)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,3, &
        & CMISS_BOUNDARY_CONDITION_FIXED,ZCoord,Err)
    ENDIF
  ENDDO

  !Apply inner surface pressure
  !NOTE: Surface pressure goes into pressure_values_set_type of the DELUDELN type
  DO NN=1,SIZE(InnerSurfaceNodes,1)
    NODE=InnerSurfaceNodes(NN)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NODE, &
        & ABS(InnerNormalXi),CMISS_BOUNDARY_CONDITION_PRESSURE,INNER_PRESSURE,Err)
    ENDIF
  ENDDO

  !Apply outer surface pressure
  DO NN=1,SIZE(OuterSurfaceNodes,1)
    NODE=OuterSurfaceNodes(NN)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NODE, &
        & ABS(OuterNormalXi),CMISS_BOUNDARY_CONDITION_PRESSURE,OUTER_PRESSURE,Err)
    ENDIF
  ENDDO

  !Fix more nodes at the base to stop free body motion
  X_FIXED=.FALSE.
  Y_FIXED=.FALSE.
  DO NN=1,SIZE(TopSurfaceNodes,1)
    NODE=TopSurfaceNodes(NN)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetGetNode(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE,1, &
        & XCoord,Err)
      CALL CMISSField_ParameterSetGetNode(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE,2, &
        & YCoord,Err)
      IF(ABS(XCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,1, &
          & CMISS_BOUNDARY_CONDITION_FIXED,XCoord,Err)
        WRITE(*,*) "FIXING NODE",NODE,"IN X DIRECTION"
        X_FIXED=.TRUE.
      ENDIF
      IF(ABS(YCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,2, &
          & CMISS_BOUNDARY_CONDITION_FIXED,YCoord,Err)
        WRITE(*,*) "FIXING NODE",NODE,"IN Y DIRECTION"
        Y_FIXED=.TRUE.
    ENDIF
    ENDIF
  ENDDO
  CALL MPI_REDUCE(X_FIXED,X_OKAY,1,MPI_LOGICAL,MPI_LOR,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_REDUCE(Y_FIXED,Y_OKAY,1,MPI_LOGICAL,MPI_LOR,0,MPI_COMM_WORLD,MPI_IERROR)
  IF(ComputationalNodeNumber==0) THEN
    IF(.NOT.(X_OKAY.AND.Y_OKAY)) THEN
      WRITE(*,*) "Free body motion could not be prevented!"
      CALL CMISSFinalise(Err)
    STOP
    ENDIF
  ENDIF

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsSolid,Err)

  !
  !================================================================================================================================
  !

  !BCs Darcy
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsDarcy,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)

    !In 'generated_mesh_routines.f90/GENERATED_MESH_ELLIPSOID_SURFACE_GET' there is a bug:
    !  BASIS=>ELLIPSOID_MESH%BASES(MESH_COMPONENT)%PTR does not account for the fact that:
    !  in 'generated_mesh_routines.f90/GENERATED_MESH_ELLIPSOID_CREATE_FINISH' the following is done:
    !  CALL MESH_NUMBER_OF_COMPONENTS_SET(GENERATED_MESH%MESH,SIZE(ELLIPSOID_MESH%BASES)/2,ERR,ERROR,*999)
    !Temporary work around, until bug fix:

    !MeshComponentNumber_dummy = DarcyVelMeshComponentNumber
    MeshComponentNumber_dummy = 3 

    !  I N N E R   S U R F A C E
    CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,MeshComponentNumber_dummy,CMISS_GENERATED_MESH_ELLIPSOID_INNER_SURFACE, &
      & InnerSurfaceNodesDarcyVel,InnerNormalXi,Err)

    write(*,*)'InnerSurfaceNodesDarcyVel = ',InnerSurfaceNodesDarcyVel
    write(*,*)'InnerNormalXi = ',InnerNormalXi

    !Set all inner surface nodes impermeable
    !MIND: CMISS_FIELD_DELVDELN_VARIABLE_TYPE -> RHS invoked in DARCY_EQUATION_FINITE_ELEMENT_CALCULATE
    !      CMISS_BOUNDARY_CONDITION_IMPERMEABLE_WALL
    DO NN=1,SIZE(InnerSurfaceNodesDarcyVel,1)
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 1
!       CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1,InnerSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_IMPERMEABLE_WALL,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)
! 
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 2PLACE_EQUATION_EQUATIONS_SET_FLUID_PRESSURE_SETUP
!       CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1,InnerSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_IMPERMEABLE_WALL,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)
! 
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 3
!       CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1,InnerSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_IMPERMEABLE_WALL,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)

!       VALUE = 1.0_CMISSDP
!       COMPONENT_NUMBER = ABS(InnerNormalXi)
!       CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1,InnerSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_IMPERMEABLE_WALL,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)

      NODE_NUMBER = InnerSurfaceNodesDarcyVel(NN)
      COMPONENT_NUMBER = ABS(InnerNormalXi)
      VALUE = 1.0_CMISSDP
      CALL CMISSField_ParameterSetUpdateNode(IndependentFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
        & CMISS_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO


    !  O U T E R   S U R F A C E
    CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,MeshComponentNumber_dummy,CMISS_GENERATED_MESH_ELLIPSOID_OUTER_SURFACE, &
      & OuterSurfaceNodesDarcyVel,OuterNormalXi,Err)

    write(*,*)'OuterSurfaceNodesDarcyVel = ',OuterSurfaceNodesDarcyVel
    write(*,*)'OuterNormalXi = ',OuterNormalXi

    !Set all outer surface nodes impermeable
    DO NN=1,SIZE(OuterSurfaceNodesDarcyVel,1)
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 1
!       CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_IMPERMEABLE_WALL,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)
! 
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 2
!       CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_IMPERMEABLE_WALL,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)
! 
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 3
!       CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_IMPERMEABLE_WALL,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)

!       VALUE = 1.0_CMISSDP
!       COMPONENT_NUMBER = ABS(OuterNormalXi)
!       CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_IMPERMEABLE_WALL,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)

      NODE_NUMBER = OuterSurfaceNodesDarcyVel(NN)
      COMPONENT_NUMBER = ABS(OuterNormalXi)
      VALUE = 1.0_CMISSDP
      CALL CMISSField_ParameterSetUpdateNode(IndependentFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
        & CMISS_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO


    !  T O P   S U R F A C E
    CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,MeshComponentNumber_dummy,CMISS_GENERATED_MESH_ELLIPSOID_TOP_SURFACE, &
      & TopSurfaceNodesDarcyVel,TopNormalXi,Err)

    write(*,*)'TopSurfaceNodesDarcyVel = ',TopSurfaceNodesDarcyVel
    write(*,*)'TopNormalXi = ',TopNormalXi

    !Set all top surface nodes to Darcy inflow BC
    DO NN=1,SIZE(TopSurfaceNodesDarcyVel,1)
      VALUE = -1.0_CMISSDP
      COMPONENT_NUMBER = 3
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMISS_FIELD_V_VARIABLE_TYPE,1,1, &
        & TopSurfaceNodesDarcyVel(NN),COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING TOP DARCY BC TO NODE", TopSurfaceNodesDarcyVel(NN)
    ENDDO


  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsDarcy,Err)

  !
  !================================================================================================================================
  !

  !Solve problem
  WRITE(*,'(A)') "Solving problem..."
  CALL CMISSProblem_Solve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"

  !
  !================================================================================================================================
  !

  !Output solution  
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"QuadraticEllipsoidDrivenDarcy","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"QuadraticEllipsoidDrivenDarcy","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)

  !
  !================================================================================================================================
  !

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM QUADRATICELLIPSOIDDRIVENDARCYEXAMPLE

