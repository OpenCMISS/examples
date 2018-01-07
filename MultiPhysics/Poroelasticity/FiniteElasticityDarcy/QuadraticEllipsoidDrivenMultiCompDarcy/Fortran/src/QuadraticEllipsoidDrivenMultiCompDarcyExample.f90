!> \file
!> \author Christian Michler
!> \brief This is an example program to solve a coupled Finite Elastiticity Darcy equation on a cylindrical geometry using OpenCMISS calls.
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

PROGRAM QUADRATICELLIPSOIDDRIVENMULTICOMPDARCYEXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OpenCMISS
  USE OpenCMISS_Iron
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

  REAL(CMISSRP), PARAMETER :: PI=3.14159_CMISSRP
  REAL(CMISSRP), PARAMETER :: LONG_AXIS=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: SHORT_AXIS=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WALL_THICKNESS=0.5_CMISSRP
  REAL(CMISSRP), PARAMETER :: CUTOFF_ANGLE=1.5708_CMISSRP
  REAL(CMISSRP), PARAMETER :: FIBRE_SLOPE_INTERSECTION=1.73205_CMISSRP !Slope of fibres in base endocardium = 60 degrees
  REAL(CMISSRP), PARAMETER :: FIBRE_SLOPE_CHANGE=-3.4641_CMISSRP !Slope change of fibres from 60 to -60 degrees in transmural direction 
  REAL(CMISSRP), PARAMETER :: SHEET_SLOPE_BASE_ENDO=1.0_CMISSRP !Slope of sheet at base endocardium 

  REAL(CMISSRP), PARAMETER :: INNER_PRESSURE=2.0_CMISSRP  !Positive is compressive
  REAL(CMISSRP), PARAMETER :: OUTER_PRESSURE=0.0_CMISSRP !Positive is compressive
  REAL(CMISSRP), PARAMETER :: C1= 2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: C2= 6.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: C3=10.0_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: NumberGlobalXElements=4  ! X ==NUMBER_GLOBAL_CIRCUMFERENTIAL_ELEMENTS
  INTEGER(CMISSIntg), PARAMETER :: NumberGlobalYElements=4  ! Y ==NUMBER_GLOBAL_LONGITUDINAL_ELEMENTS
  INTEGER(CMISSIntg), PARAMETER :: NumberGlobalZElements=1  ! Z ==NUMBER_GLOBAL_TRANSMURAL_ELEMENTS
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

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentSolidUserNumber=4
  INTEGER(CMISSIntg) :: FieldDependentSolidNumberOfVariables
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentSolidNumberOfComponents=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentFluidNumberOfComponents=4

  INTEGER(CMISSIntg), PARAMETER :: EquationSetSolidUserNumber=55
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldSolidUserNumber=25
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1
  INTEGER(CMISSIntg) :: EquationsSetFieldUserNumberDarcy
  INTEGER(CMISSIntg) :: icompartment,Ncompartments,num_var,componentnum,Nparams

  INTEGER(CMISSIntg) :: MaterialsFieldUserNumberDarcy
  INTEGER(CMISSIntg) :: EquationsSetUserNumberDarcy


  INTEGER(CMISSIntg), PARAMETER :: ControlLoopSolidNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopFluidNumber=2
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopSubiterationNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverSolidIndex=1
  INTEGER(CMISSIntg), PARAMETER :: SolverDarcyIndex=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPorosity=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPermOverVis=2
  INTEGER(CMISSIntg), PARAMETER :: SolidDisplMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolidLagrMultMeshComponentNumber=2
  INTEGER(CMISSIntg), PARAMETER :: SolidGeometryMeshComponentNumber=SolidDisplMeshComponentNumber
  !Program types


  !Program variables

  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  REAL(CMISSRP) :: FibreFieldAngle(3) 
  REAL(CMISSRP) :: nu,theta,omega,XI3,XI3delta,XI2delta, zero
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
  REAL(CMISSRP) :: XCoord,YCoord,ZCoord
  LOGICAL :: X_FIXED,Y_FIXED,X_OKAY,Y_OKAY

  INTEGER(CMISSIntg) :: GeometricFieldDarcyMeshComponentNumber, DarcyVelMeshComponentNumber, DarcyMassIncreaseMeshComponentNumber
  INTEGER(CMISSIntg) :: MeshComponentNumber_dummy

  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS

  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE

  INTEGER(CMISSIntg) :: EQUATIONS_DARCY_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_DARCY_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_DARCY_OUTPUT_TYPE

  REAL(CMISSRP) :: COORD_X, COORD_Y, COORD_Z
  REAL(CMISSRP) :: DOMAIN_X1, DOMAIN_X2, DOMAIN_Y1, DOMAIN_Y2, DOMAIN_Z1, DOMAIN_Z2
  REAL(CMISSRP) :: GEOMETRY_TOLERANCE
  INTEGER(CMISSIntg) :: EDGE_COUNT
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SOLID
  REAL(CMISSRP) :: INITIAL_FIELD_DARCY(4)
  REAL(CMISSRP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSRP) :: RELATIVE_TOLERANCE
  REAL(CMISSRP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSRP) :: LINESEARCH_ALPHA
  REAL(CMISSRP) :: VALUE
  REAL(CMISSRP) :: POROSITY_PARAM_DARCY, PERM_OVER_VIS_PARAM_DARCY

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_DARCY_DIRECT_FLAG

  !CMISS variables

  TYPE(cmfe_BasisType) :: QuadraticBasis,QuadraticCollapsedBasis,LinearBasis,LinearCollapsedBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSetSolid
  TYPE(cmfe_FieldType) :: GeometricFieldSolid,FibreFieldSolid,MaterialFieldSolid
  TYPE(cmfe_FieldType) :: DependentFieldSolid,EquationsSetFieldSolid

  TYPE(cmfe_FieldType) :: GeometricFieldDarcy

  TYPE(cmfe_FieldType), ALLOCATABLE, DIMENSION(:) :: MaterialsFieldDarcy
  TYPE(cmfe_FieldType), ALLOCATABLE, DIMENSION(:) :: EquationsSetFieldDarcy
  !Boundary conditions
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsDarcy
  !Equations sets
  TYPE(cmfe_EquationsSetType), ALLOCATABLE, DIMENSION(:) :: EquationsSetDarcy
  !Equations
  TYPE(cmfe_EquationsType), ALLOCATABLE, DIMENSION(:) :: EquationsDarcy



  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: SolverSolid,LinearSolverSolid
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsSolid
  TYPE(cmfe_ControlLoopType) :: ControlLoop

  !Solvers
  TYPE(cmfe_SolverType) :: DynamicSolverDarcy
  TYPE(cmfe_SolverType) :: LinearSolverDarcy
  !Solver equations
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsDarcy

  REAL(CMISSRP) :: DYNAMIC_SOLVER_DARCY_START_TIME
  REAL(CMISSRP) :: DYNAMIC_SOLVER_DARCY_STOP_TIME
  REAL(CMISSRP) :: DYNAMIC_SOLVER_DARCY_THETA
  REAL(CMISSRP) :: DYNAMIC_SOLVER_DARCY_TIME_INCREMENT

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: Err
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:) :: VariableTypes
  REAL(CMISSRP), ALLOCATABLE, DIMENSION(:,:) :: CouplingCoeffs,ConstitutiveParams
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
  NUMBER_OF_DIMENSIONS=3_CMISSIntg
  INITIAL_FIELD_DARCY(1)=0.0_CMISSRP
  INITIAL_FIELD_DARCY(2)=0.0_CMISSRP
  INITIAL_FIELD_DARCY(3)=0.0_CMISSRP
  INITIAL_FIELD_DARCY(4)=0.0_CMISSRP
  !Set material parameters
  POROSITY_PARAM_DARCY=0.1_CMISSRP
  PERM_OVER_VIS_PARAM_DARCY=1.0_CMISSRP
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE=CMFE_SOLVER_PROGRESS_OUTPUT
  LINEAR_SOLVER_DARCY_OUTPUT_TYPE=CMFE_SOLVER_SOLVER_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DARCY_OUTPUT=CMFE_EQUATIONS_NO_OUTPUT

  !Set time parameter
  DYNAMIC_SOLVER_DARCY_START_TIME=0.0_CMISSRP
  DYNAMIC_SOLVER_DARCY_TIME_INCREMENT=1.0e-3_CMISSRP
  DYNAMIC_SOLVER_DARCY_STOP_TIME=5_CMISSIntg * DYNAMIC_SOLVER_DARCY_TIME_INCREMENT
  DYNAMIC_SOLVER_DARCY_THETA=1.0_CMISSRP !2.0_CMISSRP/3.0_CMISSRP
  !Set result output parameter
  DYNAMIC_SOLVER_DARCY_OUTPUT_FREQUENCY=1
  !Set solver parameters
  LINEAR_SOLVER_DARCY_DIRECT_FLAG=.TRUE.
  RELATIVE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-05_CMISSRP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-10_CMISSRP
  DIVERGENCE_TOLERANCE=1.0E5_CMISSRP !default: 1.0E5
  MAXIMUM_ITERATIONS=10000_CMISSIntg !default: 100000
  RESTART_VALUE=30_CMISSIntg !default: 30
  LINESEARCH_ALPHA=1.0_CMISSRP
  icompartment =1_CMISSIntg
  Ncompartments=4_CMISSIntg

  !LinearMeshComponentNumber/QuadraticMeshComponentNumber
  DarcyVelMeshComponentNumber = LinearMeshComponentNumber
  DarcyMassIncreaseMeshComponentNumber = LinearMeshComponentNumber
!   GeometricFieldDarcyMeshComponentNumber = DarcyVelMeshComponentNumber
  GeometricFieldDarcyMeshComponentNumber = QuadraticMeshComponentNumber

  ALLOCATE (EquationsSetDarcy(Ncompartments))
  ALLOCATE (EquationsSetFieldDarcy(Ncompartments))
  ALLOCATE (MaterialsFieldDarcy(Ncompartments))
  ALLOCATE (EquationsDarcy(Ncompartments))
  !
  !================================================================================================================================
  !

  !Intialise cmiss
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
  CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_FINITE_ELEMENT_CALCULATE"],Err)

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
!   CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  NumberOfDomains=NumberOfComputationalNodes

  !Create a CS - default is 3D rectangular cartesian CS with 0,0,0 as origin
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_TypeSet(CoordinateSystem,CMFE_COORDINATE_RECTANGULAR_CARTESIAN_TYPE,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL cmfe_CoordinateSystem_OriginSet(CoordinateSystem,[0.0_CMISSRP,0.0_CMISSRP,0.0_CMISSRP],Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !Create a region and assign the CS to the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !Define basis functions - tri-linear Lagrange and tri-Quadratic Lagrange, each with collapsed variant
    !Quadratic Basis
  CALL cmfe_Basis_Initialise(QuadraticBasis,Err)
  CALL cmfe_Basis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL cmfe_Basis_InterpolationXiSet(QuadraticBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
!     & [CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err) !Have to do this
  CALL cmfe_Basis_CreateFinish(QuadraticBasis,Err)

    !Collapsed Quadratic Basis
  CALL cmfe_Basis_Initialise(QuadraticCollapsedBasis,Err)
  CALL cmfe_Basis_CreateStart(QuadraticCollapsedBasisUserNumber,QuadraticCollapsedBasis,Err)
  CALL cmfe_Basis_TypeSet(QuadraticCollapsedBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(QuadraticCollapsedBasis,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(QuadraticCollapsedBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
       & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_CollapsedXiSet(QuadraticCollapsedBasis,[CMFE_BASIS_XI_COLLAPSED, &
       & CMFE_BASIS_COLLAPSED_AT_XI0,CMFE_BASIS_NOT_COLLAPSED],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadraticCollapsedBasis, &
       & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)  
!     & [CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(QuadraticCollapsedBasis,.true.,Err) !Have to do this
  CALL cmfe_Basis_CreateFinish(QuadraticCollapsedBasis,Err)

    !Linear Basis
  CALL cmfe_Basis_Initialise(LinearBasis,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis, &
    & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
!     & [CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err) !Have to do this (unused) due to field_interp setup
  CALL cmfe_Basis_CreateFinish(LinearBasis,Err)

    !Collapsed Linear Basis
  CALL cmfe_Basis_Initialise(LinearCollapsedBasis,Err)
  CALL cmfe_Basis_CreateStart(LinearCollapsedBasisUserNumber,LinearCollapsedBasis,Err)
  CALL cmfe_Basis_TypeSet(LinearCollapsedBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(LinearCollapsedBasis,NumberOfXiCoordinates,Err)
  CALL cmfe_Basis_InterpolationXiSet(LinearCollapsedBasis,[CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
       & CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_CollapsedXiSet(LinearCollapsedBasis,[CMFE_BASIS_XI_COLLAPSED,CMFE_BASIS_COLLAPSED_AT_XI0, &
    & CMFE_BASIS_NOT_COLLAPSED],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearCollapsedBasis, &
       & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
!     & [CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(LinearCollapsedBasis,.true.,Err) !Have to do this (unused) due to field_interp setup
  CALL cmfe_Basis_CreateFinish(LinearCollapsedBasis,Err)

  !
  !================================================================================================================================
  !

  !Start the creation of a generated ellipsoid mesh
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up an ellipsoid mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_ELLIPSOID_MESH_TYPE,Err)
  !Set the quadratic and linear bases
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,[QuadraticBasis,QuadraticCollapsedBasis,LinearBasis,LinearCollapsedBasis],Err)
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[LONG_AXIS,SHORT_AXIS,WALL_THICKNESS,CUTOFF_ANGLE],Err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements],Err)
  
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !
  !================================================================================================================================
  !

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !
  !================================================================================================================================
  !

  ! --- GeometricFieldSolid ---
  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(GeometricFieldSolid,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumberSolid,Region,GeometricFieldSolid,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricFieldSolid,Decomposition,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldSolid,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(GeometricFieldSolid,FieldGeometryNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_CreateFinish(GeometricFieldSolid,Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricFieldSolid,Err)

  !
  !================================================================================================================================
  !

  ! --- GeometricFieldDarcy ---
  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(GeometricFieldDarcy,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumberDarcy,Region,GeometricFieldDarcy,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricFieldDarcy,Decomposition,Err)
  CALL cmfe_Field_TypeSet(GeometricFieldDarcy,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(GeometricFieldDarcy,FieldGeometryNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,FieldGeometryNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricFieldDarcyMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,2, &
    & GeometricFieldDarcyMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,3, &
    & GeometricFieldDarcyMeshComponentNumber,Err)
  CALL cmfe_Field_CreateFinish(GeometricFieldDarcy,Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricFieldDarcy,Err)

  !
  !================================================================================================================================
  !

  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(FibreFieldSolid,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,Region,FibreFieldSolid,Err)
  CALL cmfe_Field_TypeSet(FibreFieldSolid,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreFieldSolid,Decomposition,Err)        
  CALL cmfe_Field_GeometricFieldSet(FibreFieldSolid,GeometricFieldSolid,Err)
  CALL cmfe_Field_NumberOfVariablesSet(FibreFieldSolid,FieldFibreNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(FibreFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(FibreFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_CreateFinish(FibreFieldSolid,Err)

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
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      FibreFieldAngle=[zero,zero,zero] 
      DO component_idx=1,FieldFibreNumberOfComponents
        CALL cmfe_Field_ParameterSetUpdateNode(FibreFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
          & DerivativeUserNumber,node_idx,component_idx,FibreFieldAngle(component_idx),Err)
      ENDDO
    ENDIF
    theta=atan(FIBRE_SLOPE_CHANGE*XI3+FIBRE_SLOPE_INTERSECTION)
    DO j=2, TOTAL_NUMBER_NODES_XI(2) 
      nu=PI-XI2delta*(j-1)
      omega=PI/2+cos(2*nu)*atan(SHEET_SLOPE_BASE_ENDO)*(-2*XI3+1)
      DO i=1, TOTAL_NUMBER_NODES_XI(1)
        node_idx=node_idx+1
        CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          FibreFieldAngle=[theta,zero,omega]
          DO component_idx=1,FieldFibreNumberOfComponents
            CALL cmfe_Field_ParameterSetUpdateNode(FibreFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
              & DerivativeUserNumber, node_idx,component_idx,FibreFieldAngle(component_idx),Err)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    XI3=XI3+XI3delta
  ENDDO

  !Create a material field and attach it to the geometric field
  CALL cmfe_Field_Initialise(MaterialFieldSolid,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumber,Region,MaterialFieldSolid,Err)
  CALL cmfe_Field_TypeSet(MaterialFieldSolid,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialFieldSolid,Decomposition,Err)        
  CALL cmfe_Field_GeometricFieldSet(MaterialFieldSolid,GeometricFieldSolid,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialFieldSolid,FieldMaterialNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_CreateFinish(MaterialFieldSolid,Err)

  !Set Mooney-Rivlin constants 
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,C1,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,C2,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,C3,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS
  DO icompartment = 1,Ncompartments
    EquationsSetFieldUserNumberDarcy = 100_CMISSIntg+icompartment
    EquationsSetUserNumberDarcy = 200_CMISSIntg+icompartment
  !Create the equations set for ALE Darcy
    CALL cmfe_Field_Initialise(EquationsSetFieldDarcy(icompartment),Err)
    CALL cmfe_EquationsSet_Initialise(EquationsSetDarcy(icompartment),Err)
    CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberDarcy,Region,GeometricFieldDarcy, &
      & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMFE_EQUATIONS_SET_DARCY_EQUATION_TYPE, &
      & CMFE_EQUATIONS_SET_INCOMPRESS_ELAST_MULTI_COMP_DARCY_SUBTYPE],EquationsSetFieldUserNumberDarcy, &
      & EquationsSetFieldDarcy(icompartment),EquationsSetDarcy(icompartment),Err)
    !Finish creating the equations set
    CALL cmfe_EquationsSet_CreateFinish(EquationsSetDarcy(icompartment),Err)
    !Set the values for the equations set field to be the current compartment number (1 - N), and the total number of compartments (N)
    CALL cmfe_Field_ParameterSetUpdateConstant(EquationsSetFieldDarcy(icompartment),CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE,1,icompartment,Err)
    CALL cmfe_Field_ParameterSetUpdateConstant(EquationsSetFieldDarcy(icompartment),CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE,2,Ncompartments,Err)
  ENDDO

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create the equations_set
  CALL cmfe_Field_Initialise(EquationsSetFieldSolid,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetSolid,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetSolidUserNumber,Region,FibreFieldSolid,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_INCOMPRESS_ELAST_MULTI_COMP_DARCY_SUBTYPE], &
    & EquationsSetFieldSolidUserNumber,EquationsSetFieldSolid,EquationsSetSolid,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetSolid,Err)
  !Set the values for the equations set field to be the current compartment number (O for the finite elasticity equations_set), and the total number of compartments (N)
  !Need to store number of compartments, as finite elasticity uses this to calculate the total mass increase for the constiutive law
  CALL cmfe_Field_ParameterSetUpdateConstant(EquationsSetFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE, &
     & CMFE_FIELD_VALUES_SET_TYPE,1,0_CMISSIntg,Err)
  CALL cmfe_Field_ParameterSetUpdateConstant(EquationsSetFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE, &
     & CMFE_FIELD_VALUES_SET_TYPE,2,Ncompartments,Err)

  !
  !================================================================================================================================
  !
  ! Solid

  !Create a dependent field with two variables and four components
  CALL cmfe_Field_Initialise(DependentFieldSolid,Err)
  !
  CALL cmfe_Field_CreateStart(FieldDependentSolidUserNumber,Region,DependentFieldSolid,Err)
  !
  CALL cmfe_Field_TypeSet(DependentFieldSolid,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentFieldSolid,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentFieldSolid,GeometricFieldSolid,Err)
  CALL cmfe_Field_DependentTypeSet(DependentFieldSolid,CMFE_FIELD_DEPENDENT_TYPE,Err)
  !Create 2N+2 number of variables - 2 for solid, 2N for N Darcy compartments
  FieldDependentSolidNumberOfVariables=2*Ncompartments+2
  CALL cmfe_Field_NumberOfVariablesSet(DependentFieldSolid,FieldDependentSolidNumberOfVariables,Err)
  !create two variables for each compartment
  ALLOCATE(VariableTypes(2*Ncompartments+2))
  DO num_var=1,Ncompartments+1
     VariableTypes(2*num_var-1)=CMFE_FIELD_U_VARIABLE_TYPE+(CMFE_FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
     VariableTypes(2*num_var)=CMFE_FIELD_DELUDELN_VARIABLE_TYPE+(CMFE_FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
  ENDDO
  CALL cmfe_Field_VariableTypesSet(DependentFieldSolid,VariableTypes,Err) 
!   CALL cmfe_Field_VariableTypesSet(DependentFieldSolid,[CMFE_FIELD_U_VARIABLE_TYPE, &
!     & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_DELVDELN_VARIABLE_TYPE],Err)
    CALL cmfe_Field_DimensionSet(DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE, &
       & CMFE_FIELD_VECTOR_DIMENSION_TYPE,Err)
    CALL cmfe_Field_DimensionSet(DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE, &
       & CMFE_FIELD_VECTOR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,FieldDependentSolidNumberOfComponents,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE, &
    & FieldDependentSolidNumberOfComponents,Err)
  DO icompartment=3,2*Ncompartments+2
    CALL cmfe_Field_NumberOfComponentsSet(DependentFieldSolid,VariableTypes(icompartment),FieldDependentFluidNumberOfComponents,Err)
  ENDDO
!   CALL cmfe_Field_ComponentInterpolationSet(DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
!   CALL cmfe_Field_ComponentInterpolationSet(DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
!   CALL cmfe_Field_ComponentInterpolationSet(DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)

  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,1,SolidDisplMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,2,SolidDisplMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,3,SolidDisplMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentInterpolationSet(DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,4, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION, &
    & Err)
!   CALL cmfe_Field_ComponentInterpolationSet(DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
!   CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,4,SolidMeshComponenetNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,4,SolidLagrMultMeshComponentNumber,Err)
!   CALL cmfe_Field_ComponentInterpolationSet(DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1, &
!     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
!   CALL cmfe_Field_ComponentInterpolationSet(DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2, &
!     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
!   CALL cmfe_Field_ComponentInterpolationSet(DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3, &
!     & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1, &
    & SolidDisplMeshComponentNumber, &
    & Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2, &
    & SolidDisplMeshComponentNumber, &
    & Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3, &
    & SolidDisplMeshComponentNumber, &
    & Err)
  CALL cmfe_Field_ComponentInterpolationSet(DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
!   CALL cmfe_Field_ComponentInterpolationSet(DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
!     & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
!   CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,SolidMeshComponenetNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
    & SolidLagrMultMeshComponentNumber, &
    & Err)
  !loop over the number of compartments
  DO icompartment=3,2*Ncompartments+2
!     CALL cmfe_Field_DimensionSet(DependentFieldSolid,VariableTypes(icompartment), &
!        & CMFE_FIELD_VECTOR_DIMENSION_TYPE,Err)
    !CALL cmfe_Field_NumberOfComponentsSet(DependentFieldSolid,VariableTypes(icompartment),FieldDependentFluidNumberOfComponents,Err)
    DO componentnum=1,FieldDependentFluidNumberOfComponents-1
    !set dimension type
!     CALL cmfe_Field_DimensionSet(DependentField,VariableTypes(icompartment), &
!        & CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
      CALL cmfe_Field_ComponentInterpolationSet(DependentFieldSolid,VariableTypes(icompartment),componentnum, &
       & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
      CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,VariableTypes(icompartment),componentnum, & 
         & DarcyVelMeshComponentNumber,Err)
    ENDDO
      CALL cmfe_Field_ComponentInterpolationSet(DependentFieldSolid,VariableTypes(icompartment), &
       & FieldDependentFluidNumberOfComponents, &
       & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
!     CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,VariableTypes(icompartment), &
!       & FieldDependentFluidNumberOfComponents,MESH_COMPONENT_NUMBER_PRESSURE,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,VariableTypes(icompartment), &
      & FieldDependentFluidNumberOfComponents,DarcyMassIncreaseMeshComponentNumber,Err)
    
  ENDDO

!   CALL cmfe_Field_NumberOfComponentsSet(DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,FieldDependentFluidNumberOfComponents,Err)
!   CALL cmfe_Field_NumberOfComponentsSet(DependentFieldSolid,CMFE_FIELD_DELVDELN_VARIABLE_TYPE,FieldDependentFluidNumberOfComponents,Err)
!   !For this equation type, MESH_COMPONENT_NUMBER_PRESSURE is actually the mass increase component as the pressure is taken from the solid equations
!   CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,3,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,4,MESH_COMPONENT_NUMBER_PRESSURE,Err)
!   CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_DELVDELN_VARIABLE_TYPE,1,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_DELVDELN_VARIABLE_TYPE,2,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_DELVDELN_VARIABLE_TYPE,3,MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldSolid,CMFE_FIELD_DELVDELN_VARIABLE_TYPE,4,MESH_COMPONENT_NUMBER_PRESSURE,Err)

  !
  CALL cmfe_Field_ScalingTypeSet(DependentFieldSolid,CMFE_FIELD_UNIT_SCALING,Err)

  CALL cmfe_Field_CreateFinish(DependentFieldSolid,Err)

  !
  !================================================================================================================================
  !

  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetSolid,FieldDependentSolidUserNumber,DependentFieldSolid,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetSolid,Err)

  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetSolid,FieldMaterialUserNumber,MaterialFieldSolid,Err)  
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !
  DO icompartment = 1,Ncompartments
    CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetDarcy(icompartment),FieldDependentSolidUserNumber,&
      & DependentFieldSolid,Err)
    CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetDarcy(icompartment),Err)
  ENDDO

  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS+1
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldSolid,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldSolid,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldSolid,CMFE_FIELD_U3_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
  ENDDO
  !
  !================================================================================================================================
  !
  ALLOCATE(CouplingCoeffs(Ncompartments,Ncompartments))
  IF(Ncompartments==2)THEN
    CouplingCoeffs(1,1)=0.0E-01_CMISSRP
!     CouplingCoeffs(1,2)=-1.0E-04_CMISSRP
!     CouplingCoeffs(2,1)=-1.0E-04_CMISSRP
    CouplingCoeffs(1,2)=0.0E-01_CMISSRP
    CouplingCoeffs(2,1)=0.0E-01_CMISSRP
    CouplingCoeffs(2,2)=0.0E-01_CMISSRP
  ELSE IF(Ncompartments==3)THEN
    CouplingCoeffs(1,1)=1.0E-02_CMISSRP
    CouplingCoeffs(1,2)=1.0E-02_CMISSRP
    CouplingCoeffs(1,3)=0.0E-02_CMISSRP
    CouplingCoeffs(2,1)=1.0E-02_CMISSRP
    CouplingCoeffs(2,2)=2.0E-02_CMISSRP
    CouplingCoeffs(2,3)=1.0E-02_CMISSRP
    CouplingCoeffs(3,1)=0.0E-02_CMISSRP
    CouplingCoeffs(3,2)=1.0E-02_CMISSRP
    CouplingCoeffs(3,3)=1.0E-02_CMISSRP
  ELSE IF(Ncompartments==4)THEN
    CouplingCoeffs(1,1)=0.0E-02_CMISSRP
    CouplingCoeffs(1,2)=0.0E-02_CMISSRP
    CouplingCoeffs(1,3)=0.0E-02_CMISSRP
    CouplingCoeffs(1,4)=0.0E-02_CMISSRP
    CouplingCoeffs(2,1)=0.0E-02_CMISSRP
    CouplingCoeffs(2,2)=0.0E-02_CMISSRP
    CouplingCoeffs(2,3)=0.0E-02_CMISSRP
    CouplingCoeffs(2,4)=0.0E-02_CMISSRP
    CouplingCoeffs(3,1)=0.0E-02_CMISSRP
    CouplingCoeffs(3,2)=0.0E-02_CMISSRP
    CouplingCoeffs(3,3)=0.0E-02_CMISSRP
    CouplingCoeffs(3,4)=0.0E-02_CMISSRP
    CouplingCoeffs(4,1)=0.0E-02_CMISSRP
    CouplingCoeffs(4,2)=0.0E-02_CMISSRP
    CouplingCoeffs(4,3)=0.0E-02_CMISSRP
    CouplingCoeffs(4,4)=0.0E-02_CMISSRP
  ELSE IF(Ncompartments==5)THEN
    CouplingCoeffs(1,1)=0.0E-02_CMISSRP
    CouplingCoeffs(1,2)=0.0E-02_CMISSRP
    CouplingCoeffs(1,3)=0.0E-02_CMISSRP
    CouplingCoeffs(1,4)=0.0E-02_CMISSRP
    CouplingCoeffs(1,5)=0.0E-02_CMISSRP
    CouplingCoeffs(2,1)=0.0E-02_CMISSRP
    CouplingCoeffs(2,2)=0.0E-02_CMISSRP
    CouplingCoeffs(2,3)=0.0E-02_CMISSRP
    CouplingCoeffs(2,4)=0.0E-02_CMISSRP
    CouplingCoeffs(2,5)=0.0E-02_CMISSRP
    CouplingCoeffs(3,1)=0.0E-02_CMISSRP
    CouplingCoeffs(3,2)=0.0E-02_CMISSRP
    CouplingCoeffs(3,3)=0.0E-02_CMISSRP
    CouplingCoeffs(3,4)=0.0E-02_CMISSRP
    CouplingCoeffs(3,5)=0.0E-02_CMISSRP
    CouplingCoeffs(4,1)=0.0E-02_CMISSRP
    CouplingCoeffs(4,2)=0.0E-02_CMISSRP
    CouplingCoeffs(4,3)=0.0E-02_CMISSRP
    CouplingCoeffs(4,4)=0.0E-02_CMISSRP
    CouplingCoeffs(4,5)=0.0E-02_CMISSRP
    CouplingCoeffs(5,1)=0.0E-02_CMISSRP
    CouplingCoeffs(5,2)=0.0E-02_CMISSRP
    CouplingCoeffs(5,3)=0.0E-02_CMISSRP
    CouplingCoeffs(5,4)=0.0E-02_CMISSRP
    CouplingCoeffs(5,5)=0.0E-02_CMISSRP
  ELSE
    write(*,*) "Can't initialise coupling coefficients array."
  ENDIF
  !Define the material parameters for each compartments' constitutive law (for determining pressure)
  Nparams=3
  ALLOCATE(ConstitutiveParams(Ncompartments,Nparams))
  IF(Ncompartments==2)THEN
    ConstitutiveParams(1,1)=10.0E-01_CMISSRP
    ConstitutiveParams(1,2)=10.0E-01_CMISSRP
    ConstitutiveParams(1,3)=10.0E-01_CMISSRP
    ConstitutiveParams(2,1)=10.0E-01_CMISSRP
    ConstitutiveParams(2,2)=10.0E-01_CMISSRP
    ConstitutiveParams(2,3)=10.0E-01_CMISSRP
  ELSE IF(Ncompartments==3)THEN
    ConstitutiveParams(1,1)=1.0E-02_CMISSRP
    ConstitutiveParams(1,2)=1.0E-02_CMISSRP
    ConstitutiveParams(1,3)=0.0E-02_CMISSRP
    ConstitutiveParams(2,1)=1.0E-02_CMISSRP
    ConstitutiveParams(2,2)=2.0E-02_CMISSRP
    ConstitutiveParams(2,3)=1.0E-02_CMISSRP
    ConstitutiveParams(3,1)=0.0E-02_CMISSRP
    ConstitutiveParams(3,2)=1.0E-02_CMISSRP
    ConstitutiveParams(3,3)=1.0E-02_CMISSRP
  ELSE IF(Ncompartments==4)THEN
    ConstitutiveParams(1,1)=0.0E-02_CMISSRP
    ConstitutiveParams(1,2)=0.0E-02_CMISSRP
    ConstitutiveParams(1,3)=0.0E-02_CMISSRP
    ConstitutiveParams(2,1)=0.0E-02_CMISSRP
    ConstitutiveParams(2,2)=0.0E-02_CMISSRP
    ConstitutiveParams(2,3)=0.0E-02_CMISSRP
    ConstitutiveParams(3,1)=0.0E-02_CMISSRP
    ConstitutiveParams(3,2)=0.0E-02_CMISSRP
    ConstitutiveParams(3,3)=0.0E-02_CMISSRP
    ConstitutiveParams(4,1)=0.0E-02_CMISSRP
    ConstitutiveParams(4,2)=0.0E-02_CMISSRP
    ConstitutiveParams(4,3)=0.0E-02_CMISSRP
  ELSE IF(Ncompartments==5)THEN
    ConstitutiveParams(1,1)=0.0E-02_CMISSRP
    ConstitutiveParams(1,2)=0.0E-02_CMISSRP
    ConstitutiveParams(1,3)=0.0E-02_CMISSRP
    ConstitutiveParams(2,1)=0.0E-02_CMISSRP
    ConstitutiveParams(2,2)=0.0E-02_CMISSRP
    ConstitutiveParams(2,3)=0.0E-02_CMISSRP
    ConstitutiveParams(3,1)=0.0E-02_CMISSRP
    ConstitutiveParams(3,2)=0.0E-02_CMISSRP
    ConstitutiveParams(3,3)=0.0E-02_CMISSRP
    ConstitutiveParams(4,1)=0.0E-02_CMISSRP
    ConstitutiveParams(4,2)=0.0E-02_CMISSRP
    ConstitutiveParams(4,3)=0.0E-02_CMISSRP
    ConstitutiveParams(5,1)=0.0E-02_CMISSRP
    ConstitutiveParams(5,2)=0.0E-02_CMISSRP
    ConstitutiveParams(5,3)=0.0E-02_CMISSRP
  ELSE
    write(*,*) "Can't initialise constitutive parameters array."
  ENDIF
  !MATERIALS FIELDS
  !Auto-created field contains a U variable type to store the diffusion coefficient(s)
  !It also contains a V variable type to store the coupling coefficients 
  DO icompartment = 1,Ncompartments
    MaterialsFieldUserNumberDarcy = 400+icompartment
    CALL cmfe_Field_Initialise(MaterialsFieldDarcy(icompartment),Err)
    CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetDarcy(icompartment),MaterialsFieldUserNumberDarcy,&
         & MaterialsFieldDarcy(icompartment),Err)
    CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetDarcy(icompartment),Err)
  END DO
  DO icompartment = 1,Ncompartments
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE, &
      & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
    CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE, &
      & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)
  END DO
  DO icompartment = 1, Ncompartments
    DO COMPONENT_NUMBER=1, Ncompartments
      CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMFE_FIELD_V_VARIABLE_TYPE, &
         & CMFE_FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
!         CALL cmfe_Field_ParameterSetUpdateConstant(MaterialsFieldDarcy(icompartment),CMFE_FIELD_V_VARIABLE_TYPE, &
!           & CMFE_FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
    END DO
  END DO
  DO icompartment = 1, Ncompartments
    DO COMPONENT_NUMBER=1,Nparams
      CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMFE_FIELD_U1_VARIABLE_TYPE, &
         & CMFE_FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,ConstitutiveParams(icompartment,COMPONENT_NUMBER),Err)
!         CALL cmfe_Field_ParameterSetUpdateConstant(MaterialsFieldDarcy(icompartment),CMFE_FIELD_V_VARIABLE_TYPE, &
!           & CMFE_FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
    END DO
  END DO

  !
  !================================================================================================================================
  !

  !EQUATIONS SET EQUATIONS

  !Darcy
  DO icompartment=1,Ncompartments
    !Create the equations set equations
    CALL cmfe_Equations_Initialise(EquationsDarcy(icompartment),Err)
    CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetDarcy(icompartment),EquationsDarcy(icompartment),Err)
    !Set the equations matrices sparsity type
    CALL cmfe_Equations_SparsityTypeSet(EquationsDarcy(icompartment),CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !   !Set the equations lumping type
  !   CALL cmfe_Equations_LumpingTypeSet(EquationsDarcy,CMFE_EQUATIONS_UNLUMPED_MATRICES,Err)
    !Set the equations set output
    CALL cmfe_Equations_OutputTypeSet(EquationsDarcy(icompartment),EQUATIONS_DARCY_OUTPUT,Err)
  !Finish the equations set equations
    CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetDarcy(icompartment),Err)
  ENDDO

  !Solid
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetSolid,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetSolid,Err)   

  !
  !================================================================================================================================
  !

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
!   CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,-14.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
    & 0.0_CMISSRP,Err)

  !
  !================================================================================================================================
  !
  !PROBLEMS

  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_DARCY_TYPE, &
    & CMFE_PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
!   CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoop,2,Err)
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,DYNAMIC_SOLVER_DARCY_START_TIME,DYNAMIC_SOLVER_DARCY_STOP_TIME, &
    & DYNAMIC_SOLVER_DARCY_TIME_INCREMENT,Err)
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,DYNAMIC_SOLVER_DARCY_OUTPUT_FREQUENCY,Err)
!   CALL cmfe_ControlLoop_OutputTypeSet(ControlLoop,CMFE_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(SolverSolid,Err)
  CALL cmfe_Solver_Initialise(LinearSolverSolid,Err)
  CALL cmfe_Solver_Initialise(DynamicSolverDarcy,Err)
  CALL cmfe_Solver_Initialise(LinearSolverDarcy,Err)

  CALL cmfe_Problem_SolversCreateStart(Problem,Err)

  ! Solid
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopSubiterationNumber,ControlLoopSolidNumber,CMFE_CONTROL_LOOP_NODE], &
    & SolverSolidIndex,SolverSolid,Err)
  CALL cmfe_Solver_OutputTypeSet(SolverSolid,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
!   CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(SolverSolid,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(SolverSolid,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)

!   CALL cmfe_SolverNonLinearTypeSet(SolverSolid,CMFE_SOLVER_NONLINEAR_NEWTON,Err)
!   CALL cmfe_Solver_LibraryTypeSet(SolverSolid,CMFE_SOLVER_PETSC_LIBRARY,Err)

  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(SolverSolid,ABSOLUTE_TOLERANCE,Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(SolverSolid,RELATIVE_TOLERANCE,Err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(SolverSolid,MAXIMUM_ITERATIONS,Err)

  CALL cmfe_Solver_NewtonLinearSolverGet(SolverSolid,LinearSolverSolid,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolverSolid,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)

  !Get the Darcy solver
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMFE_CONTROL_LOOP_NODE], &
    & SolverDarcyIndex,DynamicSolverDarcy,Err)
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(DynamicSolverDarcy,DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE,Err)
  !Set theta
  CALL cmfe_Solver_DynamicThetaSet(DynamicSolverDarcy,DYNAMIC_SOLVER_DARCY_THETA,Err)
  !CALL cmfe_SolverDynamicDynamicSet(DynamicSolverDarcy,.TRUE.,Err)
  !Get the dynamic linear solver
  CALL cmfe_Solver_DynamicLinearSolverGet(DynamicSolverDarcy,LinearSolverDarcy,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_DARCY_DIRECT_FLAG) THEN
    CALL cmfe_Solver_LinearTypeSet(LinearSolverDarcy,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LibraryTypeSet(LinearSolverDarcy,CMFE_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(LinearSolverDarcy,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolverDarcy,MAXIMUM_ITERATIONS,Err)
    CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(LinearSolverDarcy,DIVERGENCE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolverDarcy,RELATIVE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolverDarcy,ABSOLUTE_TOLERANCE,Err)
    CALL cmfe_Solver_LinearIterativeGMRESRestartSet(LinearSolverDarcy,RESTART_VALUE,Err)
  ENDIF

  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)


  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(SolverSolid,Err)
  CALL cmfe_Solver_Initialise(LinearSolverDarcy,Err)

  CALL cmfe_SolverEquations_Initialise(SolverEquationsSolid,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsDarcy,Err)

  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !
  !Get the finite elasticity solver equations
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopSubiterationNumber,ControlLoopSolidNumber,CMFE_CONTROL_LOOP_NODE], &
    & SolverSolidIndex,SolverSolid,Err)
  CALL cmfe_Solver_SolverEquationsGet(SolverSolid,SolverEquationsSolid,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsSolid,CMFE_SOLVER_SPARSE_MATRICES,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsSolid,EquationsSetSolid,EquationsSetIndex,Err)
  !
  !Get the Darcy solver equations
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMFE_CONTROL_LOOP_NODE], &
    & SolverDarcyIndex,LinearSolverDarcy,Err)
  CALL cmfe_Solver_SolverEquationsGet(LinearSolverDarcy,SolverEquationsDarcy,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsDarcy,CMFE_SOLVER_SPARSE_MATRICES,Err)
  DO icompartment=1,Ncompartments
    CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsDarcy,EquationsSetDarcy(icompartment),EquationsSetIndex,Err)
  ENDDO
  !
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsSolid,BoundaryConditions,Err)

  !Grab the list of nodes on inner, outer and top surfaces
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_ELLIPSOID_TOP_SURFACE,TopSurfaceNodes,TopNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_ELLIPSOID_INNER_SURFACE,InnerSurfaceNodes,InnerNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_ELLIPSOID_OUTER_SURFACE,OuterSurfaceNodes,OuterNormalXi,Err)

  write(*,*)'TopSurfaceNodes = ',TopSurfaceNodes

  ! ASSIGN BOUNDARY CONDITIONS
  !Fix base of the ellipsoid in z direction
  DO NN=1,SIZE(TopSurfaceNodes,1)
    NODE=TopSurfaceNodes(NN)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NODE,3, &
        & ZCoord,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE,3, &
        & CMFE_BOUNDARY_CONDITION_FIXED,ZCoord,Err)
    ENDIF
  ENDDO

  !Apply inner surface pressure
  !NOTE: Surface pressure goes into pressure_values_set_type of the DELUDELN type
  DO NN=1,SIZE(InnerSurfaceNodes,1)
    NODE=InnerSurfaceNodes(NN)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NODE, &
        & ABS(InnerNormalXi), &
        & CMFE_BOUNDARY_CONDITION_PRESSURE,INNER_PRESSURE,Err)
    ENDIF
  ENDDO

  !Apply outer surface pressure
  DO NN=1,SIZE(OuterSurfaceNodes,1)
    NODE=OuterSurfaceNodes(NN)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentFieldSolid,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,NODE, &
        & ABS(OuterNormalXi), &
        & CMFE_BOUNDARY_CONDITION_PRESSURE,OUTER_PRESSURE,Err)
    ENDIF
  ENDDO

  !Fix more nodes at the base to stop free body motion
  X_FIXED=.FALSE.
  Y_FIXED=.FALSE.
  DO NN=1,SIZE(TopSurfaceNodes,1)
    NODE=TopSurfaceNodes(NN)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NODE,1, &
        & XCoord,Err)
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NODE,2, &
        & YCoord,Err)
      IF(ABS(XCoord)<1.0E-6_CMISSRP) THEN
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE,1, &
          & CMFE_BOUNDARY_CONDITION_FIXED,XCoord,Err)
        WRITE(*,*) "FIXING NODE",NODE,"IN X DIRECTION"
        X_FIXED=.TRUE.
    ENDIF
      IF(ABS(YCoord)<1.0E-6_CMISSRP) THEN
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE,2, &
          & CMFE_BOUNDARY_CONDITION_FIXED,YCoord,Err)
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
      CALL cmfe_Finalise(Err)
    STOP
    ENDIF
  ENDIF

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsSolid,Err)

  !
  !================================================================================================================================
  !

  !BCs Darcy
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsDarcy,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)
  DO icompartment=1,Ncompartments

    !In 'generated_mesh_routines.f90/GENERATED_MESH_ELLIPSOID_SURFACE_GET' there is a bug:
    !  BASIS=>ELLIPSOID_MESH%BASES(MESH_COMPONENT)%PTR does not account for the fact that:
    !  in 'generated_mesh_routines.f90/GENERATED_MESH_ELLIPSOID_CREATE_FINISH' the following is done:
    !  CALL MESH_NUMBER_OF_COMPONENTS_SET(GENERATED_MESH%MESH,SIZE(ELLIPSOID_MESH%BASES)/2,ERR,ERROR,*999)
    !Temporary work around, until bug fix:

    !MeshComponentNumber_dummy = DarcyVelMeshComponentNumber
    MeshComponentNumber_dummy = 3

!     !inner surface
!     CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,MeshComponentNumber_dummy,CMFE_GENERATED_MESH_ELLIPSOID_INNER_SURFACE, &
!       & InnerSurfaceNodesDarcyVel,InnerNormalXi,Err)
!
!     write(*,*)'InnerSurfaceNodesDarcyVel = ',InnerSurfaceNodesDarcyVel
!
!     !Set all inner surface nodes impermeable
!     DO NN=1,SIZE(InnerSurfaceNodesDarcyVel,1)
! !       VALUE = 0.0_CMISSRP
! !       COMPONENT_NUMBER = 1
! !       CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,1,InnerSurfaceNodesDarcyVel(NN), &
! !         & COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! !       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)
! !
! !       VALUE = 0.0_CMISSRP
! !       COMPONENT_NUMBER = 2
! !       CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,1,InnerSurfaceNodesDarcyVel(NN), &
! !         & COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! !       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)
! !
!       VALUE = 0.0_CMISSRP
!       COMPONENT_NUMBER = 3
! !       CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,1,InnerSurfaceNodesDarcyVel(NN), &
! !         & COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!       CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,1,InnerSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_IMPERMEABLE_WALL,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)
!     ENDDO


!     !outer surface
!     CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,MeshComponentNumber_dummy,CMFE_GENERATED_MESH_ELLIPSOID_OUTER_SURFACE, &
!       & OuterSurfaceNodesDarcyVel,OuterNormalXi,Err)
!
!     write(*,*)'OuterSurfaceNodesDarcyVel = ',OuterSurfaceNodesDarcyVel
!
!     !Set all outer surface nodes impermeable
!     DO NN=1,SIZE(OuterSurfaceNodesDarcyVel,1)
!       VALUE = 0.0_CMISSRP
!       COMPONENT_NUMBER = 1
!       CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)
!
!       VALUE = 0.0_CMISSRP
!       COMPONENT_NUMBER = 2
!       CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)
!
!       VALUE = 0.0_CMISSRP
!       COMPONENT_NUMBER = 3
!       CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)
!     ENDDO


    !top surface
    CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,MeshComponentNumber_dummy,CMFE_GENERATED_MESH_ELLIPSOID_TOP_SURFACE, &
      & TopSurfaceNodesDarcyVel,TopNormalXi,Err)
    IF(icompartment==1) THEN
    write(*,*)'TopSurfaceNodesDarcyVel = ',TopSurfaceNodesDarcyVel

    !Set all top surface nodes to Darcy inflow BC
      DO NN=1,SIZE(TopSurfaceNodesDarcyVel,1)
        VALUE = -0.25_CMISSRP
        COMPONENT_NUMBER = 3
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMFE_FIELD_V_VARIABLE_TYPE,1,1, &
          & TopSurfaceNodesDarcyVel(NN), &
          & COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING TOP DARCY BC TO NODE", TopSurfaceNodesDarcyVel(NN)
      ENDDO
    ELSEIF(icompartment==2)THEN

    !Set all top surface nodes to Darcy inflow BC
      DO NN=1,SIZE(TopSurfaceNodesDarcyVel,1)
        VALUE = -0.25_CMISSRP
        COMPONENT_NUMBER = 3
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMFE_FIELD_U1_VARIABLE_TYPE,1,1, &
          & TopSurfaceNodesDarcyVel(NN), &
          & COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING TOP DARCY BC TO NODE", TopSurfaceNodesDarcyVel(NN)
      ENDDO
    ELSEIF(icompartment==3)THEN
    !Set all top surface nodes to Darcy inflow BC
      DO NN=1,SIZE(TopSurfaceNodesDarcyVel,1)
        VALUE = -0.25_CMISSRP
        COMPONENT_NUMBER = 3
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMFE_FIELD_U2_VARIABLE_TYPE,1,1, &
          & TopSurfaceNodesDarcyVel(NN), &
          & COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING TOP DARCY BC TO NODE", TopSurfaceNodesDarcyVel(NN)
      ENDDO
    ELSEIF(icompartment==4)THEN
    !Set all top surface nodes to Darcy inflow BC
      DO NN=1,SIZE(TopSurfaceNodesDarcyVel,1)
        VALUE = -0.25_CMISSRP
        COMPONENT_NUMBER = 3
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMFE_FIELD_U3_VARIABLE_TYPE,1,1, &
          & TopSurfaceNodesDarcyVel(NN), &
          & COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
        IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING TOP DARCY BC TO NODE", TopSurfaceNodesDarcyVel(NN)
      ENDDO
    ENDIF


  ENDDO
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsDarcy,Err)

  !
  !================================================================================================================================
  !

  !Solve problem
  WRITE(*,'(A)') "Solving problem..."
  CALL cmfe_Problem_Solve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"

  !
  !================================================================================================================================
  !

  !Output solution  
!   CALL cmfe_Fields_Initialise(Fields,Err)
!   CALL cmfe_Fields_Create(Region,Fields,Err)
!   CALL cmfe_Fields_NodesExport(Fields,"QuadraticEllipsoidDrivenMultiCompDarcy","FORTRAN",Err)
!   CALL cmfe_Fields_ElementsExport(Fields,"QuadraticEllipsoidDrivenMultiCompDarcy","FORTRAN",Err)
!   CALL cmfe_Fields_Finalise(Fields,Err)

  !
  !================================================================================================================================
  !

  CALL cmfe_Finalise(Err)

  IF (ALLOCATED(EquationsSetFieldDarcy)) DEALLOCATE(EquationsSetFieldDarcy)
  IF (ALLOCATED(EquationsSetDarcy)) DEALLOCATE(EquationsSetDarcy)
  IF (ALLOCATED(MaterialsFieldDarcy)) DEALLOCATE(MaterialsFieldDarcy)
  IF (ALLOCATED(EquationsDarcy)) DEALLOCATE(EquationsDarcy)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM QUADRATICELLIPSOIDDRIVENMULTICOMPDARCYEXAMPLE

