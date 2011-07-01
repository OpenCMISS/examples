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

PROGRAM QUADRATICELLIPSOIDDRIVENMULTICOMPDARCYEXAMPLE

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
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
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
  TYPE(CMISSFieldType) :: DependentFieldSolid,EquationsSetFieldSolid

  TYPE(CMISSFieldType) :: GeometricFieldDarcy

  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: MaterialsFieldDarcy
  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: EquationsSetFieldDarcy
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsDarcy
  !Equations sets
  TYPE(CMISSEquationsSetType), ALLOCATABLE, DIMENSION(:) :: EquationsSetDarcy
  !Equations
  TYPE(CMISSEquationsType), ALLOCATABLE, DIMENSION(:) :: EquationsDarcy



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

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: Err
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:) :: VariableTypes
  REAL(CMISSDP), ALLOCATABLE, DIMENSION(:,:) :: CouplingCoeffs,ConstitutiveParams
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
  DYNAMIC_SOLVER_DARCY_TIME_INCREMENT=1.0e-3_CMISSDP
  DYNAMIC_SOLVER_DARCY_STOP_TIME=5_CMISSIntg * DYNAMIC_SOLVER_DARCY_TIME_INCREMENT
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
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  WRITE(*,'(A)') "Program starting."

  !Set all diganostic levels on for testing
  CALL CMISSDiagnosticsSetOn(CMISSFromDiagType,(/1,2,3,4,5/),"Diagnostics",(/"PROBLEM_FINITE_ELEMENT_CALCULATE"/),Err)

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
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemTypeSet(CoordinateSystem,CMISSCoordinateRectangularCartesianType,Err)
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NumberOfSpatialCoordinates,Err)
  CALL CMISSCoordinateSystemOriginSet(CoordinateSystem,(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/),Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !Create a region and assign the CS to the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !Define basis functions - tri-linear Lagrange and tri-Quadratic Lagrange, each with collapsed variant
    !Quadratic Basis
  CALL CMISSBasisTypeInitialise(QuadraticBasis,Err)
  CALL CMISSBasisCreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasisInterpolationXiSet(QuadraticBasis,(/CMISSBasisQuadraticLagrangeInterpolation, &
    & CMISSBasisQuadraticLagrangeInterpolation,CMISSBasisQuadraticLagrangeInterpolation/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)
!     & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err) !Have to do this
  CALL CMISSBasisCreateFinish(QuadraticBasis,Err)

    !Collapsed Quadratic Basis
  CALL CMISSBasisTypeInitialise(QuadraticCollapsedBasis,Err)
  CALL CMISSBasisCreateStart(QuadraticCollapsedBasisUserNumber,QuadraticCollapsedBasis,Err)
  CALL CMISSBasisTypeSet(QuadraticCollapsedBasis,CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(QuadraticCollapsedBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasisInterpolationXiSet(QuadraticCollapsedBasis,(/CMISSBasisQuadraticLagrangeInterpolation, &
       & CMISSBasisQuadraticLagrangeInterpolation,CMISSBasisQuadraticLagrangeInterpolation/),Err)
  CALL CMISSBasisCollapsedXiSet(QuadraticCollapsedBasis,(/CMISSBasisXiCollapsed, &
       & CMISSBasisCollapsedAtXi0,CMISSBasisNotCollapsed/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(QuadraticCollapsedBasis, &
       & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)  
!     & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(QuadraticCollapsedBasis,.true.,Err) !Have to do this
  CALL CMISSBasisCreateFinish(QuadraticCollapsedBasis,Err)

    !Linear Basis
  CALL CMISSBasisTypeInitialise(LinearBasis,Err)
  CALL CMISSBasisCreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(LinearBasis, &
    & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)
!     & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err) !Have to do this (unused) due to field_interp setup
  CALL CMISSBasisCreateFinish(LinearBasis,Err)

    !Collapsed Linear Basis
  CALL CMISSBasisTypeInitialise(LinearCollapsedBasis,Err)
  CALL CMISSBasisCreateStart(LinearCollapsedBasisUserNumber,LinearCollapsedBasis,Err)
  CALL CMISSBasisTypeSet(LinearCollapsedBasis,CMISSBasisLagrangeHermiteTPType,Err)
  CALL CMISSBasisNumberOfXiSet(LinearCollapsedBasis,NumberOfXiCoordinates,Err)
  CALL CMISSBasisInterpolationXiSet(LinearCollapsedBasis,(/CMISSBasisLinearLagrangeInterpolation, &
       & CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation/),Err)
  CALL CMISSBasisCollapsedXiSet(LinearCollapsedBasis,(/CMISSBasisXiCollapsed,CMISSBasisCollapsedAtXi0,CMISSBasisNotCollapsed/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(LinearCollapsedBasis, &
       & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)
!     & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(LinearCollapsedBasis,.true.,Err) !Have to do this (unused) due to field_interp setup
  CALL CMISSBasisCreateFinish(LinearCollapsedBasis,Err)

  !
  !================================================================================================================================
  !

  !Start the creation of a generated ellipsoid mesh
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up an ellipsoid mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshEllipsoidMeshType,Err)
  !Set the quadratic and linear bases
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,[QuadraticBasis,QuadraticCollapsedBasis,LinearBasis,LinearCollapsedBasis],Err)
  !Define the mesh on the region
  CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/LONG_AXIS,SHORT_AXIS,WALL_THICKNESS,CUTOFF_ANGLE/),Err)
  CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements/),Err)
  
  !Finish the creation of a generated mesh in the region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !
  !================================================================================================================================
  !

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecompositionCalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !
  !================================================================================================================================
  !

  ! --- GeometricFieldSolid ---
  !Create a field to put the geometry (default is geometry)
  CALL CMISSFieldTypeInitialise(GeometricFieldSolid,Err)
  CALL CMISSFieldCreateStart(FieldGeometryUserNumberSolid,Region,GeometricFieldSolid,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricFieldSolid,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricFieldSolid,CMISSFieldGeometricType,Err)
  CALL CMISSFieldNumberOfVariablesSet(GeometricFieldSolid,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricFieldSolid,CMISSFieldUVariableType,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(GeometricFieldSolid,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricFieldSolid,GeneratedMesh,Err)

  !
  !================================================================================================================================
  !

  ! --- GeometricFieldDarcy ---
  !Create a field to put the geometry (default is geometry)
  CALL CMISSFieldTypeInitialise(GeometricFieldDarcy,Err)
  CALL CMISSFieldCreateStart(FieldGeometryUserNumberDarcy,Region,GeometricFieldDarcy,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricFieldDarcy,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricFieldDarcy,CMISSFieldGeometricType,Err)
  CALL CMISSFieldNumberOfVariablesSet(GeometricFieldDarcy,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricFieldDarcy,CMISSFieldUVariableType,FieldGeometryNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldDarcy,CMISSFieldUVariableType,1,GeometricFieldDarcyMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldDarcy,CMISSFieldUVariableType,2,GeometricFieldDarcyMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldDarcy,CMISSFieldUVariableType,3,GeometricFieldDarcyMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(GeometricFieldDarcy,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricFieldDarcy,GeneratedMesh,Err)

  !
  !================================================================================================================================
  !

  !Create a fibre field and attach it to the geometric field
  CALL CMISSFieldTypeInitialise(FibreFieldSolid,Err)
  CALL CMISSFieldCreateStart(FieldFibreUserNumber,Region,FibreFieldSolid,Err)
  CALL CMISSFieldTypeSet(FibreFieldSolid,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreFieldSolid,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(FibreFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSFieldNumberOfVariablesSet(FibreFieldSolid,FieldFibreNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreFieldSolid,CMISSFieldUVariableType,FieldFibreNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(FibreFieldSolid,Err)

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
    CALL CMISSDecompositionNodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      FibreFieldAngle=(/zero,zero,zero/) 
      DO component_idx=1,FieldFibreNumberOfComponents
        CALL CMISSFieldParameterSetUpdateNode(FibreFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1, &
          & DerivativeUserNumber,node_idx,component_idx,FibreFieldAngle(component_idx),Err)
      ENDDO
    ENDIF
    theta=atan(FIBRE_SLOPE_CHANGE*XI3+FIBRE_SLOPE_INTERSECTION)
    DO j=2, TOTAL_NUMBER_NODES_XI(2) 
      nu=PI-XI2delta*(j-1)
      omega=PI/2+cos(2*nu)*atan(SHEET_SLOPE_BASE_ENDO)*(-2*XI3+1)
      DO i=1, TOTAL_NUMBER_NODES_XI(1)
        node_idx=node_idx+1
        CALL CMISSDecompositionNodeDomainGet(Decomposition,node_idx,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          FibreFieldAngle=(/theta,zero,omega/)
          DO component_idx=1,FieldFibreNumberOfComponents
            CALL CMISSFieldParameterSetUpdateNode(FibreFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1, &
              & DerivativeUserNumber, node_idx,component_idx,FibreFieldAngle(component_idx),Err)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    XI3=XI3+XI3delta
  ENDDO

  !Create a material field and attach it to the geometric field
  CALL CMISSFieldTypeInitialise(MaterialFieldSolid,Err)
  CALL CMISSFieldCreateStart(FieldMaterialUserNumber,Region,MaterialFieldSolid,Err)
  CALL CMISSFieldTypeSet(MaterialFieldSolid,CMISSFieldMaterialType,Err)
  CALL CMISSFieldMeshDecompositionSet(MaterialFieldSolid,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(MaterialFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialFieldSolid,FieldMaterialNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialFieldSolid,CMISSFieldUVariableType,FieldMaterialNumberOfComponents,Err)  
  CALL CMISSFieldComponentInterpolationSet(MaterialFieldSolid,CMISSFieldUVariableType,1,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(MaterialFieldSolid,CMISSFieldUVariableType,2,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(MaterialFieldSolid,CMISSFieldUVariableType,3,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldCreateFinish(MaterialFieldSolid,Err)

  !Set Mooney-Rivlin constants 
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,C1,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,C2,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,C3,Err)

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
    CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberDarcy,Region,GeometricFieldDarcy,CMISSEquationsSetFluidMechanicsClass, & 
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

  !
  !================================================================================================================================
  !
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
  CALL CMISSFieldScalingTypeSet(DependentFieldSolid,CMISSFieldUnitScaling,Err)

  CALL CMISSFieldCreateFinish(DependentFieldSolid,Err)

  !
  !================================================================================================================================
  !

  CALL CMISSEquationsSetDependentCreateStart(EquationsSetSolid,FieldDependentSolidUserNumber,DependentFieldSolid,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetSolid,Err)

  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetSolid,FieldMaterialUserNumber,MaterialFieldSolid,Err)  
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !
  DO icompartment = 1,Ncompartments
    CALL CMISSEquationsSetDependentCreateStart(EquationsSetDarcy(icompartment),FieldDependentSolidUserNumber,&
      & DependentFieldSolid,Err)
    CALL CMISSEquationsSetDependentCreateFinish(EquationsSetDarcy(icompartment),Err)
  ENDDO

  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS+1
    CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldVVariableType,CMISSFieldValuesSetType, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
    CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldU1VariableType,CMISSFieldValuesSetType, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
    CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldU2VariableType,CMISSFieldValuesSetType, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
    CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldU3VariableType,CMISSFieldValuesSetType, &
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
  !MATERIALS FIELDS
  !Auto-created field contains a U variable type to store the diffusion coefficient(s)
  !It also contains a V variable type to store the coupling coefficients 
  DO icompartment = 1,Ncompartments
    MaterialsFieldUserNumberDarcy = 400+icompartment
    CALL CMISSFieldTypeInitialise(MaterialsFieldDarcy(icompartment),Err)
    CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetDarcy(icompartment),MaterialsFieldUserNumberDarcy,&
         & MaterialsFieldDarcy(icompartment),Err)
    CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetDarcy(icompartment),Err)
  END DO
  DO icompartment = 1,Ncompartments
    CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMISSFieldUVariableType, &
      & CMISSFieldValuesSetType, &
      & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
    CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMISSFieldUVariableType, &
      & CMISSFieldValuesSetType, &
      & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)
  END DO
  DO icompartment = 1, Ncompartments
    DO COMPONENT_NUMBER=1, Ncompartments
      CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMISSFieldVVariableType, &
         & CMISSFieldValuesSetType,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
!         CALL CMISSFieldParameterSetUpdateConstant(MaterialsFieldDarcy(icompartment),CMISSFieldVVariableType, &
!           & CMISSFieldValuesSetType,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
    END DO
  END DO
  DO icompartment = 1, Ncompartments
    DO COMPONENT_NUMBER=1,Nparams
      CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy(icompartment),CMISSFieldU1VariableType, &
         & CMISSFieldValuesSetType,COMPONENT_NUMBER,ConstitutiveParams(icompartment,COMPONENT_NUMBER),Err)
!         CALL CMISSFieldParameterSetUpdateConstant(MaterialsFieldDarcy(icompartment),CMISSFieldVVariableType, &
!           & CMISSFieldValuesSetType,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
    END DO
  END DO

  !
  !================================================================================================================================
  !

  !EQUATIONS SET EQUATIONS

  !Darcy
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

  !Solid
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetSolid,Equations,Err)
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetSolid,Err)   

  !
  !================================================================================================================================
  !

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 1,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)
!   CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,-14.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,0.0_CMISSDP,Err)

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

  !Create the problem solvers
  CALL CMISSSolverTypeInitialise(SolverSolid,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverSolid,Err)
  CALL CMISSSolverTypeInitialise(DynamicSolverDarcy,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverDarcy,Err)

  CALL CMISSProblemSolversCreateStart(Problem,Err)

  ! Solid
  CALL CMISSProblemSolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopSolidNumber,CMISSControlLoopNode/), &
    & SolverSolidIndex,SolverSolid,Err)
  CALL CMISSSolverOutputTypeSet(SolverSolid,CMISSSolverProgressOutput,Err)
!   CALL CMISSSolverNewtonJacobianCalculationTypeSet(SolverSolid,CMISSSolverNewtonJacobianFDCalculated,Err)
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(SolverSolid,CMISSSolverNewtonJacobianAnalyticCalculated,Err)

!   CALL CMISSSolverNonLinearTypeSet(SolverSolid,CMISSSolverNonlinearNewton,Err)
!   CALL CMISSSolverLibraryTypeSet(SolverSolid,CMISSSolverPETScLibrary,Err)

  CALL CMISSSolverNewtonAbsoluteToleranceSet(SolverSolid,ABSOLUTE_TOLERANCE,Err)
  CALL CMISSSolverNewtonRelativeToleranceSet(SolverSolid,RELATIVE_TOLERANCE,Err)
  CALL CMISSSolverNewtonMaximumIterationsSet(SolverSolid,MAXIMUM_ITERATIONS,Err)

  CALL CMISSSolverNewtonLinearSolverGet(SolverSolid,LinearSolverSolid,Err)
  CALL CMISSSolverLinearTypeSet(LinearSolverSolid,CMISSSolverLinearDirectSolveType,Err)

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

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsSolid,BoundaryConditions,Err)

  !Grab the list of nodes on inner, outer and top surfaces
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshEllipsoidTopSurfaceType,TopSurfaceNodes,TopNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshEllipsoidInnerSurfaceType,InnerSurfaceNodes,InnerNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshEllipsoidOuterSurfaceType,OuterSurfaceNodes,OuterNormalXi,Err)

  write(*,*)'TopSurfaceNodes = ',TopSurfaceNodes

  ! ASSIGN BOUNDARY CONDITIONS
  !Fix base of the ellipsoid in z direction
  DO NN=1,SIZE(TopSurfaceNodes,1)
    NODE=TopSurfaceNodes(NN)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetGetNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NODE,3,ZCoord,Err)
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentFieldSolid,CMISSFieldUVariableType,1,1,NODE,3, &
        & CMISSBoundaryConditionFixed,ZCoord,Err)
    ENDIF
  ENDDO

  !Apply inner surface pressure
  !NOTE: Surface pressure goes into pressure_values_set_type of the DELUDELN type
  DO NN=1,SIZE(InnerSurfaceNodes,1)
    NODE=InnerSurfaceNodes(NN)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentFieldSolid,CMISSFieldDelUDelNVariableType,1,1,NODE, &
        & ABS(InnerNormalXi), &
        & CMISSBoundaryConditionPressure,INNER_PRESSURE,Err)
    ENDIF
  ENDDO

  !Apply outer surface pressure
  DO NN=1,SIZE(OuterSurfaceNodes,1)
    NODE=OuterSurfaceNodes(NN)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentFieldSolid,CMISSFieldDelUDelNVariableType,1,1,NODE, &
        & ABS(OuterNormalXi), &
        & CMISSBoundaryConditionPressure,OUTER_PRESSURE,Err)
    ENDIF
  ENDDO

  !Fix more nodes at the base to stop free body motion
  X_FIXED=.FALSE.
  Y_FIXED=.FALSE.
  DO NN=1,SIZE(TopSurfaceNodes,1)
    NODE=TopSurfaceNodes(NN)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSFieldParameterSetGetNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NODE,1,XCoord,Err)
      CALL CMISSFieldParameterSetGetNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NODE,2,YCoord,Err)
      IF(ABS(XCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentFieldSolid,CMISSFieldUVariableType,1,1,NODE,1, &
          & CMISSBoundaryConditionFixed,XCoord,Err)
        WRITE(*,*) "FIXING NODE",NODE,"IN X DIRECTION"
        X_FIXED=.TRUE.
    ENDIF
      IF(ABS(YCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentFieldSolid,CMISSFieldUVariableType,1,1,NODE,2, &
          & CMISSBoundaryConditionFixed,YCoord,Err)
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

  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsSolid,Err)

  !
  !================================================================================================================================
  !

  !BCs Darcy
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsDarcy,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)
  DO icompartment=1,Ncompartments

    !In 'generated_mesh_routines.f90/GENERATED_MESH_ELLIPSOID_SURFACE_GET' there is a bug:
    !  BASIS=>ELLIPSOID_MESH%BASES(MESH_COMPONENT)%PTR does not account for the fact that:
    !  in 'generated_mesh_routines.f90/GENERATED_MESH_ELLIPSOID_CREATE_FINISH' the following is done:
    !  CALL MESH_NUMBER_OF_COMPONENTS_SET(GENERATED_MESH%MESH,SIZE(ELLIPSOID_MESH%BASES)/2,ERR,ERROR,*999)
    !Temporary work around, until bug fix:

    !MeshComponentNumber_dummy = DarcyVelMeshComponentNumber
    MeshComponentNumber_dummy = 3

!     !inner surface
!     CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,MeshComponentNumber_dummy,CMISSGeneratedMeshEllipsoidInnerSurfaceType, &
!       & InnerSurfaceNodesDarcyVel,InnerNormalXi,Err)
!
!     write(*,*)'InnerSurfaceNodesDarcyVel = ',InnerSurfaceNodesDarcyVel
!
!     !Set all inner surface nodes impermeable
!     DO NN=1,SIZE(InnerSurfaceNodesDarcyVel,1)
! !       VALUE = 0.0_CMISSDP
! !       COMPONENT_NUMBER = 1
! !       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,InnerSurfaceNodesDarcyVel(NN), &
! !         & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
! !       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)
! !
! !       VALUE = 0.0_CMISSDP
! !       COMPONENT_NUMBER = 2
! !       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,InnerSurfaceNodesDarcyVel(NN), &
! !         & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
! !       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)
! !
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 3
! !       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,InnerSurfaceNodesDarcyVel(NN), &
! !         & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,InnerSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionImpermeableWall,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)
!     ENDDO


!     !outer surface
!     CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,MeshComponentNumber_dummy,CMISSGeneratedMeshEllipsoidOuterSurfaceType, &
!       & OuterSurfaceNodesDarcyVel,OuterNormalXi,Err)
!
!     write(*,*)'OuterSurfaceNodesDarcyVel = ',OuterSurfaceNodesDarcyVel
!
!     !Set all outer surface nodes impermeable
!     DO NN=1,SIZE(OuterSurfaceNodesDarcyVel,1)
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 1
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)
!
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 2
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)
!
!       VALUE = 0.0_CMISSDP
!       COMPONENT_NUMBER = 3
!       CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,OuterSurfaceNodesDarcyVel(NN), &
!         & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
!       IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)
!     ENDDO


    !top surface
    CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,MeshComponentNumber_dummy,CMISSGeneratedMeshEllipsoidTopSurfaceType, &
      & TopSurfaceNodesDarcyVel,TopNormalXi,Err)
    IF(icompartment==1) THEN
    write(*,*)'TopSurfaceNodesDarcyVel = ',TopSurfaceNodesDarcyVel

    !Set all top surface nodes to Darcy inflow BC
      DO NN=1,SIZE(TopSurfaceNodesDarcyVel,1)
        VALUE = -0.25_CMISSDP
        COMPONENT_NUMBER = 3
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldVVariableType,1,1, &
          & TopSurfaceNodesDarcyVel(NN), &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING TOP DARCY BC TO NODE", TopSurfaceNodesDarcyVel(NN)
      ENDDO
    ELSEIF(icompartment==2)THEN

    !Set all top surface nodes to Darcy inflow BC
      DO NN=1,SIZE(TopSurfaceNodesDarcyVel,1)
        VALUE = -0.25_CMISSDP
        COMPONENT_NUMBER = 3
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldU1VariableType,1,1, &
          & TopSurfaceNodesDarcyVel(NN), &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING TOP DARCY BC TO NODE", TopSurfaceNodesDarcyVel(NN)
      ENDDO
    ELSEIF(icompartment==3)THEN
    !Set all top surface nodes to Darcy inflow BC
      DO NN=1,SIZE(TopSurfaceNodesDarcyVel,1)
        VALUE = -0.25_CMISSDP
        COMPONENT_NUMBER = 3
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldU2VariableType,1,1, &
          & TopSurfaceNodesDarcyVel(NN), &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING TOP DARCY BC TO NODE", TopSurfaceNodesDarcyVel(NN)
      ENDDO
    ELSEIF(icompartment==4)THEN
    !Set all top surface nodes to Darcy inflow BC
      DO NN=1,SIZE(TopSurfaceNodesDarcyVel,1)
        VALUE = -0.25_CMISSDP
        COMPONENT_NUMBER = 3
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISSFieldU3VariableType,1,1, &
          & TopSurfaceNodesDarcyVel(NN), &
          & COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)
        IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING TOP DARCY BC TO NODE", TopSurfaceNodesDarcyVel(NN)
      ENDDO
    ENDIF


  ENDDO
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsDarcy,Err)

  !
  !================================================================================================================================
  !

  !Solve problem
  WRITE(*,'(A)') "Solving problem..."
  CALL CMISSProblemSolve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"

  !
  !================================================================================================================================
  !

  !Output solution  
!   CALL CMISSFieldsTypeInitialise(Fields,Err)
!   CALL CMISSFieldsTypeCreate(Region,Fields,Err)
!   CALL CMISSFieldIONodesExport(Fields,"QuadraticEllipsoidDrivenMultiCompDarcy","FORTRAN",Err)
!   CALL CMISSFieldIOElementsExport(Fields,"QuadraticEllipsoidDrivenMultiCompDarcy","FORTRAN",Err)
!   CALL CMISSFieldsTypeFinalise(Fields,Err)

  !
  !================================================================================================================================
  !

  CALL CMISSFinalise(Err)

  IF (ALLOCATED(EquationsSetFieldDarcy)) DEALLOCATE(EquationsSetFieldDarcy)
  IF (ALLOCATED(EquationsSetDarcy)) DEALLOCATE(EquationsSetDarcy)
  IF (ALLOCATED(MaterialsFieldDarcy)) DEALLOCATE(MaterialsFieldDarcy)
  IF (ALLOCATED(BoundaryConditionsDarcy)) DEALLOCATE(BoundaryConditionsDarcy)
  IF (ALLOCATED(EquationsDarcy)) DEALLOCATE(EquationsDarcy)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM QUADRATICELLIPSOIDDRIVENMULTICOMPDARCYEXAMPLE

