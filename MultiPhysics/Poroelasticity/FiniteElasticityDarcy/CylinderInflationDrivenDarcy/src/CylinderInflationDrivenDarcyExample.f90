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

!> \example MultiPhysics/Poroelasticity/FiniteElasticityDarcy/CylinderInflationDrivenDarcy/src/CylinderInflationDrivenDarcyExample.f90
!! Example program to solve coupled FiniteElasticityDarcy equations using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/Poroelasticity/FiniteElasticityDarcy/CylinderInflationDrivenDarcy/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/Poroelasticity/FiniteElasticityDarcy/CylinderInflationDrivenDarcy/build-intel'>Linux GNU Build</a>
!!
!<

! !
! !  This example considers a coupled Finite Elasticity Darcy problem on a cylindrical geometry
! !

!> Main program

PROGRAM CYLINDERINFLATIONDRIVENDARCYEXAMPLE

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

  !\todo: don't hard code, read in + default
  REAL(CMISSDP), PARAMETER :: INNER_PRESSURE=1.0_CMISSDP !Positive is compressive
  REAL(CMISSDP), PARAMETER :: OUTER_PRESSURE=0.0_CMISSDP !Positive is compressive
  REAL(CMISSDP), PARAMETER :: LAMBDA=1.1_CMISSDP
  REAL(CMISSDP), PARAMETER :: TSI=0.0_CMISSDP    !Not yet working. Leave at 0
  REAL(CMISSDP), PARAMETER :: INNER_RAD=1.0_CMISSDP
!   REAL(CMISSDP), PARAMETER :: OUTER_RAD=1.2_CMISSDP
  REAL(CMISSDP), PARAMETER :: OUTER_RAD=1.6_CMISSDP
  REAL(CMISSDP), PARAMETER :: HEIGHT=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: deformedHeight = 1.0_CMISSDP * HEIGHT
  REAL(CMISSDP), PARAMETER :: C1=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: C2=6.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: C3=10.0_CMISSDP
  INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalXElements=2
  INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalYElements=8
  INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalZElements=2

!   !Standard test parameters (don't remove or change)
!   REAL(CMISSDP), PARAMETER :: INNER_PRESSURE=1.0_CMISSDP !Positive is compressive
!   REAL(CMISSDP), PARAMETER :: OUTER_PRESSURE=0.0_CMISSDP !Positive is compressive
!   REAL(CMISSDP), PARAMETER :: LAMBDA=1.1_CMISSDP
!   REAL(CMISSDP), PARAMETER :: TSI=0.0_CMISSDP    !Not yet used
!   REAL(CMISSDP), PARAMETER :: INNER_RAD=1.0_CMISSDP
!   REAL(CMISSDP), PARAMETER :: OUTER_RAD=1.2_CMISSDP
!   REAL(CMISSDP), PARAMETER :: HEIGHT=2.0_CMISSDP
!   REAL(CMISSDP), PARAMETER :: C1=2.0_CMISSDP
!   REAL(CMISSDP), PARAMETER :: C2=6.0_CMISSDP
!   INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalXElements=2 !\todo: don't hardcode?
!   INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalYElements=8
!   INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalZElements=2
!  Increment loop of 2

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfSpatialCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1

  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshDimensions=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfXiCoordinates=3
  INTEGER(CMISSIntg), PARAMETER :: NumberOfMeshComponents=2
  INTEGER(CMISSIntg), PARAMETER :: QuadraticMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: LinearMeshComponentNumber=2

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

!   INTEGER(CMISSIntg), PARAMETER :: FieldAnalyticUserNumber=1337

  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumberSolid=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberSolid=5
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1


  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcy=8
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberDarcy=12
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberDarcy=22

  INTEGER(CMISSIntg), PARAMETER :: ControlLoopSolidNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopFluidNumber=2
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopSubiterationNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverSolidIndex=1
  INTEGER(CMISSIntg), PARAMETER :: SolverDarcyIndex=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPorosity=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPermOverVis=2

  !Program types


  !Program variables
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: EquationsSetIndex  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  INTEGER(CMISSIntg) :: DarcyVelMeshComponentNumber, DarcyMassIncreaseMeshComponentNumber

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

  TYPE(CMISSBasisType) :: QuadraticBasis, LinearBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSetSolid
  TYPE(CMISSFieldType) :: GeometricFieldSolid,FibreField,MaterialField
  TYPE(CMISSFieldType) :: DependentField,EquationsSetFieldSolid !,AnalyticField

  TYPE(CMISSFieldType) :: GeometricFieldDarcy

  TYPE(CMISSFieldType) :: MaterialsFieldDarcy
  TYPE(CMISSFieldType) :: EquationsSetFieldDarcy

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

  !Other variables
  INTEGER(CMISSIntg) :: NN,NE,E

  INTEGER(CMISSIntg),ALLOCATABLE :: BottomSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: TopSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: InnerSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: OuterSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: BottomSurfaceNodesDarcyVel(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: TopSurfaceNodesDarcyVel(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: InnerSurfaceNodesDarcyVel(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: OuterSurfaceNodesDarcyVel(:)
  INTEGER(CMISSIntg) :: BottomNormalXi,TopNormalXi,InnerNormalXi,OuterNormalXi
  REAL(CMISSDP) :: xValue,yValue, InitialPressure
  LOGICAL :: X_FIXED,Y_FIXED
  
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
  DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE=CMISSSolverProgressOutput
  LINEAR_SOLVER_DARCY_OUTPUT_TYPE=CMISSSolverSolverOutput
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

  write(*,*) "NumberOfDomains=",NumberOfComputationalNodes
  NumberOfDomains=NumberOfComputationalNodes !1

  !Broadcast the number of elements in the X,Y and Z directions and the number of partitions to the other computational nodes
!   CALL MPI_BCAST(NumberGlobalXElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR) ! NOW A PARAMETER
!   CALL MPI_BCAST(NumberGlobalYElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR) ! NOW A PARAMETER
!   CALL MPI_BCAST(NumberGlobalZElements,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR) ! NOW A PARAMETER 
  CALL MPI_BCAST(NumberOfDomains,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !
  !================================================================================================================================
  !

  NUMBER_OF_DIMENSIONS = 3

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

  !Define basis functions - tri-linear Lagrange and tri-Quadratic Lagrange
  CALL CMISSBasisTypeInitialise(LinearBasis,Err)
  CALL CMISSBasisCreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(LinearBasis, &
!     & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)
    & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err) !Have to do this (unused) due to field_interp setup
  CALL CMISSBasisCreateFinish(LinearBasis,Err)

  CALL CMISSBasisTypeInitialise(QuadraticBasis,Err)
  CALL CMISSBasisCreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasisInterpolationXiSet(QuadraticBasis,(/CMISSBasisQuadraticLagrangeInterpolation, &
    & CMISSBasisQuadraticLagrangeInterpolation,CMISSBasisQuadraticLagrangeInterpolation/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(QuadraticBasis, &
!     & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)
    & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err) !Enable 3D interpolation on faces
  CALL CMISSBasisCreateFinish(QuadraticBasis,Err)

  !
  !================================================================================================================================
  !

  !Start the creation of a generated cylinder mesh
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshCylinderMeshType,Err)
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,[QuadraticBasis,LinearBasis],Err)
  CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/INNER_RAD, OUTER_RAD, HEIGHT/),Err)
  CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements/),Err)
  
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !
  !================================================================================================================================
  !

  !Create a decomposition
  CALL CMISSRandomSeedsSet(0_CMISSIntg,Err) !To keep the automatic decomposition same each time
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Automatic decomposition
!   CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
!   CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  !Manual decomposition
  IF(NumberOfDomains>1) THEN
    CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionUserDefinedType,Err)
    !Set all elements but last one to first domain
    CALL CMISSMeshNumberOfElementsGet(Mesh,NE,Err)
    do E=1,NE/2
      CALL CMISSDecompositionElementDomainSet(Decomposition,E,0,Err)
    enddo
    do E=NE/2+1,NE
      CALL CMISSDecompositionElementDomainSet(Decomposition,E,1,Err)
    enddo
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  ENDIF
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
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldDarcy,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldDarcy,CMISSFieldUVariableType,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldDarcy,CMISSFieldUVariableType,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(GeometricFieldDarcy,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricFieldDarcy,GeneratedMesh,Err)

  !
  !================================================================================================================================
  !

  !Create a fibre field and attach it to the geometric field
  CALL CMISSFieldTypeInitialise(FibreField,Err)
  CALL CMISSFieldCreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSFieldTypeSet(FibreField,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(FibreField,GeometricFieldSolid,Err)
  CALL CMISSFieldNumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreField,CMISSFieldUVariableType,FieldFibreNumberOfComponents,Err)  
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,1,LinearMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,2,LinearMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreField,CMISSFieldUVariableType,3,LinearMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(FibreField,Err)

  !
  !================================================================================================================================
  !

  !Create a material field and attach it to the geometric field
  CALL CMISSFieldTypeInitialise(MaterialField,Err)
  CALL CMISSFieldCreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL CMISSFieldTypeSet(MaterialField,CMISSFieldMaterialType,Err)
  CALL CMISSFieldMeshDecompositionSet(MaterialField,Decomposition,Err)        
  CALL CMISSFieldGeometricFieldSet(MaterialField,GeometricFieldSolid,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialField,CMISSFieldUVariableType,FieldMaterialNumberOfComponents,Err)  
  CALL CMISSFieldComponentInterpolationSet(MaterialField,CMISSFieldUVariableType,1,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(MaterialField,CMISSFieldUVariableType,2,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(MaterialField,CMISSFieldUVariableType,3,CMISSFieldConstantInterpolation,Err)
  CALL CMISSFieldCreateFinish(MaterialField,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,C1,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,C2,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,C3,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for ALE Darcy
  CALL CMISSFieldTypeInitialise(EquationsSetFieldDarcy,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSetDarcy,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberDarcy,Region,GeometricFieldDarcy,CMISSEquationsSetFluidMechanicsClass, &
    & CMISSEquationsSetDarcyEquationType,CMISSEquationsSetIncompressibleElasticityDrivenDarcySubtype,&
    & EquationsSetFieldUserNumberDarcy,EquationsSetFieldDarcy,EquationsSetDarcy,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSetDarcy,Err)

  !Create the equations set for the solid
  CALL CMISSFieldTypeInitialise(EquationsSetFieldSolid,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSetSolid,Err)
  CALL CMISSEquationsSetCreateStart(EquationSetUserNumberSolid,Region,FibreField,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetIncompressibleElasticityDrivenDarcySubtype, &
    & EquationsSetFieldUserNumberSolid,EquationsSetFieldSolid,EquationsSetSolid,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create a dependent field with four variables (U, DelUDelN = solid, V, DelVDelN = Darcy) and four components
  !Solid: The U, DelUDelN variables have 4 components (3 displacement, 1 pressure)
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSFieldCreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL CMISSFieldTypeSet(DependentField,CMISSFieldGeneralType,Err)
  CALL CMISSFieldMeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(DependentField,GeometricFieldSolid,Err)
  CALL CMISSFieldDependentTypeSet(DependentField,CMISSFieldDependentType,Err)
  CALL CMISSFieldNumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)

  CALL CMISSFieldVariableTypesSet(DependentField,(/CMISSFieldUVariableType, &
    & CMISSFieldDelUDelNVariableType,CMISSFieldVVariableType,CMISSFieldDelVDelNVariableType/),Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldUVariableType,FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldDelUDelNVariableType,FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldVVariableType,FieldDependentFluidNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldDelVDelNVariableType,FieldDependentFluidNumberOfComponents,Err)

  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,4,LinearMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,1,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,2,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,3,QuadraticMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,4,LinearMeshComponentNumber,Err)

  DarcyVelMeshComponentNumber = LinearMeshComponentNumber
  DarcyMassIncreaseMeshComponentNumber = LinearMeshComponentNumber

  !Darcy: The V, DelVDelN variables have 4 components (3 velocities, 1 mass increase)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldVVariableType,1,DarcyVelMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldVVariableType,2,DarcyVelMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldVVariableType,3,DarcyVelMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldVVariableType,4,DarcyMassIncreaseMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelVDelNVariableType,1,DarcyVelMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelVDelNVariableType,2,DarcyVelMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelVDelNVariableType,3,DarcyVelMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelVDelNVariableType,4,DarcyMassIncreaseMeshComponentNumber,Err)

  CALL CMISSFieldScalingTypeSet(DependentField,CMISSFieldUnitScaling,Err)

  CALL CMISSFieldCreateFinish(DependentField,Err)

  CALL CMISSEquationsSetDependentCreateStart(EquationsSetSolid,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetSolid,Err)

  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetSolid,FieldMaterialUserNumber,MaterialField,Err)  
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetSolid,Err)


  CALL CMISSEquationsSetDependentCreateStart(EquationsSetDarcy,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetDarcy,Err)

  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS+1
    CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldVVariableType,CMISSFieldValuesSetType, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for ALE Darcy
  CALL CMISSFieldTypeInitialise(MaterialsFieldDarcy,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetDarcy,MaterialsFieldUserNumberDarcy, &
    & MaterialsFieldDarcy,Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetDarcy,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)

  !
  !================================================================================================================================
  !
!   !Set up analytic field
!   CALL CMISSFieldTypeInitialise(AnalyticField,Err)
!   CALL CMISSEquationsSetAnalyticCreateStart(EquationsSetSolid,CMISSEquationsSetFiniteElasticityCylinder, &
!     & FieldAnalyticUserNumber,AnalyticField,Err)
!   !Finish the equations set analytic field variables
!   CALL CMISSEquationsSetAnalyticCreateFinish(EquationsSetSolid,Err)
! 
!   !Set the analytic parameters
!   CALL CMISSEquationsSetAnalyticUserParamSet(EquationsSetSolid,CMISSFiniteElasticityAnalyticCylinderParamPinIdx,INNER_PRESSURE,Err)
!   CALL CMISSEquationsSetAnalyticUserParamSet(EquationsSetSolid,CMISSFiniteElasticityAnalyticCylinderParamPoutIdx,OUTER_PRESSURE,Err)
!   CALL CMISSEquationsSetAnalyticUserParamSet(EquationsSetSolid,CMISSFiniteElasticityAnalyticCylinderParamLambdaIdx,LAMBDA,Err)
!   CALL CMISSEquationsSetAnalyticUserParamSet(EquationsSetSolid,CMISSFiniteElasticityAnalyticCylinderParamTsiIdx,TSI,Err)
!   CALL CMISSEquationsSetAnalyticUserParamSet(EquationsSetSolid,CMISSFiniteElasticityAnalyticCylinderParamRinIdx,INNER_RAD,Err)
!   CALL CMISSEquationsSetAnalyticUserParamSet(EquationsSetSolid,CMISSFiniteElasticityAnalyticCylinderParamRoutIdx,OUTER_RAD,Err)
!   CALL CMISSEquationsSetAnalyticUserParamSet(EquationsSetSolid,CMISSFiniteElasticityAnalyticCylinderParamC1Idx,C1,Err)
!   CALL CMISSEquationsSetAnalyticUserParamSet(EquationsSetSolid,CMISSFiniteElasticityAnalyticCylinderParamC2Idx,C2,Err)

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
    & 1,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)
!   CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,-14.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,0.0_CMISSDP,Err)

  !
  !================================================================================================================================
  !

  !Define the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemMultiPhysicsClass,CMISSProblemFiniteElasticityDarcyType, &
    & CMISSProblemQuasistaticElasticityTransientDarcySubtype,Err)
  CALL CMISSProblemCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !Create the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
!   CALL CMISSControlLoopMaximumIterationsSet(ControlLoop,1,Err)  ! this one sets the increment loop counter
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

  !Create the problem solver equations
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

  !
  !================================================================================================================================
  !

! MANUAL BC ASSIGNMENT - REPLACED BY THE ANALYTIC BC ROUTINE BELOW
  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsSolid,BoundaryConditions,Err)

  !Get surfaces - will fix two nodes on bottom face, pressure conditions inside
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshCylinderBottomSurfaceType,BottomSurfaceNodes,BottomNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshCylinderTopSurfaceType,TopSurfaceNodes,TopNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshCylinderInnerSurfaceType,InnerSurfaceNodes,InnerNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,CMISSGeneratedMeshCylinderOuterSurfaceType,OuterSurfaceNodes,OuterNormalXi,Err)

  !Set all inner surface nodes to inner pressure
  DO NN=1,SIZE(InnerSurfaceNodes,1)
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldDelUDelNVariableType,1,1, &
        & InnerSurfaceNodes(NN),abs(InnerNormalXi),CMISSBoundaryConditionPressureIncremented,INNER_PRESSURE,Err)   ! INNER_PRESSURE
    IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER PRESSURE TO NODE", InnerSurfaceNodes(NN)
  ENDDO
        
  !Set all outer surface nodes to outer pressure
  DO NN=1,SIZE(OuterSurfaceNodes,1)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldDelUDelNVariableType,1,1, &
      & OuterSurfaceNodes(NN),abs(OuterNormalXi),CMISSBoundaryConditionPressureIncremented,OUTER_PRESSURE,Err)
    IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER PRESSURE TO NODE", OuterSurfaceNodes(NN)
  ENDDO

  !Set all top nodes fixed in z plane at the set height
  DO NN=1,SIZE(TopSurfaceNodes,1)
    CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,TopSurfaceNodes(NN), &
      & 3,CMISSBoundaryConditionFixed,deformedHeight,Err)
    IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING FIXED CONDITION TO NODE", TopSurfaceNodes(NN)
  ENDDO

!   !Set all bottom nodes fixed in z plane
!   DO NN=1,SIZE(BottomSurfaceNodes,1)
!     CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,BottomSurfaceNodes(NN), &
!       & 3,CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
!     IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING FIXED CONDITION TO NODE", BottomSurfaceNodes(NN)
!   ENDDO
! 
!   !Set two nodes on the bottom surface to axial displacement only
!   X_FIXED=.FALSE.
!   Y_FIXED=.FALSE.
!   DO NN=1,SIZE(BottomSurfaceNodes,1)
!     IF (.NOT.X_FIXED) THEN
!       CALL CMISSFieldParameterSetGetNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
!         & 1,BottomSurfaceNodes(NN),1,xValue,Err)
!       IF(abs(xValue)<1e-5_CMISSDP) THEN
!         !Constrain it in x direction
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,BottomSurfaceNodes(NN),1, &
!           & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
!         X_FIXED=.TRUE.
!         WRITE(*,*) "CyliderInflationExample: SUCCESSFULLY CONSTRAINED IN X DIRECTION NODE",BottomSurfaceNodes(NN)
!     ENDIF
!     ENDIF
!     IF(.NOT.Y_FIXED) THEN
!       CALL CMISSFieldParameterSetGetNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
!         & 1,BottomSurfaceNodes(NN),2,yValue,Err)
!       IF(abs(yValue)<1e-5_CMISSDP) THEN
!         !Constrain it in y direction
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,BottomSurfaceNodes(NN),2, &
!           & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
!         Y_FIXED=.TRUE.
!         WRITE(*,*) "CyliderInflationExample: SUCCESSFULLY CONSTRAINED IN Y DIRECTION NODE",BottomSurfaceNodes(NN)
!     ENDIF
!     ENDIF
!     IF (X_FIXED.AND.Y_FIXED) EXIT
!   ENDDO
!   !Check
!   IF(.NOT.X_FIXED .OR. .NOT.Y_FIXED) THEN
!     Write(*,*) "Couldn't fix bottom surface. No node lies on x or y axis, try changing number of elements"// &
!       & " in theta coordinate"
!     STOP
!     ENDIF

  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsSolid,Err)
! END OF MANUAL BC ASSIGNMENT

!   !Set the bc using the analytic solution routine
!   CALL CMISSEquationsSetBoundaryConditionsAnalytic(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  !BCs Darcy
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsDarcy,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)

    !top surface
    CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISSGeneratedMeshCylinderTopSurfaceType, &
      & TopSurfaceNodesDarcyVel,TopNormalXi,Err)

    !Set all top surface nodes to Darcy inflow BC
    DO NN=1,SIZE(TopSurfaceNodesDarcyVel,1)
      VALUE = -1.0_CMISSDP
      COMPONENT_NUMBER = 3
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,1,1, &
        & TopSurfaceNodesDarcyVel(NN),COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)

      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING TOP DARCY BC TO NODE", TopSurfaceNodesDarcyVel(NN)
    ENDDO


    !bottom surface
    CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISSGeneratedMeshCylinderBottomSurfaceType, &
      & BottomSurfaceNodesDarcyVel,BottomNormalXi,Err)

    !Set all bottom surface nodes impermeable
    DO NN=1,SIZE(BottomSurfaceNodesDarcyVel,1)
      VALUE = 0.0_CMISSDP
      COMPONENT_NUMBER = 3
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,1,1, &
        & BottomSurfaceNodesDarcyVel(NN),COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)

      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING BOTTOM DARCY BC TO NODE", BottomSurfaceNodesDarcyVel(NN)
    ENDDO


    !inner surface
    CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISSGeneratedMeshCylinderInnerSurfaceType, &
      & InnerSurfaceNodesDarcyVel,InnerNormalXi,Err)

    !Set all inner surface nodes impermeable
    DO NN=1,SIZE(InnerSurfaceNodesDarcyVel,1)
      VALUE = 0.0_CMISSDP
      COMPONENT_NUMBER = 1
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,1,1, &
        & InnerSurfaceNodesDarcyVel(NN),COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)

      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)

      VALUE = 0.0_CMISSDP
      COMPONENT_NUMBER = 2
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,1,1, &
        & InnerSurfaceNodesDarcyVel(NN),COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)

      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)
    ENDDO


    !outer surface
    CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISSGeneratedMeshCylinderOuterSurfaceType, &
      & OuterSurfaceNodesDarcyVel,OuterNormalXi,Err)

    !Set all outer surface nodes impermeable
    DO NN=1,SIZE(OuterSurfaceNodesDarcyVel,1)
      VALUE = 0.0_CMISSDP
      COMPONENT_NUMBER = 1
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,1,1, &
        & OuterSurfaceNodesDarcyVel(NN),COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)

      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)

      VALUE = 0.0_CMISSDP
      COMPONENT_NUMBER = 2
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,1,1, &
        & OuterSurfaceNodesDarcyVel(NN),COMPONENT_NUMBER,CMISSBoundaryConditionFixed,VALUE,Err)

      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)
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

!   !Output Analytic analysis
!   Call CMISSAnalyticAnalysisOutput(DependentField,"outputs/CylinderInflation",Err)

  !
  !================================================================================================================================
  !

  !Output solution  
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"outputs/CylinderInflation","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"outputs/CylinderInflation","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)

  !
  !================================================================================================================================
  !

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM CYLINDERINFLATIONDRIVENDARCYEXAMPLE

