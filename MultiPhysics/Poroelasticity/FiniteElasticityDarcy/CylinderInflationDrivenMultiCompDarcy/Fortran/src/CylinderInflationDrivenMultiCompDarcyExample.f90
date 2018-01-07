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

  !\todo: don't hard code, read in + default
  REAL(CMISSRP), PARAMETER :: INNER_PRESSURE=1.0_CMISSRP !Positive is compressive
  REAL(CMISSRP), PARAMETER :: OUTER_PRESSURE=0.0_CMISSRP !Positive is compressive
  REAL(CMISSRP), PARAMETER :: LAMBDA=1.1_CMISSRP
  REAL(CMISSRP), PARAMETER :: TSI=0.0_CMISSRP    !Not yet working. Leave at 0
  REAL(CMISSRP), PARAMETER :: INNER_RAD=1.0_CMISSRP
!   REAL(CMISSRP), PARAMETER :: OUTER_RAD=1.2_CMISSRP
  REAL(CMISSRP), PARAMETER :: OUTER_RAD=1.6_CMISSRP
  REAL(CMISSRP), PARAMETER :: HEIGHT=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: deformedHeight = 1.0_CMISSRP * HEIGHT
  REAL(CMISSRP), PARAMETER :: C1=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: C2=6.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: C3=10.0_CMISSRP
  INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalXElements=2
  INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalYElements=8
  INTEGER(CMISSIntg), PARAMETER ::   NumberGlobalZElements=2

!   !Standard test parameters (don't remove or change)
!   REAL(CMISSRP), PARAMETER :: INNER_PRESSURE=1.0_CMISSRP !Positive is compressive
!   REAL(CMISSRP), PARAMETER :: OUTER_PRESSURE=0.0_CMISSRP !Positive is compressive
!   REAL(CMISSRP), PARAMETER :: LAMBDA=1.1_CMISSRP
!   REAL(CMISSRP), PARAMETER :: TSI=0.0_CMISSRP    !Not yet used
!   REAL(CMISSRP), PARAMETER :: INNER_RAD=1.0_CMISSRP
!   REAL(CMISSRP), PARAMETER :: OUTER_RAD=1.2_CMISSRP
!   REAL(CMISSRP), PARAMETER :: HEIGHT=2.0_CMISSRP
!   REAL(CMISSRP), PARAMETER :: C1=2.0_CMISSRP
!   REAL(CMISSRP), PARAMETER :: C2=6.0_CMISSRP
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

  TYPE(cmfe_BasisType) :: QuadraticBasis, LinearBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSetSolid
  TYPE(cmfe_FieldType) :: GeometricFieldSolid,FibreField,MaterialField
  TYPE(cmfe_FieldType) :: DependentField,EquationsSetFieldSolid !,AnalyticField

  TYPE(cmfe_FieldType) :: GeometricFieldDarcy

  TYPE(cmfe_FieldType) :: MaterialsFieldDarcy
  TYPE(cmfe_FieldType) :: EquationsSetFieldDarcy

  !Boundary conditions
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsDarcy
  !Equations sets
  TYPE(cmfe_EquationsSetType) :: EquationsSetDarcy
  !Equations
  TYPE(cmfe_EquationsType) :: EquationsDarcy

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
  REAL(CMISSRP) :: xValue,yValue, InitialPressure
  LOGICAL :: X_FIXED,Y_FIXED
  
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
  DYNAMIC_SOLVER_DARCY_STOP_TIME=2_CMISSIntg * DYNAMIC_SOLVER_DARCY_TIME_INCREMENT
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

  !Define basis functions - tri-linear Lagrange and tri-Quadratic Lagrange
  CALL cmfe_Basis_Initialise(LinearBasis,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis, &
!     & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
    & [CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err) !Have to do this (unused) due to field_interp setup
  CALL cmfe_Basis_CreateFinish(LinearBasis,Err)

  CALL cmfe_Basis_Initialise(QuadraticBasis,Err)
  CALL cmfe_Basis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL cmfe_Basis_InterpolationXiSet(QuadraticBasis,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
!     & [CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME,CMFE_BASIS_MID_QUADRATURE_SCHEME],Err)
    & [CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME,CMFE_BASIS_HIGH_QUADRATURE_SCHEME],Err)
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err) !Enable 3D interpolation on faces
  CALL cmfe_Basis_CreateFinish(QuadraticBasis,Err)

  !
  !================================================================================================================================
  !

  !Start the creation of a generated cylinder mesh
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_CYLINDER_MESH_TYPE,Err)
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,[QuadraticBasis,LinearBasis],Err)
  CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[INNER_RAD, OUTER_RAD, HEIGHT],Err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements],Err)
  
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !
  !================================================================================================================================
  !

  !Create a decomposition
  CALL cmfe_RandomSeedsSet(0_CMISSIntg,Err) !To keep the automatic decomposition same each time
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Automatic decomposition
!   CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
!   CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  !Manual decomposition
  IF(NumberOfDomains>1) THEN
    CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_USER_DEFINED_TYPE,Err)
    !Set all elements but last one to first domain
    CALL cmfe_Mesh_NumberOfElementsGet(Mesh,NE,Err)
    do E=1,NE/2
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,E,0,Err)
    enddo
    do E=NE/2+1,NE
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,E,1,Err)
    enddo
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  ENDIF
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
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_CreateFinish(GeometricFieldDarcy,Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricFieldDarcy,Err)

  !
  !================================================================================================================================
  !

  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)        
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricFieldSolid,Err)
  CALL cmfe_Field_NumberOfVariablesSet(FibreField,FieldFibreNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,FieldFibreNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,3,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !
  !================================================================================================================================
  !

  !Create a material field and attach it to the geometric field
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL cmfe_Field_TypeSet(MaterialField,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialField,Decomposition,Err)        
  CALL cmfe_Field_GeometricFieldSet(MaterialField,GeometricFieldSolid,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialField,FieldMaterialNumberOfVariables,Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,FieldMaterialNumberOfComponents,Err)  
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,2,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_ComponentInterpolationSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,3,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  CALL cmfe_Field_CreateFinish(MaterialField,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,C1,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,C2,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,C3,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for ALE Darcy
  CALL cmfe_Field_Initialise(EquationsSetFieldDarcy,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetDarcy,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberDarcy,Region,GeometricFieldDarcy, &
    & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMFE_EQUATIONS_SET_DARCY_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_INCOMPRESS_ELASTICITY_DRIVEN_DARCY_SUBTYPE],EquationsSetFieldUserNumberDarcy,EquationsSetFieldDarcy, &
    & EquationsSetDarcy,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetDarcy,Err)

  !Create the equations set for the solid
  CALL cmfe_Field_Initialise(EquationsSetFieldSolid,Err)
  CALL cmfe_EquationsSet_Initialise(EquationsSetSolid,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumberSolid,Region,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_INCOMPRESS_ELASTICITY_DRIVEN_DARCY_SUBTYPE], &
    & EquationsSetFieldUserNumberSolid,EquationsSetFieldSolid,EquationsSetSolid,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create a dependent field with four variables (U, DelUDelN = solid, V, DelVDelN = Darcy) and four components
  !Solid: The U, DelUDelN variables have 4 components (3 displacement, 1 pressure)
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL cmfe_Field_TypeSet(DependentField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentField,GeometricFieldSolid,Err)
  CALL cmfe_Field_DependentTypeSet(DependentField,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentField,FieldDependentNumberOfVariables,Err)

  CALL cmfe_Field_VariableTypesSet(DependentField,[CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_DELVDELN_VARIABLE_TYPE],Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,FieldDependentSolidNumberOfComponents,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,FieldDependentSolidNumberOfComponents,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_V_VARIABLE_TYPE,FieldDependentFluidNumberOfComponents,Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_DELVDELN_VARIABLE_TYPE,FieldDependentFluidNumberOfComponents,Err)

  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,QuadraticMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,LinearMeshComponentNumber,Err)

  DarcyVelMeshComponentNumber = LinearMeshComponentNumber
  DarcyMassIncreaseMeshComponentNumber = LinearMeshComponentNumber

  !Darcy: The V, DelVDelN variables have 4 components (3 velocities, 1 mass increase)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_V_VARIABLE_TYPE,1,DarcyVelMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_V_VARIABLE_TYPE,2,DarcyVelMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_V_VARIABLE_TYPE,3,DarcyVelMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_V_VARIABLE_TYPE,4,DarcyMassIncreaseMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELVDELN_VARIABLE_TYPE,1,DarcyVelMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELVDELN_VARIABLE_TYPE,2,DarcyVelMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELVDELN_VARIABLE_TYPE,3,DarcyVelMeshComponentNumber,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELVDELN_VARIABLE_TYPE,4, &
    & DarcyMassIncreaseMeshComponentNumber,Err)

  CALL cmfe_Field_ScalingTypeSet(DependentField,CMFE_FIELD_UNIT_SCALING,Err)

  CALL cmfe_Field_CreateFinish(DependentField,Err)

  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetSolid,FieldDependentUserNumber,DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetSolid,Err)

  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetSolid,FieldMaterialUserNumber,MaterialField,Err)  
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetSolid,Err)


  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetDarcy,FieldDependentUserNumber,DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetDarcy,Err)

  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS+1
    CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for ALE Darcy
  CALL cmfe_Field_Initialise(MaterialsFieldDarcy,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetDarcy,MaterialsFieldUserNumberDarcy, &
    & MaterialsFieldDarcy,Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetDarcy,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDarcy,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)

  !
  !================================================================================================================================
  !
!   !Set up analytic field
!   CALL cmfe_Field_Initialise(AnalyticField,Err)
!   CALL cmfe_EquationsSet_AnalyticCreateStart(EquationsSetSolid,CMFE_EQUATIONS_SET_FINITE_ELASTICITY_CYLINDER, &
!     & FieldAnalyticUserNumber,AnalyticField,Err)
!   !Finish the equations set analytic field variables
!   CALL cmfe_EquationsSet_AnalyticCreateFinish(EquationsSetSolid,Err)
! 
!   !Set the analytic parameters
!   CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSetSolid,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_PIN_IDX,INNER_PRESSURE,Err)
!   CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSetSolid,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_POUT_IDX,OUTER_PRESSURE,Err)
!   CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSetSolid,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_LAMBDA_IDX,LAMBDA,Err)
!   CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSetSolid,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_TSI_IDX,TSI,Err)
!   CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSetSolid,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_RIN_IDX,INNER_RAD,Err)
!   CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSetSolid,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_ROUT_IDX,OUTER_RAD,Err)
!   CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSetSolid,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C1_IDX,C1,Err)
!   CALL cmfe_EquationsSet_AnalyticUserParamSet(EquationsSetSolid,CMFE_FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C2_IDX,C2,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SET EQUATIONS

  !Darcy
  CALL cmfe_Equations_Initialise(EquationsDarcy,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetDarcy,EquationsDarcy,Err)
  CALL cmfe_Equations_SparsityTypeSet(EquationsDarcy,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(EquationsDarcy,EQUATIONS_DARCY_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetDarcy,Err)

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
    & 1,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
!   CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,-14.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,0.0_CMISSRP, &
    & Err)

  !
  !================================================================================================================================
  !

  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_DARCY_TYPE, &
    & CMFE_PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
!   CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoop,1,Err)  ! this one sets the increment loop counter
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

  !Darcy
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMFE_CONTROL_LOOP_NODE], &
    & SolverDarcyIndex,DynamicSolverDarcy,Err)
  CALL cmfe_Solver_OutputTypeSet(DynamicSolverDarcy,DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE,Err)
  CALL cmfe_Solver_DynamicThetaSet(DynamicSolverDarcy,DYNAMIC_SOLVER_DARCY_THETA,Err)
!   CALL cmfe_SolverDynamicDynamicSet(DynamicSolverDarcy,.TRUE.,Err)
  CALL cmfe_Solver_DynamicLinearSolverGet(DynamicSolverDarcy,LinearSolverDarcy,Err)
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

CALL cmfe_Problem_SolversCreateFinish(Problem,Err)


  !
  !================================================================================================================================
  !

  !Create the problem solver equations
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
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsDarcy,EquationsSetDarcy,EquationsSetIndex,Err)
  !
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

! MANUAL BC ASSIGNMENT - REPLACED BY THE ANALYTIC BC ROUTINE BELOW
  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsSolid,BoundaryConditions,Err)

  !Get surfaces - will fix two nodes on bottom face, pressure conditions inside
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_CYLINDER_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi, &
    & Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_CYLINDER_TOP_SURFACE,TopSurfaceNodes,TopNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_CYLINDER_INNER_SURFACE,InnerSurfaceNodes,InnerNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_CYLINDER_OUTER_SURFACE,OuterSurfaceNodes,OuterNormalXi,Err)

!   !Set all inner surface nodes to inner pressure
!   DO NN=1,SIZE(InnerSurfaceNodes,1)
!       CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,InnerSurfaceNodes(NN), &
!       & abs(InnerNormalXi),CMFE_BOUNDARY_CONDITION_PRESSURE_INCREMENTED,INNER_PRESSURE,Err)   ! INNER_PRESSURE
!     IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER PRESSURE TO NODE", InnerSurfaceNodes(NN)
!   ENDDO
!         
!   !Set all outer surface nodes to outer pressure
!   DO NN=1,SIZE(OuterSurfaceNodes,1)
!     CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,OuterSurfaceNodes(NN), &
!       & abs(OuterNormalXi),CMFE_BOUNDARY_CONDITION_PRESSURE_INCREMENTED,OUTER_PRESSURE,Err)
!     IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER PRESSURE TO NODE", OuterSurfaceNodes(NN)
!   ENDDO

!   !Set all top nodes fixed in z plane at the set height
!   DO NN=1,SIZE(TopSurfaceNodes,1)
!     CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,TopSurfaceNodes(NN), &
!       & 3,CMFE_BOUNDARY_CONDITION_FIXED,deformedHeight,Err)
!     IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING FIXED CONDITION TO NODE", TopSurfaceNodes(NN)
!   ENDDO

  !Set all bottom nodes fixed in z plane
  DO NN=1,SIZE(BottomSurfaceNodes,1)
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,BottomSurfaceNodes(NN), &
      & 3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING FIXED CONDITION TO NODE", BottomSurfaceNodes(NN)
  ENDDO

  !Set two nodes on the bottom surface to axial displacement only
  X_FIXED=.FALSE.
  Y_FIXED=.FALSE.
  DO NN=1,SIZE(BottomSurfaceNodes,1)
    IF (.NOT.X_FIXED) THEN
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & 1,1,BottomSurfaceNodes(NN),1,xValue,Err)
      IF(abs(xValue)<1e-5_CMISSRP) THEN
        !Constrain it in x direction
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
          & BottomSurfaceNodes(NN),1,&
          & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
        X_FIXED=.TRUE.
        WRITE(*,*) "CyliderInflationExample: SUCCESSFULLY CONSTRAINED IN X DIRECTION NODE",BottomSurfaceNodes(NN)
    ENDIF
    ENDIF
    IF(.NOT.Y_FIXED) THEN
      CALL cmfe_Field_ParameterSetGetNode(GeometricFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & 1,1,BottomSurfaceNodes(NN),2,yValue,Err)
      IF(abs(yValue)<1e-5_CMISSRP) THEN
        !Constrain it in y direction
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
          & BottomSurfaceNodes(NN),2,&
          & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
        Y_FIXED=.TRUE.
        WRITE(*,*) "CyliderInflationExample: SUCCESSFULLY CONSTRAINED IN Y DIRECTION NODE",BottomSurfaceNodes(NN)
    ENDIF
    ENDIF
    IF (X_FIXED.AND.Y_FIXED) EXIT
  ENDDO
  !Check
  IF(.NOT.X_FIXED .OR. .NOT.Y_FIXED) THEN
    Write(*,*) "Couldn't fix bottom surface. No node lies on x or y axis, try changing number of elements"// &
      & " in theta coordinate"
    STOP
    ENDIF

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsSolid,Err)
! END OF MANUAL BC ASSIGNMENT

!   !Set the bc using the analytic solution routine
!   CALL cmfe_EquationsSetBoundaryConditionsAnalytic(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  !BCs Darcy
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsDarcy,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)

    !top surface
    CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMFE_GENERATED_MESH_CYLINDER_TOP_SURFACE, &
      & TopSurfaceNodesDarcyVel,TopNormalXi,Err)

    !Set all top surface nodes to Darcy inflow BC
    DO NN=1,SIZE(TopSurfaceNodesDarcyVel,1)
      VALUE = -1.0_CMISSRP
      COMPONENT_NUMBER = 3
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMFE_FIELD_V_VARIABLE_TYPE,1,1, &
        & TopSurfaceNodesDarcyVel(NN),COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)

      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING TOP DARCY BC TO NODE", TopSurfaceNodesDarcyVel(NN)
    ENDDO


    !bottom surface
    CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMFE_GENERATED_MESH_CYLINDER_BOTTOM_SURFACE, &
      & BottomSurfaceNodesDarcyVel,BottomNormalXi,Err)

    !Set all bottom surface nodes impermeable
    DO NN=1,SIZE(BottomSurfaceNodesDarcyVel,1)
      VALUE = 0.0_CMISSRP
      COMPONENT_NUMBER = 3
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMFE_FIELD_V_VARIABLE_TYPE,1,1, &
        & BottomSurfaceNodesDarcyVel(NN),COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)

      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING BOTTOM DARCY BC TO NODE", BottomSurfaceNodesDarcyVel(NN)
    ENDDO


    !inner surface
    CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMFE_GENERATED_MESH_CYLINDER_INNER_SURFACE, &
      & InnerSurfaceNodesDarcyVel,InnerNormalXi,Err)

    !Set all inner surface nodes impermeable
    DO NN=1,SIZE(InnerSurfaceNodesDarcyVel,1)
      VALUE = 0.0_CMISSRP
      COMPONENT_NUMBER = 1
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMFE_FIELD_V_VARIABLE_TYPE,1,1, &
        & InnerSurfaceNodesDarcyVel(NN),COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)

      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)

      VALUE = 0.0_CMISSRP
      COMPONENT_NUMBER = 2
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMFE_FIELD_V_VARIABLE_TYPE,1,1, &
        & InnerSurfaceNodesDarcyVel(NN),COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)

      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING INNER DARCY BC TO NODE", InnerSurfaceNodesDarcyVel(NN)
    ENDDO


    !outer surface
    CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMFE_GENERATED_MESH_CYLINDER_OUTER_SURFACE, &
      & OuterSurfaceNodesDarcyVel,OuterNormalXi,Err)

    !Set all outer surface nodes impermeable
    DO NN=1,SIZE(OuterSurfaceNodesDarcyVel,1)
      VALUE = 0.0_CMISSRP
      COMPONENT_NUMBER = 1
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMFE_FIELD_V_VARIABLE_TYPE,1,1, &
        & OuterSurfaceNodesDarcyVel(NN),COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)

      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)

      VALUE = 0.0_CMISSRP
      COMPONENT_NUMBER = 2
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentField,CMFE_FIELD_V_VARIABLE_TYPE,1,1, &
        & OuterSurfaceNodesDarcyVel(NN),COMPONENT_NUMBER,CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)

      IF(Err/=0) WRITE(*,*) "ERROR WHILE ASSIGNING OUTER DARCY BC TO NODE", OuterSurfaceNodesDarcyVel(NN)
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

!   !Output Analytic analysis
!   Call cmfe_AnalyticAnalysisOutput(DependentField,"outputs/CylinderInflation",Err)

  !
  !================================================================================================================================
  !

  !Output solution  
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"outputs/CylinderInflation","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"outputs/CylinderInflation","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)

  !
  !================================================================================================================================
  !

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM CYLINDERINFLATIONDRIVENDARCYEXAMPLE

