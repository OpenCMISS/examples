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

!> \example MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDrivenDarcyMatSolve/src/IncompressibleElasticityDrivenDarcyMatSolveExample.f90
!! Example program to solve coupled FiniteElasticityDarcy equations using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDrivenDarcyMatSolve/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDrivenDarcyMatSolve/build-intel'>Linux GNU Build</a>
!!
!<

! !
! !  This example considers a coupled Finite Elasticity Darcy problem
! !

!> Main program

PROGRAM FINITEELASTICITYDARCYEXAMPLE

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
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberMatProperties=42
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcy=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMatProperties=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberDarcy=12
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberMatProperties=13
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=14
!   INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberSolid=15
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberDarcy=22
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberMatProperties=23

  INTEGER(CMISSIntg), PARAMETER :: ControlLoopSolidNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopFluidNumber=2
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopSubiterationNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverSolidIndex=1
  INTEGER(CMISSIntg), PARAMETER :: SolverMatPropertiesIndex=1
  INTEGER(CMISSIntg), PARAMETER :: SolverDarcyIndex=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPorosity=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPermOverVis=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMatPropertiesPorosity=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMatPropertiesPermOverVis=2

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=3

  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS

  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS

  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE

  INTEGER(CMISSIntg) :: EQUATIONS_DARCY_OUTPUT
  INTEGER(CMISSIntg) :: EQUATIONS_MAT_PROPERTIES_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_DARCY_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_DARCY_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_MAT_PROPERTIES_OUTPUT_TYPE

  REAL(CMISSDP) :: COORD_X, COORD_Y, COORD_Z
  REAL(CMISSDP) :: DOMAIN_X1, DOMAIN_X2, DOMAIN_Y1, DOMAIN_Y2, DOMAIN_Z1, DOMAIN_Z2
  REAL(CMISSDP) :: GEOMETRY_TOLERANCE
  INTEGER(CMISSIntg) :: EDGE_COUNT
  INTEGER(CMISSIntg) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SOLID
  REAL(CMISSDP) :: INITIAL_FIELD_DARCY(4)
  REAL(CMISSDP) :: INITIAL_FIELD_MAT_PROPERTIES(3)
  REAL(CMISSDP) :: INITIAL_FIELD_SOLID(4)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE
  REAL(CMISSDP) :: POROSITY_PARAM_MAT_PROPERTIES, PERM_OVER_VIS_PARAM_MAT_PROPERTIES
  REAL(CMISSDP) :: POROSITY_PARAM_DARCY, PERM_OVER_VIS_PARAM_DARCY

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_DARCY_DIRECT_FLAG
  LOGICAL :: LINEAR_SOLVER_MAT_PROPERTIES_DIRECT_FLAG

  !CMISS variables

  !Regions
  TYPE(CMISSRegionType) :: Region
  TYPE(CMISSRegionType) :: WorldRegion
  !Coordinate systems
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  !Basis
  TYPE(CMISSBasisType) :: BasisGeometry
  TYPE(CMISSBasisType) :: CubicBasis, QuadraticBasis, LinearBasis, Bases(2)
  !Nodes
  TYPE(CMISSNodesType) :: Nodes
  !Meshes
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh

  !Decompositions
  TYPE(CMISSDecompositionType) :: Decomposition
  !Fields
  TYPE(CMISSFieldsType) :: Fields
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: DependentFieldMatProperties
  TYPE(CMISSFieldType) :: MaterialsFieldDarcy
  TYPE(CMISSFieldType) :: MaterialsFieldMatProperties
  TYPE(CMISSFieldType) :: EquationsSetFieldDarcy
  TYPE(CMISSFieldType) :: EquationsSetFieldMatProperties
!   TYPE(CMISSFieldType) :: IndependentFieldSolid
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsDarcy
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsMatProperties
  !Equations sets
  TYPE(CMISSEquationsSetType) :: EquationsSetDarcy
  TYPE(CMISSEquationsSetType) :: EquationsSetMatProperties
  !Equations
  TYPE(CMISSEquationsType) :: EquationsDarcy
  TYPE(CMISSEquationsType) :: EquationsMatProperties
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: DynamicSolverDarcy
  TYPE(CMISSSolverType) :: LinearSolverDarcy
  TYPE(CMISSSolverType) :: LinearSolverMatProperties
!   TYPE(CMISSSolverType) :: LinearSolverSolid
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsDarcy
  TYPE(CMISSSolverEquationsType) :: SolverEquationsMatProperties

  !Other variables
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face1Nodes(:),Face2Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face3Nodes(:),Face4Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face5Nodes(:),Face6Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face7Nodes(:),Face8Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face9Nodes(:),Face10Nodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE,TARGET :: Face11Nodes(:),Face12Nodes(:)
  INTEGER(CMISSIntg),POINTER :: BoundaryNodes(:)
  INTEGER(CMISSIntg) :: FaceXi(6)
  INTEGER(CMISSIntg) :: I,VariableType
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

  INTEGER(CMISSIntg) :: BASIS_NUMBER_SOLID

  INTEGER(CMISSIntg) :: TotalNumberOfSolidNodes
  INTEGER(CMISSIntg) :: SolidMeshComponenetNumber
  INTEGER(CMISSIntg) :: SolidPressureMeshComponenetNumber

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

  TYPE(CMISSBasisType) :: BasisSolid
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsSolid
  TYPE(CMISSEquationsType) :: EquationsSolid
  TYPE(CMISSEquationsSetType) :: EquationsSetSolid
  TYPE(CMISSFieldType) :: GeometricFieldSolid,FibreFieldSolid,MaterialFieldSolid
  TYPE(CMISSFieldType) :: DependentFieldSolid,EquationsSetFieldSolid
  TYPE(CMISSSolverType) :: SolverSolid
  TYPE(CMISSSolverEquationsType) :: SolverEquationsSolid
  TYPE(CMISSMeshElementsType) :: MeshElementsSolid

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

  !Import cmHeart mesh information
  !CALL FLUID_MECHANICS_IO_READ_CMHEART(CM,Err)
!   BASIS_XI_INTERPOLATION_SOLID=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  BASIS_XI_INTERPOLATION_SOLID=CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
  !Set geometric tolerance
  GEOMETRY_TOLERANCE = 1.0E-12_CMISSDP
  !Set initial values
  INITIAL_FIELD_DARCY(1)=0.0_CMISSDP
  INITIAL_FIELD_DARCY(2)=0.0_CMISSDP
  INITIAL_FIELD_DARCY(3)=0.0_CMISSDP
  INITIAL_FIELD_DARCY(4)=0.0_CMISSDP
  INITIAL_FIELD_MAT_PROPERTIES(1)=0.0_CMISSDP
  INITIAL_FIELD_MAT_PROPERTIES(2)=0.0_CMISSDP
  INITIAL_FIELD_MAT_PROPERTIES(3)=0.0_CMISSDP
!   INITIAL_FIELD_SOLID(1)=1.0_CMISSDP
!   INITIAL_FIELD_SOLID(2)=1.0_CMISSDP
!   INITIAL_FIELD_SOLID(3)=1.0_CMISSDP
!   INITIAL_FIELD_SOLID(4)=1.0_CMISSDP
  !Set material parameters
  POROSITY_PARAM_DARCY=0.1_CMISSDP
!   PERM_OVER_VIS_PARAM_DARCY=1.0e-1_CMISSDP
  PERM_OVER_VIS_PARAM_DARCY=1.0_CMISSDP
!   PERM_OVER_VIS_PARAM_DARCY=0.1_CMISSDP
  POROSITY_PARAM_MAT_PROPERTIES=POROSITY_PARAM_DARCY
  PERM_OVER_VIS_PARAM_MAT_PROPERTIES=PERM_OVER_VIS_PARAM_DARCY
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_MAT_PROPERTIES_OUTPUT_TYPE=CMISS_SOLVER_PROGRESS_OUTPUT
  DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE=CMISS_SOLVER_PROGRESS_OUTPUT
  LINEAR_SOLVER_DARCY_OUTPUT_TYPE=CMISS_SOLVER_SOLVER_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DARCY_OUTPUT=CMISS_EQUATIONS_NO_OUTPUT
  EQUATIONS_MAT_PROPERTIES_OUTPUT=CMISS_EQUATIONS_NO_OUTPUT

  !Set time parameter
  DYNAMIC_SOLVER_DARCY_START_TIME=0.0_CMISSDP
!   DYNAMIC_SOLVER_DARCY_STOP_TIME=0.03_CMISSDP
  DYNAMIC_SOLVER_DARCY_TIME_INCREMENT=1.0e-3_CMISSDP
  DYNAMIC_SOLVER_DARCY_STOP_TIME=2_CMISSIntg * DYNAMIC_SOLVER_DARCY_TIME_INCREMENT
  DYNAMIC_SOLVER_DARCY_THETA=1.0_CMISSDP !2.0_CMISSDP/3.0_CMISSDP
  !Set result output parameter
  DYNAMIC_SOLVER_DARCY_OUTPUT_FREQUENCY=1
  !Set solver parameters
  LINEAR_SOLVER_MAT_PROPERTIES_DIRECT_FLAG=.TRUE.
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

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  !
  !================================================================================================================================
  !

  !Set diagnostics

  DIAG_LEVEL_LIST(1)=1
  DIAG_LEVEL_LIST(2)=2
  DIAG_LEVEL_LIST(3)=3
  DIAG_LEVEL_LIST(4)=4
  DIAG_LEVEL_LIST(5)=5

!   DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_FINITE_ELEMENT_CALCULATE"
!   DIAG_ROUTINE_LIST(2)="DARCY_EQUATION_PRE_SOLVE_STORE_REFERENCE_DATA"
!   DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_PRE_SOLVE_GET_SOLID_DISPLACEMENT"
!   DIAG_ROUTINE_LIST(2)="DARCY_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH"
!   DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_PRE_SOLVE_GET_SOLID_DISPLACEMENT"
!   DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS"
!   DIAG_ROUTINE_LIST(1)="FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE"
!   DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_POST_SOLVE_ADD_MASS_CORRECTION"
  DIAG_ROUTINE_LIST(1)="WRITE_IP_INFO"
!   DIAG_ROUTINE_LIST(2)="FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR"
!   DIAG_ROUTINE_LIST(3)="EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION"
!   DIAG_ROUTINE_LIST(5)="DARCY_EQUATION_PRE_SOLVE_MAT_PROPERTIES"
!   DIAG_ROUTINE_LIST(6)="FITTING_FINITE_ELEMENT_CALCULATE"
!   DIAG_ROUTINE_LIST(7)="FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE"
!   DIAG_ROUTINE_LIST(8)="FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE"
!   DIAG_ROUTINE_LIST(1)="PROBLEM_SOLVER_EQUATIONS_SOLVE"
!   DIAG_ROUTINE_LIST(1)="SOLVER_NEWTON_SOLVE"
!   DIAG_ROUTINE_LIST(2)="SOLVER_NEWTON_LINESEARCH_SOLVE"
!   DIAG_ROUTINE_LIST(1)="SOLVER_SOLUTION_UPDATE"
!   DIAG_ROUTINE_LIST(1)="FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE"

  !CMISS_ALL_DIAG_TYPE/CMISS_IN_DIAG_TYPE/CMISS_FROM_DIAG_TYPE
  CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISS_ALL_TIMING_TYPE/CMISS_IN_TIMING_TYPE/CMISS_FROM_TIMING_TYPE
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

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !REGION
  !For a volume-coupled problem, solid and fluid are based in the same region

  !Start the creation of a new region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system as defined above
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !BASES
  !Define basis functions
  CALL CMISSBasis_Initialise(LinearBasis,Err)
  CALL CMISSBasis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(LinearBasis, &
    & (/CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME/),Err)
  !CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err)
  CALL CMISSBasis_CreateFinish(LinearBasis,Err)

  CALL CMISSBasis_Initialise(QuadraticBasis,Err)
  CALL CMISSBasis_CreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasis_InterpolationXiSet(QuadraticBasis,(/CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION, &
    & CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & (/CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME/),Err)
  !CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err)
  CALL CMISSBasis_CreateFinish(QuadraticBasis,Err)

  CALL CMISSBasis_Initialise(CubicBasis,Err)
  CALL CMISSBasis_CreateStart(CubicBasisUserNumber,CubicBasis,Err)
  CALL CMISSBasis_InterpolationXiSet(CubicBasis,(/CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION, &
    & CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(CubicBasis, &
    & (/CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME,CMISS_BASIS_HIGH_QUADRATURE_SCHEME/),Err)
  !CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(CubicBasis,.true.,Err) !Enable 3D interpolation on faces
  CALL CMISSBasis_CreateFinish(CubicBasis,Err)

!   Bases(1)=QuadraticBasis
!   Bases(2)=LinearBasis
  Bases(1)=CubicBasis
  Bases(2)=QuadraticBasis

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  !CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,BasisGeometry,Err)
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,Bases,Err)
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/X_DIM,Y_DIM/),Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/),Err)
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/X_DIM,Y_DIM,Z_DIM/),Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS/),Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

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

  !Start to create a default (geometric) field on the region
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
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
!   !Set the scaling to use
!   CALL CMISSField_ScalingTypeSet(GeometricField,CMISS_FIELD_NO_SCALING,Err)

!   !Update the geometric field parameters
!   DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
!     DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!       VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
!       CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
!         & CMISS_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
!     ENDDO
!   ENDDO
!   CALL CMISSField_ParameterSetUpdateStart(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
!   CALL CMISSField_ParameterSetUpdateFinish(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

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
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricFieldSolid,Err)
!   !Set the scaling to use
!   CALL CMISSField_ScalingTypeSet(GeometricFieldSolid,CMISS_FIELD_NO_SCALING,Err)

! !---
!   !Update the geometric field parameters
!   DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
!     DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!       VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
!       CALL CMISSField_ParameterSetUpdateNode(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
!         & CMISS_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
!     ENDDO
!   ENDDO
!   CALL CMISSField_ParameterSetUpdateStart(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
!   CALL CMISSField_ParameterSetUpdateFinish(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
!---

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

  !Create the equations set for ALE Darcy
  CALL CMISSField_Initialise(EquationsSetFieldDarcy,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSetDarcy,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumberDarcy,Region,GeometricField,CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & CMISS_EQUATIONS_SET_DARCY_EQUATION_TYPE,CMISS_EQUATIONS_SET_INCOMPRESS_ELASTICITY_DRIVEN_DARCY_SUBTYPE,&
    & EquationsSetFieldUserNumberDarcy,EquationsSetFieldDarcy,EquationsSetDarcy,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSetDarcy,Err)

  !Create the equations set for deformation-dependent material properties
  CALL CMISSField_Initialise(EquationsSetFieldMatProperties,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSetMatProperties,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumberMatProperties,Region,GeometricField,CMISS_EQUATIONS_SET_FITTING_CLASS,&
    & CMISS_EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE,CMISS_EQUATIONS_SET_MAT_PROP_INRIA_MODEL_DATA_FITTING_SUBTYPE,&
    & EquationsSetFieldUserNumberMatProperties,EquationsSetFieldMatProperties,EquationsSetMatProperties,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSetMatProperties,Err)

  !Create the equations set for the solid
  CALL CMISSField_Initialise(EquationsSetFieldSolid,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSetSolid,Err)
  CALL CMISSEquationsSet_CreateStart(EquationSetSolidUserNumber,Region,FibreFieldSolid,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_INCOMPRESS_ELASTICITY_DRIVEN_DARCY_SUBTYPE,&
    & EquationsSetFieldSolidUserNumber,EquationsSetFieldSolid,EquationsSetSolid,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSetSolid,Err)

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


  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for deformation-dependent material properties
  CALL CMISSField_Initialise(DependentFieldMatProperties,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetMatProperties,DependentFieldUserNumberMatProperties, &
    & DependentFieldMatProperties,Err)
  NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES = 2
  DO COMPONENT_NUMBER=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES
    ! In 'FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP' the MESH_COMPONENT for the DEPENDENT_FIELD 
    ! of EQUATIONS_SET_MAT_PROPERTIES is currently being defaulted to GEOMETRIC_MESH_COMPONENT. Therfore:
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldMatProperties,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, &
      & SolidGeometryMeshComponentNumber,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldMatProperties,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, &
      & SolidGeometryMeshComponentNumber,Err)
  ENDDO
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetMatProperties,Err)

  !Initialise dependent field
  DO COMPONENT_NUMBER=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES
    CALL CMISSField_ComponentValuesInitialise(DependentFieldMatProperties,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & COMPONENT_NUMBER,INITIAL_FIELD_MAT_PROPERTIES(COMPONENT_NUMBER),Err)
  ENDDO

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create a dependent field with four variables (U, DelUDelN = solid, V, DelVDelN = Darcy) and four components
  CALL CMISSField_Initialise(DependentFieldSolid,Err)
  !
  CALL CMISSField_CreateStart(FieldDependentSolidUserNumber,Region,DependentFieldSolid,Err)
  !
  CALL CMISSField_TypeSet(DependentFieldSolid,CMISS_FIELD_GENERAL_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(DependentFieldSolid,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(DependentFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSField_DependentTypeSet(DependentFieldSolid,CMISS_FIELD_DEPENDENT_TYPE,Err)
  CALL CMISSField_NumberOfVariablesSet(DependentFieldSolid,FieldDependentSolidNumberOfVariables,Err)
  CALL CMISSField_VariableTypesSet(DependentFieldSolid,(/CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_DELVDELN_VARIABLE_TYPE/),Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldSolid,CMISS_FIELD_DELUDELN_VARIABLE_TYPE, &
    & FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,FieldDependentFluidNumberOfComponents,Err)
  CALL CMISSField_NumberOfComponentsSet(DependentFieldSolid,CMISS_FIELD_DELVDELN_VARIABLE_TYPE, &
    & FieldDependentFluidNumberOfComponents,Err)
  !
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,SolidDisplMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,2,SolidDisplMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,3,SolidDisplMeshComponentNumber,Err)
  CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,4, &
    & CMISS_FIELD_NODE_BASED_INTERPOLATION, &
    & Err)
!   CALL CMISSField_ComponentInterpolationSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
!   CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,4,SolidMeshComponenetNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,4,SolidLagrMultMeshComponentNumber,Err)
  !
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

  !For this equation type, MESH_COMPONENT_NUMBER_PRESSURE is actually the mass increase component as the pressure is taken from the solid equations
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,DarcyVelMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,2,DarcyVelMeshComponentNumber,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,3,DarcyVelMeshComponentNumber,Err)
!   CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,4,MESH_COMPONENT_NUMBER_PRESSURE,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,4, &
    & DarcyMassIncreaseMeshComponentNumber, &
    & Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1,DarcyVelMeshComponentNumber, &
    & Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,2,DarcyVelMeshComponentNumber, &
    & Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,3,DarcyVelMeshComponentNumber, &
    & Err)
!   CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,4,MESH_COMPONENT_NUMBER_PRESSURE,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldSolid,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,4, &
    & DarcyMassIncreaseMeshComponentNumber,Err)

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


  !Create the equations set dependent field variables for ALE Darcy
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetDarcy,FieldDependentSolidUserNumber, & ! ??? UserNumber ???
    & DependentFieldSolid,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetDarcy,Err)

  !Initialise dependent field (velocity components,pressure,mass increase)
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS+1
    CALL CMISSField_ComponentValuesInitialise(DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
  ENDDO


  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for ALE Darcy
  CALL CMISSField_Initialise(MaterialsFieldDarcy,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetDarcy,MaterialsFieldUserNumberDarcy, &
    & MaterialsFieldDarcy,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetDarcy,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDarcy,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)
  !Create the equations set materials field variables for deformation-dependent material properties
  CALL CMISSField_Initialise(MaterialsFieldMatProperties,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetMatProperties,MaterialsFieldUserNumberMatProperties, &
    & MaterialsFieldMatProperties,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetMatProperties,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldMatProperties,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberMatPropertiesPorosity,POROSITY_PARAM_MAT_PROPERTIES,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldMatProperties,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberMatPropertiesPermOverVis,PERM_OVER_VIS_PARAM_MAT_PROPERTIES,Err)


  !
  !================================================================================================================================
  !

!   !INDEPENDENT FIELDS
! 
!   !Create the equations set independent field variables for the solid
!   CALL CMISSField_Initialise(IndependentFieldSolid,Err)
!   CALL CMISSEquationsSet_IndependentCreateStart(EquationsSetSolid,IndependentFieldUserNumberSolid, &
!     & IndependentFieldSolid,Err)
!   !Set the mesh component to be used by the field components.
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!     CALL CMISSField_ComponentMeshComponentSet(IndependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, &
!       & SolidDisplMeshComponentNumber,Err)
!   ENDDO
!   !Finish the equations set independent field variables
!   CALL CMISSEquationsSet_IndependentCreateFinish(EquationsSetSolid,Err)

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

  !MatProperties
  CALL CMISSEquations_Initialise(EquationsMatProperties,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetMatProperties,EquationsMatProperties,Err)
  CALL CMISSEquations_SparsityTypeSet(EquationsMatProperties,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(EquationsMatProperties,EQUATIONS_MAT_PROPERTIES_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetMatProperties,Err)

  !Solid
  CALL CMISSEquations_Initialise(EquationsSolid,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetSolid,EquationsSolid,Err)
  CALL CMISSEquations_SparsityTypeSet(EquationsSolid,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(EquationsSolid,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetSolid,Err)

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
    & CMISS_PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE,Err)
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

  CALL CMISSSolver_Initialise(SolverSolid,Err)
  CALL CMISSSolver_Initialise(LinearSolverMatProperties,Err)
  CALL CMISSSolver_Initialise(DynamicSolverDarcy,Err)
  CALL CMISSSolver_Initialise(LinearSolverDarcy,Err)

  CALL CMISSProblem_SolversCreateStart(Problem,Err)

  ! Solid
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


  !MatProperties
  CALL CMISSProblem_SolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMISS_CONTROL_LOOP_NODE/), &
    & SolverMatPropertiesIndex,LinearSolverMatProperties,Err)
  CALL CMISSSolver_OutputTypeSet(LinearSolverMatProperties,LINEAR_SOLVER_MAT_PROPERTIES_OUTPUT_TYPE,Err)
  IF(LINEAR_SOLVER_MAT_PROPERTIES_DIRECT_FLAG) THEN
    CALL CMISSSolver_LinearTypeSet(LinearSolverMatProperties,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL CMISSSolver_LibraryTypeSet(LinearSolverMatProperties,CMISS_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL CMISSSolver_LinearTypeSet(LinearSolverMatProperties,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolverMatProperties,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolver_LinearIterativeDivergenceToleranceSet(LinearSolverMatProperties,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeRelativeToleranceSet(LinearSolverMatProperties,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(LinearSolverMatProperties,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeGMRESRestartSet(LinearSolverMatProperties,RESTART_VALUE,Err)
  ENDIF

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

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(SolverSolid,Err)
  CALL CMISSSolver_Initialise(LinearSolverMatProperties,Err)
  CALL CMISSSolver_Initialise(LinearSolverDarcy,Err)

  CALL CMISSSolverEquations_Initialise(SolverEquationsSolid,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsMatProperties,Err)
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
  !Get the deformation-dependent material properties solver equations
  CALL CMISSProblem_SolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMISS_CONTROL_LOOP_NODE/), &
    & SolverMatPropertiesIndex,LinearSolverMatProperties,Err)
  CALL CMISSSolver_SolverEquationsGet(LinearSolverMatProperties,SolverEquationsMatProperties,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsMatProperties,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsMatProperties,EquationsSetMatProperties,EquationsSetIndex,Err)
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

  !BOUNDARY CONDITIONS
  !Start the creation of the equations set boundary conditions for Darcy
!   CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsDarcy,Err)
!   CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  !Solid is computed in absolute position, rather than displacement. Thus BCs for absolute position
!   CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsSolid,Err)
!   CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsSolid,BoundaryConditionsSolid,Err)

  !--- BCs on normal velocity only
  CONDITION = CMISS_BOUNDARY_CONDITION_MOVED_WALL

!   IF( CM%D==2_CMISSIntg ) THEN !CM%D = number of dimensions, ie 2D
!     DO NODE_NUMBER=1_CMISSIntg,NUMBER_OF_NODES_GEOMETRY
!       COORD_X = CM%N(NODE_NUMBER,1_CMISSIntg)
!       COORD_Y = CM%N(NODE_NUMBER,2_CMISSIntg)
! 
!       IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) ) THEN
!         !x-velocity
!         VALUE = 1.0_CMISSDP
!         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
!           & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
!         !x-velocity
!         VALUE = 1.0_CMISSDP
!         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
!           & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
!         !y-velocity
!         VALUE = 2.0_CMISSDP
!         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
!           & NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
!         !y-velocity
!         VALUE = 2.0_CMISSDP
!         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
!           & NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
!       END IF
!     END DO
!   ELSE IF( CM%D==3_CMISSIntg ) THEN ! 3D geometry
!     DO NODE_NUMBER=1_CMISSIntg,NUMBER_OF_NODES_GEOMETRY  !What if different number of nodes geometry and velocity ?
!       COORD_X = CM%N(NODE_NUMBER,1_CMISSIntg)
!       COORD_Y = CM%N(NODE_NUMBER,2_CMISSIntg)
!       COORD_Z = CM%N(NODE_NUMBER,3_CMISSIntg)
! 
!       IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) ) THEN
! !         !x-velocity: F L U I D ( V Variable type )
!         VALUE = 0.0_CMISSDP
!         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
! ! !           & NODE_NUMBER,4_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err) !BC on pressure component
! ! !           & NODE_NUMBER,1_CMISSIntg,CMISS_BOUNDARY_CONDITION_MOVED_WALL,VALUE,Err) !inflow
!           & NODE_NUMBER,1_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err) !time-dependent inflow
! 
! !         !x-position: S O L I D ( U Variable Type)
! ! !         VALUE = 1.0_CMISSDP * DOMAIN_X1
! !         VALUE = COORD_X
! !         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
! !           & NODE_NUMBER,1_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! 
! !         EDGE_COUNT = 0_CMISSIntg
! !         IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) EDGE_COUNT = EDGE_COUNT + 1_CMISSIntg
! !         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) EDGE_COUNT = EDGE_COUNT + 1_CMISSIntg
! !         IF( (ABS(COORD_Z-DOMAIN_Z1) < GEOMETRY_TOLERANCE) ) EDGE_COUNT = EDGE_COUNT + 1_CMISSIntg
! !         IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) ) EDGE_COUNT = EDGE_COUNT + 1_CMISSIntg
! ! 
! !         IF(EDGE_COUNT == 2_CMISSIntg) THEN !it is a corner node
! !           VALUE = COORD_Y
! !           CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
! !             & NODE_NUMBER,2_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! ! 
! !           VALUE = COORD_Z
! !           CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
! !             & NODE_NUMBER,3_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! !       END IF
!       END IF
!       !
!       IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
! !         !x-velocity: F L U I D
!         VALUE = 0.0_CMISSDP
!         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
! ! !           & NODE_NUMBER,4_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err) !BC on pressure component
! ! !           & NODE_NUMBER,1_CMISSIntg,CMISS_BOUNDARY_CONDITION_MOVED_WALL,VALUE,Err) !impermeable wall, zero flux
!           & NODE_NUMBER,1_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err) !impermeable wall, zero flux
! 
! !         !x-position: S O L I D
! !         VALUE = COORD_X
! !         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
! !           & NODE_NUMBER,1_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! !         
! !         !Fix point 1
! !         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
! !           IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) ) THEN
! !             VALUE = COORD_Y
! !             CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
! !               & NODE_NUMBER,2_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! ! 
! !             VALUE = COORD_Z
! !             CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
! !               & NODE_NUMBER,3_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! !           END IF
! !         END IF
! !         !(Fix) point 2
! !         IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
! !           IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) ) THEN
! !             VALUE = COORD_Z
! !             CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
! !               & NODE_NUMBER,3_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! !           END IF
! !         END IF
! !         !(Fix) point 3
! !         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
! !           IF( (ABS(COORD_Z-DOMAIN_Z1) < GEOMETRY_TOLERANCE) ) THEN
! !             VALUE = COORD_Y
! !             CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
! !               & NODE_NUMBER,2_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! !           END IF
! !         END IF
! 
! 
!       END IF
!       !
!       IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
!         !y-velocity: F L U I D
!         VALUE = 0.0_CMISSDP
!         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
!           & NODE_NUMBER,2_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! ! 
! !         !y-position: S O L I D
! !         VALUE = 1.0_CMISSDP * DOMAIN_Y1
! !         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
! ! !           & NODE_NUMBER,2_CMISSIntg,CMISS_BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED,VALUE,Err)
! !           & NODE_NUMBER,2_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
!         !y-velocity: F L U I D
!         VALUE = 0.0_CMISSDP
!         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
!           & NODE_NUMBER,2_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! ! 
! ! !         !y-position: S O L I D
! ! !         VALUE = 1.1_CMISSDP * DOMAIN_Y2
! ! !         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
! ! !           & NODE_NUMBER,2_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Z-DOMAIN_Z1) < GEOMETRY_TOLERANCE) ) THEN
!         !z-velocity: F L U I D
!         !mass-correction: F L U I D
!         VALUE = 10.0_CMISSDP
!         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
!           & NODE_NUMBER,3_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! !           & NODE_NUMBER,3_CMISSIntg,CMISS_BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE,VALUE,Err)
! !           & NODE_NUMBER,4_CMISSIntg,CMISS_BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE,VALUE,Err)
! !           & NODE_NUMBER,4_CMISSIntg,CMISS_BOUNDARY_CONDITION_FREE,VALUE,Err)
! 
!         !z-position: S O L I D
!         VALUE = COORD_Z
!         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
!           & NODE_NUMBER,3_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) ) THEN
!         !z-velocity: F L U I D
!         VALUE = 0.0_CMISSDP
!         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
!           & NODE_NUMBER,3_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! 
! !         VALUE = COORD_X
! !         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
! !           & NODE_NUMBER,1_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! ! 
! !         VALUE = COORD_Y
! !         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
! !           & NODE_NUMBER,2_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! ! 
! !         VALUE = COORD_Z
! !         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
! !           & NODE_NUMBER,3_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! 
! !         !x-position: S O L I D
! !         VALUE = COORD_Z
! !         CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
! !           & NODE_NUMBER,3_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!         
!         !Fix point 1
!         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
!           IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
!             VALUE = COORD_Y
!             CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
!               & NODE_NUMBER,2_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
! 
!             VALUE = COORD_X
!             CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
!               & NODE_NUMBER,1_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!           END IF
!         END IF
!         !(Fix) point 2
!         IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
!           IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
!             VALUE = COORD_X
!             CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
!               & NODE_NUMBER,1_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!           END IF
!         END IF
!         !(Fix) point 3
!         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
!           IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) ) THEN
!             VALUE = COORD_Y
!             CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV, &
!               & NODE_NUMBER,2_CMISSIntg,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
!           END IF
!         END IF
! 
! 
!       END IF
!     END DO
!   END IF

  !Finish the creation of the equations set boundary conditions for Darcy
!   CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsDarcy,Err)
!   !Finish the creation of the equations set boundary conditions for the solid
!   CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsSolid,Err)

  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsMatProperties,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsMatProperties,BoundaryConditionsMatProperties,Err)
  !(No boundary conditions requrired for deformation-dependent material properties)
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsMatProperties,Err)







  !------------------------------------
  ! ASSIGN BOUNDARY CONDITIONS - SOLID
  ! Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsSolid,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsSolid,BoundaryConditionsSolid,Err)

  !Get surfaces 
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_FRONT_SURFACE, &
    & Face1Nodes,FaceXi(1),Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_BACK_SURFACE, &
    & Face2Nodes,FaceXi(2),Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_RIGHT_SURFACE, &
    & Face3Nodes,FaceXi(3),Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE, &
    & Face4Nodes,FaceXi(4),Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_TOP_SURFACE, &
    & Face5Nodes,FaceXi(5),Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,SolidDisplMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_BOTTOM_SURFACE, &
    & Face6Nodes,FaceXi(6),Err)

  ! Fix the bottom in z direction
  DO NN=1,SIZE(Face6Nodes,1)
    NODE=Face6Nodes(NN)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetGetNode(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE,3, &
        & ZCoord,Err)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,3, &
        & CMISS_BOUNDARY_CONDITION_FIXED,ZCoord,Err)
      WRITE(*,*) "FIXING NODE",NODE,"AT BOTTOM IN Z DIRECTION"
    ENDIF
  ENDDO

  ! Fix the top in z direction
  DO NN=1,SIZE(Face5Nodes,1)
    NODE=Face5Nodes(NN)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetGetNode(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE,3, &
        & ZCoord,Err)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,3, &
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
      CALL CMISSField_ParameterSetGetNode(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE,1, &
        & XCoord,Err)
      CALL CMISSField_ParameterSetGetNode(GeometricFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NODE,2, &
        & YCoord,Err)

      !Fix Origin displacement in x and y (z already fixed)
      IF(ABS(XCoord)<1.0E-6_CMISSDP.AND.ABS(YCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,1, &
          & CMISS_BOUNDARY_CONDITION_FIXED,XCoord,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,2, &
          & CMISS_BOUNDARY_CONDITION_FIXED,YCoord,Err)
        WRITE(*,*) "FIXING ORIGIN NODE",NODE,"IN X AND Y DIRECTION"
        X_FIXED=.TRUE.
        Y_FIXED=.TRUE.
      ENDIF

      !Fix nodal displacements at (X_DIM,0) in y
      IF(ABS(XCoord - X_DIM)<1.0E-6_CMISSDP .AND. ABS(YCoord)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,2, &
          & CMISS_BOUNDARY_CONDITION_FIXED,YCoord,Err)
        WRITE(*,*) "FIXING NODES",NODE,"AT (X_DIM,0) IN Y DIRECTION"
        Y_FIXED=.TRUE.
      ENDIF

      !Fix nodal displacements at (0,Y_DIM) in x
      IF(ABS(XCoord)<1.0E-6_CMISSDP .AND. ABS(YCoord - Y_DIM)<1.0E-6_CMISSDP) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsSolid,DependentFieldSolid,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NODE,1, &
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
  !------------------------------------








  !------------------------------------
  ! ASSIGN BOUNDARY CONDITIONS - FLUID
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsDarcy,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)

  !Get surfaces 
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_FRONT_SURFACE, &
    & Face7Nodes,FaceXi(1),Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_BACK_SURFACE, &
    & Face8Nodes,FaceXi(2),Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_RIGHT_SURFACE, &
    & Face9Nodes,FaceXi(3),Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE, &
    & Face10Nodes,FaceXi(4),Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_TOP_SURFACE, &
    & Face11Nodes,FaceXi(5),Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,DarcyVelMeshComponentNumber,CMISS_GENERATED_MESH_REGULAR_BOTTOM_SURFACE, &
    & Face12Nodes,FaceXi(6),Err)



  ! At the top impose Darcy velocity in z direction
  DO NN=1,SIZE(Face11Nodes,1)
    NODE=Face11Nodes(NN)
!     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
!     IF(NodeDomain==ComputationalNodeNumber) THEN
      VALUE = -2.0_CMISSDP
      COMPONENT_NUMBER = 3
      write(*,*)'Marker 0'
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,1,NODE, &
        & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      WRITE(*,*) "SPECIFIED INFLOW AT NODE",NODE,"IN Z DIRECTION"

!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,1,XCoord,Err)
!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,2,YCoord,Err)
!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,3,ZCoord,Err)
!       WRITE(*,*) "XCoord, YCoord, ZCoord = ",XCoord, YCoord, ZCoord
!     ENDIF
  ENDDO

  !All other faces are impermeable
  DO NN=1,SIZE(Face7Nodes,1)
    NODE=Face7Nodes(NN)
!     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
!     IF(NodeDomain==ComputationalNodeNumber) THEN
      VALUE = 0.0_CMISSDP
      COMPONENT_NUMBER = 1
      write(*,*)'Marker 1'
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,1,NODE, &
        & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN X DIRECTION"

!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,1,XCoord,Err)
!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,2,YCoord,Err)
!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,3,ZCoord,Err)
!       WRITE(*,*) "XCoord, YCoord, ZCoord = ",XCoord, YCoord, ZCoord
!     ENDIF
  ENDDO

  DO NN=1,SIZE(Face8Nodes,1)
    NODE=Face8Nodes(NN)
!     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
!     IF(NodeDomain==ComputationalNodeNumber) THEN
      VALUE = 0.0_CMISSDP
      COMPONENT_NUMBER = 1
      write(*,*)'Marker 2'
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,1,NODE, &
        & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN X DIRECTION"

!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,1,XCoord,Err)
!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,2,YCoord,Err)
!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,3,ZCoord,Err)
!       WRITE(*,*) "XCoord, YCoord, ZCoord = ",XCoord, YCoord, ZCoord
!     ENDIF
  ENDDO

  DO NN=1,SIZE(Face9Nodes,1)
    NODE=Face9Nodes(NN)
!     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
!     IF(NodeDomain==ComputationalNodeNumber) THEN
      VALUE = 0.0_CMISSDP
      COMPONENT_NUMBER = 2
      write(*,*)'Marker 3'
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,1,NODE, &
        & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Y DIRECTION"

!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,1,XCoord,Err)
!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,2,YCoord,Err)
!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,3,ZCoord,Err)
!       WRITE(*,*) "XCoord, YCoord, ZCoord = ",XCoord, YCoord, ZCoord
!     ENDIF
  ENDDO

  DO NN=1,SIZE(Face10Nodes,1)
    NODE=Face10Nodes(NN)
!     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
!     IF(NodeDomain==ComputationalNodeNumber) THEN
      VALUE = 0.0_CMISSDP
      COMPONENT_NUMBER = 2
      write(*,*)'Marker 4'
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,1,NODE, &
        & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Y DIRECTION"

!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,1,XCoord,Err)
!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,2,YCoord,Err)
!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,3,ZCoord,Err)
!       WRITE(*,*) "XCoord, YCoord, ZCoord = ",XCoord, YCoord, ZCoord
!     ENDIF
  ENDDO

  DO NN=1,SIZE(Face12Nodes,1)
    NODE=Face12Nodes(NN)
!     CALL CMISSDecomposition_NodeDomainGet(Decomposition,NODE,1,NodeDomain,Err)
!     IF(NodeDomain==ComputationalNodeNumber) THEN
      VALUE = 0.0_CMISSDP
      COMPONENT_NUMBER = 3
      write(*,*)'Marker 5'
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDarcy,DependentFieldSolid,CMISS_FIELD_V_VARIABLE_TYPE,1,1,NODE, &
        & COMPONENT_NUMBER,CMISS_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      WRITE(*,*) "SPECIFIED IMPERMEABLE WALL AT NODE",NODE,"IN Z DIRECTION"

!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,1,XCoord,Err)
!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,2,YCoord,Err)
!       CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,NODE,3,ZCoord,Err)
!       WRITE(*,*) "XCoord, YCoord, ZCoord = ",XCoord, YCoord, ZCoord
!     ENDIF
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

  EXPORT_FIELD_IO=.FALSE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL CMISSFields_Initialise(Fields,Err)
    CALL CMISSFields_Create(Region,Fields,Err)
    CALL CMISSFields_NodesExport(Fields,"FiniteElasticityDarcy","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields,"FiniteElasticityDarcy","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF


  !Finialise CMISS
!   CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM FINITEELASTICITYDARCYEXAMPLE
