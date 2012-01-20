!> \file
!> \authors Andrew Cookson
!> \brief This is an example program to solve a coupled diffusion & advection-diffusion equation using openCMISS calls.
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

!> \example MultiPhysics/CoupledDiffusionDiffusion/src/CoupledDiffusionDiffusionExample.f90
!! Example program to solve coupled CoupledSourceDiffusionDiffusion equations using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/CoupledDiffusionDiffusion/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/CoupledDiffusionDiffusion/build-intel'>Linux GNU Build</a>
!!
!<

! ! 
! !  This example considers a volume coupled diffusion & advection-diffusion problem - using a *shared* dependent field.
! !  The equations are solved in a partitioned manner currently.

!> Main program

PROGRAM COUPLEDDIFFUSIONADVECTIONDIFFUSIONEXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OPENCMISS
  USE FLUID_MECHANICS_IO_ROUTINES
  USE FIELDML_OUTPUT_ROUTINES
  USE FIELDML_UTIL_ROUTINES
  USE FIELDML_API
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
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=6   !combined diffusion & advection-diffusion
!  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberDiffusion=7
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberAdvectionDiffusion=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDiffusion=9
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumberAdvectionDiffusion=10
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumberDiffusion=11
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberAdvectionDiffusion=12
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberDiffusion=13
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberAdvectionDiffusion=15
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberAdvectionDiffusion=22
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberDiffusion=23

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfUserDomains=1
  INTEGER(CMISSIntg), PARAMETER :: SolverAdvectionDiffusionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverDiffusionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPorosity=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPermOverVis=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMatPropertiesPorosity=1     !??? 3 ???
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMatPropertiesPermOverVis=2     !??? 4 ???

  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS
  
  INTEGER(CMISSIntg) :: BASIS_TYPE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_CONC_ONE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_CONC_TWO
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_CONC_ONE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_CONC_TWO
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_CONC_ONE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_CONC_TWO
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS,MESH_NUMBER_OF_ALL_COMPONENTS
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_GEOMETRY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_CONC_ONE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_CONC_TWO
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_GEOMETRY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_CONC_ONE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_CONC_TWO
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_GEOMETRY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_CONC_ONE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_CONC_TWO
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ALL_NODES
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES_ADVECTION_DIFFUSION
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_WALL_NODES_ADVECTION_DIFFUSION
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES_DIFFUSION
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_WALL_NODES_DIFFUSION
  INTEGER(CMISSIntg) :: EQUATIONS_ADVECTION_DIFFUSION_OUTPUT
  INTEGER(CMISSIntg) :: EQUATIONS_DIFFUSION_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_COUNTER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_DIFFUSION_OUTPUT_TYPE


  REAL(CMISSDP) :: DOMAIN_X1, DOMAIN_X2, DOMAIN_Y1, DOMAIN_Y2, DOMAIN_Z1, DOMAIN_Z2
  REAL(CMISSDP) :: GEOMETRY_TOLERANCE

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_ADVECTION_DIFFUSION
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_ADVECTION_DIFFUSION
  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_DIFFUSION
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_DIFFUSION

  REAL(CMISSDP) :: INITIAL_FIELD_ADVECTION_DIFFUSION
  REAL(CMISSDP) :: INITIAL_FIELD_DIFFUSION
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_ADVECTION_DIFFUSION
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_DIFFUSION
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE

  REAL(CMISSDP) :: LINEAR_SOLVER_ADVECTION_DIFFUSION_START_TIME
  REAL(CMISSDP) :: LINEAR_SOLVER_ADVECTION_DIFFUSION_STOP_TIME
  REAL(CMISSDP) :: LINEAR_SOLVER_ADVECTION_DIFFUSION_TIME_INCREMENT

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_ADVECTION_DIFFUSION_DIRECT_FLAG
  LOGICAL :: LINEAR_SOLVER_DIFFUSION_DIRECT_FLAG
  LOGICAL :: INLET_WALL_NODES_ADVECTION_DIFFUSION_FLAG
  LOGICAL :: INLET_WALL_NODES_DIFFUSION_FLAG
  !CMISS variables

  !Regions
  TYPE(CMISSRegionType) :: Region
  TYPE(CMISSRegionType) :: WorldRegion
  !Coordinate systems
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  !Basis
  TYPE(CMISSBasisType) :: BasisGeometry
  TYPE(CMISSBasisType) :: BasisConcOne
  TYPE(CMISSBasisType) :: BasisConcTwo
  !Nodes
  TYPE(CMISSNodesType) :: Nodes
  !Elements
  TYPE(CMISSMeshElementsType) :: MeshElementsGeometry
  TYPE(CMISSMeshElementsType) :: MeshElementsConcOne
  TYPE(CMISSMeshElementsType) :: MeshElementsConcTwo
  !Meshes
  TYPE(CMISSMeshType) :: Mesh
  !Decompositions
  TYPE(CMISSDecompositionType) :: Decomposition
  !Fields
  TYPE(CMISSFieldsType) :: Fields
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: DependentField
!  TYPE(CMISSFieldType) :: DependentFieldDiffusion
  TYPE(CMISSFieldType) :: MaterialsFieldAdvectionDiffusion
  TYPE(CMISSFieldType) :: MaterialsFieldDiffusion
  TYPE(CMISSFieldType) :: SourceFieldAdvectionDiffusion
  TYPE(CMISSFieldType) :: SourceFieldDiffusion
  TYPE(CMISSFieldType) :: InDependentFieldAdvectionDiffusion
  TYPE(CMISSFieldType) :: EquationsSetFieldAdvectionDiffusion
  TYPE(CMISSFieldType) :: EquationsSetFieldDiffusion
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsAdvectionDiffusion
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsDiffusion
  !Equations sets
  TYPE(CMISSEquationsSetType) :: EquationsSetAdvectionDiffusion
  TYPE(CMISSEquationsSetType) :: EquationsSetDiffusion
  !Equations
  TYPE(CMISSEquationsType) :: EquationsAdvectionDiffusion
  TYPE(CMISSEquationsType) :: EquationsDiffusion
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: SolverAdvectionDiffusion
  TYPE(CMISSSolverType) :: SolverDiffusion
  TYPE(CMISSSolverType) :: LinearSolverAdvectionDiffusion
  TYPE(CMISSSolverType) :: LinearSolverDiffusion
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsAdvectionDiffusion
  TYPE(CMISSSolverEquationsType) :: SolverEquationsDiffusion

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

  !FieldML variables
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: basename = "advection_diffusion_advection"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputDirectory = ""
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputFilename = "AdvectionDiffusionAdvection.xml"

  INTEGER(C_INT) :: equationsSetDomain

  TYPE(FieldmlInfoType) :: fieldmlInfo
  
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !

  !Program variables and types (finite elasticity part)

  !Test program parameters

! 
  INTEGER(CMISSIntg) :: TotalNumberOfSolidNodes
!   INTEGER(CMISSIntg) :: NumberOfSolidMeshComponents

! 

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

  !INITIALISE OPENCMISS

  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  !
  !================================================================================================================================
  !

  !PROBLEM CONTROL PANEL

  !Import cmHeart mesh information
  CALL FLUID_MECHANICS_IO_READ_CMHEART(CM,Err)  
  BASIS_NUMBER_GEOMETRY=CM%ID_M
  BASIS_NUMBER_CONC_ONE=CM%ID_V !USE THE V COMPONENT OF CMHEART INPUT FOR CONCENTRATION ONE
  BASIS_NUMBER_CONC_TWO=CM%ID_P !USE THE p COMPONENT OF CMHEART INPUT FOR CONCENTRATION TWO
  NUMBER_OF_DIMENSIONS=CM%D
  BASIS_TYPE=CM%IT_T
  BASIS_XI_INTERPOLATION_GEOMETRY=CM%IT_M
  BASIS_XI_INTERPOLATION_CONC_ONE=CM%IT_V
  BASIS_XI_INTERPOLATION_CONC_TWO=CM%IT_P
  NUMBER_OF_NODES_GEOMETRY=CM%N_M
  NUMBER_OF_NODES_CONC_ONE=CM%N_V
  NUMBER_OF_NODES_CONC_TWO=CM%N_P
  TOTAL_NUMBER_OF_NODES=CM%N_T
  TOTAL_NUMBER_OF_ELEMENTS=CM%E_T
  NUMBER_OF_ELEMENT_NODES_GEOMETRY=CM%EN_M
  NUMBER_OF_ELEMENT_NODES_CONC_ONE=CM%EN_V
  NUMBER_OF_ELEMENT_NODES_CONC_TWO=CM%EN_P
!   !Set domain dimensions
!   DOMAIN_X1 = -5.0_CMISSDP
!   DOMAIN_X2 =  5.0_CMISSDP
!   DOMAIN_Y1 = -5.0_CMISSDP
!   DOMAIN_Y2 =  5.0_CMISSDP
!   DOMAIN_Z1 = -5.0_CMISSDP
!   DOMAIN_Z2 =  5.0_CMISSDP
  !Set domain dimensions
  DOMAIN_X1 =  0.0_CMISSDP
  DOMAIN_X2 =  1.0_CMISSDP
  DOMAIN_Y1 =  0.0_CMISSDP
  DOMAIN_Y2 =  1.0_CMISSDP
  DOMAIN_Z1 =  0.0_CMISSDP
  DOMAIN_Z2 =  1.0_CMISSDP
  !Set geometric tolerance
  GEOMETRY_TOLERANCE = 1.0E-12_CMISSDP
  !Set initial values
  INITIAL_FIELD_ADVECTION_DIFFUSION=0.5_CMISSDP
  INITIAL_FIELD_DIFFUSION=1.0_CMISSDP
  !Set initial boundary conditions
  INLET_WALL_NODES_ADVECTION_DIFFUSION_FLAG=.TRUE.
  INLET_WALL_NODES_DIFFUSION_FLAG=.TRUE.
  IF(INLET_WALL_NODES_ADVECTION_DIFFUSION_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_ADVECTION_DIFFUSION=36
    ALLOCATE(INLET_WALL_NODES_ADVECTION_DIFFUSION(NUMBER_OF_INLET_WALL_NODES_ADVECTION_DIFFUSION))
    INLET_WALL_NODES_ADVECTION_DIFFUSION=(/1,5,73,109,145,181,3,7,75,111,147,183,&
     & 25,27,85,121,157,193,37,39,91,127,163,199,49,51,97,133,169,205,61,63,103,139,175,211/)
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_ADVECTION_DIFFUSION=1.0_CMISSDP
  ENDIF  !Set material parameters
  IF(INLET_WALL_NODES_DIFFUSION_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_DIFFUSION=36
    ALLOCATE(INLET_WALL_NODES_DIFFUSION(NUMBER_OF_INLET_WALL_NODES_DIFFUSION))
    INLET_WALL_NODES_DIFFUSION=(/191,155,119,83,23,21,192,156,120,84,24,22,&
     & 198,162,126,90,36,35,204,168,132,96,48,47,210,174,138,102,60,59,216,180,144,108,72,71/)
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_DIFFUSION=1.0_CMISSDP
  ENDIF  !Set material parameters

  !Set material parameters
  !Set number of Gauss points (Mind that also material field may be interpolated)
  BASIS_XI_GAUSS_GEOMETRY=3 !4
  BASIS_XI_GAUSS_CONC_ONE=3 !4
  BASIS_XI_GAUSS_CONC_TWO=3 !4
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_TYPE=CMISS_SOLVER_PROGRESS_OUTPUT
  LINEAR_SOLVER_DIFFUSION_OUTPUT_TYPE=CMISS_SOLVER_SOLVER_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_ADVECTION_DIFFUSION_OUTPUT=CMISS_EQUATIONS_NO_OUTPUT
  EQUATIONS_DIFFUSION_OUTPUT=CMISS_EQUATIONS_NO_OUTPUT
  !Set time parameter
  LINEAR_SOLVER_ADVECTION_DIFFUSION_START_TIME=0.0_CMISSDP
  LINEAR_SOLVER_ADVECTION_DIFFUSION_STOP_TIME=1.0001_CMISSDP 
  LINEAR_SOLVER_ADVECTION_DIFFUSION_TIME_INCREMENT=0.125_CMISSDP
  !Set result output parameter
  LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_FREQUENCY=1
  !Set solver parameters
  LINEAR_SOLVER_ADVECTION_DIFFUSION_DIRECT_FLAG=.FALSE.
  LINEAR_SOLVER_DIFFUSION_DIRECT_FLAG=.FALSE.

  RELATIVE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E5_CMISSDP !default: 1.0E5
  MAXIMUM_ITERATIONS=10000_CMISSIntg !default: 100000
  RESTART_VALUE=30_CMISSIntg !default: 30
  LINESEARCH_ALPHA=1.0_CMISSDP

  !
  !================================================================================================================================
  !

  !Set diagnostics

  DIAG_LEVEL_LIST(1)=1
  DIAG_LEVEL_LIST(2)=2
  DIAG_LEVEL_LIST(3)=3
  DIAG_LEVEL_LIST(4)=4
  DIAG_LEVEL_LIST(5)=5

!   DIAG_ROUTINE_LIST(1)="DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE"
!   DIAG_ROUTINE_LIST(2)="DARCY_EQUATION_PRE_SOLVE_STORE_REFERENCE_DATA"
!   DIAG_ROUTINE_LIST(3)="DARCY_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH"
!   DIAG_ROUTINE_LIST(4)="DARCY_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS"
!   DIAG_ROUTINE_LIST(5)="DARCY_EQUATION_PRE_SOLVE_MAT_PROPERTIES"
!   DIAG_ROUTINE_LIST(6)="FITTING_FINITE_ELEMENT_CALCULATE"
!   DIAG_ROUTINE_LIST(7)="FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE"
!   DIAG_ROUTINE_LIST(8)="FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE"
!   DIAG_ROUTINE_LIST(1)="PROBLEM_SOLVER_EQUATIONS_SOLVE"
!   DIAG_ROUTINE_LIST(1)="SOLVER_NEWTON_SOLVE"
!   DIAG_ROUTINE_LIST(2)="SOLVER_NEWTON_LINESEARCH_SOLVE"
!   DIAG_ROUTINE_LIST(1)="SOLVER_SOLUTION_UPDATE"
!  DIAG_ROUTINE_LIST(1)="FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE"

  !CMISS_ALL_DIAG_TYPE/CMISS_IN_DIAG_TYPE/CMISS_FROM_DIAG_TYPE
!   CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISS_ALL_TIMING_TYPE/CMISS_IN_TIMING_TYPE/CMISS_FROM_TIMING_TYPE
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

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
  !For a volume-coupled problem, both concentrations are based in the same region

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

  !Start the creation of new bases: Geometry
  MESH_NUMBER_OF_COMPONENTS=1
  CALL CMISSBasis_Initialise(BasisGeometry,Err)
  CALL CMISSBasis_CreateStart(BASIS_NUMBER_GEOMETRY,BasisGeometry,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasis_TypeSet(BasisGeometry,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL CMISSBasis_NumberOfXiSet(BasisGeometry,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  IF(NUMBER_OF_DIMENSIONS==2) THEN
    CALL CMISSBasis_InterpolationXiSet(BasisGeometry,(/BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY/),Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisGeometry,(/BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY/),Err)
  ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
    CALL CMISSBasis_InterpolationXiSet(BasisGeometry,(/BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY, & 
      & BASIS_XI_INTERPOLATION_GEOMETRY/),Err)                         
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisGeometry,(/BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY, &
      & BASIS_XI_GAUSS_GEOMETRY/),Err)
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(BasisGeometry,Err)
  !
  !Start the creation of another basis: Concentration_One
  IF(BASIS_XI_INTERPOLATION_CONC_ONE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisConcOne=BasisGeometry
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new velocity basis
    CALL CMISSBasis_Initialise(BasisConcOne,Err)
    !Start the creation of a basis
    CALL CMISSBasis_CreateStart(BASIS_NUMBER_CONC_ONE,BasisConcOne,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasis_TypeSet(BasisConcOne,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL CMISSBasis_NumberOfXiSet(BasisConcOne,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisConcOne,(/BASIS_XI_INTERPOLATION_CONC_ONE,BASIS_XI_INTERPOLATION_CONC_ONE/),Err)
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisConcOne,(/BASIS_XI_GAUSS_CONC_ONE,BASIS_XI_GAUSS_CONC_ONE/),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisConcOne,(/BASIS_XI_INTERPOLATION_CONC_ONE,BASIS_XI_INTERPOLATION_CONC_ONE, & 
        & BASIS_XI_INTERPOLATION_CONC_ONE/),Err)                         
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisConcOne,(/BASIS_XI_GAUSS_CONC_ONE,BASIS_XI_GAUSS_CONC_ONE, & 
        & BASIS_XI_GAUSS_CONC_ONE/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasis_CreateFinish(BasisConcOne,Err)
  ENDIF
  !
  !Start the creation of another basis: Concentration_Two
  IF(BASIS_XI_INTERPOLATION_CONC_TWO==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisConcTwo=BasisGeometry
  ELSE IF(BASIS_XI_INTERPOLATION_CONC_TWO==BASIS_XI_INTERPOLATION_CONC_ONE) THEN
    BasisConcTwo=BasisConcOne
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new concentration basis
    CALL CMISSBasis_Initialise(BasisConcTwo,Err)
    !Start the creation of a basis
    CALL CMISSBasis_CreateStart(BASIS_NUMBER_CONC_TWO,BasisConcTwo,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasis_TypeSet(BasisConcTwo,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL CMISSBasis_NumberOfXiSet(BasisConcTwo,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisConcTwo,(/BASIS_XI_INTERPOLATION_CONC_TWO,BASIS_XI_INTERPOLATION_CONC_TWO/),Err)
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisConcTwo,(/BASIS_XI_GAUSS_CONC_TWO,BASIS_XI_GAUSS_CONC_TWO/),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisConcTwo,(/BASIS_XI_INTERPOLATION_CONC_TWO,BASIS_XI_INTERPOLATION_CONC_TWO, & 
        & BASIS_XI_INTERPOLATION_CONC_TWO/),Err)                         
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisConcTwo,(/BASIS_XI_GAUSS_CONC_TWO,BASIS_XI_GAUSS_CONC_TWO, & 
        & BASIS_XI_GAUSS_CONC_TWO/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasis_CreateFinish(BasisConcTwo,Err)
  ENDIF

  !
  !================================================================================================================================
  !

  !MESH
  !All types of physics utilize the same "mesh", but may be represented on individual mesh components.

  TotalNumberOfSolidNodes = NUMBER_OF_NODES_GEOMETRY
  TOTAL_NUMBER_OF_ALL_NODES = TOTAL_NUMBER_OF_NODES + TotalNumberOfSolidNodes

!   NumberOfSolidMeshComponents = 1
!   MESH_NUMBER_OF_ALL_COMPONENTS = MESH_NUMBER_OF_COMPONENTS + NumberOfSolidMeshComponents
  MESH_NUMBER_OF_ALL_COMPONENTS = MESH_NUMBER_OF_COMPONENTS

  !Start the creation of mesh nodes
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSNodes_CreateStart(Region,TOTAL_NUMBER_OF_ALL_NODES,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL CMISSMesh_NumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL CMISSMesh_NumberOfComponentsSet(Mesh,MESH_NUMBER_OF_ALL_COMPONENTS,Err)
  !
  CALL CMISSMeshElements_Initialise(MeshElementsGeometry,Err)
  CALL CMISSMeshElements_Initialise(MeshElementsConcOne,Err)
  CALL CMISSMeshElements_Initialise(MeshElementsConcTwo,Err)
  MESH_COMPONENT_NUMBER_GEOMETRY=1
  MESH_COMPONENT_NUMBER_CONC_ONE=1
  MESH_COMPONENT_NUMBER_CONC_TWO=1
  !Specify spatial mesh component
  CALL CMISSMeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_GEOMETRY,BasisGeometry,MeshElementsGeometry,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL CMISSMeshElements_NodesSet(MeshElementsGeometry,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_GEOMETRY),Err)
  ENDDO
  CALL CMISSMeshElements_CreateFinish(MeshElementsGeometry,Err)
  !Specify concentration one mesh component
  IF(BASIS_XI_INTERPOLATION_CONC_ONE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    MeshElementsConcOne=MeshElementsGeometry
  ELSE
    MESH_COMPONENT_NUMBER_CONC_ONE=MESH_COMPONENT_NUMBER_GEOMETRY+1
    CALL CMISSMeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_ONE,BasisConcOne,MeshElementsConcOne,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL CMISSMeshElements_NodesSet(MeshElementsConcOne,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_CONC_ONE),Err)
    ENDDO
    CALL CMISSMeshElements_CreateFinish(MeshElementsConcOne,Err)
  ENDIF
  !Specify concentration two mesh component
  IF(BASIS_XI_INTERPOLATION_CONC_TWO==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    MeshElementsConcTwo=MeshElementsGeometry
    MESH_COMPONENT_NUMBER_CONC_TWO=MESH_COMPONENT_NUMBER_GEOMETRY
  ELSE IF(BASIS_XI_INTERPOLATION_CONC_TWO==BASIS_XI_INTERPOLATION_CONC_ONE) THEN
    MeshElementsConcTwo=MeshElementsConcOne
    MESH_COMPONENT_NUMBER_CONC_TWO=MESH_COMPONENT_NUMBER_CONC_ONE
  ELSE
    MESH_COMPONENT_NUMBER_CONC_TWO=MESH_COMPONENT_NUMBER_CONC_ONE+1
    CALL CMISSMeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_TWO,BasisConcTwo,MeshElementsConcTwo,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL CMISSMeshElements_NodesSet(MeshElementsConcTwo,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_CONC_TWO),Err)
    ENDDO
    CALL CMISSMeshElements_CreateFinish(MeshElementsConcTwo,Err)
  ENDIF

  !Finish the creation of the mesh
  CALL CMISSMesh_CreateFinish(Mesh,Err)

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  !Create a decomposition:
  !All mesh components share the same decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,DomainUserNumber,Err)
  ! ??? Above, this should be: 'NumberOfDomains' rather than 'DomainUserNumber' ???
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the field type
  CALL CMISSField_TypeSet(GeometricField,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL CMISSField_ScalingTypeSet(GeometricField,CMISS_FIELD_NO_SCALING,Err)

  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
  ENDDO
   WRITE(*,'(A)') "now finishing creating geometric field"
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, & 
        & CMISS_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSField_ParameterSetUpdateStart(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  WRITE(*,'(A)') "geometric field made"
  !
  !================================================================================================================================
  !
 !
  !================================================================================================================================
  !
  !EQUATIONS SETS

  !Create the equations set for diffusion_one
  CALL CMISSField_Initialise(EquationsSetFieldAdvectionDiffusion,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSetAdvectionDiffusion,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumberAdvectionDiffusion,Region,GeometricField,&
    & CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE,CMISS_EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE,&
    & EquationsSetFieldUserNumberAdvectionDiffusion,EquationsSetFieldAdvectionDiffusion,EquationsSetAdvectionDiffusion,Err)
  !Set the equations set to be a linear source diffusion problem
!   CALL CMISSEquationsSet_SpecificationSet(EquationsSetAdvectionDiffusion,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
!     & CMISS_EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE,CMISS_EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSetAdvectionDiffusion,Err)

  !Create the equations set for diffusion_two
  CALL CMISSField_Initialise(EquationsSetFieldDiffusion,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSetDiffusion,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumberDiffusion,Region,GeometricField,&
    & CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMISS_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,CMISS_EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE,&
    & EquationsSetFieldUserNumberDiffusion,EquationsSetFieldDiffusion,EquationsSetDiffusion,Err)
  !Set the equations set to be a constant source diffusion problem
!   CALL CMISSEquationsSet_SpecificationSet(EquationsSetDiffusion,CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
!     & CMISS_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,CMISS_EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSetDiffusion,Err)
  !-------------------------------------------------------------------------------------
  ! DEPENDENT FIELD: Shared. We're gonna create the dependent field manually for now
  ! U & DEL_U_DEL_N: Advection-diffusion
  ! V & DEL_V_DEL_N: Diffusion
  CALL CMISSField_Initialise(DependentField,Err)
    ! -- Advection-diffusion
    CALL CMISSField_CreateStart(DependentFieldUserNumber,Region,DependentField,Err)
    CALL CMISSField_TypeSet(DependentField,CMISS_FIELD_GENERAL_TYPE,Err)  
    CALL CMISSField_MeshDecompositionSet(DependentField,Decomposition,Err)
    CALL CMISSField_GeometricFieldSet(DependentField,GeometricField,Err) 
    CALL CMISSField_DependentTypeSet(DependentField,CMISS_FIELD_DEPENDENT_TYPE,Err) 
     CALL CMISSField_NumberOfVariablesSet(DependentField,4,Err)
     CALL CMISSField_VariableTypesSet(DependentField,&
      & (/CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,CMISS_FIELD_V_VARIABLE_TYPE, &
        & CMISS_FIELD_DELVDELN_VARIABLE_TYPE/),Err)
   WRITE(*,'(A)') "set number of variables"
!     CALL CMISSField_NumberOfVariablesSet(DependentField,2,Err)
!     CALL CMISSField_VariableTypesSet(DependentField,&
!      & (/CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_DELVDELN_VARIABLE_TYPE/),Err)
     CALL CMISSField_DimensionSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE, &
                   & CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
     CALL CMISSField_DimensionSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE, &
                   & CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
     CALL CMISSField_DimensionSet(DependentField,CMISS_FIELD_V_VARIABLE_TYPE, &
                   & CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
     CALL CMISSField_DimensionSet(DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE, &
                   & CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
     CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,Err)
     CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1, & 
        & MESH_COMPONENT_NUMBER_CONC_ONE,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,1, & 
        & MESH_COMPONENT_NUMBER_CONC_ONE,Err)
   WRITE(*,'(A)') "advec-diff vars set"
!     ! -- Diffusion
     CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_V_VARIABLE_TYPE,1,Err)
     CALL CMISSField_NumberOfComponentsSet(DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_V_VARIABLE_TYPE,1, & 
        & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
      CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1, & 
        & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
  WRITE(*,'(A)') "diff vars set"
  CALL CMISSField_CreateFinish(DependentField,Err)
  WRITE(*,'(A)') "field created"
 


  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for ALE Darcy
!  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetAdvectionDiffusion,DependentFieldUserNumber, & 
    & DependentField,Err)
  !Set the mesh component to be used by the field components.
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDiffusionOne,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_CONC_ONE,Err)
!     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDiffusionOne,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_CONC_ONE,Err)
!   ENDDO
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetAdvectionDiffusion,Err)

!  CALL CMISSField_Initialise(DependentFieldDiffusion,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetDiffusion,DependentFieldUserNumber, & 
    & DependentField,Err)
  !Set the mesh component to be used by the field components.
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
!     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
!   ENDDO
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetDiffusion,Err)



  !-------------------------------------------------------------------------------------
  ! INITIALISE DEPENDENT FIELDS

!   !Initialise dependent field (concentration one components)
!   CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
!     & 1,INITIAL_FIELD_ADVECTION_DIFFUSION,Err)
! 
!   !Initialise dependent field (concentration two components)
!   CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
!     & 1,INITIAL_FIELD_DIFFUSION,Err)



  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  CALL CMISSField_Initialise(MaterialsFieldAdvectionDiffusion,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetAdvectionDiffusion, &
    & MaterialsFieldUserNumberAdvectionDiffusion,MaterialsFieldAdvectionDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetAdvectionDiffusion,Err)

  CALL CMISSField_Initialise(MaterialsFieldDiffusion,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetDiffusion,&
    & MaterialsFieldUserNumberDiffusion,MaterialsFieldDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetDiffusion,Err)

!   CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDiffusionOne,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
!     & MaterialsFieldUserNumberDiffusionOne,POROSITY_PARAM_MAT_PROPERTIES,Err)
!   CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDiffusionTwo,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
!     & MaterialsFieldUserNumberMatPropertiesPermOverVis,PERM_OVER_VIS_PARAM_MAT_PROPERTIES,Err)
  !
  !================================================================================================================================
  !

  !
  !================================================================================================================================
  !

  !Create the equations set independent field variables
  CALL CMISSField_Initialise(IndependentFieldAdvectionDiffusion,Err)
  CALL CMISSEquationsSet_IndependentCreateStart(EquationsSetAdvectionDiffusion,IndependentFieldUserNumberAdvectionDiffusion,&
    & IndependentFieldAdvectionDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_IndependentCreateFinish(EquationsSetAdvectionDiffusion,Err)

  !SOURCE FIELDS

   !create the equations set source field variables for both equations sets 
  !Create the equations set source field variables
  CALL CMISSField_Initialise(SourceFieldAdvectionDiffusion,Err)
  CALL CMISSEquationsSet_SourceCreateStart(EquationsSetAdvectionDiffusion,SourceFieldUserNumberAdvectionDiffusion,&
    & SourceFieldAdvectionDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_SourceCreateFinish(EquationsSetAdvectionDiffusion,Err)
  !Create the equations set source field variables
  CALL CMISSField_Initialise(SourceFieldDiffusion,Err)
  CALL CMISSEquationsSet_SourceCreateStart(EquationsSetDiffusion,SourceFieldUserNumberDiffusion,SourceFieldDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_SourceCreateFinish(EquationsSetDiffusion,Err)
  !
  !================================================================================================================================
  !

  !EQUATIONS

  WRITE(*,'(A)') "Creating equations set equations."

  !Create the equations set equations
  CALL CMISSEquations_Initialise(EquationsAdvectionDiffusion,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetAdvectionDiffusion,EquationsAdvectionDiffusion,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(EquationsAdvectionDiffusion,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
!   !Set the equations lumping type
!   CALL CMISSEquations_LumpingTypeSet(EquationsDarcy,CMISS_EQUATIONS_UNLUMPED_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(EquationsAdvectionDiffusion,EQUATIONS_ADVECTION_DIFFUSION_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetAdvectionDiffusion,Err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(EquationsDiffusion,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetDiffusion,EquationsDiffusion,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(EquationsDiffusion,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
!   !Set the equations lumping type
!   CALL CMISSEquations_LumpingTypeSet(EquationsDarcy,CMISS_EQUATIONS_UNLUMPED_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(EquationsDiffusion,EQUATIONS_DIFFUSION_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetDiffusion,Err)

  !
  !
  !================================================================================================================================
  !

  WRITE(*,'(A)') "start creation of a problem"
  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a coupled diffusion-diffusion problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_MULTI_PHYSICS_CLASS,CMISS_PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE, &
    & CMISS_PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoop_TimesSet(ControlLoop,LINEAR_SOLVER_ADVECTION_DIFFUSION_START_TIME,&
    & LINEAR_SOLVER_ADVECTION_DIFFUSION_STOP_TIME,LINEAR_SOLVER_ADVECTION_DIFFUSION_TIME_INCREMENT,Err)
  !Set the output timing
  CALL CMISSControlLoop_TimeOutputSet(ControlLoop,LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)


  !
  !================================================================================================================================
  !


  WRITE(*,'(A)') "Start creation of problem solvers."
  !SOLVERS

  !Start the creation of the problem solvers
  CALL CMISSSolver_Initialise(SolverAdvectionDiffusion,Err)
  CALL CMISSSolver_Initialise(SolverDiffusion,Err)
  CALL CMISSSolver_Initialise(LinearSolverAdvectionDiffusion,Err)
  CALL CMISSSolver_Initialise(LinearSolverDiffusion,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  !Get the deformation-dependent material properties solver
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverAdvectionDiffusionUserNumber,SolverAdvectionDiffusion,Err)
  WRITE(*,'(A)') "Solver one got."
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(SolverAdvectionDiffusion,LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_TYPE,Err)
  !Set the solver settings
!  IF(LINEAR_SOLVER_DIFFUSION_ONE_DIRECT_FLAG) THEN
!    CALL CMISSSolver_LinearTypeSet(LinearSolverDiffusionOne,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
!    CALL CMISSSolver_LibraryTypeSet(LinearSolverDiffusionOne,CMISS_SOLVER_MUMPS_LIBRARY,Err)
!  ELSE
!    CALL CMISSSolver_LinearTypeSet(LinearSolverDiffusionOne,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
   CALL CMISSSolver_DynamicLinearSolverGet(SolverAdvectionDiffusion,LinearSolverAdvectionDiffusion,Err)
    CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolverAdvectionDiffusion,MAXIMUM_ITERATIONS,Err)
!    CALL CMISSSolver_LinearIterativeDivergenceToleranceSet(LinearSolverDiffusionOne,DIVERGENCE_TOLERANCE,Err)
!    CALL CMISSSolver_LinearIterativeRelativeToleranceSet(LinearSolverDiffusionOne,RELATIVE_TOLERANCE,Err)
!    CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(LinearSolverDiffusionOne,ABSOLUTE_TOLERANCE,Err)
!    CALL CMISSSolver_LinearIterativeGMRESRestartSet(LinearSolverDiffusionOne,RESTART_VALUE,Err)
!  ENDIF
  !Get the Darcy solver
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverDiffusionUserNumber,SolverDiffusion,Err)
  WRITE(*,'(A)') "Solver two got."
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(SolverDiffusion,LINEAR_SOLVER_DIFFUSION_OUTPUT_TYPE,Err)
  !Set the solver settings
!  IF(LINEAR_SOLVER_DIFFUSION_ONE_DIRECT_FLAG) THEN
!    CALL CMISSSolver_LinearTypeSet(LinearSolverDiffusionOne,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
!    CALL CMISSSolver_LibraryTypeSet(LinearSolverDiffusionOne,CMISS_SOLVER_MUMPS_LIBRARY,Err)
!  ELSE
!    CALL CMISSSolver_LinearTypeSet(LinearSolverDiffusionOne,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
   CALL CMISSSolver_DynamicLinearSolverGet(SolverDiffusion,LinearSolverDiffusion,Err)
    CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolverDiffusion,MAXIMUM_ITERATIONS,Err)
!    CALL CMISSSolver_LinearIterativeDivergenceToleranceSet(LinearSolverDiffusionOne,DIVERGENCE_TOLERANCE,Err)
!    CALL CMISSSolver_LinearIterativeRelativeToleranceSet(LinearSolverDiffusionOne,RELATIVE_TOLERANCE,Err)
!    CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(LinearSolverDiffusionOne,ABSOLUTE_TOLERANCE,Err)
!    CALL CMISSSolver_LinearIterativeGMRESRestartSet(LinearSolverDiffusionOne,RESTART_VALUE,Err)
!  ENDIF


  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !


  WRITE(*,'(A)') "Start creation of the problem solver equations."
  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(SolverAdvectionDiffusion,Err)
  CALL CMISSSolver_Initialise(SolverDiffusion,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsAdvectionDiffusion,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsDiffusion,Err)


  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  !
  !Get the diffusion_one solver equations
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverAdvectionDiffusionUserNumber,SolverAdvectionDiffusion,Err)
  CALL CMISSSolver_SolverEquationsGet(SolverAdvectionDiffusion,SolverEquationsAdvectionDiffusion,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsAdvectionDiffusion,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsAdvectionDiffusion,EquationsSetAdvectionDiffusion,EquationsSetIndex,Err)
  WRITE(*,'(A)') "Solver one equations got."
  !Get the diffusion_two equations
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverDiffusionUserNumber,SolverDiffusion,Err)
  CALL CMISSSolver_SolverEquationsGet(SolverDiffusion,SolverEquationsDiffusion,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsDiffusion,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsDiffusion,EquationsSetDiffusion,EquationsSetIndex,Err)
WRITE(*,'(A)') "Solver two equations got."
 !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS
  !Start the creation of the equations set boundary conditions for advection-diffusion
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsAdvectionDiffusion,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsAdvectionDiffusion,BoundaryConditionsAdvectionDiffusion, &
    & Err)
  IF(INLET_WALL_NODES_ADVECTION_DIFFUSION_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_ADVECTION_DIFFUSION
      NODE_NUMBER=INLET_WALL_NODES_ADVECTION_DIFFUSION(NODE_COUNTER)
      CONDITION=CMISS_BOUNDARY_CONDITION_FIXED
!       DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=0.1_CMISSDP
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsAdvectionDiffusion,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
         & CMISS_NO_GLOBAL_DERIV,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
!       ENDDO
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions for advection-diffusion
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsAdvectionDiffusion,Err)
  !Start the creation of the equations set boundary conditions for diffusion
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsDiffusion,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsDiffusion,BoundaryConditionsDiffusion,Err)
  IF(INLET_WALL_NODES_DIFFUSION_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION
      NODE_NUMBER=INLET_WALL_NODES_DIFFUSION(NODE_COUNTER)
      CONDITION=CMISS_BOUNDARY_CONDITION_FIXED
!       DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=0.2_CMISSDP
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsDiffusion,DependentField,CMISS_FIELD_V_VARIABLE_TYPE,1, &
          & CMISS_NO_GLOBAL_DERIV,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_TWO,CONDITION,VALUE,Err)
!       ENDDO
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions for diffusion
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsDiffusion,Err)


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

  EXPORT_FIELD_IO=.TRUE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
!    CALL CMISSFields_Initialise(Fields,Err)
!    CALL CMISSFields_Create(Region,Fields,Err)
!    CALL CMISSFields_NodesExport(Fields,"CoupledSourceDiffusionAdvectionDiffusion","FORTRAN",Err)
!    CALL CMISSFields_ElementsExport(Fields,"CoupledSourceDiffusionAdvectionDiffusion","FORTRAN",Err)
!    CALL CMISSFields_Finalise(Fields,Err)

    CALL FieldmlOutput_InitializeInfo( Region, Mesh, NUMBER_OF_DIMENSIONS, outputDirectory, basename, fieldmlInfo, err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".geometric", region, mesh, GeometricField, &
      & CMISS_FIELD_U_VARIABLE_TYPE, err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".dependent", region, mesh, DependentField, &
      & CMISS_FIELD_U_VARIABLE_TYPE, err )

    CALL FieldmlOutput_CreateEnsembleDomain( fieldmlInfo, baseName//".equation_set", 3, equationsSetDomain, err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".equation_set.advection_diffusion", mesh, &
      & EquationsSetFieldDiffusion, CMISS_FIELD_U_VARIABLE_TYPE, equationsSetDomain, err )

    CALL FieldmlOutput_Write( fieldmlInfo, outputFilename, err )

    CALL FieldmlUtil_FinalizeInfo( fieldmlInfo )


    WRITE(*,'(A)') "Field exported!"
  ENDIF


  !Finialise CMISS
!   CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM COUPLEDDIFFUSIONADVECTIONDIFFUSIONEXAMPLE
