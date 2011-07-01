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
! !  This example considers a volume coupled diffusion-diffusion problem
! ! 

!> Main program

PROGRAM COUPLEDDIFFUSIONADVECTIONDIFFUSIONEXAMPLE

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

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberAdvectionDiffusion=6
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberDiffusion=7
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberAdvectionDiffusion=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDiffusion=9
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumberAdvectionDiffusion=10
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumberDiffusion=11
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberAdvectionDiffusion=12
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberDiffusion=13
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberAdvectionDiffusion=15
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=14

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
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
  TYPE(CMISSFieldType) :: DependentFieldAdvectionDiffusion
  TYPE(CMISSFieldType) :: DependentFieldDiffusion
  TYPE(CMISSFieldType) :: MaterialsFieldAdvectionDiffusion
  TYPE(CMISSFieldType) :: MaterialsFieldDiffusion
  TYPE(CMISSFieldType) :: SourceFieldAdvectionDiffusion
  TYPE(CMISSFieldType) :: SourceFieldDiffusion
  TYPE(CMISSFieldType) :: InDependentFieldAdvectionDiffusion
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
  LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_TYPE=CMISSSolverProgressOutput
  LINEAR_SOLVER_DIFFUSION_OUTPUT_TYPE=CMISSSolverSolverOutput
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_ADVECTION_DIFFUSION_OUTPUT=CMISSEquationsNoOutput
  EQUATIONS_DIFFUSION_OUTPUT=CMISSEquationsNoOutput
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

  !CMISSAllDiagType/CMISSInDiagType/CMISSFromDiagType
!   CALL CMISSDiagnosticsSetOn(CMISSInDiagType,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISSAllTimingType/CMISSInTimingType/CMISSFromTimingType
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

  !
  !================================================================================================================================
  !

  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !REGION
  !For a volume-coupled problem, both concentrations are based in the same region

  !Start the creation of a new region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system as defined above
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !BASES

  !Start the creation of new bases: Geometry
  MESH_NUMBER_OF_COMPONENTS=1
  CALL CMISSBasisTypeInitialise(BasisGeometry,Err)
  CALL CMISSBasisCreateStart(BASIS_NUMBER_GEOMETRY,BasisGeometry,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasisTypeSet(BasisGeometry,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL CMISSBasisNumberOfXiSet(BasisGeometry,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  IF(NUMBER_OF_DIMENSIONS==2) THEN
    CALL CMISSBasisInterpolationXiSet(BasisGeometry,(/BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY/),Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisGeometry,(/BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY/),Err)
  ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
    CALL CMISSBasisInterpolationXiSet(BasisGeometry,(/BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY, & 
      & BASIS_XI_INTERPOLATION_GEOMETRY/),Err)                         
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisGeometry,(/BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY, &
      & BASIS_XI_GAUSS_GEOMETRY/),Err)
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(BasisGeometry,Err)
  !
  !Start the creation of another basis: Concentration_One
  IF(BASIS_XI_INTERPOLATION_CONC_ONE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisConcOne=BasisGeometry
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new velocity basis
    CALL CMISSBasisTypeInitialise(BasisConcOne,Err)
    !Start the creation of a basis
    CALL CMISSBasisCreateStart(BASIS_NUMBER_CONC_ONE,BasisConcOne,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasisTypeSet(BasisConcOne,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL CMISSBasisNumberOfXiSet(BasisConcOne,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL CMISSBasisInterpolationXiSet(BasisConcOne,(/BASIS_XI_INTERPOLATION_CONC_ONE,BASIS_XI_INTERPOLATION_CONC_ONE/),Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisConcOne,(/BASIS_XI_GAUSS_CONC_ONE,BASIS_XI_GAUSS_CONC_ONE/),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL CMISSBasisInterpolationXiSet(BasisConcOne,(/BASIS_XI_INTERPOLATION_CONC_ONE,BASIS_XI_INTERPOLATION_CONC_ONE, & 
        & BASIS_XI_INTERPOLATION_CONC_ONE/),Err)                         
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisConcOne,(/BASIS_XI_GAUSS_CONC_ONE,BASIS_XI_GAUSS_CONC_ONE, & 
        & BASIS_XI_GAUSS_CONC_ONE/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(BasisConcOne,Err)
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
    CALL CMISSBasisTypeInitialise(BasisConcTwo,Err)
    !Start the creation of a basis
    CALL CMISSBasisCreateStart(BASIS_NUMBER_CONC_TWO,BasisConcTwo,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasisTypeSet(BasisConcTwo,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL CMISSBasisNumberOfXiSet(BasisConcTwo,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL CMISSBasisInterpolationXiSet(BasisConcTwo,(/BASIS_XI_INTERPOLATION_CONC_TWO,BASIS_XI_INTERPOLATION_CONC_TWO/),Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisConcTwo,(/BASIS_XI_GAUSS_CONC_TWO,BASIS_XI_GAUSS_CONC_TWO/),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL CMISSBasisInterpolationXiSet(BasisConcTwo,(/BASIS_XI_INTERPOLATION_CONC_TWO,BASIS_XI_INTERPOLATION_CONC_TWO, & 
        & BASIS_XI_INTERPOLATION_CONC_TWO/),Err)                         
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisConcTwo,(/BASIS_XI_GAUSS_CONC_TWO,BASIS_XI_GAUSS_CONC_TWO, & 
        & BASIS_XI_GAUSS_CONC_TWO/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(BasisConcTwo,Err)
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
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSNodesCreateStart(Region,TOTAL_NUMBER_OF_ALL_NODES,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL CMISSMeshNumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL CMISSMeshNumberOfComponentsSet(Mesh,MESH_NUMBER_OF_ALL_COMPONENTS,Err)
  !
  CALL CMISSMeshElementsTypeInitialise(MeshElementsGeometry,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsConcOne,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsConcTwo,Err)
  MESH_COMPONENT_NUMBER_GEOMETRY=1
  MESH_COMPONENT_NUMBER_CONC_ONE=1
  MESH_COMPONENT_NUMBER_CONC_TWO=1
  !Specify spatial mesh component
  CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_GEOMETRY,BasisGeometry,MeshElementsGeometry,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL CMISSMeshElementsNodesSet(MeshElementsGeometry,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_GEOMETRY),Err)
  ENDDO
  CALL CMISSMeshElementsCreateFinish(MeshElementsGeometry,Err)
  !Specify concentration one mesh component
  IF(BASIS_XI_INTERPOLATION_CONC_ONE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    MeshElementsConcOne=MeshElementsGeometry
  ELSE
    MESH_COMPONENT_NUMBER_CONC_ONE=MESH_COMPONENT_NUMBER_GEOMETRY+1
    CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_ONE,BasisConcOne,MeshElementsConcOne,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL CMISSMeshElementsNodesSet(MeshElementsConcOne,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_CONC_ONE),Err)
    ENDDO
    CALL CMISSMeshElementsCreateFinish(MeshElementsConcOne,Err)
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
    CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_TWO,BasisConcTwo,MeshElementsConcTwo,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL CMISSMeshElementsNodesSet(MeshElementsConcTwo,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_CONC_TWO),Err)
    ENDDO
    CALL CMISSMeshElementsCreateFinish(MeshElementsConcTwo,Err)
  ENDIF

  !Finish the creation of the mesh
  CALL CMISSMeshCreateFinish(Mesh,Err)

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  !Create a decomposition:
  !All mesh components share the same decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,DomainUserNumber,Err)
  ! ??? Above, this should be: 'NumberOfDomains' rather than 'DomainUserNumber' ???
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the field type
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL CMISSFieldScalingTypeSet(GeometricField,CMISSFieldNoScaling,Err)
  !Set the mesh component to be used by the field components.

  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
  ENDDO

  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
        & CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSFieldParameterSetUpdateStart(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for diffusion_one
  CALL CMISSEquationsSetTypeInitialise(EquationsSetAdvectionDiffusion,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberAdvectionDiffusion,Region,GeometricField,&
    & EquationsSetAdvectionDiffusion,Err)
  !Set the equations set to be a linear source diffusion problem
  CALL CMISSEquationsSetSpecificationSet(EquationsSetAdvectionDiffusion,CMISSEquationsSetClassicalFieldClass, &
    & CMISSEquationsSetAdvectionDiffusionEquationType,CMISSEquationsSetLinearSourceAdvectionDiffusionSubtype,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSetAdvectionDiffusion,Err)

  !Create the equations set for diffusion_two
  CALL CMISSEquationsSetTypeInitialise(EquationsSetDiffusion,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberDiffusion,Region,GeometricField,EquationsSetDiffusion,Err)
  !Set the equations set to be a constant source diffusion problem
  CALL CMISSEquationsSetSpecificationSet(EquationsSetDiffusion,CMISSEquationsSetClassicalFieldClass, &
    & CMISSEquationsSetDiffusionEquationType,CMISSEquationsSetConstantSourceDiffusionSubtype,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSetDiffusion,Err)

! 

  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for ALE Darcy
  CALL CMISSFieldTypeInitialise(DependentFieldAdvectionDiffusion,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetAdvectionDiffusion,DependentFieldUserNumberAdvectionDiffusion, & 
    & DependentFieldAdvectionDiffusion,Err)
  !Set the mesh component to be used by the field components.
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDiffusionOne,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_CONC_ONE,Err)
!     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDiffusionOne,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_CONC_ONE,Err)
!   ENDDO
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetAdvectionDiffusion,Err)

  !Initialise dependent field (concentration one components)
  CALL CMISSFieldComponentValuesInitialise(DependentFieldAdvectionDiffusion,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,INITIAL_FIELD_ADVECTION_DIFFUSION,Err)


  CALL CMISSFieldTypeInitialise(DependentFieldDiffusion,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetDiffusion,DependentFieldUserNumberDiffusion, & 
    & DependentFieldDiffusion,Err)
  !Set the mesh component to be used by the field components.
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
!     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
!   ENDDO
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetDiffusion,Err)

  !Initialise dependent field (concentration one components)
  CALL CMISSFieldComponentValuesInitialise(DependentFieldDiffusioN,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,INITIAL_FIELD_DIFFUSION,Err)

  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  CALL CMISSFieldTypeInitialise(MaterialsFieldAdvectionDiffusion,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetAdvectionDiffusion, &
    & MaterialsFieldUserNumberAdvectionDiffusion,MaterialsFieldAdvectionDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetAdvectionDiffusion,Err)

  CALL CMISSFieldTypeInitialise(MaterialsFieldDiffusion,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetDiffusion,&
    & MaterialsFieldUserNumberDiffusion,MaterialsFieldDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetDiffusion,Err)

!   CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDiffusionOne,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
!     & MaterialsFieldUserNumberDiffusionOne,POROSITY_PARAM_MAT_PROPERTIES,Err)
!   CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDiffusionTwo,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
!     & MaterialsFieldUserNumberMatPropertiesPermOverVis,PERM_OVER_VIS_PARAM_MAT_PROPERTIES,Err)
  !
  !================================================================================================================================
  !

  !INDEPENDENT FIELDS
! 
!   !Create the equations set independent field variables for the solid
!   CALL CMISSFieldTypeInitialise(IndependentFieldSolid,Err)
!   CALL CMISSEquationsSetIndependentCreateStart(EquationsSetSolid,IndependentFieldUserNumberSolid, & 
!     & IndependentFieldSolid,Err)
!   !Set the mesh component to be used by the field components.
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!     CALL CMISSFieldComponentMeshComponentSet(IndependentFieldSolid,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
!   ENDDO
!   !Finish the equations set independent field variables
!   CALL CMISSEquationsSetIndependentCreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  !Create the equations set independent field variables
  CALL CMISSFieldTypeInitialise(IndependentFieldAdvectionDiffusion,Err)
  CALL CMISSEquationsSetIndependentCreateStart(EquationsSetAdvectionDiffusion,IndependentFieldUserNumberAdvectionDiffusion,&
    & IndependentFieldAdvectionDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetIndependentCreateFinish(EquationsSetAdvectionDiffusion,Err)

  !SOURCE FIELDS

   !create the equations set source field variables for both equations sets 
  !Create the equations set source field variables
  CALL CMISSFieldTypeInitialise(SourceFieldAdvectionDiffusion,Err)
  CALL CMISSEquationsSetSourceCreateStart(EquationsSetAdvectionDiffusion,SourceFieldUserNumberAdvectionDiffusion,&
    & SourceFieldAdvectionDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetSourceCreateFinish(EquationsSetAdvectionDiffusion,Err)
  !Create the equations set source field variables
  CALL CMISSFieldTypeInitialise(SourceFieldDiffusion,Err)
  CALL CMISSEquationsSetSourceCreateStart(EquationsSetDiffusion,SourceFieldUserNumberDiffusion,SourceFieldDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetSourceCreateFinish(EquationsSetDiffusion,Err)
  !
  !================================================================================================================================
  !

  !EQUATIONS

  WRITE(*,'(A)') "Creating equations set equations."

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(EquationsAdvectionDiffusion,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetAdvectionDiffusion,EquationsAdvectionDiffusion,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(EquationsAdvectionDiffusion,CMISSEquationsSparseMatrices,Err)
!   !Set the equations lumping type
!   CALL CMISSEquationsLumpingTypeSet(EquationsDarcy,CMISSEquationsUnlumpedMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(EquationsAdvectionDiffusion,EQUATIONS_ADVECTION_DIFFUSION_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetAdvectionDiffusion,Err)

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(EquationsDiffusion,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetDiffusion,EquationsDiffusion,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(EquationsDiffusion,CMISSEquationsSparseMatrices,Err)
!   !Set the equations lumping type
!   CALL CMISSEquationsLumpingTypeSet(EquationsDarcy,CMISSEquationsUnlumpedMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(EquationsDiffusion,EQUATIONS_DIFFUSION_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetDiffusion,Err)

  !
  !
  !================================================================================================================================
  !

  WRITE(*,'(A)') "start creation of a problem"
  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a coupled diffusion-diffusion problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemMultiPhysicsClass,CMISSProblemDiffusionAdvectionDiffusionType, &
    & CMISSProblemCoupledSourceDiffusionAdvecDiffusionSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,LINEAR_SOLVER_ADVECTION_DIFFUSION_START_TIME,&
    & LINEAR_SOLVER_ADVECTION_DIFFUSION_STOP_TIME,LINEAR_SOLVER_ADVECTION_DIFFUSION_TIME_INCREMENT,Err)
  !Set the output timing
  CALL CMISSControlLoopTimeOutputSet(ControlLoop,LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)


  !
  !================================================================================================================================
  !


  WRITE(*,'(A)') "Start creation of problem solvers."
  !SOLVERS

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(SolverAdvectionDiffusion,Err)
  CALL CMISSSolverTypeInitialise(SolverDiffusion,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverAdvectionDiffusion,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverDiffusion,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  !Get the deformation-dependent material properties solver
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverAdvectionDiffusionUserNumber,SolverAdvectionDiffusion,Err)
  WRITE(*,'(A)') "Solver one got."
  !Set the output type
  CALL CMISSSolverOutputTypeSet(SolverAdvectionDiffusion,LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_TYPE,Err)
  !Set the solver settings
!  IF(LINEAR_SOLVER_DIFFUSION_ONE_DIRECT_FLAG) THEN
!    CALL CMISSSolverLinearTypeSet(LinearSolverDiffusionOne,CMISSSolverLinearDirectSolveType,Err)
!    CALL CMISSSolverLibraryTypeSet(LinearSolverDiffusionOne,CMISSSolverMUMPSLibrary,Err)
!  ELSE
!    CALL CMISSSolverLinearTypeSet(LinearSolverDiffusionOne,CMISSSolverLinearIterativeSolveType,Err)
   CALL CMISSSolverDynamicLinearSolverGet(SolverAdvectionDiffusion,LinearSolverAdvectionDiffusion,Err)
    CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolverAdvectionDiffusion,MAXIMUM_ITERATIONS,Err)
!    CALL CMISSSolverLinearIterativeDivergenceToleranceSet(LinearSolverDiffusionOne,DIVERGENCE_TOLERANCE,Err)
!    CALL CMISSSolverLinearIterativeRelativeToleranceSet(LinearSolverDiffusionOne,RELATIVE_TOLERANCE,Err)
!    CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(LinearSolverDiffusionOne,ABSOLUTE_TOLERANCE,Err)
!    CALL CMISSSolverLinearIterativeGMRESRestartSet(LinearSolverDiffusionOne,RESTART_VALUE,Err)
!  ENDIF
  !Get the Darcy solver
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverDiffusionUserNumber,SolverDiffusion,Err)
  WRITE(*,'(A)') "Solver two got."
  !Set the output type
  CALL CMISSSolverOutputTypeSet(SolverDiffusion,LINEAR_SOLVER_DIFFUSION_OUTPUT_TYPE,Err)
  !Set the solver settings
!  IF(LINEAR_SOLVER_DIFFUSION_ONE_DIRECT_FLAG) THEN
!    CALL CMISSSolverLinearTypeSet(LinearSolverDiffusionOne,CMISSSolverLinearDirectSolveType,Err)
!    CALL CMISSSolverLibraryTypeSet(LinearSolverDiffusionOne,CMISSSolverMUMPSLibrary,Err)
!  ELSE
!    CALL CMISSSolverLinearTypeSet(LinearSolverDiffusionOne,CMISSSolverLinearIterativeSolveType,Err)
   CALL CMISSSolverDynamicLinearSolverGet(SolverDiffusion,LinearSolverDiffusion,Err)
    CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolverDiffusion,MAXIMUM_ITERATIONS,Err)
!    CALL CMISSSolverLinearIterativeDivergenceToleranceSet(LinearSolverDiffusionOne,DIVERGENCE_TOLERANCE,Err)
!    CALL CMISSSolverLinearIterativeRelativeToleranceSet(LinearSolverDiffusionOne,RELATIVE_TOLERANCE,Err)
!    CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(LinearSolverDiffusionOne,ABSOLUTE_TOLERANCE,Err)
!    CALL CMISSSolverLinearIterativeGMRESRestartSet(LinearSolverDiffusionOne,RESTART_VALUE,Err)
!  ENDIF


  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !


  WRITE(*,'(A)') "Start creation of the problem solver equations."
  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL CMISSSolverTypeInitialise(SolverAdvectionDiffusion,Err)
  CALL CMISSSolverTypeInitialise(SolverDiffusion,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsAdvectionDiffusion,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsDiffusion,Err)


  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !
  !Get the diffusion_one solver equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverAdvectionDiffusionUserNumber,SolverAdvectionDiffusion,Err)
  CALL CMISSSolverSolverEquationsGet(SolverAdvectionDiffusion,SolverEquationsAdvectionDiffusion,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsAdvectionDiffusion,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsAdvectionDiffusion,EquationsSetAdvectionDiffusion,EquationsSetIndex,Err)
  WRITE(*,'(A)') "Solver one equations got."
  !Get the diffusion_two equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverDiffusionUserNumber,SolverDiffusion,Err)
  CALL CMISSSolverSolverEquationsGet(SolverDiffusion,SolverEquationsDiffusion,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsDiffusion,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsDiffusion,EquationsSetDiffusion,EquationsSetIndex,Err)
WRITE(*,'(A)') "Solver two equations got."
 !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS
  !Start the creation of the equations set boundary conditions for diffusion_one
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsAdvectionDiffusion,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsAdvectionDiffusion,BoundaryConditionsAdvectionDiffusion,Err)
  IF(INLET_WALL_NODES_ADVECTION_DIFFUSION_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_ADVECTION_DIFFUSION
      NODE_NUMBER=INLET_WALL_NODES_ADVECTION_DIFFUSION(NODE_COUNTER)
      CONDITION=CMISSBoundaryConditionFixed
!       DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=0.1_CMISSDP
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsAdvectionDiffusion,DependentFieldAdvectionDiffusion, &
          & CMISSFieldUVariableType,CMISSNoGlobalDerivative,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
!       ENDDO
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions for diffusion_one
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsAdvectionDiffusion,Err)
  !Start the creation of the equations set boundary conditions for diffusion_two
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsDiffusion,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsDiffusion,BoundaryConditionsDiffusion,Err)
  IF(INLET_WALL_NODES_DIFFUSION_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION
      NODE_NUMBER=INLET_WALL_NODES_DIFFUSION(NODE_COUNTER)
      CONDITION=CMISSBoundaryConditionFixed
!       DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=0.2_CMISSDP
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDiffusion,DependentFieldDiffusion,CMISSFieldUVariableType, &
          & CMISSNoGlobalDerivative,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_TWO,CONDITION,VALUE,Err)
!       ENDDO
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions for diffusion_two
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsDiffusion,Err)

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

  EXPORT_FIELD_IO=.TRUE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"CoupledSourceDiffusionAdvectionDiffusion","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"CoupledSourceDiffusionAdvectionDiffusion","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF


  !Finialise CMISS
!   CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM COUPLEDDIFFUSIONADVECTIONDIFFUSIONEXAMPLE
