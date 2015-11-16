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

  INTEGER(CMFEIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMFEIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMFEIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMFEIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMFEIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMFEIntg), PARAMETER :: DependentFieldUserNumberAdvectionDiffusion=6
  INTEGER(CMFEIntg), PARAMETER :: DependentFieldUserNumberDiffusion=7
  INTEGER(CMFEIntg), PARAMETER :: MaterialsFieldUserNumberAdvectionDiffusion=8
  INTEGER(CMFEIntg), PARAMETER :: MaterialsFieldUserNumberDiffusion=9
  INTEGER(CMFEIntg), PARAMETER :: SourceFieldUserNumberAdvectionDiffusion=10
  INTEGER(CMFEIntg), PARAMETER :: SourceFieldUserNumberDiffusion=11
  INTEGER(CMFEIntg), PARAMETER :: EquationsSetUserNumberAdvectionDiffusion=12
  INTEGER(CMFEIntg), PARAMETER :: EquationsSetUserNumberDiffusion=13
  INTEGER(CMFEIntg), PARAMETER :: IndependentFieldUserNumberAdvectionDiffusion=15
  INTEGER(CMFEIntg), PARAMETER :: ProblemUserNumber=14

  INTEGER(CMFEIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMFEIntg), PARAMETER :: SolverAdvectionDiffusionUserNumber=1
  INTEGER(CMFEIntg), PARAMETER :: SolverDiffusionUserNumber=2
  INTEGER(CMFEIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPorosity=1
  INTEGER(CMFEIntg), PARAMETER :: MaterialsFieldUserNumberDarcyPermOverVis=2
  INTEGER(CMFEIntg), PARAMETER :: MaterialsFieldUserNumberMatPropertiesPorosity=1     !??? 3 ???
  INTEGER(CMFEIntg), PARAMETER :: MaterialsFieldUserNumberMatPropertiesPermOverVis=2     !??? 4 ???

  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMFEIntg) :: NUMBER_OF_DIMENSIONS
  
  INTEGER(CMFEIntg) :: BASIS_TYPE
  INTEGER(CMFEIntg) :: BASIS_NUMBER_GEOMETRY
  INTEGER(CMFEIntg) :: BASIS_NUMBER_CONC_ONE
  INTEGER(CMFEIntg) :: BASIS_NUMBER_CONC_TWO
  INTEGER(CMFEIntg) :: BASIS_XI_GAUSS_GEOMETRY
  INTEGER(CMFEIntg) :: BASIS_XI_GAUSS_CONC_ONE
  INTEGER(CMFEIntg) :: BASIS_XI_GAUSS_CONC_TWO
  INTEGER(CMFEIntg) :: BASIS_XI_INTERPOLATION_GEOMETRY
  INTEGER(CMFEIntg) :: BASIS_XI_INTERPOLATION_CONC_ONE
  INTEGER(CMFEIntg) :: BASIS_XI_INTERPOLATION_CONC_TWO
  INTEGER(CMFEIntg) :: MESH_NUMBER_OF_COMPONENTS,MESH_NUMBER_OF_ALL_COMPONENTS
  INTEGER(CMFEIntg) :: MESH_COMPONENT_NUMBER_GEOMETRY
  INTEGER(CMFEIntg) :: MESH_COMPONENT_NUMBER_CONC_ONE
  INTEGER(CMFEIntg) :: MESH_COMPONENT_NUMBER_CONC_TWO
  INTEGER(CMFEIntg) :: NUMBER_OF_NODES_GEOMETRY
  INTEGER(CMFEIntg) :: NUMBER_OF_NODES_CONC_ONE
  INTEGER(CMFEIntg) :: NUMBER_OF_NODES_CONC_TWO
  INTEGER(CMFEIntg) :: NUMBER_OF_ELEMENT_NODES_GEOMETRY
  INTEGER(CMFEIntg) :: NUMBER_OF_ELEMENT_NODES_CONC_ONE
  INTEGER(CMFEIntg) :: NUMBER_OF_ELEMENT_NODES_CONC_TWO
  INTEGER(CMFEIntg) :: TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ALL_NODES
  INTEGER(CMFEIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMFEIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMFEIntg) :: RESTART_VALUE
  INTEGER(CMFEIntg) :: NUMBER_OF_FIXED_WALL_NODES_ADVECTION_DIFFUSION
  INTEGER(CMFEIntg) :: NUMBER_OF_INLET_WALL_NODES_ADVECTION_DIFFUSION
  INTEGER(CMFEIntg) :: NUMBER_OF_FIXED_WALL_NODES_DIFFUSION
  INTEGER(CMFEIntg) :: NUMBER_OF_INLET_WALL_NODES_DIFFUSION
  INTEGER(CMFEIntg) :: EQUATIONS_ADVECTION_DIFFUSION_OUTPUT
  INTEGER(CMFEIntg) :: EQUATIONS_DIFFUSION_OUTPUT
  INTEGER(CMFEIntg) :: COMPONENT_NUMBER
  INTEGER(CMFEIntg) :: NODE_NUMBER
  INTEGER(CMFEIntg) :: ELEMENT_NUMBER
  INTEGER(CMFEIntg) :: NODE_COUNTER
  INTEGER(CMFEIntg) :: CONDITION

  INTEGER(CMFEIntg) :: LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_FREQUENCY
  INTEGER(CMFEIntg) :: LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_TYPE
  INTEGER(CMFEIntg) :: LINEAR_SOLVER_DIFFUSION_OUTPUT_TYPE


  REAL(CMFEDP) :: DOMAIN_X1, DOMAIN_X2, DOMAIN_Y1, DOMAIN_Y2, DOMAIN_Z1, DOMAIN_Z2
  REAL(CMFEDP) :: GEOMETRY_TOLERANCE

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_ADVECTION_DIFFUSION
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_ADVECTION_DIFFUSION
  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_DIFFUSION
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_DIFFUSION

  REAL(CMFEDP) :: INITIAL_FIELD_ADVECTION_DIFFUSION
  REAL(CMFEDP) :: INITIAL_FIELD_DIFFUSION
  REAL(CMFEDP) :: BOUNDARY_CONDITIONS_ADVECTION_DIFFUSION
  REAL(CMFEDP) :: BOUNDARY_CONDITIONS_DIFFUSION
  REAL(CMFEDP) :: DIVERGENCE_TOLERANCE
  REAL(CMFEDP) :: RELATIVE_TOLERANCE
  REAL(CMFEDP) :: ABSOLUTE_TOLERANCE
  REAL(CMFEDP) :: LINESEARCH_ALPHA
  REAL(CMFEDP) :: VALUE

  REAL(CMFEDP) :: LINEAR_SOLVER_ADVECTION_DIFFUSION_START_TIME
  REAL(CMFEDP) :: LINEAR_SOLVER_ADVECTION_DIFFUSION_STOP_TIME
  REAL(CMFEDP) :: LINEAR_SOLVER_ADVECTION_DIFFUSION_TIME_INCREMENT

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_ADVECTION_DIFFUSION_DIRECT_FLAG
  LOGICAL :: LINEAR_SOLVER_DIFFUSION_DIRECT_FLAG
  LOGICAL :: INLET_WALL_NODES_ADVECTION_DIFFUSION_FLAG
  LOGICAL :: INLET_WALL_NODES_DIFFUSION_FLAG
  !CMISS variables

  !Regions
  TYPE(cmfe_RegionType) :: Region
  TYPE(cmfe_RegionType) :: WorldRegion
  !Coordinate systems
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem
  TYPE(cmfe_CoordinateSystemType) :: WorldCoordinateSystem
  !Basis
  TYPE(cmfe_BasisType) :: BasisGeometry
  TYPE(cmfe_BasisType) :: BasisConcOne
  TYPE(cmfe_BasisType) :: BasisConcTwo
  !Nodes
  TYPE(cmfe_NodesType) :: Nodes
  !Elements
  TYPE(cmfe_MeshElementsType) :: MeshElementsGeometry
  TYPE(cmfe_MeshElementsType) :: MeshElementsConcOne
  TYPE(cmfe_MeshElementsType) :: MeshElementsConcTwo
  !Meshes
  TYPE(cmfe_MeshType) :: Mesh
  !Decompositions
  TYPE(cmfe_DecompositionType) :: Decomposition
  !Fields
  TYPE(cmfe_FieldsType) :: Fields
  !Field types
  TYPE(cmfe_FieldType) :: GeometricField
  TYPE(cmfe_FieldType) :: DependentFieldAdvectionDiffusion
  TYPE(cmfe_FieldType) :: DependentFieldDiffusion
  TYPE(cmfe_FieldType) :: MaterialsFieldAdvectionDiffusion
  TYPE(cmfe_FieldType) :: MaterialsFieldDiffusion
  TYPE(cmfe_FieldType) :: SourceFieldAdvectionDiffusion
  TYPE(cmfe_FieldType) :: SourceFieldDiffusion
  TYPE(cmfe_FieldType) :: InDependentFieldAdvectionDiffusion
  !Boundary conditions
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsAdvectionDiffusion
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsDiffusion
  !Equations sets
  TYPE(cmfe_EquationsSetType) :: EquationsSetAdvectionDiffusion
  TYPE(cmfe_EquationsSetType) :: EquationsSetDiffusion
  !Equations
  TYPE(cmfe_EquationsType) :: EquationsAdvectionDiffusion
  TYPE(cmfe_EquationsType) :: EquationsDiffusion
  !Problems
  TYPE(cmfe_ProblemType) :: Problem
  !Control loops
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  !Solvers
  TYPE(cmfe_SolverType) :: SolverAdvectionDiffusion
  TYPE(cmfe_SolverType) :: SolverDiffusion
  TYPE(cmfe_SolverType) :: LinearSolverAdvectionDiffusion
  TYPE(cmfe_SolverType) :: LinearSolverDiffusion
  !Solver equations
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsAdvectionDiffusion
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsDiffusion

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMFEIntg) :: EquationsSetIndex
  INTEGER(CMFEIntg) :: Err


  INTEGER(CMFEIntg) :: DIAG_LEVEL_LIST(5)
!   CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(8) !,TIMING_ROUTINE_LIST(1)
  CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(1) !,TIMING_ROUTINE_LIST(1)

  
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !

  !Program variables and types (finite elasticity part)

  !Test program parameters

! 
  INTEGER(CMFEIntg) :: TotalNumberOfSolidNodes
!   INTEGER(CMFEIntg) :: NumberOfSolidMeshComponents

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

  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

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
!   DOMAIN_X1 = -5.0_CMFEDP
!   DOMAIN_X2 =  5.0_CMFEDP
!   DOMAIN_Y1 = -5.0_CMFEDP
!   DOMAIN_Y2 =  5.0_CMFEDP
!   DOMAIN_Z1 = -5.0_CMFEDP
!   DOMAIN_Z2 =  5.0_CMFEDP
  !Set domain dimensions
  DOMAIN_X1 =  0.0_CMFEDP
  DOMAIN_X2 =  1.0_CMFEDP
  DOMAIN_Y1 =  0.0_CMFEDP
  DOMAIN_Y2 =  1.0_CMFEDP
  DOMAIN_Z1 =  0.0_CMFEDP
  DOMAIN_Z2 =  1.0_CMFEDP
  !Set geometric tolerance
  GEOMETRY_TOLERANCE = 1.0E-12_CMFEDP
  !Set initial values
  INITIAL_FIELD_ADVECTION_DIFFUSION=0.5_CMFEDP
  INITIAL_FIELD_DIFFUSION=1.0_CMFEDP
  !Set initial boundary conditions
  INLET_WALL_NODES_ADVECTION_DIFFUSION_FLAG=.TRUE.
  INLET_WALL_NODES_DIFFUSION_FLAG=.TRUE.
  IF(INLET_WALL_NODES_ADVECTION_DIFFUSION_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_ADVECTION_DIFFUSION=36
    ALLOCATE(INLET_WALL_NODES_ADVECTION_DIFFUSION(NUMBER_OF_INLET_WALL_NODES_ADVECTION_DIFFUSION))
    INLET_WALL_NODES_ADVECTION_DIFFUSION=(/1,5,73,109,145,181,3,7,75,111,147,183,&
     & 25,27,85,121,157,193,37,39,91,127,163,199,49,51,97,133,169,205,61,63,103,139,175,211/)
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_ADVECTION_DIFFUSION=1.0_CMFEDP
  ENDIF  !Set material parameters
  IF(INLET_WALL_NODES_DIFFUSION_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_DIFFUSION=36
    ALLOCATE(INLET_WALL_NODES_DIFFUSION(NUMBER_OF_INLET_WALL_NODES_DIFFUSION))
    INLET_WALL_NODES_DIFFUSION=(/191,155,119,83,23,21,192,156,120,84,24,22,&
     & 198,162,126,90,36,35,204,168,132,96,48,47,210,174,138,102,60,59,216,180,144,108,72,71/)
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_DIFFUSION=1.0_CMFEDP
  ENDIF  !Set material parameters

  !Set material parameters
  !Set number of Gauss points (Mind that also material field may be interpolated)
  BASIS_XI_GAUSS_GEOMETRY=3 !4
  BASIS_XI_GAUSS_CONC_ONE=3 !4
  BASIS_XI_GAUSS_CONC_TWO=3 !4
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_TYPE=CMFE_SOLVER_PROGRESS_OUTPUT
  LINEAR_SOLVER_DIFFUSION_OUTPUT_TYPE=CMFE_SOLVER_SOLVER_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_ADVECTION_DIFFUSION_OUTPUT=CMFE_EQUATIONS_NO_OUTPUT
  EQUATIONS_DIFFUSION_OUTPUT=CMFE_EQUATIONS_NO_OUTPUT
  !Set time parameter
  LINEAR_SOLVER_ADVECTION_DIFFUSION_START_TIME=0.0_CMFEDP
  LINEAR_SOLVER_ADVECTION_DIFFUSION_STOP_TIME=1.0001_CMFEDP 
  LINEAR_SOLVER_ADVECTION_DIFFUSION_TIME_INCREMENT=0.125_CMFEDP
  !Set result output parameter
  LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_FREQUENCY=1
  !Set solver parameters
  LINEAR_SOLVER_ADVECTION_DIFFUSION_DIRECT_FLAG=.FALSE.
  LINEAR_SOLVER_DIFFUSION_DIRECT_FLAG=.FALSE.

  RELATIVE_TOLERANCE=1.0E-10_CMFEDP !default: 1.0E-05_CMFEDP
  ABSOLUTE_TOLERANCE=1.0E-10_CMFEDP !default: 1.0E-10_CMFEDP
  DIVERGENCE_TOLERANCE=1.0E5_CMFEDP !default: 1.0E5
  MAXIMUM_ITERATIONS=10000_CMFEIntg !default: 100000
  RESTART_VALUE=30_CMFEIntg !default: 30
  LINESEARCH_ALPHA=1.0_CMFEDP


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

  !CMFE_ALL_DIAG_TYPE/CMFE_IN_DIAG_TYPE/CMFE_FROM_DIAG_TYPE
!   CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMFE_ALL_TIMING_TYPE/CMFE_IN_TIMING_TYPE/CMFE_FROM_TIMING_TYPE
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

  !
  !================================================================================================================================
  !

  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,NUMBER_OF_DIMENSIONS,Err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !REGION
  !For a volume-coupled problem, both concentrations are based in the same region

  !Start the creation of a new region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system as defined above
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !BASES

  !Start the creation of new bases: Geometry
  MESH_NUMBER_OF_COMPONENTS=1
  CALL cmfe_Basis_Initialise(BasisGeometry,Err)
  CALL cmfe_Basis_CreateStart(BASIS_NUMBER_GEOMETRY,BasisGeometry,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL cmfe_Basis_TypeSet(BasisGeometry,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL cmfe_Basis_NumberOfXiSet(BasisGeometry,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  IF(NUMBER_OF_DIMENSIONS==2) THEN
    CALL cmfe_Basis_InterpolationXiSet(BasisGeometry,(/BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY/),Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisGeometry,(/BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY/),Err)
  ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
    CALL cmfe_Basis_InterpolationXiSet(BasisGeometry,(/BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY, & 
      & BASIS_XI_INTERPOLATION_GEOMETRY/),Err)                         
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisGeometry,(/BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY, &
      & BASIS_XI_GAUSS_GEOMETRY/),Err)
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(BasisGeometry,Err)
  !
  !Start the creation of another basis: Concentration_One
  IF(BASIS_XI_INTERPOLATION_CONC_ONE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisConcOne=BasisGeometry
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new velocity basis
    CALL cmfe_Basis_Initialise(BasisConcOne,Err)
    !Start the creation of a basis
    CALL cmfe_Basis_CreateStart(BASIS_NUMBER_CONC_ONE,BasisConcOne,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL cmfe_Basis_TypeSet(BasisConcOne,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL cmfe_Basis_NumberOfXiSet(BasisConcOne,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisConcOne,(/BASIS_XI_INTERPOLATION_CONC_ONE,BASIS_XI_INTERPOLATION_CONC_ONE/),Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisConcOne,(/BASIS_XI_GAUSS_CONC_ONE,BASIS_XI_GAUSS_CONC_ONE/),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisConcOne,(/BASIS_XI_INTERPOLATION_CONC_ONE,BASIS_XI_INTERPOLATION_CONC_ONE, & 
        & BASIS_XI_INTERPOLATION_CONC_ONE/),Err)                         
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisConcOne,(/BASIS_XI_GAUSS_CONC_ONE,BASIS_XI_GAUSS_CONC_ONE, & 
        & BASIS_XI_GAUSS_CONC_ONE/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL cmfe_Basis_CreateFinish(BasisConcOne,Err)
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
    CALL cmfe_Basis_Initialise(BasisConcTwo,Err)
    !Start the creation of a basis
    CALL cmfe_Basis_CreateStart(BASIS_NUMBER_CONC_TWO,BasisConcTwo,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL cmfe_Basis_TypeSet(BasisConcTwo,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL cmfe_Basis_NumberOfXiSet(BasisConcTwo,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisConcTwo,(/BASIS_XI_INTERPOLATION_CONC_TWO,BASIS_XI_INTERPOLATION_CONC_TWO/),Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisConcTwo,(/BASIS_XI_GAUSS_CONC_TWO,BASIS_XI_GAUSS_CONC_TWO/),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisConcTwo,(/BASIS_XI_INTERPOLATION_CONC_TWO,BASIS_XI_INTERPOLATION_CONC_TWO, & 
        & BASIS_XI_INTERPOLATION_CONC_TWO/),Err)                         
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisConcTwo,(/BASIS_XI_GAUSS_CONC_TWO,BASIS_XI_GAUSS_CONC_TWO, & 
        & BASIS_XI_GAUSS_CONC_TWO/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL cmfe_Basis_CreateFinish(BasisConcTwo,Err)
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
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,TOTAL_NUMBER_OF_ALL_NODES,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,MESH_NUMBER_OF_ALL_COMPONENTS,Err)
  !
  CALL cmfe_MeshElements_Initialise(MeshElementsGeometry,Err)
  CALL cmfe_MeshElements_Initialise(MeshElementsConcOne,Err)
  CALL cmfe_MeshElements_Initialise(MeshElementsConcTwo,Err)
  MESH_COMPONENT_NUMBER_GEOMETRY=1
  MESH_COMPONENT_NUMBER_CONC_ONE=1
  MESH_COMPONENT_NUMBER_CONC_TWO=1
  !Specify spatial mesh component
  CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_GEOMETRY,BasisGeometry,MeshElementsGeometry,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL cmfe_MeshElements_NodesSet(MeshElementsGeometry,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_GEOMETRY),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(MeshElementsGeometry,Err)
  !Specify concentration one mesh component
  IF(BASIS_XI_INTERPOLATION_CONC_ONE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    MeshElementsConcOne=MeshElementsGeometry
  ELSE
    MESH_COMPONENT_NUMBER_CONC_ONE=MESH_COMPONENT_NUMBER_GEOMETRY+1
    CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_ONE,BasisConcOne,MeshElementsConcOne,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL cmfe_MeshElements_NodesSet(MeshElementsConcOne,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_CONC_ONE),Err)
    ENDDO
    CALL cmfe_MeshElements_CreateFinish(MeshElementsConcOne,Err)
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
    CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_TWO,BasisConcTwo,MeshElementsConcTwo,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL cmfe_MeshElements_NodesSet(MeshElementsConcTwo,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_CONC_TWO),Err)
    ENDDO
    CALL cmfe_MeshElements_CreateFinish(MeshElementsConcTwo,Err)
  ENDIF

  !Finish the creation of the mesh
  CALL cmfe_Mesh_CreateFinish(Mesh,Err)

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  !Create a decomposition:
  !All mesh components share the same decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,DomainUserNumber,Err)
  ! ??? Above, this should be: 'NumberOfDomains' rather than 'DomainUserNumber' ???
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the field type
  CALL cmfe_Field_TypeSet(GeometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL cmfe_Field_ScalingTypeSet(GeometricField,CMFE_FIELD_NO_SCALING,Err)
  !Set the mesh component to be used by the field components.

  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
  ENDDO

  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL cmfe_Field_ParameterSetUpdateStart(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for diffusion_one
  CALL cmfe_EquationsSet_Initialise(EquationsSetAdvectionDiffusion,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberAdvectionDiffusion,Region,GeometricField, &
    & [EquationsSetAdvectionDiffusion,Err])
  !Set the equations set to be a linear source diffusion problem
  CALL cmfe_EquationsSet_SpecificationSet(EquationsSetAdvectionDiffusion,CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetAdvectionDiffusion,Err)

  !Create the equations set for diffusion_two
  CALL cmfe_EquationsSet_Initialise(EquationsSetDiffusion,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberDiffusion,Region,GeometricField,[EquationsSetDiffusion,Err])
  !Set the equations set to be a constant source diffusion problem
  CALL cmfe_EquationsSet_SpecificationSet(EquationsSetDiffusion,CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSetDiffusion,Err)

! 

  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for ALE Darcy
  CALL cmfe_Field_Initialise(DependentFieldAdvectionDiffusion,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetAdvectionDiffusion,DependentFieldUserNumberAdvectionDiffusion, & 
    & DependentFieldAdvectionDiffusion,Err)
  !Set the mesh component to be used by the field components.
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!     CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDiffusionOne,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_CONC_ONE,Err)
!     CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDiffusionOne,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_CONC_ONE,Err)
!   ENDDO
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetAdvectionDiffusion,Err)

  !Initialise dependent field (concentration one components)
  CALL cmfe_Field_ComponentValuesInitialise(DependentFieldAdvectionDiffusion,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,INITIAL_FIELD_ADVECTION_DIFFUSION,Err)


  CALL cmfe_Field_Initialise(DependentFieldDiffusion,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetDiffusion,DependentFieldUserNumberDiffusion, & 
    & DependentFieldDiffusion,Err)
  !Set the mesh component to be used by the field components.
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!     CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
!     CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
!   ENDDO
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetDiffusion,Err)

  !Initialise dependent field (concentration one components)
  CALL cmfe_Field_ComponentValuesInitialise(DependentFieldDiffusioN,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & 1,INITIAL_FIELD_DIFFUSION,Err)

  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  CALL cmfe_Field_Initialise(MaterialsFieldAdvectionDiffusion,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetAdvectionDiffusion, &
    & MaterialsFieldUserNumberAdvectionDiffusion,MaterialsFieldAdvectionDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetAdvectionDiffusion,Err)

  CALL cmfe_Field_Initialise(MaterialsFieldDiffusion,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetDiffusion,&
    & MaterialsFieldUserNumberDiffusion,MaterialsFieldDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetDiffusion,Err)

!   CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDiffusionOne,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
!     & MaterialsFieldUserNumberDiffusionOne,POROSITY_PARAM_MAT_PROPERTIES,Err)
!   CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDiffusionTwo,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
!     & MaterialsFieldUserNumberMatPropertiesPermOverVis,PERM_OVER_VIS_PARAM_MAT_PROPERTIES,Err)
  !
  !================================================================================================================================
  !

  !INDEPENDENT FIELDS
! 
!   !Create the equations set independent field variables for the solid
!   CALL cmfe_Field_Initialise(IndependentFieldSolid,Err)
!   CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetSolid,IndependentFieldUserNumberSolid, & 
!     & IndependentFieldSolid,Err)
!   !Set the mesh component to be used by the field components.
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!     CALL cmfe_Field_ComponentMeshComponentSet(IndependentFieldSolid,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
!       & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
!   ENDDO
!   !Finish the equations set independent field variables
!   CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  !Create the equations set independent field variables
  CALL cmfe_Field_Initialise(IndependentFieldAdvectionDiffusion,Err)
  CALL cmfe_EquationsSet_IndependentCreateStart(EquationsSetAdvectionDiffusion,IndependentFieldUserNumberAdvectionDiffusion,&
    & IndependentFieldAdvectionDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_IndependentCreateFinish(EquationsSetAdvectionDiffusion,Err)

  !SOURCE FIELDS

   !create the equations set source field variables for both equations sets 
  !Create the equations set source field variables
  CALL cmfe_Field_Initialise(SourceFieldAdvectionDiffusion,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(EquationsSetAdvectionDiffusion,SourceFieldUserNumberAdvectionDiffusion,&
    & SourceFieldAdvectionDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(EquationsSetAdvectionDiffusion,Err)
  !Create the equations set source field variables
  CALL cmfe_Field_Initialise(SourceFieldDiffusion,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(EquationsSetDiffusion,SourceFieldUserNumberDiffusion,SourceFieldDiffusion,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(EquationsSetDiffusion,Err)
  !
  !================================================================================================================================
  !

  !EQUATIONS

  WRITE(*,'(A)') "Creating equations set equations."

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(EquationsAdvectionDiffusion,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetAdvectionDiffusion,EquationsAdvectionDiffusion,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(EquationsAdvectionDiffusion,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
!   !Set the equations lumping type
!   CALL cmfe_Equations_LumpingTypeSet(EquationsDarcy,CMFE_EQUATIONS_UNLUMPED_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(EquationsAdvectionDiffusion,EQUATIONS_ADVECTION_DIFFUSION_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetAdvectionDiffusion,Err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(EquationsDiffusion,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetDiffusion,EquationsDiffusion,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(EquationsDiffusion,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
!   !Set the equations lumping type
!   CALL cmfe_Equations_LumpingTypeSet(EquationsDarcy,CMFE_EQUATIONS_UNLUMPED_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(EquationsDiffusion,EQUATIONS_DIFFUSION_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetDiffusion,Err)

  !
  !
  !================================================================================================================================
  !

  WRITE(*,'(A)') "start creation of a problem"
  !PROBLEMS

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
    & CMFE_PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE,CMFE_PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,LINEAR_SOLVER_ADVECTION_DIFFUSION_START_TIME,&
    & LINEAR_SOLVER_ADVECTION_DIFFUSION_STOP_TIME,LINEAR_SOLVER_ADVECTION_DIFFUSION_TIME_INCREMENT,Err)
  !Set the output timing
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)


  !
  !================================================================================================================================
  !


  WRITE(*,'(A)') "Start creation of problem solvers."
  !SOLVERS

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(SolverAdvectionDiffusion,Err)
  CALL cmfe_Solver_Initialise(SolverDiffusion,Err)
  CALL cmfe_Solver_Initialise(LinearSolverAdvectionDiffusion,Err)
  CALL cmfe_Solver_Initialise(LinearSolverDiffusion,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  !Get the deformation-dependent material properties solver
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverAdvectionDiffusionUserNumber,SolverAdvectionDiffusion,Err)
  WRITE(*,'(A)') "Solver one got."
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(SolverAdvectionDiffusion,LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_TYPE,Err)
  !Set the solver settings
!  IF(LINEAR_SOLVER_DIFFUSION_ONE_DIRECT_FLAG) THEN
!    CALL cmfe_Solver_LinearTypeSet(LinearSolverDiffusionOne,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
!    CALL cmfe_Solver_LibraryTypeSet(LinearSolverDiffusionOne,CMFE_SOLVER_MUMPS_LIBRARY,Err)
!  ELSE
!    CALL cmfe_Solver_LinearTypeSet(LinearSolverDiffusionOne,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
   CALL cmfe_Solver_DynamicLinearSolverGet(SolverAdvectionDiffusion,LinearSolverAdvectionDiffusion,Err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolverAdvectionDiffusion,MAXIMUM_ITERATIONS,Err)
!    CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(LinearSolverDiffusionOne,DIVERGENCE_TOLERANCE,Err)
!    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolverDiffusionOne,RELATIVE_TOLERANCE,Err)
!    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolverDiffusionOne,ABSOLUTE_TOLERANCE,Err)
!    CALL cmfe_Solver_LinearIterativeGMRESRestartSet(LinearSolverDiffusionOne,RESTART_VALUE,Err)
!  ENDIF
  !Get the Darcy solver
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverDiffusionUserNumber,SolverDiffusion,Err)
  WRITE(*,'(A)') "Solver two got."
  !Set the output type
  CALL cmfe_Solver_OutputTypeSet(SolverDiffusion,LINEAR_SOLVER_DIFFUSION_OUTPUT_TYPE,Err)
  !Set the solver settings
!  IF(LINEAR_SOLVER_DIFFUSION_ONE_DIRECT_FLAG) THEN
!    CALL cmfe_Solver_LinearTypeSet(LinearSolverDiffusionOne,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
!    CALL cmfe_Solver_LibraryTypeSet(LinearSolverDiffusionOne,CMFE_SOLVER_MUMPS_LIBRARY,Err)
!  ELSE
!    CALL cmfe_Solver_LinearTypeSet(LinearSolverDiffusionOne,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
   CALL cmfe_Solver_DynamicLinearSolverGet(SolverDiffusion,LinearSolverDiffusion,Err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolverDiffusion,MAXIMUM_ITERATIONS,Err)
!    CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(LinearSolverDiffusionOne,DIVERGENCE_TOLERANCE,Err)
!    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolverDiffusionOne,RELATIVE_TOLERANCE,Err)
!    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolverDiffusionOne,ABSOLUTE_TOLERANCE,Err)
!    CALL cmfe_Solver_LinearIterativeGMRESRestartSet(LinearSolverDiffusionOne,RESTART_VALUE,Err)
!  ENDIF


  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !


  WRITE(*,'(A)') "Start creation of the problem solver equations."
  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(SolverAdvectionDiffusion,Err)
  CALL cmfe_Solver_Initialise(SolverDiffusion,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsAdvectionDiffusion,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsDiffusion,Err)


  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !
  !Get the diffusion_one solver equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverAdvectionDiffusionUserNumber,SolverAdvectionDiffusion,Err)
  CALL cmfe_Solver_SolverEquationsGet(SolverAdvectionDiffusion,SolverEquationsAdvectionDiffusion,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsAdvectionDiffusion,CMFE_SOLVER_SPARSE_MATRICES,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsAdvectionDiffusion,EquationsSetAdvectionDiffusion,EquationsSetIndex,Err)
  WRITE(*,'(A)') "Solver one equations got."
  !Get the diffusion_two equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverDiffusionUserNumber,SolverDiffusion,Err)
  CALL cmfe_Solver_SolverEquationsGet(SolverDiffusion,SolverEquationsDiffusion,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsDiffusion,CMFE_SOLVER_SPARSE_MATRICES,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsDiffusion,EquationsSetDiffusion,EquationsSetIndex,Err)
WRITE(*,'(A)') "Solver two equations got."
 !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS
  !Start the creation of the equations set boundary conditions for diffusion_one
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsAdvectionDiffusion,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsAdvectionDiffusion,BoundaryConditionsAdvectionDiffusion, &
    & Err)
  IF(INLET_WALL_NODES_ADVECTION_DIFFUSION_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_ADVECTION_DIFFUSION
      NODE_NUMBER=INLET_WALL_NODES_ADVECTION_DIFFUSION(NODE_COUNTER)
      CONDITION=CMFE_BOUNDARY_CONDITION_FIXED
!       DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=0.1_CMFEDP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsAdvectionDiffusion,DependentFieldAdvectionDiffusion, &
          & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
!       ENDDO
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions for diffusion_one
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsAdvectionDiffusion,Err)
  !Start the creation of the equations set boundary conditions for diffusion_two
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsDiffusion,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsDiffusion,BoundaryConditionsDiffusion,Err)
  IF(INLET_WALL_NODES_DIFFUSION_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION
      NODE_NUMBER=INLET_WALL_NODES_DIFFUSION(NODE_COUNTER)
      CONDITION=CMFE_BOUNDARY_CONDITION_FIXED
!       DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=0.2_CMFEDP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDiffusion,DependentFieldDiffusion,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_TWO,CONDITION,VALUE,Err)
!       ENDDO
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions for diffusion_two
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsDiffusion,Err)

  !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn of PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  WRITE(*,'(A)') "Solving problem..."
  CALL cmfe_Problem_Solve(Problem,Err)
  WRITE(*,'(A)') "Problem solved!"


  !
  !================================================================================================================================
  !

  !OUTPUT

  EXPORT_FIELD_IO=.TRUE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"CoupledSourceDiffusionAdvectionDiffusion","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"CoupledSourceDiffusionAdvectionDiffusion","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF


  !Finialise CMISS
!   CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM COUPLEDDIFFUSIONADVECTIONDIFFUSIONEXAMPLE
