!> \file
!> \authors Andrew Cookson
!> \brief This is an example program to solve coupled multi-compartment diffusion equations in monolithic scheme using OpenCMISS calls.
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

!> \example MultiPhysics/CoupledDiffusionDiffusion/MonolithicSchemeTest/src/MonolithicSchemeTestExample.f90
!! Example program to solve coupled MonolithicSchemeTest equations using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/CoupledDiffusionDiffusion/MonolithicSchemeTest/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/CoupledDiffusionDiffusion/MonolithicSchemeTest/build-intel'>Linux GNU Build</a>
!!
!<

! ! 
! !  This example considers a volume coupled multi-compartment diffusion model, as a means of testing the monolithic assembly of
! !  such a system of equations, for use in more complicated problem in future.
! !  This example will initially couple together three diffusion equations, with the transfer between each equation being proportional to
! !  the concentration difference between each equation.

!> Main program

PROGRAM MONOLITHICSCHEMETESTEXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OpenCMISS
  USE OpenCMISS_Iron
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

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=6
  INTEGER(CMISSIntg) :: MaterialsFieldUserNumberDiffusion
  INTEGER(CMISSIntg) :: SourceFieldUserNumberDiffusion
  INTEGER(CMISSIntg) :: EquationsSetUserNumberDiffusion
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=20
  INTEGER(CMISSIntg) :: EquationsSetFieldUserNumberDiffusion
  INTEGER(CMISSIntg) :: icompartment,Ncompartments,num_var
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfUserDomains=1
  INTEGER(CMISSIntg), PARAMETER :: SolverDiffusionUserNumber=1
  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS

  INTEGER(CMISSIntg) :: NUMBER_OF_COMPARTMENTS
  
  INTEGER(CMISSIntg) :: BASIS_TYPE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_CONC_ONE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_CONC_TWO
  INTEGER(CMISSIntg) :: BASIS_NUMBER_CONC_THREE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_CONC_ONE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_CONC_TWO
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_CONC_THREE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_CONC_ONE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_CONC_TWO
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_CONC_THREE
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS,MESH_NUMBER_OF_ALL_COMPONENTS
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_GEOMETRY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_CONC_ONE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_CONC_TWO
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_CONC_THREE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_GEOMETRY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_CONC_ONE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_CONC_TWO
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_CONC_THREE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_GEOMETRY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_CONC_ONE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_CONC_TWO
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_CONC_THREE
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ALL_NODES
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES_DIFFUSION_ONE
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_WALL_NODES_DIFFUSION_ONE
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES_DIFFUSION_TWO
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_WALL_NODES_DIFFUSION_TWO
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES_DIFFUSION_THREE
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_WALL_NODES_DIFFUSION_THREE
  INTEGER(CMISSIntg) :: EQUATIONS_DIFFUSION_OUTPUT
  INTEGER(CMISSIntg) :: EQUATIONS_DIFFUSION_TWO_OUTPUT
  INTEGER(CMISSIntg) :: EQUATIONS_DIFFUSION_THREE_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_COUNTER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: LINEAR_SOLVER_DIFFUSION_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_DIFFUSION_OUTPUT_TYPE


  REAL(CMISSRP) :: DOMAIN_X1, DOMAIN_X2, DOMAIN_Y1, DOMAIN_Y2, DOMAIN_Z1, DOMAIN_Z2
  REAL(CMISSRP) :: GEOMETRY_TOLERANCE

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_DIFFUSION_ONE
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_DIFFUSION_ONE
  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_DIFFUSION_TWO
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_DIFFUSION_TWO
  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_DIFFUSION_THREE
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_DIFFUSION_THREE

  REAL(CMISSRP) :: INITIAL_FIELD_DIFFUSION_ONE
  REAL(CMISSRP) :: INITIAL_FIELD_DIFFUSION_TWO
  REAL(CMISSRP) :: INITIAL_FIELD_DIFFUSION_THREE
  REAL(CMISSRP) :: BOUNDARY_CONDITIONS_DIFFUSION_ONE
  REAL(CMISSRP) :: BOUNDARY_CONDITIONS_DIFFUSION_TWO
  REAL(CMISSRP) :: BOUNDARY_CONDITIONS_DIFFUSION_THREE
  REAL(CMISSRP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSRP) :: RELATIVE_TOLERANCE
  REAL(CMISSRP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSRP) :: LINESEARCH_ALPHA
  REAL(CMISSRP) :: VALUE

  REAL(CMISSRP) :: LINEAR_SOLVER_DIFFUSION_START_TIME
  REAL(CMISSRP) :: LINEAR_SOLVER_DIFFUSION_STOP_TIME
  REAL(CMISSRP) :: LINEAR_SOLVER_DIFFUSION_TIME_INCREMENT

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_DIFFUSION_DIRECT_FLAG
  LOGICAL :: INLET_WALL_NODES_DIFFUSION_ONE_FLAG
  LOGICAL :: INLET_WALL_NODES_DIFFUSION_TWO_FLAG
  LOGICAL :: INLET_WALL_NODES_DIFFUSION_THREE_FLAG
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
  TYPE(cmfe_BasisType) :: BasisConcThree
  
  !TYPE(cmfe_BasisType), ALLOCATABLE, DIMENSION(:) :: BasicConc

  !Nodes
  TYPE(cmfe_NodesType) :: Nodes
  !Elements
  TYPE(cmfe_MeshElementsType) :: MeshElementsGeometry
  TYPE(cmfe_MeshElementsType) :: MeshElementsConcOne
  TYPE(cmfe_MeshElementsType) :: MeshElementsConcTwo
  TYPE(cmfe_MeshElementsType) :: MeshElementsConcThree

  TYPE(cmfe_MeshElementsType), ALLOCATABLE, DIMENSION(:) :: MeshElementsConc
  !Meshes
  TYPE(cmfe_MeshType) :: Mesh
  !Decompositions
  TYPE(cmfe_DecompositionType) :: Decomposition
  !Fields
  TYPE(cmfe_FieldsType) :: Fields
  !Field types
  TYPE(cmfe_FieldType) :: GeometricField
  TYPE(cmfe_FieldType) :: DependentField
  TYPE(cmfe_FieldType), ALLOCATABLE, DIMENSION(:) :: MaterialsFieldDiffusion
  TYPE(cmfe_FieldType), ALLOCATABLE, DIMENSION(:) :: SourceFieldDiffusion
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditionsDiffusion
  TYPE(cmfe_EquationsSetType), ALLOCATABLE, DIMENSION(:) :: EquationsSetDiffusion
  TYPE(cmfe_EquationsType), ALLOCATABLE, DIMENSION(:) :: EquationsDiffusion
  TYPE(cmfe_FieldType), ALLOCATABLE, DIMENSION(:) :: EquationsSetFieldDiffusion
  !Problems
  TYPE(cmfe_ProblemType) :: Problem
  !Control loops
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  !Solvers
  TYPE(cmfe_SolverType) :: SolverDiffusion
  TYPE(cmfe_SolverType) :: LinearSolverDiffusion
  !Solver equations
  TYPE(cmfe_SolverEquationsType) :: SolverEquationsDiffusion

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: EquationsSetIndexTwo
  INTEGER(CMISSIntg) :: EquationsSetIndexThree
  INTEGER(CMISSIntg) :: Err


  !Array containing the field variable types that will be used (for ease of incorporating inside a loop)
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:) :: VariableTypes
  REAL(CMISSRP), ALLOCATABLE, DIMENSION(:,:) :: CouplingCoeffs

  INTEGER(CMISSIntg) :: DIAG_LEVEL_LIST(5)
!   CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(8) !,TIMING_ROUTINE_LIST(1)
  CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(1) !,TIMING_ROUTINE_LIST(1)
 
  INTEGER(CMISSIntg) :: TotalNumberOfSolidNodes
!   INTEGER(CMISSIntg) :: NumberOfSolidMeshComponents
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
  !CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
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
  !Set domain dimensions
  DOMAIN_X1 =  0.0_CMISSRP
  DOMAIN_X2 =  1.0_CMISSRP
  DOMAIN_Y1 =  0.0_CMISSRP
  DOMAIN_Y2 =  1.0_CMISSRP
  DOMAIN_Z1 =  0.0_CMISSRP
  DOMAIN_Z2 =  1.0_CMISSRP
  !Set geometric tolerance
  GEOMETRY_TOLERANCE = 1.0E-12_CMISSRP
  !Set initial values
  INITIAL_FIELD_DIFFUSION_ONE=1.0_CMISSRP
  INITIAL_FIELD_DIFFUSION_TWO=1.0_CMISSRP
  INITIAL_FIELD_DIFFUSION_THREE=1.0_CMISSRP
  !Set initial boundary conditions
  INLET_WALL_NODES_DIFFUSION_ONE_FLAG=.TRUE.
  IF(INLET_WALL_NODES_DIFFUSION_ONE_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_DIFFUSION_ONE=36
    ALLOCATE(INLET_WALL_NODES_DIFFUSION_ONE(NUMBER_OF_INLET_WALL_NODES_DIFFUSION_ONE))
    INLET_WALL_NODES_DIFFUSION_ONE=[191,155,119,83,23,21,192,156,120,84,24,22,&
     & 198,162,126,90,36,35,204,168,132,96,48,47,210,174,138,102,60,59,216,180,144,108,72,71]
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_DIFFUSION_ONE=1.0_CMISSRP
  ENDIF  !Set material parameters

  INLET_WALL_NODES_DIFFUSION_TWO_FLAG=.TRUE.
  IF(INLET_WALL_NODES_DIFFUSION_TWO_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_DIFFUSION_TWO=36
    ALLOCATE(INLET_WALL_NODES_DIFFUSION_TWO(NUMBER_OF_INLET_WALL_NODES_DIFFUSION_TWO))
    INLET_WALL_NODES_DIFFUSION_TWO=[191,155,119,83,23,21,192,156,120,84,24,22,&
     & 198,162,126,90,36,35,204,168,132,96,48,47,210,174,138,102,60,59,216,180,144,108,72,71]
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_DIFFUSION_TWO=1.0_CMISSRP
  ENDIF  !Set material parameters

  INLET_WALL_NODES_DIFFUSION_THREE_FLAG=.TRUE.
  IF(INLET_WALL_NODES_DIFFUSION_THREE_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_DIFFUSION_THREE=36
    ALLOCATE(INLET_WALL_NODES_DIFFUSION_THREE(NUMBER_OF_INLET_WALL_NODES_DIFFUSION_THREE))
    INLET_WALL_NODES_DIFFUSION_THREE=[191,155,119,83,23,21,192,156,120,84,24,22,&
     & 198,162,126,90,36,35,204,168,132,96,48,47,210,174,138,102,60,59,216,180,144,108,72,71]
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_DIFFUSION_THREE=1.0_CMISSRP
  ENDIF  !Set material parameters



  !Set material parameters
  !Set number of Gauss points (Mind that also material field may be interpolated)
  BASIS_XI_GAUSS_GEOMETRY=3 !4
  BASIS_XI_GAUSS_CONC_ONE=3 !4
  BASIS_XI_GAUSS_CONC_TWO=3 !4
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_DIFFUSION_OUTPUT_TYPE=CMFE_SOLVER_SOLVER_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DIFFUSION_OUTPUT=CMFE_EQUATIONS_MATRIX_OUTPUT
  EQUATIONS_DIFFUSION_TWO_OUTPUT=CMFE_EQUATIONS_NO_OUTPUT
  EQUATIONS_DIFFUSION_THREE_OUTPUT=CMFE_EQUATIONS_NO_OUTPUT
  !Set time parameter
  LINEAR_SOLVER_DIFFUSION_START_TIME=0.0_CMISSRP
  LINEAR_SOLVER_DIFFUSION_STOP_TIME=0.25_CMISSRP 
  LINEAR_SOLVER_DIFFUSION_TIME_INCREMENT=0.125_CMISSRP
  !Set result output parameter
  LINEAR_SOLVER_DIFFUSION_OUTPUT_FREQUENCY=1
  !Set solver parameters
  LINEAR_SOLVER_DIFFUSION_DIRECT_FLAG=.FALSE.

  RELATIVE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-05_CMISSRP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSRP !default: 1.0E-10_CMISSRP
  DIVERGENCE_TOLERANCE=1.0E5_CMISSRP !default: 1.0E5
  MAXIMUM_ITERATIONS=10000_CMISSIntg !default: 100000
  RESTART_VALUE=30_CMISSIntg !default: 30
  LINESEARCH_ALPHA=1.0_CMISSRP


  icompartment =1_CMISSIntg
  Ncompartments=4_CMISSIntg
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
!   DIAG_ROUTINE_LIST(2)="PROBLEM_SOLVER_EQUATIONS_SOLVE"
!   DIAG_ROUTINE_LIST(3)="SOLVER_SOLUTION_UPDATE"

  !CMFE_ALL_DIAG_TYPE/CMFE_IN_DIAG_TYPE/CMFE_FROM_DIAG_TYPE
!   CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMFE_ALL_TIMING_TYPE/CMFE_IN_TIMING_TYPE/CMFE_FROM_TIMING_TYPE
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

  !ALLOCATE THE ARRAYS
  ALLOCATE (EquationsSetDiffusion(Ncompartments))
  ALLOCATE (EquationsSetFieldDiffusion(Ncompartments))
  ALLOCATE (MaterialsFieldDiffusion(Ncompartments))
  ALLOCATE (SourceFieldDiffusion(Ncompartments))
  ALLOCATE (EquationsDiffusion(Ncompartments))

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
    CALL cmfe_Basis_InterpolationXiSet(BasisGeometry,[BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisGeometry,[BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY],Err)
  ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
    CALL cmfe_Basis_InterpolationXiSet(BasisGeometry,[BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY, & 
      & BASIS_XI_INTERPOLATION_GEOMETRY],Err)                         
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisGeometry,[BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY, &
      & BASIS_XI_GAUSS_GEOMETRY],Err)
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(BasisGeometry,Err)
  !
  !Start the creation of another basis: Concentration_One
!  IF(BASIS_XI_INTERPOLATION_CONC_ONE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisConcOne=BasisGeometry
!  ELSE
!     MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
!     !Initialise a new velocity basis
!     CALL cmfe_Basis_Initialise(BasisConcOne,Err)
!     !Start the creation of a basis
!     CALL cmfe_Basis_CreateStart(BASIS_NUMBER_CONC_ONE,BasisConcOne,Err)
!     !Set the basis type (Lagrange/Simplex)
!     CALL cmfe_Basis_TypeSet(BasisConcOne,BASIS_TYPE,Err)
!     !Set the basis xi number
!     CALL cmfe_Basis_NumberOfXiSet(BasisConcOne,NUMBER_OF_DIMENSIONS,Err)
!     !Set the basis xi interpolation and number of Gauss points
!     IF(NUMBER_OF_DIMENSIONS==2) THEN
!       CALL cmfe_Basis_InterpolationXiSet(BasisConcOne,[BASIS_XI_INTERPOLATION_CONC_ONE,BASIS_XI_INTERPOLATION_CONC_ONE],Err)
!       CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisConcOne,[BASIS_XI_GAUSS_CONC_ONE,BASIS_XI_GAUSS_CONC_ONE],Err)
!     ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
!       CALL cmfe_Basis_InterpolationXiSet(BasisConcOne,[BASIS_XI_INTERPOLATION_CONC_ONE,BASIS_XI_INTERPOLATION_CONC_ONE, & 
!         & BASIS_XI_INTERPOLATION_CONC_ONE],Err)                         
!       CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisConcOne,[BASIS_XI_GAUSS_CONC_ONE,BASIS_XI_GAUSS_CONC_ONE, & 
!         & BASIS_XI_GAUSS_CONC_ONE],Err)
!     ENDIF
!     !Finish the creation of the basis
!     CALL cmfe_Basis_CreateFinish(BasisConcOne,Err)
!   ENDIF
  !
  !Start the creation of another basis: Concentration_Two
!  IF(BASIS_XI_INTERPOLATION_CONC_TWO==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisConcTwo=BasisGeometry
!  ELSE IF(BASIS_XI_INTERPOLATION_CONC_TWO==BASIS_XI_INTERPOLATION_CONC_ONE) THEN
!     BasisConcTwo=BasisConcOne
!   ELSE
!     MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
!     !Initialise a new concentration basis
!     CALL cmfe_Basis_Initialise(BasisConcTwo,Err)
!     !Start the creation of a basis
!     CALL cmfe_Basis_CreateStart(BASIS_NUMBER_CONC_TWO,BasisConcTwo,Err)
!     !Set the basis type (Lagrange/Simplex)
!     CALL cmfe_Basis_TypeSet(BasisConcTwo,BASIS_TYPE,Err)
!     !Set the basis xi number
!     CALL cmfe_Basis_NumberOfXiSet(BasisConcTwo,NUMBER_OF_DIMENSIONS,Err)
!     !Set the basis xi interpolation and number of Gauss points
!     IF(NUMBER_OF_DIMENSIONS==2) THEN
!       CALL cmfe_Basis_InterpolationXiSet(BasisConcTwo,[BASIS_XI_INTERPOLATION_CONC_TWO,BASIS_XI_INTERPOLATION_CONC_TWO],Err)
!       CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisConcTwo,[BASIS_XI_GAUSS_CONC_TWO,BASIS_XI_GAUSS_CONC_TWO],Err)
!     ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
!       CALL cmfe_Basis_InterpolationXiSet(BasisConcTwo,[BASIS_XI_INTERPOLATION_CONC_TWO,BASIS_XI_INTERPOLATION_CONC_TWO, & 
!         & BASIS_XI_INTERPOLATION_CONC_TWO],Err)                         
!       CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisConcTwo,[BASIS_XI_GAUSS_CONC_TWO,BASIS_XI_GAUSS_CONC_TWO, & 
!         & BASIS_XI_GAUSS_CONC_TWO],Err)
!     ENDIF
!     !Finish the creation of the basis
!     CALL cmfe_Basis_CreateFinish(BasisConcTwo,Err)
!   ENDIF

    BasisConcThree=BasisGeometry
  !
  !================================================================================================================================
  !
  !MESH
  !All types of physics utilize the same "mesh", but may be represented on individual mesh components.
  TotalNumberOfSolidNodes = NUMBER_OF_NODES_GEOMETRY
  TOTAL_NUMBER_OF_ALL_NODES = TOTAL_NUMBER_OF_NODES + TotalNumberOfSolidNodes
  MESH_NUMBER_OF_ALL_COMPONENTS = MESH_NUMBER_OF_COMPONENTS

  !Start the creation of mesh nodes
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,TOTAL_NUMBER_OF_ALL_NODES,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,MESH_NUMBER_OF_ALL_COMPONENTS,Err)
  !
  CALL cmfe_MeshElements_Initialise(MeshElementsGeometry,Err)
  CALL cmfe_MeshElements_Initialise(MeshElementsConcOne,Err)
  CALL cmfe_MeshElements_Initialise(MeshElementsConcTwo,Err)
  CALL cmfe_MeshElements_Initialise(MeshElementsConcThree,Err)
  MESH_COMPONENT_NUMBER_GEOMETRY=1
  MESH_COMPONENT_NUMBER_CONC_ONE=1
  MESH_COMPONENT_NUMBER_CONC_TWO=1
  MESH_COMPONENT_NUMBER_CONC_THREE=1
  !Specify spatial mesh component
  CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_GEOMETRY,BasisGeometry,MeshElementsGeometry,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL cmfe_MeshElements_NodesSet(MeshElementsGeometry,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_GEOMETRY),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(MeshElementsGeometry,Err)
  !Specify concentration one mesh component
! !  IF(BASIS_XI_INTERPOLATION_CONC_ONE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
     MeshElementsConcOne=MeshElementsGeometry
! !   ELSE
! !     MESH_COMPONENT_NUMBER_CONC_ONE=MESH_COMPONENT_NUMBER_GEOMETRY+1
!      CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_ONE,BasisConcOne,MeshElementsConcOne,Err)
! !     DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
! !       CALL cmfe_MeshElements_NodesSet(MeshElementsConcOne,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
! !         & 1:NUMBER_OF_ELEMENT_NODES_CONC_ONE),Err)
! !     ENDDO
 !    CALL cmfe_MeshElements_CreateFinish(MeshElementsConcOne,Err)
! !  ENDIF
!   !Specify concentration two mesh component
! !  IF(BASIS_XI_INTERPOLATION_CONC_TWO==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
     MeshElementsConcTwo=MeshElementsGeometry
! !    MESH_COMPONENT_NUMBER_CONC_TWO=MESH_COMPONENT_NUMBER_GEOMETRY
! !  ELSE IF(BASIS_XI_INTERPOLATION_CONC_TWO==BASIS_XI_INTERPOLATION_CONC_ONE) THEN
! !    MeshElementsConcTwo=MeshElementsConcOne
! !    MESH_COMPONENT_NUMBER_CONC_TWO=MESH_COMPONENT_NUMBER_CONC_ONE
! !  ELSE
! !    MESH_COMPONENT_NUMBER_CONC_TWO=MESH_COMPONENT_NUMBER_CONC_ONE+1
!     CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_TWO,BasisConcTwo,MeshElementsConcTwo,Err)
! !    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
! !      CALL cmfe_MeshElements_NodesSet(MeshElementsConcTwo,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
! !        & 1:NUMBER_OF_ELEMENT_NODES_CONC_TWO),Err)
! !    ENDDO
!     CALL cmfe_MeshElements_CreateFinish(MeshElementsConcTwo,Err)
!  ! ENDIF
!   !Specify concentration three mesh component
! !  IF(BASIS_XI_INTERPOLATION_CONC_THREE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
     MeshElementsConcThree=MeshElementsGeometry
! !    MESH_COMPONENT_NUMBER_CONC_THREE=MESH_COMPONENT_NUMBER_GEOMETRY
! !   ELSE IF(BASIS_XI_INTERPOLATION_CONC_THREE==BASIS_XI_INTERPOLATION_CONC_ONE) THEN
! !     MeshElementsConcThree=MeshElementsConcOne
! !    MESH_COMPONENT_NUMBER_CONC_THREE=MESH_COMPONENT_NUMBER_CONC_ONE
! !  ELSE
! !    MESH_COMPONENT_NUMBER_CONC_TWO=MESH_COMPONENT_NUMBER_CONC_ONE+1
!     CALL cmfe_MeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_THREE,BasisConcThree,MeshElementsConcThree,Err)
! !    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
! !      CALL cmfe_MeshElements_NodesSet(MeshElementsConcTwo,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
! !        & 1:NUMBER_OF_ELEMENT_NODES_CONC_TWO),Err)
! !    ENDDO
!     CALL cmfe_MeshElements_CreateFinish(MeshElementsConcThree,Err)
! !  ENDIF

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
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfUserDomains,Err)
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
        & 1,CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL cmfe_Field_ParameterSetUpdateStart(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  !
  !================================================================================================================================
  !
  !EQUATIONS SETS - USING NEW ARGUMENTS TO ALLOW FOR MULTI-COMPARTMENT MODELS
  DO icompartment = 1,Ncompartments

    EquationsSetFieldUserNumberDiffusion = 100_CMISSIntg+icompartment
    EquationsSetUserNumberDiffusion = 200_CMISSIntg+icompartment
    CALL cmfe_Field_Initialise(EquationsSetFieldDiffusion(icompartment),Err)
    CALL cmfe_EquationsSet_Initialise(EquationsSetDiffusion(icompartment),Err)
    CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumberDiffusion,Region,GeometricField, &
      & [CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS,CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE, &
      & CMFE_EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE],EquationsSetFieldUserNumberDiffusion, &
      & EquationsSetFieldDiffusion(icompartment),EquationsSetDiffusion(icompartment),Err)
    !Finish creating the equations set
    CALL cmfe_EquationsSet_CreateFinish(EquationsSetDiffusion(icompartment),Err)
    CALL cmfe_Field_ParameterSetUpdateConstant(RegionUserNumber,EquationsSetFielduserNumberDiffusion,CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE,1,icompartment,Err)
    CALL cmfe_Field_ParameterSetUpdateConstant(EquationsSetFieldDiffusion(icompartment),CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE,2,Ncompartments,Err)
  END DO 

  !-------------------------------------------------------------------------------------
  ! DEPENDENT FIELD: Shared
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_Field_CreateStart(DependentFieldUserNumber,Region,DependentField,Err)
  CALL cmfe_Field_TypeSet(DependentField,CMFE_FIELD_GENERAL_TYPE,Err)  
  CALL cmfe_Field_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentField,GeometricField,Err) 
  CALL cmfe_Field_DependentTypeSet(DependentField,CMFE_FIELD_DEPENDENT_TYPE,Err) 
  !Create 2N number of variables
  CALL cmfe_Field_NumberOfVariablesSet(DependentField,2*Ncompartments,Err) 
  !create two variables for each compartment
  ALLOCATE(VariableTypes(2*Ncompartments))
  DO num_var=1,Ncompartments
     VariableTypes(2*num_var-1)=CMFE_FIELD_U_VARIABLE_TYPE+(CMFE_FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
     VariableTypes(2*num_var)=CMFE_FIELD_DELUDELN_VARIABLE_TYPE+(CMFE_FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
  ENDDO
  CALL cmfe_Field_VariableTypesSet(DependentField,VariableTypes,Err) 
  !loop over the number of compartments
  DO icompartment=1,2*Ncompartments
    !set dimension type
    CALL cmfe_Field_DimensionSet(DependentField,VariableTypes(icompartment), &
       & CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
    CALL cmfe_Field_NumberOfComponentsSet(DependentField,VariableTypes(icompartment),1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,VariableTypes(icompartment),1, & 
       & MESH_COMPONENT_NUMBER_CONC_ONE,Err)
  ENDDO
  CALL cmfe_Field_CreateFinish(DependentField,Err)
  
  DO icompartment = 1,Ncompartments
    CALL cmfe_EquationsSet_DependentCreateStart(EquationsSetDiffusion(icompartment),DependentFieldUserNumber,&
      & DependentField,Err)
    CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSetDiffusion(icompartment),Err)
  ENDDO

!   !Set the mesh component to be used by the field components.
! !   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
! !     CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !     CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !   ENDDO
!   !Set the mesh component to be used by the field components.
! !   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
! !     CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !     CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !   ENDDO
!   !Set the mesh component to be used by the field components.
! !   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
! !     CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMFE_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !     CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !   ENDDO
!   !Finish the equations set dependent field variables

  !-------------------------------------------------------------------------------------
  ! INITIALISE DEPENDENT FIELDS
!   !Initialise dependent field (concentration one components)
!   CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
!     & 1,INITIAL_FIELD_DIFFUSION_ONE,Err)

!   !Initialise dependent field (concentration two components)
!   CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U2_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
!     & 1,INITIAL_FIELD_DIFFUSION_TWO,Err)

!   !Initialise dependent field (concentration three components)
!   CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U3_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
!     & 1,INITIAL_FIELD_DIFFUSION_THREE,Err)


  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS - create the materials field
  !Auto-created field contains a U variable type to store the diffusion coefficient(s)
  !It also contains a V variable type to store the coupling coefficients 
  DO icompartment = 1,Ncompartments
    MaterialsFieldUserNumberDiffusion = 400+icompartment
    CALL cmfe_Field_Initialise(MaterialsFieldDiffusion(icompartment),Err)
    CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSetDiffusion(icompartment),MaterialsFieldUserNumberDiffusion,&
         & MaterialsFieldDiffusion(icompartment),Err)
    CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSetDiffusion(icompartment),Err)
  END DO
  !Initialise the coupling coefficients
  !Need to devise a neater way of specifying these components - e.g. specify only the upper diagonal components, and then automatically fill out the rest
  ALLOCATE(CouplingCoeffs(Ncompartments,Ncompartments))
  CouplingCoeffs(1,1)=1.0_CMISSRP
  CouplingCoeffs(1,2)=1.0_CMISSRP
  CouplingCoeffs(1,3)=1.0_CMISSRP
  CouplingCoeffs(2,1)=1.0_CMISSRP
  CouplingCoeffs(2,2)=1.0_CMISSRP
  CouplingCoeffs(2,3)=1.0_CMISSRP
  CouplingCoeffs(3,1)=1.0_CMISSRP
  CouplingCoeffs(3,2)=1.0_CMISSRP
  CouplingCoeffs(3,3)=1.0_CMISSRP

  DO icompartment = 1, Ncompartments
    DO COMPONENT_NUMBER=1, Ncompartments
      CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDiffusion(icompartment),CMFE_FIELD_V_VARIABLE_TYPE, &
         & CMFE_FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
    END DO
  END DO

!   CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDiffusionOne,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
!     & MaterialsFieldUserNumberDiffusionOne,POROSITY_PARAM_MAT_PROPERTIES,Err)
!   CALL cmfe_Field_ComponentValuesInitialise(MaterialsFieldDiffusionTwo,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
!     & MaterialsFieldUserNumberMatPropertiesPermOverVis,PERM_OVER_VIS_PARAM_MAT_PROPERTIES,Err)

  !
  !================================================================================================================================
  !

  !SOURCE FIELDS

   !create the equations set source field variables for both equations sets 
  DO icompartment = 1,Ncompartments
    SourceFieldUserNumberDiffusion=700_CMISSIntg+icompartment
    CALL cmfe_Field_Initialise(SourceFieldDiffusion(icompartment),Err)
    CALL cmfe_EquationsSet_SourceCreateStart(EquationsSetDiffusion(icompartment),SourceFieldUserNumberDiffusion,&
         & SourceFieldDiffusion(icompartment),Err)
    CALL cmfe_EquationsSet_SourceCreateFinish(EquationsSetDiffusion(icompartment),Err)
  END DO 

  !
  !================================================================================================================================
  !

  !EQUATIONS
  DO icompartment=1,Ncompartments
    CALL cmfe_Equations_Initialise(EquationsDiffusion(icompartment),Err)
    CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSetDiffusion(icompartment),EquationsDiffusion(icompartment),Err)
    CALL cmfe_Equations_SparsityTypeSet(EquationsDiffusion(icompartment),CMFE_EQUATIONS_SPARSE_MATRICES,Err)
    CALL cmfe_Equations_OutputTypeSet(EquationsDiffusion(icompartment),EQUATIONS_DIFFUSION_OUTPUT,Err)
    CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSetDiffusion(icompartment),Err)
  ENDDO
  !
  !================================================================================================================================
  !
  !PROBLEMS

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
    & CMFE_PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE,CMFE_PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,LINEAR_SOLVER_DIFFUSION_START_TIME,&
    & LINEAR_SOLVER_DIFFUSION_STOP_TIME,LINEAR_SOLVER_DIFFUSION_TIME_INCREMENT,Err)
  !Set the output timing
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,LINEAR_SOLVER_DIFFUSION_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !
  !SOLVERS

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(SolverDiffusion,Err)
  CALL cmfe_Solver_Initialise(LinearSolverDiffusion,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)

  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverDiffusionUserNumber,SolverDiffusion,Err)
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
  !SOLVER EQUATIONS
  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(SolverDiffusion,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsDiffusion,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,SolverDiffusionUserNumber,SolverDiffusion,Err)
  CALL cmfe_Solver_SolverEquationsGet(SolverDiffusion,SolverEquationsDiffusion,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsDiffusion,CMFE_SOLVER_SPARSE_MATRICES,Err)

  DO icompartment=1,Ncompartments
    CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsDiffusion,EquationsSetDiffusion(icompartment),&
         & EquationsSetIndex,Err)
  ENDDO

  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsDiffusion,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquationsDiffusion,BoundaryConditionsDiffusion,Err)
  DO icompartment=1,Ncompartments
!     IF(INLET_WALL_NODES_DIFFUSION_FLAG(icompartment)) THEN
!       DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION(icompartment)
!         NODE_NUMBER=INLET_WALL_NODES_DIFFUSION(icompartment,NODE_COUNTER)
!         CONDITION=CMFE_BOUNDARY_CONDITION_FIXED
!           VALUE=0.2_CMISSRP
!           CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDiffusion(icompartment),CMFE_FIELD_U_VARIABLE_TYPE,CMFE_NO_GLOBAL_DERIV, & 
!             & NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
!       ENDDO
!     ENDIF
   IF(icompartment==1)THEN
    IF(INLET_WALL_NODES_DIFFUSION_ONE_FLAG) THEN
      DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION_ONE
        NODE_NUMBER=INLET_WALL_NODES_DIFFUSION_ONE(NODE_COUNTER)
        CONDITION=CMFE_BOUNDARY_CONDITION_FIXED
          VALUE=0.2_CMISSRP
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDiffusion,DependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
            & 1,CMFE_NO_GLOBAL_DERIV, & 
            & NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
      ENDDO
    ENDIF
   ENDIF

   IF(icompartment==2)THEN
    IF(INLET_WALL_NODES_DIFFUSION_TWO_FLAG) THEN
      DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION_TWO
        NODE_NUMBER=INLET_WALL_NODES_DIFFUSION_TWO(NODE_COUNTER)
        CONDITION=CMFE_BOUNDARY_CONDITION_FIXED
          VALUE=0.4_CMISSRP
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDiffusion,DependentField,CMFE_FIELD_V_VARIABLE_TYPE,1, &
            & CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
      ENDDO
    ENDIF
   ENDIF

   IF(icompartment==3)THEN
    IF(INLET_WALL_NODES_DIFFUSION_THREE_FLAG) THEN
      DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION_THREE
        NODE_NUMBER=INLET_WALL_NODES_DIFFUSION_THREE(NODE_COUNTER)
        CONDITION=CMFE_BOUNDARY_CONDITION_FIXED
          VALUE=0.6_CMISSRP
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDiffusion,DependentField,CMFE_FIELD_U1_VARIABLE_TYPE,1, &
            & CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
      ENDDO
    ENDIF
   ENDIF

   IF(icompartment==4)THEN
    IF(INLET_WALL_NODES_DIFFUSION_THREE_FLAG) THEN
      DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION_THREE
        NODE_NUMBER=INLET_WALL_NODES_DIFFUSION_THREE(NODE_COUNTER)
        CONDITION=CMFE_BOUNDARY_CONDITION_FIXED
          VALUE=0.8_CMISSRP
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDiffusion,DependentField,CMFE_FIELD_U2_VARIABLE_TYPE,1, &
            & CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
      ENDDO
    ENDIF
   ENDIF

   IF(icompartment==5)THEN
    IF(INLET_WALL_NODES_DIFFUSION_THREE_FLAG) THEN
      DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION_THREE
        NODE_NUMBER=INLET_WALL_NODES_DIFFUSION_THREE(NODE_COUNTER)
        CONDITION=CMFE_BOUNDARY_CONDITION_FIXED
          VALUE=1.0_CMISSRP
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsDiffusion,DependentField,CMFE_FIELD_U3_VARIABLE_TYPE,1, &
            & CMFE_NO_GLOBAL_DERIV,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
      ENDDO
    ENDIF
   ENDIF

  ENDDO
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquationsDiffusion,Err)

  !
  !================================================================================================================================
  !
  !RUN SOLVERS

  !Turn off PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  CALL cmfe_Problem_Solve(Problem,Err)
  !
  !================================================================================================================================
  !
  !OUTPUT

  EXPORT_FIELD_IO=.TRUE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"MonolithicMultiCompDiffusionTest","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"MonolithicMultiCompDiffusionTest","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF

  !Finialise CMISS
!   CALL cmfe_Finalise(Err)
   WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM MONOLITHICSCHEMETESTEXAMPLE
