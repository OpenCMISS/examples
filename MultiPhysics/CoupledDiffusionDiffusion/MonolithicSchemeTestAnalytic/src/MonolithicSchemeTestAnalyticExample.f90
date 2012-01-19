!> \file
!> \authors Andrew Cookson
!> \brief This is an example program to solve coupled multi-compartment diffusion equations in monolithic scheme using openCMISS calls.
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

PROGRAM MONOLITHICSCHEMETESTFIELDMLEXAMPLE

  !
  !================================================================================================================================
  !

  !PROGRAM LIBRARIES

  USE OPENCMISS
  USE FIELDML_OUTPUT_ROUTINES
  USE FIELDML_UTIL_ROUTINES
  USE FIELDML_API
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

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=3.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=9
  INTEGER(CMISSIntg) :: MaterialsFieldUserNumberDiffusion
!   INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDiffusionTwo=10
!   INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDiffusionThree=11
  INTEGER(CMISSIntg) :: SourceFieldUserNumberDiffusion
!   INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumberDiffusionTwo=16
!   INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumberDiffusionThree=17
  INTEGER(CMISSIntg) :: EquationsSetUserNumberDiffusion
  INTEGER(CMISSIntg) :: AnalyticFieldUserNumberDiffusion
!   INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberDiffusionTwo=13
!   INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberDiffusionThree=14
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=20
!!!!!!!!!!!!!!!!!!!EQUATIONS FIELD PARAMETERS
  INTEGER(CMISSIntg), PARAMETER :: EquationsFieldDiffusionOne=33
  INTEGER(CMISSIntg), PARAMETER :: EquationsFieldDiffusionTwo=34
  INTEGER(CMISSIntg), PARAMETER :: EquationsFieldDiffusionThree=35





  INTEGER(CMISSIntg) :: EquationsSetFieldUserNumberDiffusion
  INTEGER(CMISSIntg) :: icompartment,Ncompartments,num_var
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: NumberOfUserDomains=1
  INTEGER(CMISSIntg), PARAMETER :: SolverDiffusionUserNumber=1
  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables
  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS

  INTEGER(CMISSIntg) :: MPI_IERROR

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


  REAL(CMISSDP) :: DOMAIN_X1, DOMAIN_X2, DOMAIN_Y1, DOMAIN_Y2, DOMAIN_Z1, DOMAIN_Z2
  REAL(CMISSDP) :: GEOMETRY_TOLERANCE

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_DIFFUSION_ONE
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_DIFFUSION_ONE
  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_DIFFUSION_TWO
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_DIFFUSION_TWO
  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_DIFFUSION_THREE
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_DIFFUSION_THREE

  REAL(CMISSDP) :: INITIAL_FIELD_DIFFUSION_ONE
  REAL(CMISSDP) :: INITIAL_FIELD_DIFFUSION_TWO
  REAL(CMISSDP) :: INITIAL_FIELD_DIFFUSION_THREE
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_DIFFUSION_ONE
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_DIFFUSION_TWO
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_DIFFUSION_THREE
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE

  REAL(CMISSDP) :: LINEAR_SOLVER_DIFFUSION_START_TIME
  REAL(CMISSDP) :: LINEAR_SOLVER_DIFFUSION_STOP_TIME
  REAL(CMISSDP) :: LINEAR_SOLVER_DIFFUSION_TIME_INCREMENT

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_DIFFUSION_DIRECT_FLAG
  LOGICAL :: INLET_WALL_NODES_DIFFUSION_ONE_FLAG
  LOGICAL :: INLET_WALL_NODES_DIFFUSION_TWO_FLAG
  LOGICAL :: INLET_WALL_NODES_DIFFUSION_THREE_FLAG
  !CMISS variables
  CHARACTER(C_CHAR), PARAMETER :: NUL=C_NULL_CHAR
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
  TYPE(CMISSBasisType) :: BasisConcThree
  
  !TYPE(CMISSBasisType), ALLOCATABLE, DIMENSION(:) :: BasicConc

  !Nodes
  TYPE(CMISSNodesType) :: Nodes
  !Elements
  TYPE(CMISSMeshElementsType) :: MeshElementsGeometry
  TYPE(CMISSMeshElementsType) :: MeshElementsConcOne
  TYPE(CMISSMeshElementsType) :: MeshElementsConcTwo
  TYPE(CMISSMeshElementsType) :: MeshElementsConcThree

  TYPE(CMISSMeshElementsType), ALLOCATABLE, DIMENSION(:) :: MeshElementsConc
  !Meshes
  TYPE(CMISSMeshType) :: Mesh
  !Decompositions
  TYPE(CMISSDecompositionType) :: Decomposition
  !Fields
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: DependentField
  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: MaterialsFieldDiffusion
  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: SourceFieldDiffusion
  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: AnalyticFieldDiffusion
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsDiffusion
  TYPE(CMISSEquationsSetType), ALLOCATABLE, DIMENSION(:) :: EquationsSetDiffusion
  TYPE(CMISSEquationsType), ALLOCATABLE, DIMENSION(:) :: EquationsDiffusion
  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: EquationsSetFieldDiffusion
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: SolverDiffusion
  TYPE(CMISSSolverType) :: LinearSolverDiffusion
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsDiffusion

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: EquationsSetIndexTwo
  INTEGER(CMISSIntg) :: EquationsSetIndexThree
  INTEGER(CMISSIntg) :: Err


  !Array containing the field variable types that will be used (for ease of incorporating inside a loop)
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:) :: VariableTypes
  REAL(CMISSDP), ALLOCATABLE, DIMENSION(:,:) :: CouplingCoeffs

  INTEGER(CMISSIntg) :: DIAG_LEVEL_LIST(5)
!   CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(8) !,TIMING_ROUTINE_LIST(1)
  CHARACTER(LEN=255) :: DIAG_ROUTINE_LIST(1) !,TIMING_ROUTINE_LIST(1)
 
  INTEGER(CMISSIntg) :: TotalNumberOfSolidNodes
!   INTEGER(CMISSIntg) :: NumberOfSolidMeshComponents
  !
  !--------------------------------------------------------------------------------------------------------------------------------
  !


  !FieldML variables
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputDirectory = ""
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: outputFilename = "MonolithicMultiCompDiffusion.xml"
  CHARACTER(KIND=C_CHAR,LEN=*), PARAMETER :: basename = "monolithic_multicomp_diffusion"

  TYPE(FieldmlInfoType) :: fieldmlInfo




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
  !CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
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
!   BASIS_XI_INTERPOLATION_GEOMETRY=2
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
!   DOMAIN_X1 =  0.0_CMISSDP
!   DOMAIN_X2 =  1.0_CMISSDP
!   DOMAIN_Y1 =  0.0_CMISSDP
!   DOMAIN_Y2 =  1.0_CMISSDP
!   DOMAIN_Z1 =  0.0_CMISSDP
!   DOMAIN_Z2 =  1.0_CMISSDP
!   !Set geometric tolerance
!   GEOMETRY_TOLERANCE = 1.0E-12_CMISSDP
!   !Set initial values
!   INITIAL_FIELD_DIFFUSION_ONE=1.0_CMISSDP
!   INITIAL_FIELD_DIFFUSION_TWO=1.0_CMISSDP
!   INITIAL_FIELD_DIFFUSION_THREE=1.0_CMISSDP
!   !Set initial boundary conditions
!   INLET_WALL_NODES_DIFFUSION_ONE_FLAG=.TRUE.
!   IF(INLET_WALL_NODES_DIFFUSION_ONE_FLAG) THEN
!     NUMBER_OF_INLET_WALL_NODES_DIFFUSION_ONE=36
!     ALLOCATE(INLET_WALL_NODES_DIFFUSION_ONE(NUMBER_OF_INLET_WALL_NODES_DIFFUSION_ONE))
!     INLET_WALL_NODES_DIFFUSION_ONE=(/191,155,119,83,23,21,192,156,120,84,24,22,&
!      & 198,162,126,90,36,35,204,168,132,96,48,47,210,174,138,102,60,59,216,180,144,108,72,71/)
!     !Set initial boundary conditions
!     BOUNDARY_CONDITIONS_DIFFUSION_ONE=1.0_CMISSDP
!   ENDIF  !Set material parameters
! 
!   INLET_WALL_NODES_DIFFUSION_TWO_FLAG=.TRUE.
!   IF(INLET_WALL_NODES_DIFFUSION_TWO_FLAG) THEN
!     NUMBER_OF_INLET_WALL_NODES_DIFFUSION_TWO=36
!     ALLOCATE(INLET_WALL_NODES_DIFFUSION_TWO(NUMBER_OF_INLET_WALL_NODES_DIFFUSION_TWO))
!     INLET_WALL_NODES_DIFFUSION_TWO=(/191,155,119,83,23,21,192,156,120,84,24,22,&
!      & 198,162,126,90,36,35,204,168,132,96,48,47,210,174,138,102,60,59,216,180,144,108,72,71/)
!     !Set initial boundary conditions
!     BOUNDARY_CONDITIONS_DIFFUSION_TWO=1.0_CMISSDP
!   ENDIF  !Set material parameters
! 
!   INLET_WALL_NODES_DIFFUSION_THREE_FLAG=.TRUE.
!   IF(INLET_WALL_NODES_DIFFUSION_THREE_FLAG) THEN
!     NUMBER_OF_INLET_WALL_NODES_DIFFUSION_THREE=36
!     ALLOCATE(INLET_WALL_NODES_DIFFUSION_THREE(NUMBER_OF_INLET_WALL_NODES_DIFFUSION_THREE))
!     INLET_WALL_NODES_DIFFUSION_THREE=(/191,155,119,83,23,21,192,156,120,84,24,22,&
!      & 198,162,126,90,36,35,204,168,132,96,48,47,210,174,138,102,60,59,216,180,144,108,72,71/)
!     !Set initial boundary conditions
!     BOUNDARY_CONDITIONS_DIFFUSION_THREE=1.0_CMISSDP
!   ENDIF  !Set material parameters

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NUMBER_GLOBAL_X_ELEMENTS=40
  NUMBER_GLOBAL_Y_ELEMENTS=40
  NUMBER_GLOBAL_Z_ELEMENTS=0

  IF(NUMBER_GLOBAL_Z_ELEMENTS==0)THEN
    NUMBER_OF_DIMENSIONS=2
  ELSE
    NUMBER_OF_DIMENSIONS=3
  ENDIF
  NUMBER_OF_DOMAINS=NumberOfComputationalNodes

  !Set all diganostic levels on for testing

  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)




  !Set material parameters
  !Set number of Gauss points (Mind that also material field may be interpolated)
  BASIS_XI_GAUSS_GEOMETRY=3 !4
  BASIS_XI_GAUSS_CONC_ONE=3 !4
  BASIS_XI_GAUSS_CONC_TWO=3 !4
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_DIFFUSION_OUTPUT_TYPE=CMISS_SOLVER_MATRIX_OUTPUT
  !LINEAR_SOLVER_DIFFUSION_OUTPUT_TYPE=CMISS_SOLVER_NO_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DIFFUSION_OUTPUT=CMISS_EQUATIONS_MATRIX_OUTPUT
  !EQUATIONS_DIFFUSION_OUTPUT=CMISS_EQUATIONS_NO_OUTPUT
  EQUATIONS_DIFFUSION_TWO_OUTPUT=CMISS_EQUATIONS_NO_OUTPUT
  EQUATIONS_DIFFUSION_THREE_OUTPUT=CMISS_EQUATIONS_NO_OUTPUT
  !Set time parameter
  LINEAR_SOLVER_DIFFUSION_START_TIME=0.0_CMISSDP
  LINEAR_SOLVER_DIFFUSION_STOP_TIME=0.3000001_CMISSDP 
  LINEAR_SOLVER_DIFFUSION_TIME_INCREMENT=0.005_CMISSDP
  !Set result output parameter
  LINEAR_SOLVER_DIFFUSION_OUTPUT_FREQUENCY=1
  !Set solver parameters
  LINEAR_SOLVER_DIFFUSION_DIRECT_FLAG=.FALSE.

  RELATIVE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E5_CMISSDP !default: 1.0E5
  MAXIMUM_ITERATIONS=10000_CMISSIntg !default: 100000
  RESTART_VALUE=30_CMISSIntg !default: 30
  LINESEARCH_ALPHA=1.0_CMISSDP


  icompartment =1_CMISSIntg
  Ncompartments=2_CMISSIntg
  !
  !================================================================================================================================
  !

  !Set diagnostics

   !DIAG_LEVEL_LIST(1)=1
!   DIAG_LEVEL_LIST(2)=2
!   DIAG_LEVEL_LIST(3)=3
!   DIAG_LEVEL_LIST(4)=4
!   DIAG_LEVEL_LIST(5)=5

!   DIAG_ROUTINE_LIST(1)="DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE"
!   DIAG_ROUTINE_LIST(2)="PROBLEM_SOLVER_EQUATIONS_SOLVE"
   !DIAG_ROUTINE_LIST(1)="SOLVER_DYNAMIC_CREATE_FINISH"

  !CMISS_ALL_DIAG_TYPE/CMISS_IN_DIAG_TYPE/CMISS_FROM_DIAG_TYPE
   !CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISS_ALL_TIMING_TYPE/CMISS_IN_TIMING_TYPE/CMISS_FROM_TIMING_TYPE
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

  !ALLOCATE THE ARRAYS
  ALLOCATE (EquationsSetDiffusion(Ncompartments))
  ALLOCATE (EquationsSetFieldDiffusion(Ncompartments))
  ALLOCATE (MaterialsFieldDiffusion(Ncompartments))
  ALLOCATE (SourceFieldDiffusion(Ncompartments))
  ALLOCATE (EquationsDiffusion(Ncompartments))
  ALLOCATE (AnalyticFieldDiffusion(Ncompartments))

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
!  IF(BASIS_XI_INTERPOLATION_CONC_ONE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisConcOne=BasisGeometry
!  ELSE
!     MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
!     !Initialise a new velocity basis
!     CALL CMISSBasis_Initialise(BasisConcOne,Err)
!     !Start the creation of a basis
!     CALL CMISSBasis_CreateStart(BASIS_NUMBER_CONC_ONE,BasisConcOne,Err)
!     !Set the basis type (Lagrange/Simplex)
!     CALL CMISSBasis_TypeSet(BasisConcOne,BASIS_TYPE,Err)
!     !Set the basis xi number
!     CALL CMISSBasis_NumberOfXiSet(BasisConcOne,NUMBER_OF_DIMENSIONS,Err)
!     !Set the basis xi interpolation and number of Gauss points
!     IF(NUMBER_OF_DIMENSIONS==2) THEN
!       CALL CMISSBasis_InterpolationXiSet(BasisConcOne,(/BASIS_XI_INTERPOLATION_CONC_ONE,BASIS_XI_INTERPOLATION_CONC_ONE/),Err)
!       CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisConcOne,(/BASIS_XI_GAUSS_CONC_ONE,BASIS_XI_GAUSS_CONC_ONE/),Err)
!     ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
!       CALL CMISSBasis_InterpolationXiSet(BasisConcOne,(/BASIS_XI_INTERPOLATION_CONC_ONE,BASIS_XI_INTERPOLATION_CONC_ONE, & 
!         & BASIS_XI_INTERPOLATION_CONC_ONE/),Err)                         
!       CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisConcOne,(/BASIS_XI_GAUSS_CONC_ONE,BASIS_XI_GAUSS_CONC_ONE, & 
!         & BASIS_XI_GAUSS_CONC_ONE/),Err)
!     ENDIF
!     !Finish the creation of the basis
!     CALL CMISSBasis_CreateFinish(BasisConcOne,Err)
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
!     CALL CMISSBasis_Initialise(BasisConcTwo,Err)
!     !Start the creation of a basis
!     CALL CMISSBasis_CreateStart(BASIS_NUMBER_CONC_TWO,BasisConcTwo,Err)
!     !Set the basis type (Lagrange/Simplex)
!     CALL CMISSBasis_TypeSet(BasisConcTwo,BASIS_TYPE,Err)
!     !Set the basis xi number
!     CALL CMISSBasis_NumberOfXiSet(BasisConcTwo,NUMBER_OF_DIMENSIONS,Err)
!     !Set the basis xi interpolation and number of Gauss points
!     IF(NUMBER_OF_DIMENSIONS==2) THEN
!       CALL CMISSBasis_InterpolationXiSet(BasisConcTwo,(/BASIS_XI_INTERPOLATION_CONC_TWO,BASIS_XI_INTERPOLATION_CONC_TWO/),Err)
!       CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisConcTwo,(/BASIS_XI_GAUSS_CONC_TWO,BASIS_XI_GAUSS_CONC_TWO/),Err)
!     ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
!       CALL CMISSBasis_InterpolationXiSet(BasisConcTwo,(/BASIS_XI_INTERPOLATION_CONC_TWO,BASIS_XI_INTERPOLATION_CONC_TWO, & 
!         & BASIS_XI_INTERPOLATION_CONC_TWO/),Err)                         
!       CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisConcTwo,(/BASIS_XI_GAUSS_CONC_TWO,BASIS_XI_GAUSS_CONC_TWO, & 
!         & BASIS_XI_GAUSS_CONC_TWO/),Err)
!     ENDIF
!     !Finish the creation of the basis
!     CALL CMISSBasis_CreateFinish(BasisConcTwo,Err)
!   ENDIF

    BasisConcThree=BasisGeometry
  !
  !================================================================================================================================
  !
!   !MESH
!   !All types of physics utilize the same "mesh", but may be represented on individual mesh components.
!   TotalNumberOfSolidNodes = NUMBER_OF_NODES_GEOMETRY
!   TOTAL_NUMBER_OF_ALL_NODES = TOTAL_NUMBER_OF_NODES + TotalNumberOfSolidNodes
!   MESH_NUMBER_OF_ALL_COMPONENTS = MESH_NUMBER_OF_COMPONENTS
! 
!   !Start the creation of mesh nodes
!   CALL CMISSNodes_Initialise(Nodes,Err)
!   CALL CMISSNodes_CreateStart(Region,TOTAL_NUMBER_OF_ALL_NODES,Nodes,Err)
!   CALL CMISSNodes_CreateFinish(Nodes,Err)
!   !Start the creation of the mesh
!   CALL CMISSMesh_Initialise(Mesh,Err)
!   CALL CMISSMesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
!   !Set number of mesh elements
!   CALL CMISSMesh_NumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
!   !Set number of mesh components
!   CALL CMISSMesh_NumberOfComponentsSet(Mesh,MESH_NUMBER_OF_ALL_COMPONENTS,Err)
!   !
!   CALL CMISSMeshElements_Initialise(MeshElementsGeometry,Err)
!   CALL CMISSMeshElements_Initialise(MeshElementsConcOne,Err)
!   CALL CMISSMeshElements_Initialise(MeshElementsConcTwo,Err)
!   CALL CMISSMeshElements_Initialise(MeshElementsConcThree,Err)
  MESH_COMPONENT_NUMBER_GEOMETRY=1
  MESH_COMPONENT_NUMBER_CONC_ONE=1
  MESH_COMPONENT_NUMBER_CONC_TWO=1
  MESH_COMPONENT_NUMBER_CONC_THREE=1
!   !Specify spatial mesh component
!   CALL CMISSMeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_GEOMETRY,BasisGeometry,MeshElementsGeometry,Err)
!   DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
!     CALL CMISSMeshElements_NodesSet(MeshElementsGeometry,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_GEOMETRY),Err)
!   ENDDO
!   CALL CMISSMeshElements_CreateFinish(MeshElementsGeometry,Err)
!   !Specify concentration one mesh component
! ! !  IF(BASIS_XI_INTERPOLATION_CONC_ONE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
!      MeshElementsConcOne=MeshElementsGeometry
! ! !   ELSE
! ! !     MESH_COMPONENT_NUMBER_CONC_ONE=MESH_COMPONENT_NUMBER_GEOMETRY+1
! !      CALL CMISSMeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_ONE,BasisConcOne,MeshElementsConcOne,Err)
! ! !     DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
! ! !       CALL CMISSMeshElements_NodesSet(MeshElementsConcOne,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
! ! !         & 1:NUMBER_OF_ELEMENT_NODES_CONC_ONE),Err)
! ! !     ENDDO
!  !    CALL CMISSMeshElements_CreateFinish(MeshElementsConcOne,Err)
! ! !  ENDIF
! !   !Specify concentration two mesh component
! ! !  IF(BASIS_XI_INTERPOLATION_CONC_TWO==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
!      MeshElementsConcTwo=MeshElementsGeometry
! ! !    MESH_COMPONENT_NUMBER_CONC_TWO=MESH_COMPONENT_NUMBER_GEOMETRY
! ! !  ELSE IF(BASIS_XI_INTERPOLATION_CONC_TWO==BASIS_XI_INTERPOLATION_CONC_ONE) THEN
! ! !    MeshElementsConcTwo=MeshElementsConcOne
! ! !    MESH_COMPONENT_NUMBER_CONC_TWO=MESH_COMPONENT_NUMBER_CONC_ONE
! ! !  ELSE
! ! !    MESH_COMPONENT_NUMBER_CONC_TWO=MESH_COMPONENT_NUMBER_CONC_ONE+1
! !     CALL CMISSMeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_TWO,BasisConcTwo,MeshElementsConcTwo,Err)
! ! !    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
! ! !      CALL CMISSMeshElements_NodesSet(MeshElementsConcTwo,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
! ! !        & 1:NUMBER_OF_ELEMENT_NODES_CONC_TWO),Err)
! ! !    ENDDO
! !     CALL CMISSMeshElements_CreateFinish(MeshElementsConcTwo,Err)
! !  ! ENDIF
! !   !Specify concentration three mesh component
! ! !  IF(BASIS_XI_INTERPOLATION_CONC_THREE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
!      MeshElementsConcThree=MeshElementsGeometry
! ! !    MESH_COMPONENT_NUMBER_CONC_THREE=MESH_COMPONENT_NUMBER_GEOMETRY
! ! !   ELSE IF(BASIS_XI_INTERPOLATION_CONC_THREE==BASIS_XI_INTERPOLATION_CONC_ONE) THEN
! ! !     MeshElementsConcThree=MeshElementsConcOne
! ! !    MESH_COMPONENT_NUMBER_CONC_THREE=MESH_COMPONENT_NUMBER_CONC_ONE
! ! !  ELSE
! ! !    MESH_COMPONENT_NUMBER_CONC_TWO=MESH_COMPONENT_NUMBER_CONC_ONE+1
! !     CALL CMISSMeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_THREE,BasisConcThree,MeshElementsConcThree,Err)
! ! !    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
! ! !      CALL CMISSMeshElements_NodesSet(MeshElementsConcTwo,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
! ! !        & 1:NUMBER_OF_ELEMENT_NODES_CONC_TWO),Err)
! ! !    ENDDO
! !     CALL CMISSMeshElements_CreateFinish(MeshElementsConcThree,Err)
! ! !  ENDIF
! 
!   !Finish the creation of the mesh
!   CALL CMISSMesh_CreateFinish(Mesh,Err)

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,BasisGeometry,Err)   
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/WIDTH,HEIGHT/),Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/),Err)
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,(/WIDTH,HEIGHT,LENGTH/),Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS/),Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

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
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfUserDomains,Err)
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
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
  !Update the geometric field parameters
!   DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
!     DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!       VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
!       CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
!         & CMISS_NO_GLOBAL_DERIV,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
!     ENDDO
!   ENDDO
!   CALL CMISSField_ParameterSetUpdateStart(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
!   CALL CMISSField_ParameterSetUpdateFinish(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !
  !================================================================================================================================
  !
  !EQUATIONS SETS - USING NEW ARGUMENTS TO ALLOW FOR MULTI-COMPARTMENT MODELS
  DO icompartment = 1,Ncompartments

    EquationsSetFieldUserNumberDiffusion = 100_CMISSIntg+icompartment
    EquationsSetUserNumberDiffusion = 200_CMISSIntg+icompartment
    CALL CMISSField_Initialise(EquationsSetFieldDiffusion(icompartment),Err)
    CALL CMISSEquationsSet_Initialise(EquationsSetDiffusion(icompartment),Err)
    CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumberDiffusion,Region,GeometricField,&
         & CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
         & CMISS_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,CMISS_EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE,&
         & EquationsSetFieldUserNumberDiffusion,EquationsSetFieldDiffusion(icompartment),EquationsSetDiffusion(icompartment),&
         & Err)
    !Finish creating the equations set
    CALL CMISSEquationsSet_CreateFinish(EquationsSetDiffusion(icompartment),Err)
    CALL CMISSField_ParameterSetUpdateConstant(RegionUserNumber,EquationsSetFielduserNumberDiffusion,CMISS_FIELD_U_VARIABLE_TYPE, &
      & CMISS_FIELD_VALUES_SET_TYPE,1,icompartment,Err)
    CALL CMISSField_ParameterSetUpdateConstant(EquationsSetFieldDiffusion(icompartment),CMISS_FIELD_U_VARIABLE_TYPE, &
      & CMISS_FIELD_VALUES_SET_TYPE,2,Ncompartments,Err)
  END DO 

  !-------------------------------------------------------------------------------------
  ! DEPENDENT FIELD: Shared
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSField_CreateStart(DependentFieldUserNumber,Region,DependentField,Err)
  CALL CMISSField_TypeSet(DependentField,CMISS_FIELD_GENERAL_TYPE,Err)  
  CALL CMISSField_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(DependentField,GeometricField,Err) 
  CALL CMISSField_DependentTypeSet(DependentField,CMISS_FIELD_DEPENDENT_TYPE,Err) 
  !Create 2N number of variables
  CALL CMISSField_NumberOfVariablesSet(DependentField,2*Ncompartments,Err) 
  !create two variables for each compartment
  ALLOCATE(VariableTypes(2*Ncompartments))
  DO num_var=1,Ncompartments
     VariableTypes(2*num_var-1)=CMISS_FIELD_U_VARIABLE_TYPE+(CMISS_FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
     VariableTypes(2*num_var)=CMISS_FIELD_DELUDELN_VARIABLE_TYPE+(CMISS_FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
  ENDDO
  CALL CMISSField_VariableTypesSet(DependentField,VariableTypes,Err) 
  !loop over the number of compartments
  DO icompartment=1,2*Ncompartments
    !set dimension type
    CALL CMISSField_DimensionSet(DependentField,VariableTypes(icompartment), &
       & CMISS_FIELD_SCALAR_DIMENSION_TYPE,Err)
    CALL CMISSField_NumberOfComponentsSet(DependentField,VariableTypes(icompartment),1,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,VariableTypes(icompartment),1, & 
       & MESH_COMPONENT_NUMBER_CONC_ONE,Err)
  ENDDO
  CALL CMISSField_CreateFinish(DependentField,Err)
  
  DO icompartment = 1,Ncompartments
    CALL CMISSEquationsSet_DependentCreateStart(EquationsSetDiffusion(icompartment),DependentFieldUserNumber,&
      & DependentField,Err)
    CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetDiffusion(icompartment),Err)
  ENDDO

!   !Set the mesh component to be used by the field components.
! !   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
! !     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !   ENDDO
!   !Set the mesh component to be used by the field components.
! !   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
! !     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !   ENDDO
!   !Set the mesh component to be used by the field components.
! !   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
! !     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !     CALL CMISSField_ComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !   ENDDO
!   !Finish the equations set dependent field variables

  !-------------------------------------------------------------------------------------
  ! INITIALISE DEPENDENT FIELDS
!   !Initialise dependent field (concentration one components)
!   CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
!     & 1,1.0_CMISSDP,Err)

!   !Initialise dependent field (concentration two components)
!   CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U2_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
!     & 1,INITIAL_FIELD_DIFFUSION_TWO,Err)

!   !Initialise dependent field (concentration three components)
!   CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U3_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
!     & 1,INITIAL_FIELD_DIFFUSION_THREE,Err)


  !
  !================================================================================================================================
  !
  ALLOCATE(CouplingCoeffs(Ncompartments,Ncompartments))
  IF(Ncompartments==2)THEN
!     CouplingCoeffs(1,1)=1.0E-01_CMISSDP
! !     CouplingCoeffs(1,2)=-1.0E-04_CMISSDP
! !     CouplingCoeffs(2,1)=-1.0E-04_CMISSDP
!     CouplingCoeffs(1,2)=1.0E-01_CMISSDP
!     CouplingCoeffs(2,1)=1.0E-01_CMISSDP
!     CouplingCoeffs(2,2)=1.0E-01_CMISSDP
    CouplingCoeffs(1,1)=1.0E-01_CMISSDP
    CouplingCoeffs(1,2)=1.0E-01_CMISSDP
    CouplingCoeffs(2,1)=1.0E-01_CMISSDP
    CouplingCoeffs(2,2)=1.0E-01_CMISSDP
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
  !MATERIALS FIELDS - create the materials field
  !Auto-created field contains a U variable type to store the diffusion coefficient(s)
  !It also contains a V variable type to store the coupling coefficients 
  DO icompartment = 1,Ncompartments
    MaterialsFieldUserNumberDiffusion = 400+icompartment
    CALL CMISSField_Initialise(MaterialsFieldDiffusion(icompartment),Err)
    CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetDiffusion(icompartment),MaterialsFieldUserNumberDiffusion,&
         & MaterialsFieldDiffusion(icompartment),Err)
    CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetDiffusion(icompartment),Err)
  END DO
  !Initialise the coupling coefficients
  !Need to devise a neater way of specifying these components - e.g. specify only the upper diagonal components, and then automatically fill out the rest


  DO icompartment = 1, Ncompartments
    DO COMPONENT_NUMBER=1, Ncompartments
!       CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDiffusion(icompartment),CMISS_FIELD_V_VARIABLE_TYPE, &
!          & CMISS_FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
        CALL CMISSField_ParameterSetUpdateConstant(MaterialsFieldDiffusion(icompartment),CMISS_FIELD_V_VARIABLE_TYPE, &
          & CMISS_FIELD_VALUES_SET_TYPE,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
    END DO
  END DO

!   CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDiffusionOne,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
!     & MaterialsFieldUserNumberDiffusionOne,POROSITY_PARAM_MAT_PROPERTIES,Err)
!   CALL CMISSField_ComponentValuesInitialise(MaterialsFieldDiffusionTwo,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
!     & MaterialsFieldUserNumberMatPropertiesPermOverVis,PERM_OVER_VIS_PARAM_MAT_PROPERTIES,Err)

  !
  !================================================================================================================================
  !

  !SOURCE FIELDS

   !create the equations set source field variables for both equations sets 
  DO icompartment = 1,Ncompartments
    SourceFieldUserNumberDiffusion=700_CMISSIntg+icompartment
    CALL CMISSField_Initialise(SourceFieldDiffusion(icompartment),Err)
    CALL CMISSEquationsSet_SourceCreateStart(EquationsSetDiffusion(icompartment),SourceFieldUserNumberDiffusion,&
         & SourceFieldDiffusion(icompartment),Err)
    CALL CMISSEquationsSet_SourceCreateFinish(EquationsSetDiffusion(icompartment),Err)
  END DO 

  !Create the equations set analytic field variables
  DO icompartment = 1,Ncompartments
    AnalyticFieldUserNumberDiffusion=900_CMISSIntg+icompartment
    CALL CMISSField_Initialise(AnalyticFieldDiffusion(icompartment),Err)

    CALL CMISSEquationsSet_AnalyticCreateStart(EquationsSetDiffusion(icompartment), &
      & CMISS_EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_TWO_DIM,&
      & AnalyticFieldUserNumberDiffusion,AnalyticFieldDiffusion(icompartment),Err)
  
    !Finish the equations set analytic field variables
    CALL CMISSEquationsSet_AnalyticCreateFinish(EquationsSetDiffusion(icompartment),Err)
  END DO 

  !
  !================================================================================================================================
  !

  !EQUATIONS
  DO icompartment=1,Ncompartments
    CALL CMISSEquations_Initialise(EquationsDiffusion(icompartment),Err)
    CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetDiffusion(icompartment),EquationsDiffusion(icompartment),Err)
    CALL CMISSEquations_SparsityTypeSet(EquationsDiffusion(icompartment),CMISS_EQUATIONS_SPARSE_MATRICES,Err)
    CALL CMISSEquations_OutputTypeSet(EquationsDiffusion(icompartment),EQUATIONS_DIFFUSION_OUTPUT,Err)
    CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetDiffusion(icompartment),Err)
  ENDDO
  !
  !================================================================================================================================
  !
  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a coupled diffusion-diffusion problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_MULTI_PHYSICS_CLASS,CMISS_PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE, &
    & CMISS_PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoop_TimesSet(ControlLoop,LINEAR_SOLVER_DIFFUSION_START_TIME,&
    & LINEAR_SOLVER_DIFFUSION_STOP_TIME,LINEAR_SOLVER_DIFFUSION_TIME_INCREMENT,Err)
  !Set the output timing
  CALL CMISSControlLoop_TimeOutputSet(ControlLoop,LINEAR_SOLVER_DIFFUSION_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !
  !SOLVERS

  !Start the creation of the problem solvers
  CALL CMISSSolver_Initialise(SolverDiffusion,Err)
  CALL CMISSSolver_Initialise(LinearSolverDiffusion,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)

  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverDiffusionUserNumber,SolverDiffusion,Err)
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
  !SOLVER EQUATIONS
  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(SolverDiffusion,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsDiffusion,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  !
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverDiffusionUserNumber,SolverDiffusion,Err)
  CALL CMISSSolver_SolverEquationsGet(SolverDiffusion,SolverEquationsDiffusion,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsDiffusion,CMISS_SOLVER_SPARSE_MATRICES,Err)

  DO icompartment=1,Ncompartments
    CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsDiffusion,EquationsSetDiffusion(icompartment),&
         & EquationsSetIndex,Err)
  ENDDO

  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsDiffusion,BoundaryConditionsDiffusion,Err)
  CALL CMISSSolverEquations_BoundaryConditionsAnalytic(SolverEquationsDiffusion,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsDiffusion,Err)

  !
  !================================================================================================================================
  !
  !RUN SOLVERS

  !Turn off PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  CALL CMISSProblem_Solve(Problem,Err)
  !
  !================================================================================================================================
  !
  !OUTPUT

  Call CMISSAnalyticAnalysisOutput(DependentField,"MultiCompDiffusionAnalytics_2Comp_2D_x40_y40_L_T1",Err)


  EXPORT_FIELD_IO=.TRUE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
!     CALL CMISSFields_Initialise(Fields,Err)
!     CALL CMISSFields_Create(Region,Fields,Err)
!     CALL CMISSFields_NodesExport(Fields,"MonolithicMultiCompDiffusionTest","FORTRAN",Err)
!     CALL CMISSFields_ElementsExport(Fields,"MonolithicMultiCompDiffusionTest","FORTRAN",Err)
!     CALL CMISSFields_Finalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
    CALL FieldmlOutput_InitializeInfo( Region, Mesh, 3, outputDirectory, basename, fieldmlInfo, err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".geometric", region, mesh, GeometricField, & 
      & CMISS_FIELD_U_VARIABLE_TYPE, err )

!     CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".dependent", region, mesh, DependentField, err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".dependent_U", region, mesh, DependentField, &
      & CMISS_FIELD_U_VARIABLE_TYPE, err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".dependent_V", region, mesh, DependentField, &
      & CMISS_FIELD_V_VARIABLE_TYPE, err )
    
    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".equations_set_field_1", region, mesh, EquationsSetFieldDiffusion(1), &
      & CMISS_FIELD_U_VARIABLE_TYPE, err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".equations_set_field_2", region, mesh, EquationsSetFieldDiffusion(2), &
      & CMISS_FIELD_U_VARIABLE_TYPE, err )

!     CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".equations_set_field_3", region, mesh, EquationsSetFieldDiffusion(3), &
!       & CMISS_FIELD_U_VARIABLE_TYPE, err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".source_1", region,mesh,SourceFieldDiffusion(1), &
      & CMISS_FIELD_U_VARIABLE_TYPE,&
      & err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".source_2", region,mesh,SourceFieldDiffusion(2), &
      & CMISS_FIELD_U_VARIABLE_TYPE,&
      & err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".materials_1", region, mesh, MaterialsFieldDiffusion(1), &
      & CMISS_FIELD_U_VARIABLE_TYPE,err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".materials_2", region, mesh, MaterialsFieldDiffusion(2), &
      & CMISS_FIELD_U_VARIABLE_TYPE,err )

    !CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".analytic", region, mesh, AnalyticField, err )
    
    CALL FieldmlOutput_Write( fieldmlInfo, outputFilename, err )
    
    CALL FieldmlUtil_FinalizeInfo( fieldmlInfo )
  ENDIF

  !Finialise CMISS
!   CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM MONOLITHICSCHEMETESTFIELDMLEXAMPLE
