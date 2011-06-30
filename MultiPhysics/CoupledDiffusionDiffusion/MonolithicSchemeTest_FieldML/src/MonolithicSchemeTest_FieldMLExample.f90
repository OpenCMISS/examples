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
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: DependentField
  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: MaterialsFieldDiffusion
  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: SourceFieldDiffusion
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
  !CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)
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
  DOMAIN_X1 =  0.0_CMISSDP
  DOMAIN_X2 =  1.0_CMISSDP
  DOMAIN_Y1 =  0.0_CMISSDP
  DOMAIN_Y2 =  1.0_CMISSDP
  DOMAIN_Z1 =  0.0_CMISSDP
  DOMAIN_Z2 =  1.0_CMISSDP
  !Set geometric tolerance
  GEOMETRY_TOLERANCE = 1.0E-12_CMISSDP
  !Set initial values
  INITIAL_FIELD_DIFFUSION_ONE=1.0_CMISSDP
  INITIAL_FIELD_DIFFUSION_TWO=1.0_CMISSDP
  INITIAL_FIELD_DIFFUSION_THREE=1.0_CMISSDP
  !Set initial boundary conditions
  INLET_WALL_NODES_DIFFUSION_ONE_FLAG=.TRUE.
  IF(INLET_WALL_NODES_DIFFUSION_ONE_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_DIFFUSION_ONE=36
    ALLOCATE(INLET_WALL_NODES_DIFFUSION_ONE(NUMBER_OF_INLET_WALL_NODES_DIFFUSION_ONE))
    INLET_WALL_NODES_DIFFUSION_ONE=(/191,155,119,83,23,21,192,156,120,84,24,22,&
     & 198,162,126,90,36,35,204,168,132,96,48,47,210,174,138,102,60,59,216,180,144,108,72,71/)
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_DIFFUSION_ONE=1.0_CMISSDP
  ENDIF  !Set material parameters

  INLET_WALL_NODES_DIFFUSION_TWO_FLAG=.TRUE.
  IF(INLET_WALL_NODES_DIFFUSION_TWO_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_DIFFUSION_TWO=36
    ALLOCATE(INLET_WALL_NODES_DIFFUSION_TWO(NUMBER_OF_INLET_WALL_NODES_DIFFUSION_TWO))
    INLET_WALL_NODES_DIFFUSION_TWO=(/191,155,119,83,23,21,192,156,120,84,24,22,&
     & 198,162,126,90,36,35,204,168,132,96,48,47,210,174,138,102,60,59,216,180,144,108,72,71/)
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_DIFFUSION_TWO=1.0_CMISSDP
  ENDIF  !Set material parameters

  INLET_WALL_NODES_DIFFUSION_THREE_FLAG=.TRUE.
  IF(INLET_WALL_NODES_DIFFUSION_THREE_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_DIFFUSION_THREE=36
    ALLOCATE(INLET_WALL_NODES_DIFFUSION_THREE(NUMBER_OF_INLET_WALL_NODES_DIFFUSION_THREE))
    INLET_WALL_NODES_DIFFUSION_THREE=(/191,155,119,83,23,21,192,156,120,84,24,22,&
     & 198,162,126,90,36,35,204,168,132,96,48,47,210,174,138,102,60,59,216,180,144,108,72,71/)
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_DIFFUSION_THREE=1.0_CMISSDP
  ENDIF  !Set material parameters



  !Set material parameters
  !Set number of Gauss points (Mind that also material field may be interpolated)
  BASIS_XI_GAUSS_GEOMETRY=3 !4
  BASIS_XI_GAUSS_CONC_ONE=3 !4
  BASIS_XI_GAUSS_CONC_TWO=3 !4
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_DIFFUSION_OUTPUT_TYPE=CMISSSolverSolverMatrixOutput
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DIFFUSION_OUTPUT=CMISSEquationsMatrixOutput
  EQUATIONS_DIFFUSION_TWO_OUTPUT=CMISSEquationsNoOutput
  EQUATIONS_DIFFUSION_THREE_OUTPUT=CMISSEquationsNoOutput
  !Set time parameter
  LINEAR_SOLVER_DIFFUSION_START_TIME=0.0_CMISSDP
  LINEAR_SOLVER_DIFFUSION_STOP_TIME=2.1_CMISSDP 
  LINEAR_SOLVER_DIFFUSION_TIME_INCREMENT=0.1_CMISSDP
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

  DIAG_LEVEL_LIST(1)=1
  DIAG_LEVEL_LIST(2)=2
  DIAG_LEVEL_LIST(3)=3
  DIAG_LEVEL_LIST(4)=4
  DIAG_LEVEL_LIST(5)=5

!   DIAG_ROUTINE_LIST(1)="DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE"
!   DIAG_ROUTINE_LIST(2)="PROBLEM_SOLVER_EQUATIONS_SOLVE"
!   DIAG_ROUTINE_LIST(3)="SOLVER_SOLUTION_UPDATE"

  !CMISSAllDiagType/CMISSInDiagType/CMISSFromDiagType
!   CALL CMISSDiagnosticsSetOn(CMISSInDiagType,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISSAllTimingType/CMISSInTimingType/CMISSFromTimingType
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
!  IF(BASIS_XI_INTERPOLATION_CONC_ONE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisConcOne=BasisGeometry
!  ELSE
!     MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
!     !Initialise a new velocity basis
!     CALL CMISSBasisTypeInitialise(BasisConcOne,Err)
!     !Start the creation of a basis
!     CALL CMISSBasisCreateStart(BASIS_NUMBER_CONC_ONE,BasisConcOne,Err)
!     !Set the basis type (Lagrange/Simplex)
!     CALL CMISSBasisTypeSet(BasisConcOne,BASIS_TYPE,Err)
!     !Set the basis xi number
!     CALL CMISSBasisNumberOfXiSet(BasisConcOne,NUMBER_OF_DIMENSIONS,Err)
!     !Set the basis xi interpolation and number of Gauss points
!     IF(NUMBER_OF_DIMENSIONS==2) THEN
!       CALL CMISSBasisInterpolationXiSet(BasisConcOne,(/BASIS_XI_INTERPOLATION_CONC_ONE,BASIS_XI_INTERPOLATION_CONC_ONE/),Err)
!       CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisConcOne,(/BASIS_XI_GAUSS_CONC_ONE,BASIS_XI_GAUSS_CONC_ONE/),Err)
!     ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
!       CALL CMISSBasisInterpolationXiSet(BasisConcOne,(/BASIS_XI_INTERPOLATION_CONC_ONE,BASIS_XI_INTERPOLATION_CONC_ONE, & 
!         & BASIS_XI_INTERPOLATION_CONC_ONE/),Err)                         
!       CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisConcOne,(/BASIS_XI_GAUSS_CONC_ONE,BASIS_XI_GAUSS_CONC_ONE, & 
!         & BASIS_XI_GAUSS_CONC_ONE/),Err)
!     ENDIF
!     !Finish the creation of the basis
!     CALL CMISSBasisCreateFinish(BasisConcOne,Err)
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
!     CALL CMISSBasisTypeInitialise(BasisConcTwo,Err)
!     !Start the creation of a basis
!     CALL CMISSBasisCreateStart(BASIS_NUMBER_CONC_TWO,BasisConcTwo,Err)
!     !Set the basis type (Lagrange/Simplex)
!     CALL CMISSBasisTypeSet(BasisConcTwo,BASIS_TYPE,Err)
!     !Set the basis xi number
!     CALL CMISSBasisNumberOfXiSet(BasisConcTwo,NUMBER_OF_DIMENSIONS,Err)
!     !Set the basis xi interpolation and number of Gauss points
!     IF(NUMBER_OF_DIMENSIONS==2) THEN
!       CALL CMISSBasisInterpolationXiSet(BasisConcTwo,(/BASIS_XI_INTERPOLATION_CONC_TWO,BASIS_XI_INTERPOLATION_CONC_TWO/),Err)
!       CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisConcTwo,(/BASIS_XI_GAUSS_CONC_TWO,BASIS_XI_GAUSS_CONC_TWO/),Err)
!     ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
!       CALL CMISSBasisInterpolationXiSet(BasisConcTwo,(/BASIS_XI_INTERPOLATION_CONC_TWO,BASIS_XI_INTERPOLATION_CONC_TWO, & 
!         & BASIS_XI_INTERPOLATION_CONC_TWO/),Err)                         
!       CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisConcTwo,(/BASIS_XI_GAUSS_CONC_TWO,BASIS_XI_GAUSS_CONC_TWO, & 
!         & BASIS_XI_GAUSS_CONC_TWO/),Err)
!     ENDIF
!     !Finish the creation of the basis
!     CALL CMISSBasisCreateFinish(BasisConcTwo,Err)
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
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSNodesCreateStart(Region,TOTAL_NUMBER_OF_ALL_NODES,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL CMISSMeshNumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL CMISSMeshNumberOfComponentsSet(Mesh,MESH_NUMBER_OF_ALL_COMPONENTS,Err)
  !
  CALL CMISSMeshElementsTypeInitialise(MeshElementsGeometry,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsConcOne,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsConcTwo,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsConcThree,Err)
  MESH_COMPONENT_NUMBER_GEOMETRY=1
  MESH_COMPONENT_NUMBER_CONC_ONE=1
  MESH_COMPONENT_NUMBER_CONC_TWO=1
  MESH_COMPONENT_NUMBER_CONC_THREE=1
  !Specify spatial mesh component
  CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_GEOMETRY,BasisGeometry,MeshElementsGeometry,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL CMISSMeshElementsNodesSet(MeshElementsGeometry,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_GEOMETRY),Err)
  ENDDO
  CALL CMISSMeshElementsCreateFinish(MeshElementsGeometry,Err)
  !Specify concentration one mesh component
! !  IF(BASIS_XI_INTERPOLATION_CONC_ONE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
     MeshElementsConcOne=MeshElementsGeometry
! !   ELSE
! !     MESH_COMPONENT_NUMBER_CONC_ONE=MESH_COMPONENT_NUMBER_GEOMETRY+1
!      CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_ONE,BasisConcOne,MeshElementsConcOne,Err)
! !     DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
! !       CALL CMISSMeshElementsNodesSet(MeshElementsConcOne,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
! !         & 1:NUMBER_OF_ELEMENT_NODES_CONC_ONE),Err)
! !     ENDDO
 !    CALL CMISSMeshElementsCreateFinish(MeshElementsConcOne,Err)
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
!     CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_TWO,BasisConcTwo,MeshElementsConcTwo,Err)
! !    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
! !      CALL CMISSMeshElementsNodesSet(MeshElementsConcTwo,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
! !        & 1:NUMBER_OF_ELEMENT_NODES_CONC_TWO),Err)
! !    ENDDO
!     CALL CMISSMeshElementsCreateFinish(MeshElementsConcTwo,Err)
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
!     CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_CONC_THREE,BasisConcThree,MeshElementsConcThree,Err)
! !    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
! !      CALL CMISSMeshElementsNodesSet(MeshElementsConcTwo,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, & 
! !        & 1:NUMBER_OF_ELEMENT_NODES_CONC_TWO),Err)
! !    ENDDO
!     CALL CMISSMeshElementsCreateFinish(MeshElementsConcThree,Err)
! !  ENDIF

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
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfUserDomains,Err)
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
        & 1,CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSFieldParameterSetUpdateStart(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  !
  !================================================================================================================================
  !
  !EQUATIONS SETS - USING NEW ARGUMENTS TO ALLOW FOR MULTI-COMPARTMENT MODELS
  DO icompartment = 1,Ncompartments

    EquationsSetFieldUserNumberDiffusion = 100_CMISSIntg+icompartment
    EquationsSetUserNumberDiffusion = 200_CMISSIntg+icompartment
    CALL CMISSFieldTypeInitialise(EquationsSetFieldDiffusion(icompartment),Err)
    CALL CMISSEquationsSetTypeInitialise(EquationsSetDiffusion(icompartment),Err)
    CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberDiffusion,Region,GeometricField,&
         & CMISSEquationsSetClassicalFieldClass, &
         & CMISSEquationsSetDiffusionEquationType,CMISSEquationsSetMultiCompTransportDiffusionSubtype,&
         & EquationsSetFieldUserNumberDiffusion,EquationsSetFieldDiffusion(icompartment),EquationsSetDiffusion(icompartment),&
         & Err)
    !Finish creating the equations set
    CALL CMISSEquationsSetCreateFinish(EquationsSetDiffusion(icompartment),Err)
    CALL CMISSFieldParameterSetUpdateConstant(RegionUserNumber,EquationsSetFielduserNumberDiffusion,CMISSFieldUVariableType, &
      & CMISSFieldValuesSetType,1,icompartment,Err)
    CALL CMISSFieldParameterSetUpdateConstant(EquationsSetFieldDiffusion(icompartment),CMISSFieldUVariableType, &
      & CMISSFieldValuesSetType,2,Ncompartments,Err)
  END DO 

  !-------------------------------------------------------------------------------------
  ! DEPENDENT FIELD: Shared
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSFieldCreateStart(DependentFieldUserNumber,Region,DependentField,Err)
  CALL CMISSFieldTypeSet(DependentField,CMISSFieldGeneralType,Err)  
  CALL CMISSFieldMeshDecompositionSet(DependentField,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(DependentField,GeometricField,Err) 
  CALL CMISSFieldDependentTypeSet(DependentField,CMISSFieldDependentType,Err) 
  !Create 2N number of variables
  CALL CMISSFieldNumberOfVariablesSet(DependentField,2*Ncompartments,Err) 
  !create two variables for each compartment
  ALLOCATE(VariableTypes(2*Ncompartments))
  DO num_var=1,Ncompartments
     VariableTypes(2*num_var-1)=CMISSFieldUVariableType+(CMISSFieldNumberOfVariableSubtypes*(num_var-1))
     VariableTypes(2*num_var)=CMISSFieldDelUDelNVariableType+(CMISSFieldNumberOfVariableSubtypes*(num_var-1))
  ENDDO
  CALL CMISSFieldVariableTypesSet(DependentField,VariableTypes,Err) 
  !loop over the number of compartments
  DO icompartment=1,2*Ncompartments
    !set dimension type
    CALL CMISSFieldDimensionSet(DependentField,VariableTypes(icompartment), &
       & CMISSFieldScalarDimensionType,Err)
    CALL CMISSFieldNumberOfComponentsSet(DependentField,VariableTypes(icompartment),1,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,VariableTypes(icompartment),1, & 
       & MESH_COMPONENT_NUMBER_CONC_ONE,Err)
  ENDDO
  CALL CMISSFieldCreateFinish(DependentField,Err)
  
  DO icompartment = 1,Ncompartments
    CALL CMISSEquationsSetDependentCreateStart(EquationsSetDiffusion(icompartment),DependentFieldUserNumber,&
      & DependentField,Err)
    CALL CMISSEquationsSetDependentCreateFinish(EquationsSetDiffusion(icompartment),Err)
  ENDDO

!   !Set the mesh component to be used by the field components.
! !   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !   ENDDO
!   !Set the mesh component to be used by the field components.
! !   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !   ENDDO
!   !Set the mesh component to be used by the field components.
! !   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDiffusionTwo,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
! !       & MESH_COMPONENT_NUMBER_CONC_TWO,Err)
! !   ENDDO
!   !Finish the equations set dependent field variables

  !-------------------------------------------------------------------------------------
  ! INITIALISE DEPENDENT FIELDS
!   !Initialise dependent field (concentration one components)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,1.0_CMISSDP,Err)

!   !Initialise dependent field (concentration two components)
!   CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU2VariableType,CMISSFieldValuesSetType, & 
!     & 1,INITIAL_FIELD_DIFFUSION_TWO,Err)

!   !Initialise dependent field (concentration three components)
!   CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldU3VariableType,CMISSFieldValuesSetType, & 
!     & 1,INITIAL_FIELD_DIFFUSION_THREE,Err)


  !
  !================================================================================================================================
  !
  ALLOCATE(CouplingCoeffs(Ncompartments,Ncompartments))
  IF(Ncompartments==2)THEN
    CouplingCoeffs(1,1)=10.0E-01_CMISSDP
!     CouplingCoeffs(1,2)=-1.0E-04_CMISSDP
!     CouplingCoeffs(2,1)=-1.0E-04_CMISSDP
    CouplingCoeffs(1,2)=10.0E-01_CMISSDP
    CouplingCoeffs(2,1)=10.0E-01_CMISSDP
    CouplingCoeffs(2,2)=10.0E-01_CMISSDP
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
    CALL CMISSFieldTypeInitialise(MaterialsFieldDiffusion(icompartment),Err)
    CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetDiffusion(icompartment),MaterialsFieldUserNumberDiffusion,&
         & MaterialsFieldDiffusion(icompartment),Err)
    CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetDiffusion(icompartment),Err)
  END DO
  !Initialise the coupling coefficients
  !Need to devise a neater way of specifying these components - e.g. specify only the upper diagonal components, and then automatically fill out the rest


  DO icompartment = 1, Ncompartments
    DO COMPONENT_NUMBER=1, Ncompartments
!       CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDiffusion(icompartment),CMISSFieldVVariableType, &
!          & CMISSFieldValuesSetType,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
        CALL CMISSFieldParameterSetUpdateConstant(MaterialsFieldDiffusion(icompartment),CMISSFieldVVariableType, &
          & CMISSFieldValuesSetType,COMPONENT_NUMBER,CouplingCoeffs(icompartment,COMPONENT_NUMBER),Err)
    END DO
  END DO

!   CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDiffusionOne,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
!     & MaterialsFieldUserNumberDiffusionOne,POROSITY_PARAM_MAT_PROPERTIES,Err)
!   CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDiffusionTwo,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
!     & MaterialsFieldUserNumberMatPropertiesPermOverVis,PERM_OVER_VIS_PARAM_MAT_PROPERTIES,Err)

  !
  !================================================================================================================================
  !

  !SOURCE FIELDS

   !create the equations set source field variables for both equations sets 
  DO icompartment = 1,Ncompartments
    SourceFieldUserNumberDiffusion=700_CMISSIntg+icompartment
    CALL CMISSFieldTypeInitialise(SourceFieldDiffusion(icompartment),Err)
    CALL CMISSEquationsSetSourceCreateStart(EquationsSetDiffusion(icompartment),SourceFieldUserNumberDiffusion,&
         & SourceFieldDiffusion(icompartment),Err)
    CALL CMISSEquationsSetSourceCreateFinish(EquationsSetDiffusion(icompartment),Err)
  END DO 

  !
  !================================================================================================================================
  !

  !EQUATIONS
  DO icompartment=1,Ncompartments
    CALL CMISSEquationsTypeInitialise(EquationsDiffusion(icompartment),Err)
    CALL CMISSEquationsSetEquationsCreateStart(EquationsSetDiffusion(icompartment),EquationsDiffusion(icompartment),Err)
    CALL CMISSEquationsSparsityTypeSet(EquationsDiffusion(icompartment),CMISSEquationsSparseMatrices,Err)
    CALL CMISSEquationsOutputTypeSet(EquationsDiffusion(icompartment),EQUATIONS_DIFFUSION_OUTPUT,Err)
    CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetDiffusion(icompartment),Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a coupled diffusion-diffusion problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemMultiPhysicsClass,CMISSProblemMultiCompartmentTransportType, &
    & CMISSProblemStandardMultiCompartmentTransportSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,LINEAR_SOLVER_DIFFUSION_START_TIME,&
    & LINEAR_SOLVER_DIFFUSION_STOP_TIME,LINEAR_SOLVER_DIFFUSION_TIME_INCREMENT,Err)
  !Set the output timing
  CALL CMISSControlLoopTimeOutputSet(ControlLoop,LINEAR_SOLVER_DIFFUSION_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !
  !SOLVERS

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(SolverDiffusion,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverDiffusion,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)

  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverDiffusionUserNumber,SolverDiffusion,Err)
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
  !SOLVER EQUATIONS
  !Start the creation of the problem solver equations
  CALL CMISSSolverTypeInitialise(SolverDiffusion,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsDiffusion,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverDiffusionUserNumber,SolverDiffusion,Err)
  CALL CMISSSolverSolverEquationsGet(SolverDiffusion,SolverEquationsDiffusion,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsDiffusion,CMISSSolverEquationsSparseMatrices,Err)

  DO icompartment=1,Ncompartments
    CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsDiffusion,EquationsSetDiffusion(icompartment),&
         & EquationsSetIndex,Err)
  ENDDO

  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsDiffusion,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsDiffusion,BoundaryConditionsDiffusion,Err)
  DO icompartment=1,Ncompartments
!     IF(INLET_WALL_NODES_DIFFUSION_FLAG(icompartment)) THEN
!       DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION(icompartment)
!         NODE_NUMBER=INLET_WALL_NODES_DIFFUSION(icompartment,NODE_COUNTER)
!         CONDITION=CMISSBoundaryConditionFixed
!           VALUE=0.2_CMISSDP
!           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDiffusion,DependentField,CMISSFieldUVariableType,1,CMISSNoGlobalDerivative, & 
!             & NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
!       ENDDO
!     ENDIF
    IF(icompartment==1)THEN
      IF(INLET_WALL_NODES_DIFFUSION_ONE_FLAG) THEN
        DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION_ONE
          NODE_NUMBER=INLET_WALL_NODES_DIFFUSION_ONE(NODE_COUNTER)
          CONDITION=CMISSBoundaryConditionFixed
            VALUE=1.0_CMISSDP
            CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDiffusion,DependentField,CMISSFieldUVariableType, &
              & 1,CMISSNoGlobalDerivative, & 
              & NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDIF

    IF(icompartment==2)THEN
      IF(INLET_WALL_NODES_DIFFUSION_TWO_FLAG) THEN
        DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION_TWO
          NODE_NUMBER=INLET_WALL_NODES_DIFFUSION_TWO(NODE_COUNTER)
          CONDITION=CMISSBoundaryConditionFixed
            VALUE=0.0_CMISSDP
            CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDiffusion,DependentField,CMISSFieldDelVDelNVariableType, &
              & 1,CMISSNoGlobalDerivative, & 
              & NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDIF

!    IF(icompartment==3)THEN
!     IF(INLET_WALL_NODES_DIFFUSION_THREE_FLAG) THEN
!       DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION_THREE
!         NODE_NUMBER=INLET_WALL_NODES_DIFFUSION_THREE(NODE_COUNTER)
!         CONDITION=CMISSBoundaryConditionFixed
!           VALUE=0.0_CMISSDP
!           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDiffusion,DependentField,CMISSFieldU1VariableType,1,CMISSNoGlobalDerivative, & 
!             & NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
!       ENDDO
!     ENDIF
!    ENDIF

    IF(icompartment==4)THEN
      IF(INLET_WALL_NODES_DIFFUSION_THREE_FLAG) THEN
        DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION_THREE
          NODE_NUMBER=INLET_WALL_NODES_DIFFUSION_THREE(NODE_COUNTER)
          CONDITION=CMISSBoundaryConditionFixed
            VALUE=0.0_CMISSDP
            CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDiffusion,DependentField,CMISSFieldU2VariableType,1, &
              & CMISSNoGlobalDerivative,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDIF

    IF(icompartment==5)THEN
      IF(INLET_WALL_NODES_DIFFUSION_THREE_FLAG) THEN
        DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_DIFFUSION_THREE
          NODE_NUMBER=INLET_WALL_NODES_DIFFUSION_THREE(NODE_COUNTER)
          CONDITION=CMISSBoundaryConditionFixed
            VALUE=1.0_CMISSDP
            CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDiffusion,DependentField,CMISSFieldU3VariableType,1, &
              & CMISSNoGlobalDerivative,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONC_ONE,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDIF
  ENDDO
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsDiffusion,Err)

  !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn off PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  CALL CMISSProblemSolve(Problem,Err)
  !
  !================================================================================================================================
  !
  !OUTPUT

  EXPORT_FIELD_IO=.TRUE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
!     CALL CMISSFieldsTypeInitialise(Fields,Err)
!     CALL CMISSFieldsTypeCreate(Region,Fields,Err)
!     CALL CMISSFieldIONodesExport(Fields,"MonolithicMultiCompDiffusionTest","FORTRAN",Err)
!     CALL CMISSFieldIOElementsExport(Fields,"MonolithicMultiCompDiffusionTest","FORTRAN",Err)
!     CALL CMISSFieldsTypeFinalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
    CALL FieldmlOutput_InitializeInfo( Region, Mesh, 3, outputDirectory, basename, fieldmlInfo, err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".geometric", region, mesh, GeometricField, & 
      & CMISSFieldUVariableType, err )

!     CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".dependent", region, mesh, DependentField, err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".dependent_U", region, mesh, DependentField, &
      & CMISSFieldUVariableType, err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".dependent_V", region, mesh, DependentField, &
      & CMISSFieldVVariableType, err )

!     CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".dependent_U1", region, mesh, DependentField, &
!       & CMISSFieldVVariableType, err )
    
    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".equations_set_field_1", region, mesh, EquationsSetFieldDiffusion(1), &
      & CMISSFieldUVariableType, err )

    CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".equations_set_field_2", region, mesh, EquationsSetFieldDiffusion(2), &
      & CMISSFieldUVariableType, err )

!     CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".equations_set_field_3", region, mesh, EquationsSetFieldDiffusion(3), &
!       & CMISSFieldUVariableType, err )

    !CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".source", region, mesh, SourceField, err )

    !CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".materials", region, mesh, MaterialsField, err )

    !CALL FieldmlOutput_AddField( fieldmlInfo, baseName//".analytic", region, mesh, AnalyticField, err )
    
    CALL FieldmlOutput_Write( fieldmlInfo, outputFilename, err )
    
    CALL FieldmlUtil_FinalizeInfo( fieldmlInfo )
  ENDIF

  !Finialise CMISS
!   CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM MONOLITHICSCHEMETESTFIELDMLEXAMPLE
