!> \file
!> $Id$
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

!> \example MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDrivenDarcy/src/IncompressibleElasticityDrivenDarcyExample.f90
!! Example program to solve coupled FiniteElasticityDarcy equations using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDrivenDarcy/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDrivenDarcy/build-intel'>Linux GNU Build</a>
!!
!<

! !
! !  This example considers a coupled Finite Elasticity Multi-Compartment Darcy problem
! !

!> Main program

PROGRAM FINITEELASTICITYMULTICOMPDARCYEXAMPLE

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
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberDarcy=6
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberMatProperties=42
  INTEGER(CMISSIntg) :: MaterialsFieldUserNumberDarcy
  INTEGER(CMISSIntg) :: MaterialsFieldUserNumberMatProperties
  INTEGER(CMISSIntg) :: EquationsSetUserNumberDarcy
  INTEGER(CMISSIntg) :: EquationsSetUserNumberMatProperties
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberSolid=15
  INTEGER(CMISSIntg) :: EquationsSetFieldUserNumberDarcy
  INTEGER(CMISSIntg) :: EquationsSetFieldUserNumberMatProperties
  INTEGER(CMISSIntg) :: icompartment,Ncompartments,num_var,componentnum
  INTEGER(CMISSIntg), PARAMETER :: NumberOfDomains=1
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

  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS
  INTEGER(CMISSIntg) :: BASIS_TYPE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_GEOMETRY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SOLID
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS,MESH_NUMBER_OF_ALL_COMPONENTS
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_GEOMETRY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_GEOMETRY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_GEOMETRY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_PRESSURE
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ALL_NODES
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS
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
  TYPE(CMISSBasisType) :: BasisVelocity
  TYPE(CMISSBasisType) :: BasisPressure
  !Nodes
  TYPE(CMISSNodesType) :: Nodes
  !Elements
  TYPE(CMISSMeshElementsType) :: MeshElementsGeometry
  TYPE(CMISSMeshElementsType) :: MeshElementsVelocity
  TYPE(CMISSMeshElementsType) :: MeshElementsPressure
  !Meshes
  TYPE(CMISSMeshType) :: Mesh
  !Decompositions
  TYPE(CMISSDecompositionType) :: Decomposition
  !Fields
  TYPE(CMISSFieldsType) :: Fields
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: DependentFieldMatProperties
  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: MaterialsFieldDarcy
  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: MaterialsFieldMatProperties
  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: EquationsSetFieldDarcy
  TYPE(CMISSFieldType), ALLOCATABLE, DIMENSION(:) :: EquationsSetFieldMatProperties
  TYPE(CMISSFieldType) :: IndependentFieldSolid
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType), ALLOCATABLE, DIMENSION(:) :: BoundaryConditionsDarcy
  TYPE(CMISSBoundaryConditionsType), ALLOCATABLE, DIMENSION(:) :: BoundaryConditionsMatProperties
  !Equations sets
  TYPE(CMISSEquationsSetType), ALLOCATABLE, DIMENSION(:) :: EquationsSetDarcy
  TYPE(CMISSEquationsSetType), ALLOCATABLE, DIMENSION(:) :: EquationsSetMatProperties
  !Equations
  TYPE(CMISSEquationsType), ALLOCATABLE, DIMENSION(:) :: EquationsDarcy
  TYPE(CMISSEquationsType), ALLOCATABLE, DIMENSION(:) :: EquationsMatProperties
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: DynamicSolverDarcy
  TYPE(CMISSSolverType) :: LinearSolverDarcy
  TYPE(CMISSSolverType) :: LinearSolverMatProperties
  TYPE(CMISSSolverType) :: LinearSolverSolid
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsDarcy
  TYPE(CMISSSolverEquationsType) :: SolverEquationsMatProperties

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables

  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: Err
  !Array containing the field variable types that will be used (for ease of incorporating inside a loop)
  INTEGER(CMISSIntg), ALLOCATABLE, DIMENSION(:) :: VariableTypes
  REAL(CMISSDP), ALLOCATABLE, DIMENSION(:,:) :: CouplingCoeffs

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

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometrySolidUserNumber=51
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometrySolidNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometrySolidNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldFibreSolidUserNumber=52
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreSolidNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreSolidNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialSolidUserNumber=53
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialSolidNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialSolidNumberOfComponents=3

  INTEGER(CMISSIntg), PARAMETER :: FieldDependentSolidUserNumber=54
  INTEGER(CMISSIntg) :: FieldDependentSolidNumberOfVariables
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentSolidNumberOfComponents=4
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentFluidNumberOfComponents=4  !(u,v,w,m)

  INTEGER(CMISSIntg), PARAMETER :: EquationSetSolidUserNumber=55
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldSolidUserNumber=25
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
  TYPE(CMISSFieldType) :: GeometricFieldSolid,FibreFieldSolid,MaterialFieldSolid,DependentFieldSolid,EquationsSetFieldSolid
  TYPE(CMISSSolverType) :: SolverSolid
  TYPE(CMISSSolverEquationsType) :: SolverEquationsSolid
  TYPE(CMISSMeshElementsType) :: MeshElementsSolid

  !End - Program variables and types (finite elasticity part)

  INTEGER(CMISSIntg) :: dummy_counter

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

  !CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  !
  !================================================================================================================================
  !

  !PROBLEM CONTROL PANEL

  !Import cmHeart mesh information
  CALL FLUID_MECHANICS_IO_READ_CMHEART(CM,Err)
  BASIS_NUMBER_GEOMETRY=CM%ID_M
  BASIS_NUMBER_VELOCITY=CM%ID_V
  BASIS_NUMBER_PRESSURE=CM%ID_P
  NUMBER_OF_DIMENSIONS=CM%D
  BASIS_TYPE=CM%IT_T
  BASIS_XI_INTERPOLATION_GEOMETRY=CM%IT_M
  BASIS_XI_INTERPOLATION_VELOCITY=CM%IT_V
  BASIS_XI_INTERPOLATION_PRESSURE=CM%IT_P
  BASIS_XI_INTERPOLATION_SOLID=CMISSBasisLinearLagrangeInterpolation
  NUMBER_OF_NODES_GEOMETRY=CM%N_M
  NUMBER_OF_NODES_VELOCITY=CM%N_V
  NUMBER_OF_NODES_PRESSURE=CM%N_P
  TOTAL_NUMBER_OF_NODES=CM%N_T
  TOTAL_NUMBER_OF_ELEMENTS=CM%E_T
  NUMBER_OF_ELEMENT_NODES_GEOMETRY=CM%EN_M
  NUMBER_OF_ELEMENT_NODES_VELOCITY=CM%EN_V
  NUMBER_OF_ELEMENT_NODES_PRESSURE=CM%EN_P
  !Set domain dimensions
  DOMAIN_X1 = -5.0_CMISSDP
  DOMAIN_X2 =  5.0_CMISSDP
  DOMAIN_Y1 = -5.0_CMISSDP
  DOMAIN_Y2 =  5.0_CMISSDP
  DOMAIN_Z1 = -5.0_CMISSDP
  DOMAIN_Z2 =  5.0_CMISSDP
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
  PERM_OVER_VIS_PARAM_DARCY=1.0e-1_CMISSDP
  POROSITY_PARAM_MAT_PROPERTIES=POROSITY_PARAM_DARCY
  PERM_OVER_VIS_PARAM_MAT_PROPERTIES=PERM_OVER_VIS_PARAM_DARCY
  !Set number of Gauss points (Mind that also material field may be interpolated)
  BASIS_XI_GAUSS_GEOMETRY=3 !4
  BASIS_XI_GAUSS_VELOCITY=3 !4
  BASIS_XI_GAUSS_PRESSURE=3 !4
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_MAT_PROPERTIES_OUTPUT_TYPE=CMISSSolverProgressOutput
  DYNAMIC_SOLVER_DARCY_OUTPUT_TYPE=CMISSSolverProgressOutput
  LINEAR_SOLVER_DARCY_OUTPUT_TYPE=CMISSSolverSolverOutput
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_DARCY_OUTPUT=CMISSEquationsNoOutput
  EQUATIONS_MAT_PROPERTIES_OUTPUT=CMISSEquationsNoOutput

  !Set time parameter
  DYNAMIC_SOLVER_DARCY_START_TIME=0.0_CMISSDP
!   DYNAMIC_SOLVER_DARCY_STOP_TIME=0.03_CMISSDP
  DYNAMIC_SOLVER_DARCY_TIME_INCREMENT=1.0e-2_CMISSDP
  DYNAMIC_SOLVER_DARCY_STOP_TIME=1_CMISSIntg * DYNAMIC_SOLVER_DARCY_TIME_INCREMENT
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

!   DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_FINITE_ELEMENT_CALCULATE"
!   DIAG_ROUTINE_LIST(2)="DARCY_EQUATION_PRE_SOLVE_STORE_REFERENCE_DATA"
!   DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_PRE_SOLVE_GET_SOLID_DISPLACEMENT"
!   DIAG_ROUTINE_LIST(2)="DARCY_EQUATION_PRE_SOLVE_ALE_UPDATE_MESH"
!   DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_PRE_SOLVE_GET_SOLID_DISPLACEMENT"
!   DIAG_ROUTINE_LIST(1)="DARCY_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS"
  DIAG_ROUTINE_LIST(1)="FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE"
!   DIAG_ROUTINE_LIST(2)="FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR"
!   DIAG_ROUTINE_LIST(3)="EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION"
!   DIAG_ROUTINE_LIST(5)="DARCY_EQUATION_PRE_SOLVE_MAT_PROPERTIES"
!   DIAG_ROUTINE_LIST(6)="DATA_FITTING_FINITE_ELEMENT_CALCULATE"
!   DIAG_ROUTINE_LIST(7)="FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE"
!   DIAG_ROUTINE_LIST(8)="FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE"
!   DIAG_ROUTINE_LIST(1)="PROBLEM_SOLVER_EQUATIONS_SOLVE"
!   DIAG_ROUTINE_LIST(1)="SOLVER_NEWTON_SOLVE"
!   DIAG_ROUTINE_LIST(2)="SOLVER_NEWTON_LINESEARCH_SOLVE"
!   DIAG_ROUTINE_LIST(1)="SOLVER_SOLUTION_UPDATE"
!   DIAG_ROUTINE_LIST(1)="FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE"

  !CMISSAllDiagType/CMISSInDiagType/CMISSFromDiagType
!   CALL CMISSDiagnosticsSetOn(CMISSInDiagType,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

  !CMISSAllTimingType/CMISSInTimingType/CMISSFromTimingType
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)
  ALLOCATE (EquationsSetDarcy(Ncompartments))
  ALLOCATE (EquationsSetFieldDarcy(Ncompartments))
  ALLOCATE (MaterialsFieldDarcy(Ncompartments))
  ALLOCATE (MaterialsFieldMatProperties(Ncompartments))
  ALLOCATE (BoundaryConditionsDarcy(Ncompartments))
  ALLOCATE (BoundaryConditionsMatProperties(Ncompartments))
  ALLOCATE (EquationsDarcy(Ncompartments))
  ALLOCATE (EquationsSetFieldMatProperties(Ncompartments))
  ALLOCATE (EquationsSetMatProperties(Ncompartments))
  ALLOCATE (EquationsMatProperties(Ncompartments))
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
  !For a volume-coupled problem, solid and fluid are based in the same region

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
  !Start the creation of another basis: Velocity
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisVelocity=BasisGeometry
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new velocity basis
    CALL CMISSBasisTypeInitialise(BasisVelocity,Err)
    !Start the creation of a basis
    CALL CMISSBasisCreateStart(BASIS_NUMBER_VELOCITY,BasisVelocity,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasisTypeSet(BasisVelocity,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL CMISSBasisNumberOfXiSet(BasisVelocity,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL CMISSBasisInterpolationXiSet(BasisVelocity,(/BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY/),Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisVelocity,(/BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY/),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL CMISSBasisInterpolationXiSet(BasisVelocity,(/BASIS_XI_INTERPOLATION_VELOCITY,BASIS_XI_INTERPOLATION_VELOCITY, &
        & BASIS_XI_INTERPOLATION_VELOCITY/),Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisVelocity,(/BASIS_XI_GAUSS_VELOCITY,BASIS_XI_GAUSS_VELOCITY, &
        & BASIS_XI_GAUSS_VELOCITY/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(BasisVelocity,Err)
  ENDIF
  !
  !Start the creation of another basis: Pressure
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    BasisPressure=BasisGeometry
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_VELOCITY) THEN
    BasisPressure=BasisVelocity
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new pressure basis
    CALL CMISSBasisTypeInitialise(BasisPressure,Err)
    !Start the creation of a basis
    CALL CMISSBasisCreateStart(BASIS_NUMBER_PRESSURE,BasisPressure,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasisTypeSet(BasisPressure,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL CMISSBasisNumberOfXiSet(BasisPressure,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL CMISSBasisInterpolationXiSet(BasisPressure,(/BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE/),Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisPressure,(/BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE/),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL CMISSBasisInterpolationXiSet(BasisPressure,(/BASIS_XI_INTERPOLATION_PRESSURE,BASIS_XI_INTERPOLATION_PRESSURE, &
        & BASIS_XI_INTERPOLATION_PRESSURE/),Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisPressure,(/BASIS_XI_GAUSS_PRESSURE,BASIS_XI_GAUSS_PRESSURE, &
        & BASIS_XI_GAUSS_PRESSURE/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(BasisPressure,Err)
  ENDIF

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Start the creation of another basis: Solid (uses Lagrange Hermite !)

  BASIS_NUMBER_SOLID = BASIS_NUMBER_PRESSURE + 1

  CALL CMISSBasisTypeInitialise(BasisSolid,Err)
  CALL CMISSBasisCreateStart(BASIS_NUMBER_SOLID,BasisSolid,Err)

  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasisTypeSet(BasisSolid,BASIS_TYPE,Err)

  CALL CMISSBasisNumberOfXiSet(BasisSolid,NUMBER_OF_DIMENSIONS,Err)

!   CALL CMISSBasisInterpolationXiSet(BasisSolid,(/CMISSBasisLinearLagrangeInterpolation, &
!     & CMISSBasisLinearLagrangeInterpolation,CMISSBasisLinearLagrangeInterpolation/),Err)

  CALL CMISSBasisInterpolationXiSet(BasisSolid,(/BASIS_XI_INTERPOLATION_GEOMETRY,BASIS_XI_INTERPOLATION_GEOMETRY, &
    & BASIS_XI_INTERPOLATION_GEOMETRY/),Err)

!   CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisSolid, &
!     & (/CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme,CMISSBasisMidQuadratureScheme/),Err)

  CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisSolid,(/BASIS_XI_GAUSS_GEOMETRY,BASIS_XI_GAUSS_GEOMETRY, &
    & BASIS_XI_GAUSS_GEOMETRY/),Err)

  CALL CMISSBasisCreateFinish(BasisSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

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
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL CMISSMeshNumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL CMISSMeshNumberOfComponentsSet(Mesh,MESH_NUMBER_OF_ALL_COMPONENTS,Err)
  !
  CALL CMISSMeshElementsTypeInitialise(MeshElementsGeometry,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsVelocity,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsPressure,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsSolid,Err)
  MESH_COMPONENT_NUMBER_GEOMETRY=1
  MESH_COMPONENT_NUMBER_VELOCITY=1
  MESH_COMPONENT_NUMBER_PRESSURE=1
  SolidMeshComponenetNumber = 1
  SolidPressureMeshComponenetNumber = 1
  !Specify spatial mesh component
  CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_GEOMETRY,BasisGeometry,MeshElementsGeometry,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL CMISSMeshElementsNodesSet(MeshElementsGeometry,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_GEOMETRY),Err)
  ENDDO
  CALL CMISSMeshElementsCreateFinish(MeshElementsGeometry,Err)
  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    MeshElementsVelocity=MeshElementsGeometry
  ELSE
    MESH_COMPONENT_NUMBER_VELOCITY=MESH_COMPONENT_NUMBER_GEOMETRY+1
    CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_VELOCITY,BasisVelocity,MeshElementsVelocity,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL CMISSMeshElementsNodesSet(MeshElementsVelocity,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, &
        & 1:NUMBER_OF_ELEMENT_NODES_VELOCITY),Err)
    ENDDO
    CALL CMISSMeshElementsCreateFinish(MeshElementsVelocity,Err)
  ENDIF
  !Specify pressure mesh component
  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_GEOMETRY) THEN
    MeshElementsPressure=MeshElementsGeometry
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_GEOMETRY
  ELSE IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_VELOCITY) THEN
    MeshElementsPressure=MeshElementsVelocity
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_VELOCITY
  ELSE
    MESH_COMPONENT_NUMBER_PRESSURE=MESH_COMPONENT_NUMBER_VELOCITY+1
    CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_PRESSURE,BasisPressure,MeshElementsPressure,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL CMISSMeshElementsNodesSet(MeshElementsPressure,ELEMENT_NUMBER,CM%P(ELEMENT_NUMBER, &
        & 1:NUMBER_OF_ELEMENT_NODES_PRESSURE),Err)
    ENDDO
    CALL CMISSMeshElementsCreateFinish(MeshElementsPressure,Err)
  ENDIF

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid: Specify Solid Mesh Component

!   !Simply inherit ...
  MeshElementsSolid=MeshElementsGeometry
  SolidMeshComponenetNumber=MESH_COMPONENT_NUMBER_GEOMETRY
  SolidPressureMeshComponenetNumber=MESH_COMPONENT_NUMBER_PRESSURE


  !... or specify:
!   SolidMeshComponenetNumber=MESH_COMPONENT_NUMBER_PRESSURE+1
!
!   CALL CMISSMeshElementsCreateStart(Mesh,SolidMeshComponenetNumber,BasisSolid,MeshElementsSolid,Err)
!
!   DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
!     CALL CMISSMeshElementsNodesSet(MeshElementsSolid,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER, &
!       & 1:NUMBER_OF_ELEMENT_NODES_GEOMETRY),Err)
!   ENDDO
!
!   CALL CMISSMeshElementsCreateFinish(MeshElementsSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  !Finish the creation of the mesh
  CALL CMISSMeshCreateFinish(Mesh,Err)

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  !Create a decomposition:
  !All mesh components (associated with G.Projection / Darcy / solid) share the same decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the field type
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)
  !Set the decomposition to use

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

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create a decomposition

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSFieldTypeInitialise(GeometricFieldSolid,Err)
  CALL CMISSFieldCreateStart(FieldGeometrySolidUserNumber,Region,GeometricFieldSolid,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricFieldSolid,Decomposition,Err)
  CALL CMISSFieldTypeSet(GeometricFieldSolid,CMISSFieldGeometricType,Err)
  CALL CMISSFieldNumberOfVariablesSet(GeometricFieldSolid,FieldGeometrySolidNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricFieldSolid,CMISSFieldUVariableType,FieldGeometrySolidNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,1,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,2,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,3,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldCreateFinish(GeometricFieldSolid,Err)

!---
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
      VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
      CALL CMISSFieldParameterSetUpdateNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
        & CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSFieldParameterSetUpdateStart(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
!---

  !Create a fibre field and attach it to the geometric field
  CALL CMISSFieldTypeInitialise(FibreFieldSolid,Err)
  CALL CMISSFieldCreateStart(FieldFibreSolidUserNumber,Region,FibreFieldSolid,Err)
  CALL CMISSFieldTypeSet(FibreFieldSolid,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreFieldSolid,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(FibreFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSFieldNumberOfVariablesSet(FibreFieldSolid,FieldFibreSolidNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreFieldSolid,CMISSFieldUVariableType,FieldFibreSolidNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,1,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,2,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,3,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldCreateFinish(FibreFieldSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

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
    CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberDarcy,Region,GeometricField,CMISSEquationsSetFluidMechanicsClass, & 
      & CMISSEquationsSetDarcyEquationType,CMISSEquationsSetIncompressibleElastMultiCompDarcySubtype,&
      & EquationsSetFieldUserNumberDarcy,EquationsSetFieldDarcy(icompartment),EquationsSetDarcy(icompartment),Err)
    !Finish creating the equations set
    CALL CMISSEquationsSetCreateFinish(EquationsSetDarcy(icompartment),Err)
    CALL CMISSFieldParameterSetUpdateConstant(EquationsSetFieldDarcy(icompartment),CMISSFieldUVariableType, &
      & CMISSFieldValuesSetType,1,icompartment,Err)
    CALL CMISSFieldParameterSetUpdateConstant(EquationsSetFieldDarcy(icompartment),CMISSFieldUVariableType, &
      & CMISSFieldValuesSetType,2,Ncompartments,Err)
  ENDDO
  DO icompartment = 1,Ncompartments
    EquationsSetFieldUserNumberMatProperties = 150_CMISSIntg+icompartment
    EquationsSetUserNumberMatProperties = 250_CMISSIntg+icompartment
    !Create the equations set for deformation-dependent material properties
    CALL CMISSFieldTypeInitialise(EquationsSetFieldMatProperties(icompartment),Err)
    CALL CMISSEquationsSetTypeInitialise(EquationsSetMatProperties(icompartment),Err)
    CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberMatProperties,Region,GeometricField, & 
      & CMISSEquationsSetFittingClass,&
      & CMISSEquationsSetDataFittingEquationType,CMISSEquationsSetMatPropertiesInriaModelDataFittingSubtype,&
      & EquationsSetFieldUserNumberMatProperties,EquationsSetFieldMatProperties(icompartment), &
      & EquationsSetMatProperties(icompartment),Err)
    !Finish creating the equations set
    CALL CMISSEquationsSetCreateFinish(EquationsSetMatProperties(icompartment),Err)
!     CALL CMISSFieldParameterSetUpdateConstant(EquationsSetFieldMatProperties(icompartment),CMISSFieldUVariableType, &
!       & CMISSFieldValuesSetType,1,icompartment,Err)
!     CALL CMISSFieldParameterSetUpdateConstant(EquationsSetFieldMatProperties(icompartment),CMISSFieldUVariableType, &
!       & CMISSFieldValuesSetType,2,Ncompartments,Err)
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
    CALL CMISSFieldParameterSetUpdateConstant(EquationsSetFieldSolid,CMISSFieldUVariableType, &
      & CMISSFieldValuesSetType,1,1,Err)
    CALL CMISSFieldParameterSetUpdateConstant(EquationsSetFieldSolid,CMISSFieldUVariableType, &
      & CMISSFieldValuesSetType,2,Ncompartments,Err)


  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid Materials Field

  !Create a material field and attach it to the geometric field
  CALL CMISSFieldTypeInitialise(MaterialFieldSolid,Err)
  !
  CALL CMISSFieldCreateStart(FieldMaterialSolidUserNumber,Region,MaterialFieldSolid,Err)
  !
  CALL CMISSFieldTypeSet(MaterialFieldSolid,CMISSFieldMaterialType,Err)
  CALL CMISSFieldMeshDecompositionSet(MaterialFieldSolid,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(MaterialFieldSolid,GeometricFieldSolid,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialFieldSolid,FieldMaterialSolidNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialFieldSolid,CMISSFieldUVariableType,FieldMaterialSolidNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldSolid,CMISSFieldUVariableType,1,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldSolid,CMISSFieldUVariableType,2,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldSolid,CMISSFieldUVariableType,3,SolidMeshComponenetNumber,Err)
  !
  CALL CMISSFieldCreateFinish(MaterialFieldSolid,Err)

  !Set material parameters
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1.0_CMISSDP,Err)
!   CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,2.0e3_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,1.0_CMISSDP,Err)
!   CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,33.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,10.0_CMISSDP,Err)

  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetSolid,FieldMaterialSolidUserNumber,MaterialFieldSolid,Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------


  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for deformation-dependent material properties
  CALL CMISSFieldTypeInitialise(DependentFieldMatProperties,Err)

  DO icompartment=1,Ncompartments
    CALL CMISSEquationsSetDependentCreateStart(EquationsSetMatProperties(icompartment),DependentFieldUserNumberMatProperties, &
      & DependentFieldMatProperties,Err)
    !Set the mesh component to be used by the field components.
!NEED TO GENERALISE THE FOLLOWING COMMENTED OUT LINES OF CODE TO WORK FOR ARBITRARY VARIABLE TYPES
!     NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES = 2
!     DO COMPONENT_NUMBER=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES
!       CALL CMISSFieldComponentMeshComponentSet(DependentFieldMatProperties,CMISSFieldUVariableType,COMPONENT_NUMBER, &
!         & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
!       CALL CMISSFieldComponentMeshComponentSet(DependentFieldMatProperties,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, &
!         & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
!     ENDDO
    !Finish the equations set dependent field variables
    CALL CMISSEquationsSetDependentCreateFinish(EquationsSetMatProperties(icompartment),Err)
  ENDDO


  !Initialise dependent field
  DO COMPONENT_NUMBER=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES
    CALL CMISSFieldComponentValuesInitialise(DependentFieldMatProperties,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
      & COMPONENT_NUMBER,INITIAL_FIELD_MAT_PROPERTIES(COMPONENT_NUMBER),Err)
  ENDDO

  !--------------------------------------------------------------------------------------------------------------------------------
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
  CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,1,CMISSFieldNodeBasedInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,2,CMISSFieldNodeBasedInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,3,CMISSFieldNodeBasedInterpolation,Err)

  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,1,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,2,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,3,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,4,CMISSFieldNodeBasedInterpolation,Err)
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,4,CMISSFieldElementBasedInterpolation,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,4,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,4,SolidPressureMeshComponenetNumber,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,1, &
    & CMISSFieldNodeBasedInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,2, &
    & CMISSFieldNodeBasedInterpolation,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,3, &
    & CMISSFieldNodeBasedInterpolation,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,1,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,2,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,3,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,4, &
    & CMISSFieldNodeBasedInterpolation,Err)
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,4, &
!     & CMISSFieldElementBasedInterpolation,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,4,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,4,SolidPressureMeshComponenetNumber, &
    & Err)
  !loop over the number of compartments
  DO icompartment=3,2*Ncompartments+2
    CALL CMISSFieldDimensionSet(DependentFieldSolid,VariableTypes(icompartment), &
       & CMISSFieldVectorDimensionType,Err)
    CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,VariableTypes(icompartment),FieldDependentFluidNumberOfComponents,Err)
    DO componentnum=1,FieldDependentFluidNumberOfComponents-1
    !set dimension type
!     CALL CMISSFieldDimensionSet(DependentField,VariableTypes(icompartment), &
!        & CMISSFieldScalarDimensionType,Err)
      CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,VariableTypes(icompartment),componentnum, &
       & CMISSFieldNodeBasedInterpolation,Err)
      CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,VariableTypes(icompartment),componentnum, & 
         & MESH_COMPONENT_NUMBER_VELOCITY,Err)
    ENDDO
      CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,VariableTypes(icompartment), &
       & FieldDependentFluidNumberOfComponents, &
       & CMISSFieldNodeBasedInterpolation,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,VariableTypes(icompartment), &
      & FieldDependentFluidNumberOfComponents,MESH_COMPONENT_NUMBER_PRESSURE,Err)
!     CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,VariableTypes(icompartment), &
!       & FieldDependentFluidNumberOfComponents,MESH_COMPONENT_NUMBER_VELOCITY,Err)
    
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
  CALL CMISSFieldCreateFinish(DependentFieldSolid,Err)
  !
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetSolid,FieldDependentSolidUserNumber,DependentFieldSolid,Err)
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetSolid,Err)

!   !Initialise dependent field (solid displacement and pressure)
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS !+1
!     CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
!       & COMPONENT_NUMBER,INITIAL_FIELD_SOLID(COMPONENT_NUMBER),Err)
!   ENDDO

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  DO icompartment = 1,Ncompartments
    CALL CMISSEquationsSetDependentCreateStart(EquationsSetDarcy(icompartment),FieldDependentSolidUserNumber,&
      & DependentFieldSolid,Err)
    CALL CMISSEquationsSetDependentCreateFinish(EquationsSetDarcy(icompartment),Err)
  ENDDO


  !Create the equations set dependent field variables for ALE Darcy
!   CALL CMISSFieldTypeInitialise(DependentFieldDarcy,Err)
!   CALL CMISSEquationsSetDependentCreateStart(EquationsSetDarcy,DependentFieldUserNumberDarcy, & ! ??? UserNumber ???
!   CALL CMISSEquationsSetDependentCreateStart(EquationsSetDarcy,FieldDependentSolidUserNumber, & ! ??? UserNumber ???
!     & DependentFieldSolid,Err)
! !   !Set the mesh component to be used by the field components.
! !   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,CMISSFieldUVariableType,COMPONENT_NUMBER, &
! !       & MESH_COMPONENT_NUMBER_VELOCITY,Err)
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, &
! !       & MESH_COMPONENT_NUMBER_VELOCITY,Err)
! !   ENDDO
! !   COMPONENT_NUMBER=NUMBER_OF_DIMENSIONS+1
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,CMISSFieldUVariableType,COMPONENT_NUMBER, &
! !       & MESH_COMPONENT_NUMBER_PRESSURE,Err)
! !     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, &
! !       & MESH_COMPONENT_NUMBER_PRESSURE,Err)
!   !Finish the equations set dependent field variables
!   CALL CMISSEquationsSetDependentCreateFinish(EquationsSetDarcy,Err)

  !Initialise dependent field (velocity components,pressure,mass increase)
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS+1
    CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldVVariableType,CMISSFieldValuesSetType, &
      & COMPONENT_NUMBER,INITIAL_FIELD_DARCY(COMPONENT_NUMBER),Err)
  ENDDO
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

  !
  !================================================================================================================================
  !

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
  !Create the equations set materials field variables for ALE Darcy
!   CALL CMISSFieldTypeInitialise(MaterialsFieldDarcy,Err)
!   CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetDarcy,MaterialsFieldUserNumberDarcy, &
!     & MaterialsFieldDarcy,Err)
!   !Finish the equations set materials field variables
!   CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetDarcy,Err)
!   CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
!     & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
!   CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
!     & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)
  !Create the equations set materials field variables for deformation-dependent material properties
  DO icompartment = 1,Ncompartments
    MaterialsFieldUserNumberMatProperties = 500+icompartment
    CALL CMISSFieldTypeInitialise(MaterialsFieldMatProperties(icompartment),Err)
    CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetMatProperties(icompartment),MaterialsFieldUserNumberMatProperties, &
      & MaterialsFieldMatProperties(icompartment),Err)
    !Finish the equations set materials field variables
    CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetMatProperties(icompartment),Err)
  ENDDO
  DO icompartment = 1,Ncompartments
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldMatProperties(icompartment),CMISSFieldUVariableType, & 
    & CMISSFieldValuesSetType, &
    & MaterialsFieldUserNumberMatPropertiesPorosity,POROSITY_PARAM_MAT_PROPERTIES,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldMatProperties(icompartment),CMISSFieldUVariableType, &
    & CMISSFieldValuesSetType, &
    & MaterialsFieldUserNumberMatPropertiesPermOverVis,PERM_OVER_VIS_PARAM_MAT_PROPERTIES,Err)
  END DO

  !
  !================================================================================================================================
  !

  !INDEPENDENT FIELDS

  !Create the equations set independent field variables for the solid
  CALL CMISSFieldTypeInitialise(IndependentFieldSolid,Err)
  CALL CMISSEquationsSetIndependentCreateStart(EquationsSetSolid,IndependentFieldUserNumberSolid, &
    & IndependentFieldSolid,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
    CALL CMISSFieldComponentMeshComponentSet(IndependentFieldSolid,CMISSFieldUVariableType,COMPONENT_NUMBER, &
      & MESH_COMPONENT_NUMBER_GEOMETRY,Err)
  ENDDO
  !Finish the equations set independent field variables
  CALL CMISSEquationsSetIndependentCreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS

 
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
  write(*,*) "round and round we go!"
  ENDDO
  DO icompartment=1,Ncompartments
    !Create the equations set equations
    CALL CMISSEquationsTypeInitialise(EquationsMatProperties(icompartment),Err)
    CALL CMISSEquationsSetEquationsCreateStart(EquationsSetMatProperties(icompartment),EquationsMatProperties(icompartment),Err) 
    !Set the equations matrices sparsity type
    CALL CMISSEquationsSparsityTypeSet(EquationsMatProperties(icompartment),CMISSEquationsSparseMatrices,Err)
    !Set the equations set output
    CALL CMISSEquationsOutputTypeSet(EquationsMatProperties(icompartment),EQUATIONS_MAT_PROPERTIES_OUTPUT,Err)
    !Finish the equations set equations
    CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetMatProperties(icompartment),Err)
  ENDDO
  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(EquationsSolid,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetSolid,EquationsSolid,Err)
  CALL CMISSEquationsSparsityTypeSet(EquationsSolid,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(EquationsSolid,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  !
  !================================================================================================================================
  !

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 1,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,0.0_CMISSDP, &
    & Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  !BOUNDARY CONDITIONS
  !Start the creation of the equations set boundary conditions for Darcy
  DO icompartment=1,Ncompartments  
    CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsDarcy(icompartment),Err)
    CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSetDarcy(icompartment),BoundaryConditionsDarcy(icompartment),Err)
    CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSetDarcy(icompartment),Err)
  ENDDO




  !Prescribe boundary conditions (absolute nodal parameters)
  !Solid is computed in absolute position, rather than displacement. Thus BCs for absolute position
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsSolid,Err)
  CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSetSolid,BoundaryConditionsSolid,Err)

  !--- BCs on normal velocity only
  CONDITION = CMISSBoundaryConditionMovedWall

  dummy_counter = 0_CMISSIntg
!   IF( CM%D==2_CMISSIntg ) THEN !CM%D = number of dimensions, ie 2D
!     DO NODE_NUMBER=1_CMISSIntg,NUMBER_OF_NODES_GEOMETRY
!       COORD_X = CM%N(NODE_NUMBER,1_CMISSIntg)
!       COORD_Y = CM%N(NODE_NUMBER,2_CMISSIntg)
! 
!       IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) ) THEN
!         !x-velocity
!         VALUE = 1.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
!           & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
!         !x-velocity
!         VALUE = 1.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
!           & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
!         !y-velocity
!         VALUE = 2.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
!           & NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
!         !y-velocity
!         VALUE = 2.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
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
!         dummy_counter = dummy_counter + 1_CMISSIntg
!         !x-velocity: F L U I D ( V Variable type )
!         VALUE = 11.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
! !           & NODE_NUMBER,4_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err) !BC on pressure component
! !           & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionMovedWall,VALUE,Err) !inflow
!           & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err) !time-dependent inflow
! 
!         !x-position: S O L I D ( U Variable Type)
!         VALUE = 1.0_CMISSDP * DOMAIN_X1
!         ! CMISSNoGlobalDerivative is equal to 1
!         !CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,CMISSFieldUVariableType,1, &
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
!           & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
!         dummy_counter = dummy_counter + 1_CMISSIntg
!         !x-velocity: F L U I D
!         VALUE = 10.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
! !           & NODE_NUMBER,4_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err) !BC on pressure component
! !           & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionMovedWall,VALUE,Err) !impermeable wall, zero flux
!           & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err) !impermeable wall, zero flux
! 
!         !x-position: S O L I D
! !         VALUE = 1.0_CMISSDP * DOMAIN_X2
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,CMISSFieldUVariableType,1, &
! !           & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
!         dummy_counter = dummy_counter + 1_CMISSIntg
!         !y-velocity: F L U I D
!         VALUE = 0.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
!           & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionMovedWall,VALUE,Err)
! 
!         !y-position: S O L I D
! !         VALUE = 1.1_CMISSDP * DOMAIN_Y1
!         VALUE = 1.0_CMISSDP * DOMAIN_Y1
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
! !           & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionMovedWallIncremented,VALUE,Err)
!           & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
!         dummy_counter = dummy_counter + 1_CMISSIntg
!         !y-velocity: F L U I D
!         VALUE = 0.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
!           & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionMovedWall,VALUE,Err)
! 
! !         !y-position: S O L I D
! !         VALUE = 1.1_CMISSDP * DOMAIN_Y2
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,CMISSFieldUVariableType,1, &
! !           & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Z-DOMAIN_Z1) < GEOMETRY_TOLERANCE) ) THEN
!         dummy_counter = dummy_counter + 1_CMISSIntg
!         !z-velocity: F L U I D
!         VALUE = 0.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
!           & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionMovedWall,VALUE,Err)
! 
!         !z-position: S O L I D
!         VALUE = 1.0_CMISSDP * DOMAIN_Z1
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
!           & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) ) THEN
!         dummy_counter = dummy_counter + 1_CMISSIntg
!         !z-velocity: F L U I D
!         VALUE = 0.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
!           & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionMovedWall,VALUE,Err)
! 
! !         !z-position: S O L I D
! !         VALUE = 1.0_CMISSDP * DOMAIN_Z2
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,CMISSFieldUVariableType,1, &
! !           & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
!       END IF
!     END DO
!   END IF

  !Finish the creation of the equations set boundary conditions for Darcy

  !Finish the creation of the equations set boundary conditions for the solid
  CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSetSolid,Err)
  !
  DO icompartment=1,Ncompartments 
   !Start the creation of the equations set boundary conditions for deformation-dependent material properties
    CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsMatProperties(icompartment),Err)
    CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSetMatProperties(icompartment), &
      &  BoundaryConditionsMatProperties(icompartment),Err)
    !(No boundary conditions requrired for deformation-dependent material properties)
    !Finish the creation of the equations set boundary conditions for deformation-dependent material properties
    CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSetMatProperties(icompartment),Err)
  ENDDO

  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a ALE Darcy problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemMultiPhysicsClass,CMISSProblemFiniteElasticityDarcyType, &
    & CMISSProblemQuasistaticElasticityTransientDarcySubtype,Err) !DO WE NEED TO DEFINE A NEW SUBTYPE HERE FOR THE MULTICOMPARTMENT VERSION?
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
!   CALL CMISSControlLoopMaximumIterationsSet(ControlLoop,2,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,DYNAMIC_SOLVER_DARCY_START_TIME,DYNAMIC_SOLVER_DARCY_STOP_TIME, &
    & DYNAMIC_SOLVER_DARCY_TIME_INCREMENT,Err)
  !Set the output timing
  CALL CMISSControlLoopTimeOutputSet(ControlLoop,DYNAMIC_SOLVER_DARCY_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)


  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(SolverSolid,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverMatProperties,Err)
  CALL CMISSSolverTypeInitialise(DynamicSolverDarcy,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverDarcy,Err)

  CALL CMISSProblemSolversCreateStart(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Get the finite elasticity solver
  CALL CMISSProblemSolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopSolidNumber,CMISSControlLoopNode/), &
    & SolverSolidIndex,SolverSolid,Err)
  CALL CMISSSolverOutputTypeSet(SolverSolid,CMISSSolverProgressOutput,Err)
!   CALL CMISSSolverNewtonJacobianCalculationTypeSet(SolverSolid,CMISSSolverNewtonJacobianFDCalculated,Err)
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(SolverSolid,CMISSSolverNewtonJacobianAnalyticCalculated,Err)

!   CALL CMISSSolverNonLinearTypeSet(SolverSolid,CMISSSolverNonlinearNewton,Err)
!   CALL CMISSSolverLibraryTypeSet(SolverSolid,CMISSSolverPETScLibrary,Err)

!   CALL CMISSSolverNewtonLinearSolverGet(SolverSolid,LinearSolverSolid,Err)
!   CALL CMISSSolverLinearTypeSet(LinearSolverSolid,CMISSSolverLinearDirectSolveType,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  !Get the deformation-dependent material properties solver
  CALL CMISSProblemSolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMISSControlLoopNode/), &
    & SolverMatPropertiesIndex,LinearSolverMatProperties,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(LinearSolverMatProperties,LINEAR_SOLVER_MAT_PROPERTIES_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_MAT_PROPERTIES_DIRECT_FLAG) THEN
    CALL CMISSSolverLinearTypeSet(LinearSolverMatProperties,CMISSSolverLinearDirectSolveType,Err)
    CALL CMISSSolverLibraryTypeSet(LinearSolverMatProperties,CMISSSolverMUMPSLibrary,Err)
  ELSE
    CALL CMISSSolverLinearTypeSet(LinearSolverMatProperties,CMISSSolverLinearIterativeSolveType,Err)
    CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolverMatProperties,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolverLinearIterativeDivergenceToleranceSet(LinearSolverMatProperties,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeRelativeToleranceSet(LinearSolverMatProperties,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(LinearSolverMatProperties,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeGMRESRestartSet(LinearSolverMatProperties,RESTART_VALUE,Err)
  ENDIF

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
  CALL CMISSSolverTypeInitialise(LinearSolverMatProperties,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverDarcy,Err)

  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsSolid,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsMatProperties,Err)
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
  !Get the deformation-dependent material properties solver equations
  CALL CMISSProblemSolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMISSControlLoopNode/), &
    & SolverMatPropertiesIndex,LinearSolverMatProperties,Err)
  CALL CMISSSolverSolverEquationsGet(LinearSolverMatProperties,SolverEquationsMatProperties,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsMatProperties,CMISSSolverEquationsSparseMatrices,Err)
  DO icompartment=1,Ncompartments
    CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsMatProperties,EquationsSetMatProperties(icompartment), &
       & EquationsSetIndex,Err)
  ENDDO
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

  WRITE(*,*)'dummy_counter = ',dummy_counter

  EXPORT_FIELD_IO=.FALSE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"FiniteElasticityMultiCompDarcy","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"FiniteElasticityMultiCompDarcy","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF


  !Finialise CMISS
!   CALL CMISSFinalise(Err)

  IF (ALLOCATED(EquationsSetFieldDarcy)) DEALLOCATE(EquationsSetFieldDarcy)
  IF (ALLOCATED(EquationsSetDarcy)) DEALLOCATE(EquationsSetDarcy)
  IF (ALLOCATED(MaterialsFieldDarcy)) DEALLOCATE(MaterialsFieldDarcy)
  IF (ALLOCATED(MaterialsFieldMatProperties)) DEALLOCATE(MaterialsFieldMatProperties)
  IF (ALLOCATED(BoundaryConditionsDarcy)) DEALLOCATE(BoundaryConditionsDarcy)
  IF (ALLOCATED(BoundaryConditionsMatProperties)) DEALLOCATE(BoundaryConditionsMatProperties)
  IF (ALLOCATED(EquationsDarcy)) DEALLOCATE(EquationsDarcy)
  IF (ALLOCATED(EquationsSetFieldMatProperties)) DEALLOCATE(EquationsSetFieldMatProperties)
  IF (ALLOCATED(EquationsSetMatProperties)) DEALLOCATE(EquationsSetMatProperties)
  IF (ALLOCATED(EquationsMatProperties)) DEALLOCATE(EquationsMatProperties)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM FINITEELASTICITYMULTICOMPDARCYEXAMPLE
