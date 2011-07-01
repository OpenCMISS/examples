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

!> \example MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDrivenDarcy/src/IncompressibleElasticityDrivenDarcyExample.f90
!! Example program to solve coupled FiniteElasticityDarcy equations using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDrivenDarcy/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/MultiPhysics/Poroelasticity/FiniteElasticityDarcy/IncompressibleElasticityDrivenDarcy/build-intel'>Linux GNU Build</a>
!!
!<

! !
! !  This example considers a coupled Finite Elasticity Darcy problem
! !

!> Main program

PROGRAM INCOMPELASTDRIVENDARCYANALYTICDARCYEXAMPLE

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

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=1.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: QuadraticBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: CubicBasisUserNumber=3

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberDarcy=6
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberMatProperties=42
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberDarcy=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberMatProperties=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberDarcy=12
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberMatProperties=13
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberSolid=15
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberDarcy=22
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumberMatProperties=23
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumberDarcy=27


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

  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfVariables=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryNumberOfComponents=3

  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS

  INTEGER(CMISSIntg) :: MPI_IERROR

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
  TYPE(CMISSBasisType) :: BasisVelocity
  TYPE(CMISSBasisType) :: BasisPressure
  TYPE(CMISSBasisType) :: CubicBasis, QuadraticBasis, LinearBasis, Bases(2)
  !Nodes
  TYPE(CMISSNodesType) :: Nodes
  !Elements
  TYPE(CMISSMeshElementsType) :: MeshElementsGeometry
  TYPE(CMISSMeshElementsType) :: MeshElementsVelocity
  TYPE(CMISSMeshElementsType) :: MeshElementsPressure
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
  TYPE(CMISSFieldType) :: IndependentFieldSolid
  TYPE(CMISSFieldType) :: AnalyticFieldDarcy
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

  INTEGER(CMISSIntg), PARAMETER :: DisplacementMeshComponentNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureMeshComponentNumber=2

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
  NUMBER_GLOBAL_X_ELEMENTS=2
  NUMBER_GLOBAL_Y_ELEMENTS=2
  NUMBER_GLOBAL_Z_ELEMENTS=2

  IF(NUMBER_GLOBAL_Z_ELEMENTS==0)THEN
    NUMBER_OF_DIMENSIONS=2
  ELSE
    NUMBER_OF_DIMENSIONS=3
  ENDIF
  !PROBLEM CONTROL PANEL

  !Import cmHeart mesh information
  !CALL FLUID_MECHANICS_IO_READ_CMHEART(CM,Err)
  BASIS_XI_INTERPOLATION_SOLID=CMISSBasisLinearLagrangeInterpolation
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
!   PERM_OVER_VIS_PARAM_DARCY=0.1_CMISSDP
  POROSITY_PARAM_MAT_PROPERTIES=POROSITY_PARAM_DARCY
  PERM_OVER_VIS_PARAM_MAT_PROPERTIES=PERM_OVER_VIS_PARAM_DARCY
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

  !CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

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

  !CMISSAllDiagType/CMISSInDiagType/CMISSFromDiagType
  CALL CMISSDiagnosticsSetOn(CMISSInDiagType,DIAG_LEVEL_LIST,"Diagnostics",DIAG_ROUTINE_LIST,Err)

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
  !Define basis functions
  CALL CMISSBasisTypeInitialise(LinearBasis,Err)
  CALL CMISSBasisCreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(LinearBasis, &
    & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
  !CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(LinearBasis,.true.,Err)
  CALL CMISSBasisCreateFinish(LinearBasis,Err)

  CALL CMISSBasisTypeInitialise(QuadraticBasis,Err)
  CALL CMISSBasisCreateStart(QuadraticBasisUserNumber,QuadraticBasis,Err)
  CALL CMISSBasisInterpolationXiSet(QuadraticBasis,(/CMISSBasisQuadraticLagrangeInterpolation, &
    & CMISSBasisQuadraticLagrangeInterpolation,CMISSBasisQuadraticLagrangeInterpolation/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(QuadraticBasis, &
    & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
  !CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(QuadraticBasis,.true.,Err)
  CALL CMISSBasisCreateFinish(QuadraticBasis,Err)

  CALL CMISSBasisTypeInitialise(CubicBasis,Err)
  CALL CMISSBasisCreateStart(CubicBasisUserNumber,CubicBasis,Err)
  CALL CMISSBasisInterpolationXiSet(CubicBasis,(/CMISSBasisCubicLagrangeInterpolation, &
    & CMISSBasisCubicLagrangeInterpolation,CMISSBasisCubicLagrangeInterpolation/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(CubicBasis, &
    & (/CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme,CMISSBasisHighQuadratureScheme/),Err)
  !CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(CubicBasis,.true.,Err) !Enable 3D interpolation on faces
  CALL CMISSBasisCreateFinish(CubicBasis,Err)

  Bases(1)=CubicBasis
  Bases(2)=QuadraticBasis

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  !CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,BasisGeometry,Err)
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Bases,Err)
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/WIDTH,HEIGHT/),Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/),Err)
  ELSE
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/WIDTH,HEIGHT,LENGTH/),Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS/),Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

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
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the field type
  CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)
  !Set the decomposition to use
  CALL CMISSFieldNumberOfVariablesSet(GeometricField,FieldGeometryNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,3,Err)  
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)
  !Set the mesh component to be used by the field components.
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

!   !Update the geometric field parameters
!   DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
!     DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!       VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
!       CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
!         & CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
!     ENDDO
!   ENDDO
!   CALL CMISSFieldParameterSetUpdateStart(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
!   CALL CMISSFieldParameterSetUpdateFinish(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create a decomposition

!   !Create a field to put the geometry (defualt is geometry)
!   CALL CMISSFieldTypeInitialise(GeometricFieldSolid,Err)
!   CALL CMISSFieldCreateStart(FieldGeometrySolidUserNumber,Region,GeometricFieldSolid,Err)
!   CALL CMISSFieldMeshDecompositionSet(GeometricFieldSolid,Decomposition,Err)
!   CALL CMISSFieldTypeSet(GeometricFieldSolid,CMISSFieldGeometricType,Err)
!   CALL CMISSFieldNumberOfVariablesSet(GeometricFieldSolid,FieldGeometrySolidNumberOfVariables,Err)
!   CALL CMISSFieldNumberOfComponentsSet(GeometricFieldSolid,CMISSFieldUVariableType,FieldGeometrySolidNumberOfComponents,Err)
!   CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,1,SolidMeshComponenetNumber,Err)
!   CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,2,SolidMeshComponenetNumber,Err)
!   CALL CMISSFieldComponentMeshComponentSet(GeometricFieldSolid,CMISSFieldUVariableType,3,SolidMeshComponenetNumber,Err)
!   CALL CMISSFieldCreateFinish(GeometricFieldSolid,Err)
! 
! !---
!   !Update the geometric field parameters
!   DO NODE_NUMBER=1,NUMBER_OF_NODES_GEOMETRY
!     DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!       VALUE=CM%N(NODE_NUMBER,COMPONENT_NUMBER)
!       CALL CMISSFieldParameterSetUpdateNode(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
!         & CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
!     ENDDO
!   ENDDO
!   CALL CMISSFieldParameterSetUpdateStart(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
!   CALL CMISSFieldParameterSetUpdateFinish(GeometricFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
!---

  !Create a fibre field and attach it to the geometric field
  CALL CMISSFieldTypeInitialise(FibreFieldSolid,Err)
  CALL CMISSFieldCreateStart(FieldFibreSolidUserNumber,Region,FibreFieldSolid,Err)
  CALL CMISSFieldTypeSet(FibreFieldSolid,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreFieldSolid,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(FibreFieldSolid,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(FibreFieldSolid,FieldFibreSolidNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(FibreFieldSolid,CMISSFieldUVariableType,FieldFibreSolidNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,1,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,2,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(FibreFieldSolid,CMISSFieldUVariableType,3,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldCreateFinish(FibreFieldSolid,Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for ALE Darcy
  CALL CMISSFieldTypeInitialise(EquationsSetFieldDarcy,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSetDarcy,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberDarcy,Region,GeometricField,CMISSEquationsSetFluidMechanicsClass, &
    & CMISSEquationsSetDarcyEquationType,CMISSEquationsSetIncompressibleElasticityDrivenDarcySubtype,&
    & EquationsSetFieldUserNumberDarcy,EquationsSetFieldDarcy,EquationsSetDarcy,Err)
  !Set the equations set to be a ALE Darcy problem
!   CALL CMISSEquationsSetSpecificationSet(EquationsSetDarcy,CMISSEquationsSetFluidMechanicsClass, &
!     & CMISSEquationsSetDarcyEquationType,CMISSEquationsSetIncompressibleElasticityDrivenDarcySubtype,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSetDarcy,Err)

  !Create the equations set for deformation-dependent material properties
  CALL CMISSFieldTypeInitialise(EquationsSetFieldMatProperties,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSetMatProperties,Err)
!   CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberMatProperties,Region,GeometricField,CMISSEquationsSetFittingClass,&
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberMatProperties,Region,GeometricField,CMISSEquationsSetFittingClass,&
    & CMISSEquationsSetDataFittingEquationType,CMISSEquationsSetMatPropertiesInriaModelDataFittingSubtype,&
    & EquationsSetFieldUserNumberMatProperties,EquationsSetFieldMatProperties,EquationsSetMatProperties,Err)
  !Set the equations set to be a deformation-dependent material properties problem
!   CALL CMISSEquationsSetSpecificationSet(EquationsSetMatProperties,CMISSEquationsSetFittingClass, &
!     & CMISSEquationsSetDataFittingEquationType,CMISSEquationsSetMatPropertiesInriaModelDataFittingSubtype,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSetMatProperties,Err)


  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Create the equations_set
  CALL CMISSFieldTypeInitialise(EquationsSetFieldSolid,Err)
  CALL CMISSEquationsSetTypeInitialise(EquationsSetSolid,Err)
  CALL CMISSEquationsSetCreateStart(EquationSetSolidUserNumber,Region,FibreFieldSolid,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetIncompressibleElasticityDrivenDarcySubtype,&
    & EquationsSetFieldSolidUserNumber,EquationsSetFieldSolid,EquationsSetSolid,Err)
!   CALL CMISSEquationsSetSpecificationSet(EquationsSetSolid,CMISSEquationsSetElasticityClass, &
!     & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetIncompressibleElasticityDrivenDarcySubtype,Err)
  CALL CMISSEquationsSetCreateFinish(EquationsSetSolid,Err)

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
  CALL CMISSFieldGeometricFieldSet(MaterialFieldSolid,GeometricField,Err)
  CALL CMISSFieldNumberOfVariablesSet(MaterialFieldSolid,FieldMaterialSolidNumberOfVariables,Err)
  CALL CMISSFieldNumberOfComponentsSet(MaterialFieldSolid,CMISSFieldUVariableType,FieldMaterialSolidNumberOfComponents,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldSolid,CMISSFieldUVariableType,1,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldSolid,CMISSFieldUVariableType,2,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(MaterialFieldSolid,CMISSFieldUVariableType,3,DisplacementMeshComponentNumber,Err)
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
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetMatProperties,DependentFieldUserNumberMatProperties, &
    & DependentFieldMatProperties,Err)
  !Set the mesh component to be used by the field components.
  NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES = 2
  DO COMPONENT_NUMBER=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_MAT_PROPERTIES
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldMatProperties,CMISSFieldUVariableType,COMPONENT_NUMBER, &
      & DisplacementMeshComponentNumber,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldMatProperties,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, &
      & DisplacementMeshComponentNumber,Err)
  ENDDO
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetMatProperties,Err)

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
  CALL CMISSFieldGeometricFieldSet(DependentFieldSolid,GeometricField,Err)
  CALL CMISSFieldDependentTypeSet(DependentFieldSolid,CMISSFieldDependentType,Err)
  CALL CMISSFieldNumberOfVariablesSet(DependentFieldSolid,FieldDependentSolidNumberOfVariables,Err)
  CALL CMISSFieldVariableTypesSet(DependentFieldSolid,(/CMISSFieldUVariableType, &
    & CMISSFieldDelUDelNVariableType,CMISSFieldVVariableType,CMISSFieldDelVDelNVariableType/),Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldUVariableType,FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,FieldDependentSolidNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldVVariableType,FieldDependentFluidNumberOfComponents,Err)
  CALL CMISSFieldNumberOfComponentsSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,FieldDependentFluidNumberOfComponents,Err)
  !
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,1,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,2,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,3,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,4,CMISSFieldNodeBasedInterpolation,Err)
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldUVariableType,4,CMISSFieldElementBasedInterpolation,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,4,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldUVariableType,4,PressureMeshComponentNumber,Err)
  !
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,1,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,2,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,3,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,4, &
    & CMISSFieldNodeBasedInterpolation,Err)
!   CALL CMISSFieldComponentInterpolationSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,4, &
!     & CMISSFieldElementBasedInterpolation,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,4,SolidMeshComponenetNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelUDelNVariableType,4,PressureMeshComponentNumber, &
    & Err)

  !For this equation type, MESH_COMPONENT_NUMBER_PRESSURE is actually the mass increase component as the pressure is taken from the solid equations
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,1,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,2,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,3,DisplacementMeshComponentNumber,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,4,MESH_COMPONENT_NUMBER_PRESSURE,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldVVariableType,4,PressureMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,1,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,2,DisplacementMeshComponentNumber,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,3,DisplacementMeshComponentNumber,Err)
!   CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,4,MESH_COMPONENT_NUMBER_PRESSURE,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentFieldSolid,CMISSFieldDelVDelNVariableType,4,PressureMeshComponentNumber,Err)

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


  !Create the equations set dependent field variables for ALE Darcy
!   CALL CMISSFieldTypeInitialise(DependentFieldDarcy,Err)
!   CALL CMISSEquationsSetDependentCreateStart(EquationsSetDarcy,DependentFieldUserNumberDarcy, & ! ??? UserNumber ???
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetDarcy,FieldDependentSolidUserNumber, & ! ??? UserNumber ???
    & DependentFieldSolid,Err)
!   !Set the mesh component to be used by the field components.
!   DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,CMISSFieldUVariableType,COMPONENT_NUMBER, &
!       & MESH_COMPONENT_NUMBER_VELOCITY,Err)
!     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, &
!       & MESH_COMPONENT_NUMBER_VELOCITY,Err)
!   ENDDO
!   COMPONENT_NUMBER=NUMBER_OF_DIMENSIONS+1
!     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,CMISSFieldUVariableType,COMPONENT_NUMBER, &
!       & MESH_COMPONENT_NUMBER_PRESSURE,Err)
!     CALL CMISSFieldComponentMeshComponentSet(DependentFieldDarcy,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, &
!       & MESH_COMPONENT_NUMBER_PRESSURE,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetDarcy,Err)

  !Initialise dependent field (velocity components,pressure,mass increase)
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS+1
    CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldVVariableType,CMISSFieldValuesSetType, &
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
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetDarcy,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & MaterialsFieldUserNumberDarcyPorosity,POROSITY_PARAM_DARCY,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldDarcy,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & MaterialsFieldUserNumberDarcyPermOverVis,PERM_OVER_VIS_PARAM_DARCY,Err)
  !Create the equations set materials field variables for deformation-dependent material properties
  CALL CMISSFieldTypeInitialise(MaterialsFieldMatProperties,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetMatProperties,MaterialsFieldUserNumberMatProperties, &
    & MaterialsFieldMatProperties,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetMatProperties,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldMatProperties,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & MaterialsFieldUserNumberMatPropertiesPorosity,POROSITY_PARAM_MAT_PROPERTIES,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldMatProperties,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & MaterialsFieldUserNumberMatPropertiesPermOverVis,PERM_OVER_VIS_PARAM_MAT_PROPERTIES,Err)


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
      & DisplacementMeshComponentNumber,Err)
  ENDDO
  !Finish the equations set independent field variables
  CALL CMISSEquationsSetIndependentCreateFinish(EquationsSetSolid,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(EquationsDarcy,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetDarcy,EquationsDarcy,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(EquationsDarcy,CMISSEquationsSparseMatrices,Err)
!   !Set the equations lumping type
!   CALL CMISSEquationsLumpingTypeSet(EquationsDarcy,CMISSEquationsUnlumpedMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(EquationsDarcy,EQUATIONS_DARCY_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetDarcy,Err)

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(EquationsMatProperties,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetMatProperties,EquationsMatProperties,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(EquationsMatProperties,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(EquationsMatProperties,EQUATIONS_MAT_PROPERTIES_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetMatProperties,Err)

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
  CALL CMISSFieldTypeInitialise(AnalyticFieldDarcy,Err)
    CALL CMISSEquationsSetAnalyticCreateStart(EquationsSetDarcy,&
      & CMISSEquationsSetIncompElastDarcyAnalyticDarcy,&
      & AnalyticFieldUserNumberDarcy,AnalyticFieldDarcy,Err)
  
    !Finish the equations set analytic field variables
    CALL CMISSEquationsSetAnalyticCreateFinish(EquationsSetDarcy,Err)



  !--------------------------------------------------------------------------------------------------------------------------------
  ! Solid

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 1,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentFieldSolid,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,0.0_CMISSDP, &
    & Err)

  ! end Solid
  !--------------------------------------------------------------------------------------------------------------------------------

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
    & CMISSProblemQuasistaticElasticityTransientDarcySubtype,Err)
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
!   !Set the output type
!   CALL CMISSControlLoopOutputTypeSet(ControlLoop,CMISSControlLoopProgressOutput,Err)
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

  CALL CMISSSolverNewtonAbsoluteToleranceSet(SolverSolid,ABSOLUTE_TOLERANCE,Err)
  CALL CMISSSolverNewtonRelativeToleranceSet(SolverSolid,RELATIVE_TOLERANCE,Err)
  CALL CMISSSolverNewtonMaximumIterationsSet(SolverSolid,MAXIMUM_ITERATIONS,Err)

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
!   CALL CMISSSolverDynamicDynamicSet(DynamicSolverDarcy,.TRUE.,Err)
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
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsMatProperties,EquationsSetMatProperties,EquationsSetIndex,Err)
  !
  !Get the Darcy solver equations
  CALL CMISSProblemSolverGet(Problem,(/ControlLoopSubiterationNumber,ControlLoopFluidNumber,CMISSControlLoopNode/), &
    & SolverDarcyIndex,LinearSolverDarcy,Err)
  CALL CMISSSolverSolverEquationsGet(LinearSolverDarcy,SolverEquationsDarcy,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsDarcy,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsDarcy,EquationsSetDarcy,EquationsSetIndex,Err)
  !
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS
  !Start the creation of the equations set boundary conditions for Darcy
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsDarcy,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsDarcy,BoundaryConditionsDarcy,Err)

  CALL CMISSProblemSolverEquationsBoundaryConditionsAnalytic(SolverEquationsDarcy,Err)


  !Prescribe boundary conditions (absolute nodal parameters)
  !Solid is computed in absolute position, rather than displacement. Thus BCs for absolute position
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsSolid,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsSolid,BoundaryConditionsSolid,Err)

!   !--- BCs on normal velocity only
!   CONDITION = CMISSBoundaryConditionMovedWall
! 
!   IF( CM%D==2_CMISSIntg ) THEN !CM%D = number of dimensions, ie 2D
!     DO NODE_NUMBER=1_CMISSIntg,NUMBER_OF_NODES_GEOMETRY
!       COORD_X = CM%N(NODE_NUMBER,1_CMISSIntg)
!       COORD_Y = CM%N(NODE_NUMBER,2_CMISSIntg)
! 
!       IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) ) THEN
!         !x-velocity
!         VALUE = 1.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
!           & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
!         !x-velocity
!         VALUE = 1.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
!           & NODE_NUMBER,1_CMISSIntg,CONDITION,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
!         !y-velocity
!         VALUE = 2.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
!           & NODE_NUMBER,2_CMISSIntg,CONDITION,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
!         !y-velocity
!         VALUE = 2.0_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
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
! !         VALUE = 10.0_CMISSDP
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,,DependentFieldCMISSFieldVVariableType,CMISSNoGlobalDerivative, &
! ! !           & NODE_NUMBER,4_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err) !BC on pressure component
! ! !           & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionMovedWall,VALUE,Err) !inflow
! !           & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err) !time-dependent inflow
! 
! !         !x-position: S O L I D ( U Variable Type)
! ! !         VALUE = 1.0_CMISSDP * DOMAIN_X1
! !         VALUE = COORD_X
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
! !           & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! 
! !         EDGE_COUNT = 0_CMISSIntg
! !         IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) EDGE_COUNT = EDGE_COUNT + 1_CMISSIntg
! !         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) EDGE_COUNT = EDGE_COUNT + 1_CMISSIntg
! !         IF( (ABS(COORD_Z-DOMAIN_Z1) < GEOMETRY_TOLERANCE) ) EDGE_COUNT = EDGE_COUNT + 1_CMISSIntg
! !         IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) ) EDGE_COUNT = EDGE_COUNT + 1_CMISSIntg
! ! 
! !         IF(EDGE_COUNT == 2_CMISSIntg) THEN !it is a corner node
! !           VALUE = COORD_Y
! !           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
! !             & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! ! 
! !           VALUE = COORD_Z
! !           CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
! !             & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! !         END IF
!       END IF
!       !
!       IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
! !         !x-velocity: F L U I D
! !         VALUE = 10.0_CMISSDP
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
! ! !           & NODE_NUMBER,4_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err) !BC on pressure component
! ! !           & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionMovedWall,VALUE,Err) !impermeable wall, zero flux
! !           & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err) !impermeable wall, zero flux
! 
! !         !x-position: S O L I D
! !         VALUE = COORD_X
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,1, &
! !           & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! !         
! !         !Fix point 1
! !         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
! !           IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) ) THEN
! !             VALUE = COORD_Y
! !             CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
! !               & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! ! 
! !             VALUE = COORD_Z
! !             CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
! !               & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! !           END IF
! !         END IF
! !         !(Fix) point 2
! !         IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
! !           IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) ) THEN
! !             VALUE = COORD_Z
! !             CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
! !               & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! !           END IF
! !         END IF
! !         !(Fix) point 3
! !         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
! !           IF( (ABS(COORD_Z-DOMAIN_Z1) < GEOMETRY_TOLERANCE) ) THEN
! !             VALUE = COORD_Y
! !             CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
! !               & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! !           END IF
! !         END IF
! 
! 
!       END IF
!       !
!       IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
! !         !y-velocity: F L U I D
! !         VALUE = 0.0_CMISSDP
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
! !           & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionMovedWall,VALUE,Err)
! ! 
! !         !y-position: S O L I D
! !         VALUE = 1.0_CMISSDP * DOMAIN_Y1
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
! ! !           & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionMovedWallIncremented,VALUE,Err)
! !           & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
! !         !y-velocity: F L U I D
! !         VALUE = 0.0_CMISSDP
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
! !           & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionMovedWall,VALUE,Err)
! ! 
! ! !         !y-position: S O L I D
! ! !         VALUE = 1.1_CMISSDP * DOMAIN_Y2
! ! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,1, &
! ! !           & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Z-DOMAIN_Z1) < GEOMETRY_TOLERANCE) ) THEN
!         !z-velocity: F L U I D
!         !mass-correction: F L U I D
!         VALUE = 0.1_CMISSDP
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
! !           & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! !           & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionCorrectionMassIncrease,VALUE,Err)
! !           & NODE_NUMBER,4_CMISSIntg,CMISSBoundaryConditionCorrectionMassIncrease,VALUE,Err)
!           & NODE_NUMBER,4_CMISSIntg,CMISSBoundaryConditionFree,VALUE,Err)
! 
! !         !z-position: S O L I D
! !         VALUE = 1.0_CMISSDP * DOMAIN_Z1
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
! !           & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
!       END IF
!       !
!       IF( (ABS(COORD_Z-DOMAIN_Z2) < GEOMETRY_TOLERANCE) ) THEN
! !         !z-velocity: F L U I D
! !         VALUE = 0.0_CMISSDP
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsDarcy,DependentField,CMISSFieldVVariableType,CMISSNoGlobalDerivative, &
! !           & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! 
! !         VALUE = COORD_X
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,1, &
! !           & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! ! 
! !         VALUE = COORD_Y
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
! !           & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! ! 
! !         VALUE = COORD_Z
! !         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
! !           & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! 
!         !x-position: S O L I D
!         VALUE = COORD_Z
!         CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,1, &
!           & NODE_NUMBER,3_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
!         
!         !Fix point 1
!         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
!           IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
!             VALUE = COORD_Y
!             CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
!               & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
! 
!             VALUE = COORD_X
!             CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
!               & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
!           END IF
!         END IF
!         !(Fix) point 2
!         IF( (ABS(COORD_Y-DOMAIN_Y1) < GEOMETRY_TOLERANCE) ) THEN
!           IF( (ABS(COORD_X-DOMAIN_X2) < GEOMETRY_TOLERANCE) ) THEN
!             VALUE = COORD_X
!             CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
!               & NODE_NUMBER,1_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
!           END IF
!         END IF
!         !(Fix) point 3
!         IF( (ABS(COORD_Y-DOMAIN_Y2) < GEOMETRY_TOLERANCE) ) THEN
!           IF( (ABS(COORD_X-DOMAIN_X1) < GEOMETRY_TOLERANCE) ) THEN
!             VALUE = COORD_Y
!             CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsSolid,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
!               & NODE_NUMBER,2_CMISSIntg,CMISSBoundaryConditionFixed,VALUE,Err)
!           END IF
!         END IF
! 
! 
!       END IF
!     END DO
!   END IF

  !Finish the creation of the equations set boundary conditions for Darcy
  !CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsDarcy,Err)
  !Finish the creation of the equations set boundary conditions for the solid
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsSolid,Err)
  !
  !Start the creation of the equations set boundary conditions for deformation-dependent material properties
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsMatProperties,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsMatProperties,BoundaryConditionsMatProperties,Err)
  !(No boundary conditions requrired for deformation-dependent material properties)
  !Finish the creation of the equations set boundary conditions for deformation-dependent material properties
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsMatProperties,Err)


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

  EXPORT_FIELD_IO=.FALSE.
  IF(EXPORT_FIELD_IO) THEN
    WRITE(*,'(A)') "Exporting fields..."
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"FiniteElasticityDarcy","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"FiniteElasticityDarcy","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF


  !Finialise CMISS
!   CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM INCOMPELASTDRIVENDARCYANALYTICDARCYEXAMPLE
