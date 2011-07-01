!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve the ALE formulation of the diffusion equation using openCMISS calls.
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
!> The Original Code is openCMISS
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

!> \example ClassicalField/Diffusion/src/DiffusionIO_ALEExample.f90
!! Example program to solve a diffusion equation using openCMISS calls.
!! \htmlinclude ClassicalField/Diffusion/history.html
!<

!> Main program
PROGRAM DIFFUSIONIOALEEXAMPLE


  USE OPENCMISS
  USE FLUID_MECHANICS_IO_ROUTINES
  USE MPI


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=1337
  TYPE(CMISSFieldType) :: EquationsSetField


  !Test program parameters
  
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopNode=0
  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1

  !Program types

  TYPE(EXPORT_CONTAINER):: CM  

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS

  INTEGER(CMISSIntg) :: BASIS_TYPE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_SPACE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_CONCENTRATION
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_CONCENTRATION
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_CONCENTRATION
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_SPACE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_CONCENTRATION
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_CONCENTRATION
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_CONCENTRATION
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE
!  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NUMBER_OF_FIXED_WALL_NODES_ADVECTION_DIFFUSION
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_WALL_NODES_ADVECTION_DIFFUSION

  INTEGER(CMISSIntg) :: EQUATIONS_ADVECTION_DIFFUSION_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_COUNTER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_TYPE

  INTEGER, ALLOCATABLE, DIMENSION(:):: FIXED_WALL_NODES_ADVECTION_DIFFUSION
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_ADVECTION_DIFFUSION

  REAL(CMISSDP) :: INITIAL_FIELD_ADVECTION_DIFFUSION(3)
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_ADVECTION_DIFFUSION(3)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE
  REAL(CMISSDP) :: DIFF_COEFF_PARAM_ADVECTION_DIFFUSION(3)

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LINEAR_SOLVER_ADVECTION_DIFFUSION_DIRECT_FLAG
  LOGICAL :: FIXED_WALL_NODES_ADVECTION_DIFFUSION_FLAG
  LOGICAL :: INLET_WALL_NODES_ADVECTION_DIFFUSION_FLAG

  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  
  INTEGER(CMISSIntg) :: MPI_IERROR

  !INTEGER(INTG) :: first_global_dof,first_local_dof,first_local_rank,last_global_dof,last_local_dof,last_local_rank,rank_idx
  !INTEGER(INTG) :: EQUATIONS_SET_INDEX
  !TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOF_MAPPING
  
    !CMISS variables

  !Regions
  TYPE(CMISSRegionType) :: Region
  TYPE(CMISSRegionType) :: WorldRegion
  !Coordinate systems
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  !Basis
  TYPE(CMISSBasisType) :: BasisSpace
  TYPE(CMISSBasisType) :: BasisConcentration
  !Nodes
  TYPE(CMISSNodesType) :: Nodes
  !Elements
  TYPE(CMISSMeshElementsType) :: MeshElementsSpace
  TYPE(CMISSMeshElementsType) :: MeshElementsConcentration
  !Meshes
  TYPE(CMISSMeshType) :: Mesh
  !Decompositions
  TYPE(CMISSDecompositionType) :: Decomposition
  !Fields
  TYPE(CMISSFieldsType) :: Fields
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: DependentField
  TYPE(CMISSFieldType) :: MaterialsField
  TYPE(CMISSFieldType) :: IndependentField
  TYPE(CMISSFieldType) :: SourceField
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  !Equations sets
  TYPE(CMISSEquationsSetType) :: EquationsSet
  !Equations
  TYPE(CMISSEquationsType) :: Equations  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: Solver, LinearSolver
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquations
!   TYPE(CMISSBasisType) :: Basis
!   TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
!   TYPE(CMISSDecompositionType) :: Decomposition
!   TYPE(CMISSEquationsType) :: Equations
!   TYPE(CMISSEquationsSetType) :: EquationsSet
!   TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialsField,AnalyticField
!   TYPE(CMISSFieldsType) :: Fields
!   TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
!   TYPE(CMISSMeshType) :: Mesh
!   TYPE(CMISSProblemType) :: Problem
!   TYPE(CMISSControlLoopType) :: ControlLoop
!   TYPE(CMISSRegionType) :: Region,WorldRegion
!   TYPE(CMISSSolverType) :: Solver, LinearSolver
!   TYPE(CMISSSolverEquationsType) :: SolverEquations

  LOGICAL :: EXPORT_FIELD

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
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

  !
  !================================================================================================================================
  !

  !PROBLEM CONTROL PANEL

  !Import cmHeart mesh information
  CALL FLUID_MECHANICS_IO_READ_CMHEART(CM,Err)  
  BASIS_NUMBER_SPACE=CM%ID_M
  BASIS_NUMBER_CONCENTRATION=CM%ID_V
  NUMBER_OF_DIMENSIONS=CM%D
  BASIS_TYPE=CM%IT_T
!  BASIS_XI_INTERPOLATION_SPACE=CM%IT_M
  BASIS_XI_INTERPOLATION_SPACE=CMISSBasisLinearLagrangeInterpolation

  BASIS_XI_INTERPOLATION_CONCENTRATION=CM%IT_V
  NUMBER_OF_NODES_SPACE=CM%N_M
  NUMBER_OF_NODES_CONCENTRATION=CM%N_V
  TOTAL_NUMBER_OF_NODES=CM%N_T
  TOTAL_NUMBER_OF_ELEMENTS=CM%E_T
  NUMBER_OF_ELEMENT_NODES_SPACE=CM%EN_M
  NUMBER_OF_ELEMENT_NODES_CONCENTRATION=CM%EN_V
  !Set initial values
  INITIAL_FIELD_ADVECTION_DIFFUSION(1)=0.0_CMISSDP
  INITIAL_FIELD_ADVECTION_DIFFUSION(2)=0.0_CMISSDP
  INITIAL_FIELD_ADVECTION_DIFFUSION(3)=0.0_CMISSDP
  !Set boundary conditions - what condition should be applied?
  !Set boundary conditions
  FIXED_WALL_NODES_ADVECTION_DIFFUSION_FLAG=.FALSE.
  INLET_WALL_NODES_ADVECTION_DIFFUSION_FLAG=.TRUE.
  IF(FIXED_WALL_NODES_ADVECTION_DIFFUSION_FLAG) THEN
    NUMBER_OF_FIXED_WALL_NODES_ADVECTION_DIFFUSION=1
    ALLOCATE(FIXED_WALL_NODES_ADVECTION_DIFFUSION(NUMBER_OF_FIXED_WALL_NODES_ADVECTION_DIFFUSION))
    FIXED_WALL_NODES_ADVECTION_DIFFUSION=(/42/)
  ENDIF
  IF(INLET_WALL_NODES_ADVECTION_DIFFUSION_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_ADVECTION_DIFFUSION=105
    ALLOCATE(INLET_WALL_NODES_ADVECTION_DIFFUSION(NUMBER_OF_INLET_WALL_NODES_ADVECTION_DIFFUSION))
    INLET_WALL_NODES_ADVECTION_DIFFUSION=(/4,11,12,13,29,33,34,35,47,51,52,53,69,71,83,87,88,89,100,102,103,104,112,114,115, & 
    & 116,126,128,137,141,142,143,154,156,157,158,166,168,169,170,180,182,191,195,196,197,208,210,211, & 
    & 212,220,222,223,224,234,236,245,249,250,251,262,264,265,266,274,276,277,278,288,290,299,303,304, & 
    & 305,316,318,319,320,328,330,331,332,342,344,353,357,358,359,370,372,373,374,382,384,385,386,396, & 
    & 398,411,412,426,427,438,439,450/)

    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_ADVECTION_DIFFUSION(1)=42.0_CMISSDP
    BOUNDARY_CONDITIONS_ADVECTION_DIFFUSION(2)=0.0_CMISSDP
    BOUNDARY_CONDITIONS_ADVECTION_DIFFUSION(3)=0.0_CMISSDP
  ENDIF  !Set material parameters
  DIFF_COEFF_PARAM_ADVECTION_DIFFUSION(1)=1.0_CMISSDP
  DIFF_COEFF_PARAM_ADVECTION_DIFFUSION(2)=1.0_CMISSDP
  DIFF_COEFF_PARAM_ADVECTION_DIFFUSION(3)=1.0_CMISSDP
  !Set interpolation parameters
  BASIS_XI_GAUSS_SPACE=3
  BASIS_XI_GAUSS_CONCENTRATION=3
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LINEAR_SOLVER_ADVECTION_DIFFUSION_OUTPUT_TYPE=CMISSSolverNoOutput
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_ADVECTION_DIFFUSION_OUTPUT=CMISSEquationsNoOutput
  !Set solver parameters
  LINEAR_SOLVER_ADVECTION_DIFFUSION_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E20 !default: 1.0E5
  MAXIMUM_ITERATIONS=100000 !default: 100000
  RESTART_VALUE=3000 !default: 30
  LINESEARCH_ALPHA=1.0

  !
  !================================================================================================================================


  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

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

  !Start the creation of new bases
  MESH_NUMBER_OF_COMPONENTS=1
  CALL CMISSBasisTypeInitialise(BasisSpace,Err)
  CALL CMISSBasisCreateStart(BASIS_NUMBER_SPACE,BasisSpace,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasisTypeSet(BasisSpace,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL CMISSBasisNumberOfXiSet(BasisSpace,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  IF(NUMBER_OF_DIMENSIONS==2) THEN
    CALL CMISSBasisInterpolationXiSet(BasisSpace,(/BASIS_XI_INTERPOLATION_SPACE,BASIS_XI_INTERPOLATION_SPACE/),Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisSpace,(/BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE/),Err)
  ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
    CALL CMISSBasisInterpolationXiSet(BasisSpace,(/BASIS_XI_INTERPOLATION_SPACE,BASIS_XI_INTERPOLATION_SPACE, & 
      & BASIS_XI_INTERPOLATION_SPACE/),Err)                         
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisSpace,(/BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE,BASIS_XI_GAUSS_SPACE/),Err)
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(BasisSpace,Err)
  !Start the creation of another basis
  IF(BASIS_XI_INTERPOLATION_CONCENTRATION==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisConcentration=BasisSpace
  ELSE
    MESH_NUMBER_OF_COMPONENTS=MESH_NUMBER_OF_COMPONENTS+1
    !Initialise a new pressure basis
    CALL CMISSBasisTypeInitialise(BasisConcentration,Err)
    !Start the creation of a basis
    CALL CMISSBasisCreateStart(BASIS_NUMBER_CONCENTRATION,BasisConcentration,Err)
    !Set the basis type (Lagrange/Simplex)
    CALL CMISSBasisTypeSet(BasisConcentration,BASIS_TYPE,Err)
    !Set the basis xi number
    CALL CMISSBasisNumberOfXiSet(BasisConcentration,NUMBER_OF_DIMENSIONS,Err)
    !Set the basis xi interpolation and number of Gauss points
    IF(NUMBER_OF_DIMENSIONS==2) THEN
      CALL CMISSBasisInterpolationXiSet(BasisConcentration,(/BASIS_XI_INTERPOLATION_CONCENTRATION, &
        & BASIS_XI_INTERPOLATION_CONCENTRATION/),Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisConcentration,(/BASIS_XI_GAUSS_CONCENTRATION, &
        & BASIS_XI_GAUSS_CONCENTRATION/),Err)
    ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
      CALL CMISSBasisInterpolationXiSet(BasisConcentration,(/BASIS_XI_INTERPOLATION_CONCENTRATION, & 
        & BASIS_XI_INTERPOLATION_CONCENTRATION, BASIS_XI_INTERPOLATION_CONCENTRATION/),Err)                         
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisConcentration,(/BASIS_XI_GAUSS_CONCENTRATION, & 
        & BASIS_XI_GAUSS_CONCENTRATION, BASIS_XI_GAUSS_CONCENTRATION/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(BasisConcentration,Err)
  ENDIF
  !
  !================================================================================================================================
  !

  !MESH

  !Start the creation of mesh nodes
  CALL CMISSNodesTypeInitialise(Nodes,Err)
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSNodesCreateStart(Region,TOTAL_NUMBER_OF_NODES,Nodes,Err)
  CALL CMISSNodesCreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL CMISSMeshCreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL CMISSMeshNumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL CMISSMeshNumberOfComponentsSet(Mesh,MESH_NUMBER_OF_COMPONENTS,Err)
  !Specify spatial mesh component
  CALL CMISSMeshElementsTypeInitialise(MeshElementsSpace,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsConcentration,Err)
  MESH_COMPONENT_NUMBER_SPACE=1
  MESH_COMPONENT_NUMBER_CONCENTRATION=1
  CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_SPACE,BasisSpace,MeshElementsSpace,Err)
  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
    CALL CMISSMeshElementsNodesSet(MeshElementsSpace,ELEMENT_NUMBER,CM%M(ELEMENT_NUMBER,1:NUMBER_OF_ELEMENT_NODES_SPACE),Err)
  ENDDO
  CALL CMISSMeshElementsCreateFinish(MeshElementsSpace,Err)
  !Specify pressure mesh component
  IF(BASIS_XI_INTERPOLATION_CONCENTRATION==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsConcentration=MeshElementsSpace
    MESH_COMPONENT_NUMBER_CONCENTRATION=MESH_COMPONENT_NUMBER_SPACE
  ELSE
    MESH_COMPONENT_NUMBER_CONCENTRATION=MESH_COMPONENT_NUMBER_SPACE+1
    CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_CONCENTRATION,BasisConcentration,MeshElementsConcentration,Err)
    DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
      CALL CMISSMeshElementsNodesSet(MeshElementsConcentration,ELEMENT_NUMBER,CM%V(ELEMENT_NUMBER, & 
        & 1:NUMBER_OF_ELEMENT_NODES_CONCENTRATION),Err)
    ENDDO
    CALL CMISSMeshElementsCreateFinish(MeshElementsConcentration,Err)
  ENDIF
  !Finish the creation of the mesh
  CALL CMISSMeshCreateFinish(Mesh,Err)

  !
  !================================================================================================================================
  !

 !GEOMETRIC FIELD

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,DomainUserNumber,Err)
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
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDDO
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)
  !Update the geometric field parameters
  DO NODE_NUMBER=1,NUMBER_OF_NODES_SPACE
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
  
  !Create the equations_set
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
    CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,CMISSEquationsSetClassicalFieldClass, &
    & CMISSEquationsSetDiffusionEquationType,CMISSEquationsSetNoSourceALEDiffusionSubtype,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Set the equations set to be a standard Laplace problem
  
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)
 !
  !================================================================================================================================
  !

 
  !Create the equations set dependent field variables
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field variables
  CALL CMISSFieldTypeInitialise(MaterialsField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

  !
  !================================================================================================================================
  !


  !EQUATIONS  

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(Equations,EQUATIONS_ADVECTION_DIFFUSION_OUTPUT,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)

  !
  !================================================================================================================================
  !

  !Create the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a No Source Diffusion problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemClassicalFieldClass,CMISSProblemDiffusionEquationType, &
    & CMISSProblemNoSourceALEDiffusionSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Create the problem control
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  !Get the control loop
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,0.0_CMISSDP,3.0_CMISSDP,0.1_CMISSDP,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !


  !SOLVERS

  !Start the creation of the problem solvers

! !   !For the Direct Solver MUMPS, uncomment the below two lines and comment out the above five
! !   CALL SOLVER_LINEAR_TYPE_SET(LINEAR_SOLVER,SOLVER_LINEAR_DIRECT_SOLVE_TYPE,ERR,ERROR,*999)
! !   CALL SOLVER_LINEAR_DIRECT_TYPE_SET(LINEAR_SOLVER,SOLVER_DIRECT_MUMPS,ERR,ERROR,*999) 
! 
!   CALL PROBLEM_SOLVERS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  CALL CMISSSolverDynamicLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolver,1000,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)
  !
  !================================================================================================================================
  !
  !SOLVER EQUATIONS


  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)  
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

 !BOUNDARY CONDITIONS

  !Start the creation of the equations set boundary conditions for Poisson
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set fixed wall nodes
  IF(FIXED_WALL_NODES_ADVECTION_DIFFUSION_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_FIXED_WALL_NODES_ADVECTION_DIFFUSION
      NODE_NUMBER=FIXED_WALL_NODES_ADVECTION_DIFFUSION(NODE_COUNTER)
      CONDITION=CMISSBoundaryConditionFixed
!       DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=0.0_CMISSDP
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
          & NODE_NUMBER,MESH_COMPONENT_NUMBER_CONCENTRATION,CONDITION,VALUE,Err)
!       ENDDO
    ENDDO
  ENDIF
  !Set velocity boundary conditions
  IF(INLET_WALL_NODES_ADVECTION_DIFFUSION_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_ADVECTION_DIFFUSION
      NODE_NUMBER=INLET_WALL_NODES_ADVECTION_DIFFUSION(NODE_COUNTER)
      CONDITION=CMISSBoundaryConditionFixed
!       DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=0.1_CMISSDP
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,CMISSNoGlobalDerivative, &
          & NODE_NUMBER,MESH_COMPONENT_NUMBER_CONCENTRATION,CONDITION,VALUE,Err)
!       ENDDO
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)
 !
  !================================================================================================================================
  !

  !RUN SOLVERS

  !Turn of PETSc error handling
  !CALL PETSC_ERRORHANDLING_SET_ON(ERR,ERROR,*999)

  !Solve the problem
  CALL CMISSProblemSolve(Problem,Err)


  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"DiffusionIO_ALE","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"DiffusionIO_ALE","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)

  ENDIF
  
  !Output timing summary
  !CALL TIMING_SUMMARY_OUTPUT(ERR,ERROR,*999)

  !Calculate the stop times and write out the elapsed user and system times
!   CALL CPU_TIMER(USER_CPU,STOP_USER_TIME,ERR,ERROR,*999)
!   CALL CPU_TIMER(SYSTEM_CPU,STOP_SYSTEM_TIME,ERR,ERROR,*999)
! 
!   CALL WRITE_STRING_TWO_VALUE(GENERAL_OUTPUT_TYPE,"User time = ",STOP_USER_TIME(1)-START_USER_TIME(1),", System time = ", &
!     & STOP_SYSTEM_TIME(1)-START_SYSTEM_TIME(1),ERR,ERROR,*999)
!   
  !CALL CMISS_FINALISE(ERR,ERROR,*999)
  !CALL CMISSFinalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  

  STOP
  
END PROGRAM DIFFUSIONIOALEEXAMPLE
