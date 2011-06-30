!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve an advection-diffusion equation using openCMISS calls.
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

!> \example ClassicalField/AdvectionDiffusion/AdvectionDiffusionIO/src/AdvectionDiffusionIOExample.f90
!! Example program to solve a diffusion equation using openCMISS calls.
!!
!! \htmlinclude ClassicalField/AdvectionDiffusion/AdvectionDiffusionIO/history.html
!<

!> Main program
PROGRAM ADVECTIONDIFFUSIONIOEXAMPLE

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
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberAdvecDiff=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberAdvecDiff=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberAdvecDiff=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopNode=0
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberAdvecDiff=12
  !INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumberAdvecDiff=14
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
  TYPE(CMISSFieldType) :: DependentFieldAdvecDiff
  TYPE(CMISSFieldType) :: MaterialsFieldAdvecDiff
  TYPE(CMISSFieldType) :: IndependentFieldAdvecDiff
  TYPE(CMISSFieldType) :: SourceFieldAdvecDiff
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsAdvecDiff
  !Equations sets
  TYPE(CMISSEquationsSetType) :: EquationsSetAdvecDiff
  !Equations
  TYPE(CMISSEquationsType) :: EquationsAdvecDiff
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: SolverAdvecDiff, LinearSolverAdvecDiff
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsAdvecDiff


#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
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

  !INITIALISE OPENCMISS

  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

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
  BASIS_XI_INTERPOLATION_SPACE=CM%IT_M
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
        & 1,CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
    ENDDO
  ENDDO
  CALL CMISSFieldParameterSetUpdateStart(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS
 
  !Create the equations_set
  CALL CMISSEquationsSetTypeInitialise(EquationsSetAdvecDiff,Err)
    CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberAdvecDiff,Region,GeometricField,CMISSEquationsSetClassicalFieldClass,&
    & CMISSEquationsSetAdvectionDiffusionEquationType,CMISSEquationsSetNoSourceAdvectionDiffusionSubtype,& 
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSetAdvecDiff,Err)
  !Set the equations set to be a standard Laplace problem
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSetAdvecDiff,Err)

  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables
  CALL CMISSFieldTypeInitialise(DependentFieldAdvecDiff,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetAdvecDiff,DependentFieldUserNumberAdvecDiff,DependentFieldAdvecDiff,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetAdvecDiff,Err)

  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set material field variables
  CALL CMISSFieldTypeInitialise(MaterialsFieldAdvecDiff,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetAdvecDiff,MaterialsFieldUserNumberAdvecDiff,MaterialsFieldAdvecDiff,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetAdvecDiff,Err)

  !
  !================================================================================================================================
  !

  !SOURCE FIELDS


  !Create the equations set source field variables
!   CALL CMISSFieldTypeInitialise(SourceFieldAdvecDiff,Err)
!   CALL CMISSEquationsSetSourceCreateStart(EquationsSetAdvecDiff,SourceFieldUserNumberAdvecDiff,SourceFieldAdvecDiff,Err)
!   CALL CMISSFieldComponentInterpolationSet(SourceFieldAdvecDiff,CMISSFieldUVariableType,1,CMISSFieldNodeBasedInterpolation,Err)
!   !Finish the equations set dependent field variables
!   CALL CMISSEquationsSetSourceCreateFinish(EquationsSetAdvecDiff,Err)


  !
  !================================================================================================================================
  !

  !INDEPENDENT FIELDS

 ! CALL CMISSFieldParameterSetDataGet(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,GEOMETRIC_PARAMETERS,Err)


  !Create the equations set independent field variables
  CALL CMISSFieldTypeInitialise(IndependentFieldAdvecDiff,Err)
  CALL CMISSEquationsSetIndependentCreateStart(EquationsSetAdvecDiff,IndependentFieldUserNumberAdvecDiff, &
    & IndependentFieldAdvecDiff,Err)
!   IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
!   CALL CMISSFieldComponentInterpolationSet(IndependentField,CMISSFieldUVariableType,1,CMISSFieldNodeBasedInterpolation,Err) 
!   CALL CMISSFieldComponentInterpolationSet(IndependentField,CMISSFieldUVariableType,2,CMISSFieldNodeBasedInterpolation,Err)
!   ENDIF
  !Set the mesh component to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(IndependentFieldAdvecDiff,CMISSFieldUVariableType,MESH_COMPONENT_NUMBER_SPACE, & 
    & MESH_COMPONENT_NUMBER_CONCENTRATION,Err)

  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetIndependentCreateFinish(EquationsSetAdvecDiff,Err)
 
  !
  !================================================================================================================================
  !


  !EQUATIONS


  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(EquationsAdvecDiff,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetAdvecDiff,EquationsAdvecDiff,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(EquationsAdvecDiff,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(EquationsAdvecDiff,EQUATIONS_ADVECTION_DIFFUSION_OUTPUT,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetAdvecDiff,Err)

  !
  !================================================================================================================================
  !

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a No Source Diffusion problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemClassicalFieldClass,CMISSProblemAdvectionDiffusionEquationType, &
    & CMISSProblemNoSourceAdvectionDiffusionSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblemControlLoopGet(Problem,ControlLoopNode,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,0.0_CMISSDP,0.03_CMISSDP,0.01_CMISSDP,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !


  !SOLVERS

  !Start the creation of the problem solvers
!  
! !   !For the Direct Solver MUMPS, uncomment the below two lines and comment out the above five
! !   CALL SOLVER_LINEAR_TYPE_SET(LINEAR_SOLVER,SOLVER_LINEAR_DIRECT_SOLVE_TYPE,ERR,ERROR,*999)
! !   CALL SOLVER_LINEAR_DIRECT_TYPE_SET(LINEAR_SOLVER,SOLVER_DIRECT_MUMPS,ERR,ERROR,*999) 
! 

  CALL CMISSSolverTypeInitialise(SolverAdvecDiff,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverAdvecDiff,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,SolverAdvecDiff,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  CALL CMISSSolverOutputTypeSet(SolverAdvecDiff,CMISSSolverProgressOutput,Err)
  CALL CMISSSolverDynamicLinearSolverGet(SolverAdvecDiff,LinearSolverAdvecDiff,Err)
  CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolverAdvecDiff,300,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !
  !SOLVER EQUATIONS

  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(SolverAdvecDiff,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsAdvecDiff,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,SolverAdvecDiff,Err)
  CALL CMISSSolverSolverEquationsGet(SolverAdvecDiff,SolverEquationsAdvecDiff,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsAdvecDiff,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)  
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsAdvecDiff,EquationsSetAdvecDiff,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

 !BOUNDARY CONDITIONS

  !Start the creation of the equations set boundary conditions for Poisson
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsAdvecDiff,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsAdvecDiff,BoundaryConditionsAdvecDiff,Err)
  !Set fixed wall nodes
  IF(FIXED_WALL_NODES_ADVECTION_DIFFUSION_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_FIXED_WALL_NODES_ADVECTION_DIFFUSION
      NODE_NUMBER=FIXED_WALL_NODES_ADVECTION_DIFFUSION(NODE_COUNTER)
      CONDITION=CMISSBoundaryConditionFixed
!       DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
        VALUE=0.0_CMISSDP
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsAdvecDiff,DependentFieldAdvecDiff,CMISSFieldUVariableType, &
          & CMISSNoGlobalDerivative, &
          & 1,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONCENTRATION,CONDITION,VALUE,Err)
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
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsAdvecDiff,DependentFieldAdvecDiff,CMISSFieldUVariableType, &
          & CMISSNoGlobalDerivative, &
          & 1,NODE_NUMBER,MESH_COMPONENT_NUMBER_CONCENTRATION,CONDITION,VALUE,Err)
!       ENDDO
    ENDDO
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsAdvecDiff,Err)

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
    CALL CMISSFieldIONodesExport(Fields,"AdvectionDiffusionIO","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"AdvectionDiffusionIO","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF
  

  !Finialise CMISS
  !CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM ADVECTIONDIFFUSIONIOEXAMPLE
