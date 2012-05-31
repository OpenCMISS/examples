!> Main program

PROGRAM NAVIERSTOKES1DTRANSIENTEXAMPLE

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
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumberNavierStokes=6
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokes=7
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberBifurcation=8
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldUserNumberNavierStokes=12
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberNavierStokes=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberBifurcation=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: Basis_Lagrange_Hermite_TP_Type=1
  INTEGER(CMISSIntg), PARAMETER :: Basis_Quadratic_Lagrange_Interpolation=2

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverBifurcationUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverNavierStokesUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesMu=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesRho=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesE=3
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesH0=4
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesA0=5
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesA1=6
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesA2=7

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS
  
  INTEGER(CMISSIntg) :: BASIS_TYPE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_SPACE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_AREA
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_AREA
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_AREA
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_SPACE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_AREA
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_AREA
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_AREA
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE
  INTEGER(CMISSIntg) :: NUMBER_OF_OUTLET_WALL_NODES_NAVIER_STOKES
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES

  INTEGER(CMISSIntg) :: EQUATIONS_NAVIER_STOKES_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: NODE_COUNTER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE

  INTEGER, ALLOCATABLE, DIMENSION(:):: OUTLET_WALL_NODES_NAVIER_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_NAVIER_STOKES

  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_NAVIER_STOKES(3)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE1
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE1
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE,VALUE1
  REAL(CMISSDP) :: MU_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: RHO_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: E_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: H0_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: A0_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: A1_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: A2_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: Q0,Q1,Q2,A0,A1,A2,Q3,Q4,Q5,Q6,Q7,Q8

  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_THETA
  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT

  LOGICAL :: LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG
  LOGICAL :: OUTLET_WALL_NODES_NAVIER_STOKES_FLAG
  LOGICAL :: INLET_WALL_NODES_NAVIER_STOKES_FLAG

  !CMISS variables

  !Regions
  TYPE(CMISSRegionType) :: Region
  TYPE(CMISSRegionType) :: WorldRegion
  !Coordinate systems
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  !Basis
  TYPE(CMISSBasisType) :: BasisSpace
  TYPE(CMISSBasisType) :: BasisVelocity
  TYPE(CMISSBasisType) :: BasisArea
  !Nodes
  TYPE(CMISSNodesType) :: Nodes
  !Elements
  TYPE(CMISSMeshElementsType) :: MeshElementsSpace
  TYPE(CMISSMeshElementsType) :: MeshElementsVelocity
  TYPE(CMISSMeshElementsType) :: MeshElementsArea
  !Meshes
  TYPE(CMISSMeshType) :: Mesh
  !Decompositions
  TYPE(CMISSDecompositionType) :: Decomposition
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: DependentFieldNavierStokes
  TYPE(CMISSFieldType) :: MaterialsFieldNavierStokes
  TYPE(CMISSFieldType) :: IndependentFieldNavierStokes
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsNavierStokes
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsBifurcation
  !Equations sets
  TYPE(CMISSEquationsSetType) :: EquationsSetNavierStokes
  TYPE(CMISSEquationsSetType) :: EquationsSetBifurcation
  !Equations
  TYPE(CMISSEquationsType) :: EquationsNavierStokes
  TYPE(CMISSEquationsType) :: EquationsBifurcation
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: DynamicSolverNavierStokes
  TYPE(CMISSSolverType) :: NonlinearSolverNavierStokes
  TYPE(CMISSSolverType) :: LinearSolverNavierStokes
  TYPE(CMISSSolverType) :: NonlinearSolverBifurcation
  TYPE(CMISSSolverType) :: LinearSolverBifurcation
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsNavierStokes
  TYPE(CMISSSolverEquationsType) :: SolverEquationsBifurcation

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,BoundaryNodeDomain
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
  CALL CMISSDiagnosticsSetOn(CMISS_ALL_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",[""],Err)
  CALL CMISSOutputSetOn("Testing",Err)
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  !
  !================================================================================================================================
  !

  !CHECK COMPUTATIONAL NODE

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !
  !================================================================================================================================
  !

  !PROBLEM CONTROL PANEL

  BASIS_NUMBER_SPACE=1
  BASIS_NUMBER_VELOCITY=2
  BASIS_NUMBER_AREA=3
  NUMBER_OF_DIMENSIONS=1
  BASIS_TYPE=1
  BASIS_XI_INTERPOLATION_SPACE=2
  BASIS_XI_INTERPOLATION_VELOCITY=2
  BASIS_XI_INTERPOLATION_AREA=2
  NUMBER_OF_NODES_SPACE=7
  NUMBER_OF_NODES_VELOCITY=7
  NUMBER_OF_NODES_AREA=7
  TOTAL_NUMBER_OF_NODES=21
  TOTAL_NUMBER_OF_ELEMENTS=3
  NUMBER_OF_ELEMENT_NODES_SPACE=3
  NUMBER_OF_ELEMENT_NODES_VELOCITY=3
  NUMBER_OF_ELEMENT_NODES_AREA=3

  !Set material parameters
  MU_PARAM_NAVIER_STOKES=0.0035_CMISSDP !Mu(Pa.s)
  RHO_PARAM_NAVIER_STOKES=1050.0_CMISSDP !Rho(kg/m3)
  E_PARAM_NAVIER_STOKES=0.8E+8_CMISSDP !Elasticity(Pa)
  H0_PARAM_NAVIER_STOKES=0.5E-3_CMISSDP !Wall Thickness(m)
  A0_PARAM_NAVIER_STOKES=20.0e-6_CMISSDP !Wall Area(m2)
  A1_PARAM_NAVIER_STOKES=10.0e-6_CMISSDP
  A2_PARAM_NAVIER_STOKES=10.0e-6_CMISSDP
  Q0=25.0e-6_CMISSDP
  Q1=12.5e-6_CMISSDP
  Q2=12.5e-6_CMISSDP
  A0=20.0e-6_CMISSDP
  A1=10.0e-6_CMISSDP
  A2=10.0e-6_CMISSDP


  OUTLET_WALL_NODES_NAVIER_STOKES_FLAG=.TRUE.
  INLET_WALL_NODES_NAVIER_STOKES_FLAG=.TRUE.

  IF(INLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES=1
    ALLOCATE(INLET_WALL_NODES_NAVIER_STOKES(NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES))
    INLET_WALL_NODES_NAVIER_STOKES=[1]
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_NAVIER_STOKES(1)=Q0
  ENDIF

  IF(OUTLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    NUMBER_OF_OUTLET_WALL_NODES_NAVIER_STOKES=2
    ALLOCATE(OUTLET_WALL_NODES_NAVIER_STOKES(NUMBER_OF_OUTLET_WALL_NODES_NAVIER_STOKES))
    OUTLET_WALL_NODES_NAVIER_STOKES=[5,7]
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_NAVIER_STOKES(2)=Q1
    BOUNDARY_CONDITIONS_NAVIER_STOKES(3)=Q2
  ENDIF

  !Set interpolation parameters
  BASIS_XI_GAUSS_SPACE=3
  BASIS_XI_GAUSS_VELOCITY=3
  BASIS_XI_GAUSS_AREA=3
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISS_SOLVER_PROGRESS_OUTPUT
  NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISS_SOLVER_SOLVER_OUTPUT
  LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISS_SOLVER_PROGRESS_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_NAVIER_STOKES_OUTPUT=CMISS_EQUATIONS_NO_OUTPUT
  !Set time parameter
  DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME=0.0_CMISSDP
  DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME=0.01_CMISSDP
  DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT=0.01_CMISSDP
  DYNAMIC_SOLVER_NAVIER_STOKES_THETA=1.0_CMISSDP/2.0_CMISSDP
  !Set result output parameter
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY=1
  !Set solver parameters
  LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG=.FALSE.
  RELATIVE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE=1.0E-10_CMISSDP !default: 1.0E-10_CMISSDP
  RELATIVE_TOLERANCE1=1.0E-10_CMISSDP !default: 1.0E-05_CMISSDP
  ABSOLUTE_TOLERANCE1=1.0E-10_CMISSDP !default: 1.0E-10_CMISSDP
  DIVERGENCE_TOLERANCE=1.0E20 !default: 1.0E5
  MAXIMUM_ITERATIONS=100000 !default: 100000
  RESTART_VALUE=3000 !default: 30
  LINESEARCH_ALPHA=1.0

  !
  !================================================================================================================================
  !

  !COORDINATE SYSTEM

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system dimension
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !
  !================================================================================================================================
  !

  !REGION

  !Start the creation of a new region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the region label
  CALL CMISSRegion_LabelSet(Region,'OpenCMISS',Err)
  !Set the regions coordinate system as defined above
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL CMISSRegion_CreateFinish(Region,Err)

  !
  !================================================================================================================================
  !

  !BASES

  !Start the creation of SPACE bases
  MESH_NUMBER_OF_COMPONENTS=1
  CALL CMISSBasis_Initialise(BasisSpace,Err)
  CALL CMISSBasis_CreateStart(BASIS_NUMBER_SPACE,BasisSpace,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasis_TypeSet(BasisSpace,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL CMISSBasis_NumberOfXiSet(BasisSpace,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  CALL CMISSBasis_InterpolationXiSet(BasisSpace,(/BASIS_XI_INTERPOLATION_SPACE/),Err)
  CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisSpace,(/BASIS_XI_GAUSS_SPACE/),Err)
  !Finish the creation of the basis
  CALL CMISSBasis_CreateFinish(BasisSpace,Err)

  !Start the creation of VELOCITY basis
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisVelocity=BasisSpace
  ENDIF

  !Start the creation of AREA basis
  IF(BASIS_XI_INTERPOLATION_AREA==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisArea=BasisSpace
  ENDIF

  !
  !================================================================================================================================
  !

  !MESH

  !Start the creation of mesh nodes
  CALL CMISSNodes_Initialise(Nodes,Err)
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSNodes_CreateStart(Region,TOTAL_NUMBER_OF_NODES,Nodes,Err)
  CALL CMISSNodes_CreateFinish(Nodes,Err)
  !Start the creation of the mesh
  CALL CMISSMesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_DIMENSIONS,Mesh,Err)
  !Set number of mesh elements
  CALL CMISSMesh_NumberOfElementsSet(Mesh,TOTAL_NUMBER_OF_ELEMENTS,Err)
  !Set number of mesh components
  CALL CMISSMesh_NumberOfComponentsSet(Mesh,MESH_NUMBER_OF_COMPONENTS,Err)

  !Specify spatial mesh component
  CALL CMISSMeshElements_Initialise(MeshElementsSpace,Err)
  CALL CMISSMeshElements_Initialise(MeshElementsVelocity,Err)
  CALL CMISSMeshElements_Initialise(MeshElementsArea,Err)

  MESH_COMPONENT_NUMBER_SPACE=1
  MESH_COMPONENT_NUMBER_VELOCITY=1
  MESH_COMPONENT_NUMBER_AREA=1

  CALL CMISSMeshElements_CreateStart(Mesh,MESH_COMPONENT_NUMBER_SPACE,BasisSpace,MeshElementsSpace,Err)

  CALL CMISSMeshElements_NodesSet(MeshElementsSpace,1,[1,2,3],Err)
  CALL CMISSMeshElements_NodesSet(MeshElementsSpace,2,[3,4,5],Err)
  CALL CMISSMeshElements_NodesSet(MeshElementsSpace,3,[3,6,7],Err)
  CALL CMISSMeshElements_LocalElementNodeVersionSet(MeshElementsSpace,1,1,CMISS_NO_GLOBAL_DERIV,3,Err) 
  CALL CMISSMeshElements_LocalElementNodeVersionSet(MeshElementsSpace,2,2,CMISS_NO_GLOBAL_DERIV,1,Err) 
  CALL CMISSMeshElements_LocalElementNodeVersionSet(MeshElementsSpace,3,3,CMISS_NO_GLOBAL_DERIV,1,Err) 
  !GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex

  CALL CMISSMeshElements_CreateFinish(MeshElementsSpace,Err)

  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsVelocity=MeshElementsSpace
  ENDIF

  !Specify area mesh component
  IF(BASIS_XI_INTERPOLATION_AREA==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsArea=MeshElementsSpace
  ENDIF

  !Finish the creation of the mesh
  CALL CMISSMesh_CreateFinish(Mesh,Err)

  !
  !================================================================================================================================
  !

  !GEOMETRIC FIELD

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,DomainUserNumber,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region

  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the field label
  CALL CMISSField_VariableLabelSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,'Coordinates',Err)
  !Set the field type
  CALL CMISSField_TypeSet(GeometricField,CMISS_FIELD_GEOMETRIC_TYPE,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the scaling to use
  CALL CMISSField_ScalingTypeSet(GeometricField,CMISS_FIELD_NO_SCALING,Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,2
    CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDDO
  !Finish creating the field
  CALL CMISSField_CreateFinish(GeometricField,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,1,1,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,2,1,0.5_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,3,1,1.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 2,CMISS_NO_GLOBAL_DERIV,3,1,1.1_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 3,CMISS_NO_GLOBAL_DERIV,3,1,1.1_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,4,1,1.6_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,5,1,2.1_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,6,1,1.6_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,7,1,2.1_CMISSDP,Err)

  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,1,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,2,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,3,2,0.0_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 2,CMISS_NO_GLOBAL_DERIV,3,2,0.05_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 3,CMISS_NO_GLOBAL_DERIV,3,2,-0.05_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,4,2,0.25_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,5,2,0.6_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,6,2,-0.25_CMISSDP,Err)
  CALL CMISSField_ParameterSetUpdateNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,7,2,-0.6_CMISSDP,Err)

  !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for 1D Navier-Stokes
  CALL CMISSEquationsSet_Initialise(EquationsSetNavierStokes,Err)
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField, &
    & CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMISS_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE, &
    & CMISS_EQUATIONS_SET_1DTRANSIENT_NAVIER_STOKES_SUBTYPE,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSetNavierStokes,Err)
  !Set the equations set to be a dynamic Navier-Stokes problem
  
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSetNavierStokes,Err)

  !Create the equations set for Bifurcation
  CALL CMISSEquationsSet_Initialise(EquationsSetBifurcation,Err)
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(EquationsSetUserNumberBifurcation,Region,GeometricField, &
    & CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMISS_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE, &
    & CMISS_EQUATIONS_SET_BIFURCATION_NAVIER_STOKES_SUBTYPE,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSetBifurcation,Err)
  !Set the equations set to be a static Nonlinear problem
  
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSetBifurcation,Err)

  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for static Nonlinear
  CALL CMISSField_Initialise(DependentFieldNavierStokes,Err)

  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetBifurcation,DependentFieldUserNumberNavierStokes, & 
    & DependentFieldNavierStokes,Err)
  !Set the field label
  CALL CMISSField_VariableLabelSet(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,'General',Err)
  CALL CMISSField_VariableLabelSet(DependentFieldNavierStokes,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,'Derivatives',Err)
  !Set the mesh component to be used by the field components.
  COMPONENT_NUMBER=1
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
    & MESH_COMPONENT_NUMBER_VELOCITY,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldNavierStokes,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, & 
    & MESH_COMPONENT_NUMBER_VELOCITY,Err)
  COMPONENT_NUMBER=2
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
    & MESH_COMPONENT_NUMBER_AREA,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentFieldNavierStokes,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,COMPONENT_NUMBER, &
    & MESH_COMPONENT_NUMBER_AREA,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetBifurcation,Err)

  !Create the equations set dependent field variables for dynamic Navier-Stokes
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSetNavierStokes,DependentFieldUserNumberNavierStokes, & 
    & DependentFieldNavierStokes,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSetNavierStokes,Err)

  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,1,1,Q0,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,2,1,Q0,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,3,1,Q0,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 2,CMISS_NO_GLOBAL_DERIV,3,1,Q1,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 3,CMISS_NO_GLOBAL_DERIV,3,1,Q2,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,4,1,Q1,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,5,1,Q1,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,6,1,Q2,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,7,1,Q2,Err)

  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,1,2,A0,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,2,2,A0,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,3,2,A0,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 2,CMISS_NO_GLOBAL_DERIV,3,2,A1,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 3,CMISS_NO_GLOBAL_DERIV,3,2,A2,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,4,2,A1,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,5,2,A1,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,6,2,A2,Err)
  CALL CMISSField_ParameterSetUpdateNode(DependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & 1,CMISS_NO_GLOBAL_DERIV,7,2,A2,Err)
  !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE

  !
  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for dynamic Navier-Stokes
  CALL CMISSField_Initialise(MaterialsFieldNavierStokes,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetNavierStokes,MaterialsFieldUserNumberNavierStokes, & 
    & MaterialsFieldNavierStokes,Err)
  !Set the field label
  CALL CMISSField_VariableLabelSet(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,'Materials',Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetNavierStokes,Err)

  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSetBifurcation,MaterialsFieldUserNumberNavierStokes, & 
    & MaterialsFieldNavierStokes,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSetBifurcation,Err)

  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberNavierStokesMu,MU_PARAM_NAVIER_STOKES,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & MaterialsFieldUserNumberNavierStokesRho,RHO_PARAM_NAVIER_STOKES,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberNavierStokesE,E_PARAM_NAVIER_STOKES,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberNavierStokesH0,H0_PARAM_NAVIER_STOKES,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberNavierStokesA0,A0_PARAM_NAVIER_STOKES,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberNavierStokesA1,A1_PARAM_NAVIER_STOKES,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialsFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & MaterialsFieldUserNumberNavierStokesA2,A2_PARAM_NAVIER_STOKES,Err)

  !
  !================================================================================================================================x
  !

  !INDEPENDENT FIELDS

  !Create the equations set independent field variables for Bifurcations Solver
  CALL CMISSField_Initialise(IndependentFieldNavierStokes,Err)

  CALL CMISSEquationsSet_IndependentCreateStart(EquationsSetBifurcation,IndependentFieldUserNumberNavierStokes, & 
    & IndependentFieldNavierStokes,Err)
  CALL CMISSField_VariableLabelSet(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,'Characteristics',Err)
  !Set the mesh component to be used by the field components.
  DO COMPONENT_NUMBER=1,2
    CALL CMISSField_ComponentMeshComponentSet(IndependentFieldNavierStokes,CMISS_FIELD_U_VARIABLE_TYPE,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_SPACE,Err)
  ENDDO
  !Finish the equations set independent field variables
  CALL CMISSEquationsSet_IndependentCreateFinish(EquationsSetBifurcation,Err)

  !Create the equations set independent field variables for 1D Transient Solver
  CALL CMISSEquationsSet_IndependentCreateStart(EquationsSetNavierStokes,IndependentFieldUserNumberNavierStokes, & 
    & IndependentFieldNavierStokes,Err)
  !Finish the equations set independent field variables
  CALL CMISSEquationsSet_IndependentCreateFinish(EquationsSetNavierStokes,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS

  !Create the equations set equations
  CALL CMISSEquations_Initialise(EquationsNavierStokes,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetNavierStokes,EquationsNavierStokes,Err)
  !Set the equations matrices sparsity type
  !CALL CMISSEquations_SparsityTypeSet(EquationsNavierStokes,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_SparsityTypeSet(EquationsNavierStokes,CMISS_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations lumping type
  CALL CMISSEquations_LumpingTypeSet(EquationsNavierStokes,CMISS_EQUATIONS_UNLUMPED_MATRICES,Err)
  !Set the equations set output
!  CALL CMISSEquations_OutputTypeSet(EquationsNavierStokes,EQUATIONS_NAVIER_STOKES_OUTPUT,Err)
  CALL CMISSEquations_OutputTypeSet(EquationsNavierStokes,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetNavierStokes,Err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(EquationsBifurcation,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSetBifurcation,EquationsBifurcation,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(EquationsBifurcation,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
!  CALL CMISSEquations_OutputTypeSet(EquationsBifurcation,EQUATIONS_NAVIER_STOKES_OUTPUT,Err)
  CALL CMISSEquations_OutputTypeSet(EquationsBifurcation,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSetBifurcation,Err)

  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a dynamic Navier-Stokes problem
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_FLUID_MECHANICS_CLASS,CMISS_PROBLEM_NAVIER_STOKES_EQUATION_TYPE, &
    & CMISS_PROBLEM_1DTRANSIENT_NAVIER_STOKES_SUBTYPE,Err)
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblem_ControlLoopGet(Problem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoop_TimesSet(ControlLoop,DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME,DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME, & 
    & DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT,Err)
  !Set the output timing
  CALL CMISSControlLoop_TimeOutputSet(ControlLoop,DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL CMISSSolver_Initialise(DynamicSolverNavierStokes,Err)
  CALL CMISSSolver_Initialise(NonlinearSolverNavierStokes,Err)
  CALL CMISSSolver_Initialise(LinearSolverNavierStokes,Err)
  CALL CMISSSolver_Initialise(NonlinearSolverBifurcation,Err)
  CALL CMISSSolver_Initialise(LinearSolverBifurcation,Err)

  CALL CMISSProblem_SolversCreateStart(Problem,Err)

  !Get the dynamic dynamic solver
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(DynamicSolverNavierStokes,DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set theta
  CALL CMISSSolver_DynamicThetaSet(DynamicSolverNavierStokes,DYNAMIC_SOLVER_NAVIER_STOKES_THETA,Err)
  !Get the dynamic nonlinear solver
  CALL CMISSSolver_DynamicNonlinearSolverGet(DynamicSolverNavierStokes,NonlinearSolverNavierStokes,Err)
  !Set the nonlinear Jacobian type
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolverNavierStokes,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED, &
    & Err)
!  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolverNavierStokes,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(NonlinearSolverNavierStokes,NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(NonlinearSolverNavierStokes,ABSOLUTE_TOLERANCE,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(NonlinearSolverNavierStokes,RELATIVE_TOLERANCE,Err)
  !Get the dynamic nonlinear linear solver
  CALL CMISSSolver_NewtonLinearSolverGet(NonlinearSolverNavierStokes,LinearSolverNavierStokes,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(LinearSolverNavierStokes,LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG) THEN
    CALL CMISSSolver_LinearTypeSet(LinearSolverNavierStokes,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL CMISSSolver_LibraryTypeSet(LinearSolverNavierStokes,CMISS_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL CMISSSolver_LinearTypeSet(LinearSolverNavierStokes,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolverNavierStokes,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolver_LinearIterativeDivergenceToleranceSet(LinearSolverNavierStokes,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeRelativeToleranceSet(LinearSolverNavierStokes,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(LinearSolverNavierStokes,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeGMRESRestartSet(LinearSolverNavierStokes,RESTART_VALUE,Err)
  ENDIF

  !Get the static Nonlinear solver
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverBifurcationUserNumber,NonlinearSolverBifurcation,Err)
  !Set the nonlinear Jacobian type
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolverBifurcation,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED, &
    & Err)
!  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolverBifurcation,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(NonlinearSolverBifurcation,NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(NonlinearSolverBifurcation,ABSOLUTE_TOLERANCE1,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(NonlinearSolverBifurcation,RELATIVE_TOLERANCE1,Err)
  !Get the nonlinear linear solver
  CALL CMISSSolver_NewtonLinearSolverGet(NonlinearSolverBifurcation,LinearSolverBifurcation,Err)
  !Set the output type
  CALL CMISSSolver_OutputTypeSet(LinearSolverBifurcation,LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG) THEN
    CALL CMISSSolver_LinearTypeSet(LinearSolverBifurcation,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL CMISSSolver_LibraryTypeSet(LinearSolverBifurcation,CMISS_SOLVER_MUMPS_LIBRARY,Err)
  ELSE
    CALL CMISSSolver_LinearTypeSet(LinearSolverBifurcation,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolverBifurcation,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolver_LinearIterativeDivergenceToleranceSet(LinearSolverBifurcation,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolver_LinearIterativeRelativeToleranceSet(LinearSolverBifurcation,RELATIVE_TOLERANCE1,Err)
    CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(LinearSolverBifurcation,ABSOLUTE_TOLERANCE1,Err)
    CALL CMISSSolver_LinearIterativeGMRESRestartSet(LinearSolverBifurcation,RESTART_VALUE,Err)
  ENDIF

  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(DynamicSolverNavierStokes,Err)
  CALL CMISSSolver_Initialise(NonlinearSolverBifurcation,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsNavierStokes,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquationsBifurcation,Err)

  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)

  !Get the dynamic Navier-Stokes solver equations
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  CALL CMISSSolver_SolverEquationsGet(DynamicSolverNavierStokes,SolverEquationsNavierStokes,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsNavierStokes,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !Add in the equations set
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsNavierStokes,EquationsSetNavierStokes,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations

  !Get the static Nonlinear solver equations
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,SolverBifurcationUserNumber,NonlinearSolverBifurcation,Err)
  CALL CMISSSolver_SolverEquationsGet(NonlinearSolverBifurcation,SolverEquationsBifurcation,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquationsBifurcation,CMISS_SOLVER_SPARSE_MATRICES,Err)
  !Add in the equations set
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquationsBifurcation,EquationsSetBifurcation,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations

  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  !Start the creation of the equations set boundary conditions for Navier-Stokes
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsNavierStokes,Err)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsBifurcation,Err)

  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsBifurcation,BoundaryConditionsBifurcation,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsBifurcation,Err)

  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquationsNavierStokes,BoundaryConditionsNavierStokes,Err)

  !Set velocity boundary conditions
  IF(INLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    CONDITION=CMISS_BOUNDARY_CONDITION_FIXED_INLET
    VALUE=BOUNDARY_CONDITIONS_NAVIER_STOKES(1)
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
      & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV,1,1,CONDITION,VALUE,Err)
  ENDIF

  !Set area boundary conditions
  IF(OUTLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    CONDITION=CMISS_BOUNDARY_CONDITION_FIXED
    VALUE=BOUNDARY_CONDITIONS_NAVIER_STOKES(2)
    VALUE1=BOUNDARY_CONDITIONS_NAVIER_STOKES(3)
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
      & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV,5,1,CONDITION,VALUE,Err)
    CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes, &
      & CMISS_FIELD_U_VARIABLE_TYPE,1,CMISS_NO_GLOBAL_DERIV,7,1,CONDITION,VALUE1,Err)
  ENDIF

  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquationsNavierStokes,Err)

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

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM NAVIERSTOKES1DTRANSIENTEXAMPLE
