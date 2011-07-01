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
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumberNavierStokes=8
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: Basis_Quadratic_Lagrange_Interpolation=2
  INTEGER(CMISSIntg), PARAMETER :: Basis_Lagrange_Hermite_TP_Type=1

  INTEGER(CMISSIntg), PARAMETER :: DomainUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolverNavierStokesUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesMu=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesRho=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesE=3
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesH0=4
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesA0=5
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumberNavierStokesSigma=6

  !Program types

  TYPE(EXPORT_CONTAINER):: CM

  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_DIMENSIONS
  
  INTEGER(CMISSIntg) :: BASIS_TYPE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_SPACE
  INTEGER(CMISSIntg) :: BASIS_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_NUMBER_AREA
!  INTEGER(CMISSIntg) :: BASIS_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_AREA
!  INTEGER(CMISSIntg) :: BASIS_XI_GAUSS_PRESSURE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_SPACE
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_VELOCITY
  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_AREA
!  INTEGER(CMISSIntg) :: BASIS_XI_INTERPOLATION_PRESSURE
  INTEGER(CMISSIntg) :: MESH_NUMBER_OF_COMPONENTS
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_SPACE
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_VELOCITY
  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_AREA
!  INTEGER(CMISSIntg) :: MESH_COMPONENT_NUMBER_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_AREA
!  INTEGER(CMISSIntg) :: NUMBER_OF_NODES_PRESSURE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_SPACE
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_VELOCITY
  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_AREA
!  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENT_NODES_PRESSURE
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_NODES
  INTEGER(CMISSIntg) :: TOTAL_NUMBER_OF_ELEMENTS
  INTEGER(CMISSIntg) :: MAXIMUM_ITERATIONS
  INTEGER(CMISSIntg) :: RESTART_VALUE
! INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER(CMISSIntg) :: NUMBER_OF_OUTLET_WALL_NODES_NAVIER_STOKES
  INTEGER(CMISSIntg) :: NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES

  INTEGER(CMISSIntg) :: EQUATIONS_NAVIER_STOKES_OUTPUT
  INTEGER(CMISSIntg) :: COMPONENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_NUMBER
  INTEGER(CMISSIntg) :: ELEMENT_NUMBER
  INTEGER(CMISSIntg) :: NODE_COUNTER
  INTEGER(CMISSIntg) :: CONDITION

  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY
  INTEGER(CMISSIntg) :: DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE
  INTEGER(CMISSIntg) :: LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE

  INTEGER, ALLOCATABLE, DIMENSION(:):: OUTLET_WALL_NODES_NAVIER_STOKES
  INTEGER, ALLOCATABLE, DIMENSION(:):: INLET_WALL_NODES_NAVIER_STOKES


  REAL(CMISSDP) :: INITIAL_FIELD_NAVIER_STOKES(2)
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_NAVIER_STOKES(2)
  REAL(CMISSDP) :: DIVERGENCE_TOLERANCE
  REAL(CMISSDP) :: RELATIVE_TOLERANCE
  REAL(CMISSDP) :: ABSOLUTE_TOLERANCE
  REAL(CMISSDP) :: LINESEARCH_ALPHA
  REAL(CMISSDP) :: VALUE
  REAL(CMISSDP) :: MU_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: RHO_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: E_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: H0_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: A0_PARAM_NAVIER_STOKES
  REAL(CMISSDP) :: SIGMA_PARAM_NAVIER_STOKES

  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME
  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_THETA
  REAL(CMISSDP) :: DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT

  LOGICAL :: EXPORT_FIELD_IO
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
!  TYPE(CMISSBasisType) :: BasisPressure
  !Nodes
  TYPE(CMISSNodesType) :: Nodes
  !Elements
  TYPE(CMISSMeshElementsType) :: MeshElementsSpace
  TYPE(CMISSMeshElementsType) :: MeshElementsVelocity
  TYPE(CMISSMeshElementsType) :: MeshElementsArea
!  TYPE(CMISSMeshElementsType) :: MeshElementsPressure
  !Meshes
  TYPE(CMISSMeshType) :: Mesh
  !Decompositions
  TYPE(CMISSDecompositionType) :: Decomposition
  !Fields
  TYPE(CMISSFieldsType) :: Fields
  !Field types
  TYPE(CMISSFieldType) :: GeometricField
  TYPE(CMISSFieldType) :: DependentFieldNavierStokes
  TYPE(CMISSFieldType) :: MaterialsFieldNavierStokes
  !Boundary conditions
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditionsNavierStokes
  !Equations sets
  TYPE(CMISSEquationsSetType) :: EquationsSetNavierStokes
  !Equations
  TYPE(CMISSEquationsType) :: EquationsNavierStokes
  !Problems
  TYPE(CMISSProblemType) :: Problem
  !Control loops
  TYPE(CMISSControlLoopType) :: ControlLoop
  !Solvers
  TYPE(CMISSSolverType) :: DynamicSolverNavierStokes
  TYPE(CMISSSolverType) :: NonlinearSolverNavierStokes
  TYPE(CMISSSolverType) :: LinearSolverNavierStokes
  !Solver equations
  TYPE(CMISSSolverEquationsType) :: SolverEquationsNavierStokes

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

  !================================================================================================================================
  !

  !INITIALISE OPENCMISS

  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL CMISSDiagnosticsSetOn(CMISSInDiagType,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"],Err)
  CALL CMISSOutputSetOn("Testing",Err)
  !CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

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
!  BASIS_NUMBER_PRESSURE=4
  NUMBER_OF_DIMENSIONS=1
  BASIS_TYPE=1
  BASIS_XI_INTERPOLATION_SPACE=2
  BASIS_XI_INTERPOLATION_VELOCITY=2
  BASIS_XI_INTERPOLATION_AREA=2
!  BASIS_XI_INTERPOLATION_PRESSURE=2
  NUMBER_OF_NODES_SPACE=9
  NUMBER_OF_NODES_VELOCITY=9
  NUMBER_OF_NODES_AREA=9
!  NUMBER_OF_NODES_PRESSURE=3
  TOTAL_NUMBER_OF_NODES=27
  TOTAL_NUMBER_OF_ELEMENTS=3
  NUMBER_OF_ELEMENT_NODES_SPACE=3
  NUMBER_OF_ELEMENT_NODES_VELOCITY=3
  NUMBER_OF_ELEMENT_NODES_AREA=3
!  NUMBER_OF_ELEMENT_NODES_PRESSURE=3
  !Set initial values
  INITIAL_FIELD_NAVIER_STOKES(1)=.05_CMISSDP
  INITIAL_FIELD_NAVIER_STOKES(2)=18.1_CMISSDP
!  INITIAL_FIELD_NAVIER_STOKES(3)=1000.0_CMISSDP
  !Set boundary conditions
  OUTLET_WALL_NODES_NAVIER_STOKES_FLAG=.TRUE.
  INLET_WALL_NODES_NAVIER_STOKES_FLAG=.TRUE.

  IF(INLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES=1
    ALLOCATE(INLET_WALL_NODES_NAVIER_STOKES(NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES))
    INLET_WALL_NODES_NAVIER_STOKES=[1]
    !Set initial boundary conditions
    BOUNDARY_CONDITIONS_NAVIER_STOKES(1)=.05_CMISSDP
!    BOUNDARY_CONDITIONS_NAVIER_STOKES(3)=1000.0_CMISSDP
  ENDIF

  IF(OUTLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    NUMBER_OF_OUTLET_WALL_NODES_NAVIER_STOKES=2
    ALLOCATE(OUTLET_WALL_NODES_NAVIER_STOKES(NUMBER_OF_OUTLET_WALL_NODES_NAVIER_STOKES))
    OUTLET_WALL_NODES_NAVIER_STOKES=[5,7]
   !Set initial boundary conditions
    BOUNDARY_CONDITIONS_NAVIER_STOKES(2)=11.4_CMISSDP
  ENDIF

  !Set material parameters
  MU_PARAM_NAVIER_STOKES=0.0035_CMISSDP !g/(mm.s)=kg/(m.s)=Pa.s
  RHO_PARAM_NAVIER_STOKES=1.050_CMISSDP !g/cm3
  E_PARAM_NAVIER_STOKES=800_CMISSDP !KPa
  H0_PARAM_NAVIER_STOKES=0.5_CMISSDP !mm
  A0_PARAM_NAVIER_STOKES=18.1_CMISSDP !mm2
  SIGMA_PARAM_NAVIER_STOKES=0.5_CMISSDP
  !Set interpolation parameters
  BASIS_XI_GAUSS_SPACE=3
  BASIS_XI_GAUSS_VELOCITY=3
  BASIS_XI_GAUSS_AREA=3
!  BASIS_XI_GAUSS_PRESSURE=3
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISSSolverProgressOutput
  NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISSSolverProgressOutput
  LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE=CMISSSolverProgressOutput
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EQUATIONS_NAVIER_STOKES_OUTPUT=CMISSEquationsNoOutput
  !Set time parameter
  DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME=0.0_CMISSDP
  DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME=200.0_CMISSDP
  DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT=100.0_CMISSDP
  DYNAMIC_SOLVER_NAVIER_STOKES_THETA=1.0_CMISSDP/2.0_CMISSDP
  !Set result output parameter
  DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY=1
  !Set solver parameters
  LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG=.FALSE.
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

  !Start the creation of SPACE bases
  MESH_NUMBER_OF_COMPONENTS=1
  CALL CMISSBasisTypeInitialise(BasisSpace,Err)
  CALL CMISSBasisCreateStart(BASIS_NUMBER_SPACE,BasisSpace,Err)
  !Set the basis type (Lagrange/Simplex)
  CALL CMISSBasisTypeSet(BasisSpace,BASIS_TYPE,Err)
  !Set the basis xi number
  CALL CMISSBasisNumberOfXiSet(BasisSpace,NUMBER_OF_DIMENSIONS,Err)
  !Set the basis xi interpolation and number of Gauss points
  CALL CMISSBasisInterpolationXiSet(BasisSpace,(/BASIS_XI_INTERPOLATION_SPACE/),Err)
  CALL CMISSBasisQuadratureNumberOfGaussXiSet(BasisSpace,(/BASIS_XI_GAUSS_SPACE/),Err)
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(BasisSpace,Err)

  !Start the creation of VELOCITY basis
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisVelocity=BasisSpace
  ENDIF

  !Start the creation of AREA basis
  IF(BASIS_XI_INTERPOLATION_AREA==BASIS_XI_INTERPOLATION_SPACE) THEN
    BasisArea=BasisSpace
  ENDIF

  !Start the creation of PRESSURE basis
!  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_SPACE) THEN
!    BasisPressure=BasisSpace
!  ENDIF

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
  CALL CMISSMeshElementsTypeInitialise(MeshElementsVelocity,Err)
  CALL CMISSMeshElementsTypeInitialise(MeshElementsArea,Err)
!  CALL CMISSMeshElementsTypeInitialise(MeshElementsPressure,Err)
  MESH_COMPONENT_NUMBER_SPACE=1
  MESH_COMPONENT_NUMBER_VELOCITY=1
  MESH_COMPONENT_NUMBER_AREA=1
!  MESH_COMPONENT_NUMBER_PRESSURE=1
  CALL CMISSMeshElementsCreateStart(Mesh,MESH_COMPONENT_NUMBER_SPACE,BasisSpace,MeshElementsSpace,Err)

!  DO ELEMENT_NUMBER=1,TOTAL_NUMBER_OF_ELEMENTS
!    CALL CMISSMeshElementsNodesSet(MeshElementsSpace,ELEMENT_NUMBER,[(2*ELEMENT_NUMBER)-1, &
!     & (2*ELEMENT_NUMBER),(2*ELEMENT_NUMBER)+1],Err)
!  ENDDO

  CALL CMISSMeshElementsNodesSet(MeshElementsSpace,1,[1,2,3],Err)
  CALL CMISSMeshElementsNodesSet(MeshElementsSpace,2,[3,4,5],Err)
  CALL CMISSMeshElementsNodesSet(MeshElementsSpace,3,[3,6,7],Err)
  CALL CMISSMeshElementsLocalElementNodeVersionSet(MeshElementsSpace,1,1,CMISSNoGlobalDerivative,3,Err) 
  CALL CMISSMeshElementsLocalElementNodeVersionSet(MeshElementsSpace,2,5,CMISSNoGlobalDerivative,1,Err) 
  CALL CMISSMeshElementsLocalElementNodeVersionSet(MeshElementsSpace,3,6,CMISSNoGlobalDerivative,1,Err) 
  !GlobalElementNumber,VersionNumber,DerivativeNumber,ElementNodeIndex

  CALL CMISSMeshElementsCreateFinish(MeshElementsSpace,Err)

  !Specify velocity mesh component
  IF(BASIS_XI_INTERPOLATION_VELOCITY==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsVelocity=MeshElementsSpace
  ENDIF

  !Specify area mesh component
  IF(BASIS_XI_INTERPOLATION_AREA==BASIS_XI_INTERPOLATION_SPACE) THEN
    MeshElementsArea=MeshElementsSpace
  ENDIF

  !Specify pressure mesh component
!  IF(BASIS_XI_INTERPOLATION_PRESSURE==BASIS_XI_INTERPOLATION_SPACE) THEN
!    MeshElementsPressure=MeshElementsSpace
!  ENDIF
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
!  DO NODE_NUMBER=1,NUMBER_OF_NODES_SPACE
!    DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
!      VALUE=(1-1)*0.5
!      CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
!        & 1,CMISSNoGlobalDerivative,NODE_NUMBER,COMPONENT_NUMBER,VALUE,Err)
!    ENDDO
!  ENDDO

  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,1,1,0.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,2,1,0.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,3,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 2,CMISSNoGlobalDerivative,3,1,1.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 3,CMISSNoGlobalDerivative,3,1,1.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 4,CMISSNoGlobalDerivative,3,1,1.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 5,CMISSNoGlobalDerivative,3,1,1.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 6,CMISSNoGlobalDerivative,3,1,1.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,4,1,2.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,5,1,2.5_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,6,1,2.0_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,7,1,2.5_CMISSDP,Err)
  !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE


  CALL CMISSFieldParameterSetUpdateStart(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)

  !
  !================================================================================================================================
  !

  !EQUATIONS SETS

  !Create the equations set for dynamic Navier-Stokes
  CALL CMISSEquationsSetTypeInitialise(EquationsSetNavierStokes,Err)
  CALL CMISSFieldTypeInitialise(EquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumberNavierStokes,Region,GeometricField, &
    & CMISSEquationsSetFluidMechanicsClass,CMISSEquationsSetNavierStokesEquationType, &
    & CMISSEquationsSet1DTransientNavierStokesSubtype,EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSetNavierStokes,Err)
  !Set the equations set to be a dynamic Navier-Stokes problem
  
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSetNavierStokes,Err)


  !
  !================================================================================================================================
  !

  !DEPENDENT FIELDS

  !Create the equations set dependent field variables for dynamic Navier-Stokes
  CALL CMISSFieldTypeInitialise(DependentFieldNavierStokes,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSetNavierStokes,DependentFieldUserNumberNavierStokes, & 
    & DependentFieldNavierStokes,Err)
  !Set the mesh component to be used by the field components.
  COMPONENT_NUMBER=1
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldUVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, & 
      & MESH_COMPONENT_NUMBER_VELOCITY,Err)


  COMPONENT_NUMBER=2
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldUVariableType,COMPONENT_NUMBER, &
      & MESH_COMPONENT_NUMBER_AREA,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, &
      & MESH_COMPONENT_NUMBER_AREA,Err)


!  COMPONENT_NUMBER=3
!    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldUVariableType,COMPONENT_NUMBER, &
!      & MESH_COMPONENT_NUMBER_PRESSURE,Err)
!    CALL CMISSFieldComponentMeshComponentSet(DependentFieldNavierStokes,CMISSFieldDeludelnVariableType,COMPONENT_NUMBER, &
!      & MESH_COMPONENT_NUMBER_PRESSURE,Err)


  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSetNavierStokes,Err)

  !Initialise dependent field
  DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS+1
    CALL CMISSFieldComponentValuesInitialise(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
      & COMPONENT_NUMBER,INITIAL_FIELD_NAVIER_STOKES(COMPONENT_NUMBER),Err)
  ENDDO
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,1,1,.05_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,2,1,.05_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,3,1,.05_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 2,CMISSNoGlobalDerivative,3,1,.05_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 3,CMISSNoGlobalDerivative,3,1,.04_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 4,CMISSNoGlobalDerivative,3,1,.04_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 5,CMISSNoGlobalDerivative,3,1,.04_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 6,CMISSNoGlobalDerivative,3,1,.04_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,4,1,.04_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,5,1,.04_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,6,1,.04_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,7,1,.04_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,1,2,18.1_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,2,2,18.1_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,3,2,18.1_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 2,CMISSNoGlobalDerivative,3,2,18.1_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 3,CMISSNoGlobalDerivative,3,2,11.4_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 4,CMISSNoGlobalDerivative,3,2,11.4_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 5,CMISSNoGlobalDerivative,3,2,11.4_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 6,CMISSNoGlobalDerivative,3,2,11.4_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,4,2,11.4_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,5,2,11.4_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,6,2,11.4_CMISSDP,Err)
  CALL CMISSFieldParameterSetUpdateNode(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & 1,CMISSNoGlobalDerivative,7,2,11.4_CMISSDP,Err)
  !VERSION_NUMBER,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE

  CALL CMISSFieldParameterSetUpdateStart(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)
  CALL CMISSFieldParameterSetUpdateFinish(DependentFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType,Err)

  !================================================================================================================================
  !

  !MATERIALS FIELDS

  !Create the equations set materials field variables for dynamic Navier-Stokes
  CALL CMISSFieldTypeInitialise(MaterialsFieldNavierStokes,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSetNavierStokes,MaterialsFieldUserNumberNavierStokes, & 
    & MaterialsFieldNavierStokes,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSetNavierStokes,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & MaterialsFieldUserNumberNavierStokesMu,MU_PARAM_NAVIER_STOKES,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, & 
    & MaterialsFieldUserNumberNavierStokesRho,RHO_PARAM_NAVIER_STOKES,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & MaterialsFieldUserNumberNavierStokesE,E_PARAM_NAVIER_STOKES,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & MaterialsFieldUserNumberNavierStokesH0,H0_PARAM_NAVIER_STOKES,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & MaterialsFieldUserNumberNavierStokesA0,A0_PARAM_NAVIER_STOKES,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialsFieldNavierStokes,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & MaterialsFieldUserNumberNavierStokesSigma,SIGMA_PARAM_NAVIER_STOKES,Err)

  !================================================================================================================================x
  !

  !EQUATIONS


  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(EquationsNavierStokes,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSetNavierStokes,EquationsNavierStokes,Err)
  !Set the equations matrices sparsity type
  !CALL CMISSEquationsSparsityTypeSet(EquationsNavierStokes,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsSparsityTypeSet(EquationsNavierStokes,CMISSEquationsFullMatrices,Err)
  !Set the equations lumping type
  CALL CMISSEquationsLumpingTypeSet(EquationsNavierStokes,CMISSEquationsUnlumpedMatrices,Err)
  !Set the equations set output
  !CALL CMISSEquationsOutputTypeSet(EquationsNavierStokes,EQUATIONS_NAVIER_STOKES_OUTPUT,Err)
  !Set the equations set output
  CALL CMISSEquationsOutputTypeSet(EquationsNavierStokes,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSetNavierStokes,Err)


  !
  !================================================================================================================================
  !

  !PROBLEMS

  !Start the creation of a problem.
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a dynamic Navier-Stokes problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemFluidMechanicsClass,CMISSProblemNavierStokesEquationType, &
    & CMISSProblem1DTransientNavierStokesSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)
  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMISSProblemControlLoopGet(Problem,CMISSControlLoopNode,ControlLoop,Err)
  !Set the times
  CALL CMISSControlLoopTimesSet(ControlLoop,DYNAMIC_SOLVER_NAVIER_STOKES_START_TIME,DYNAMIC_SOLVER_NAVIER_STOKES_STOP_TIME, & 
    & DYNAMIC_SOLVER_NAVIER_STOKES_TIME_INCREMENT,Err)
  !Set the output timing
  CALL CMISSControlLoopTimeOutputSet(ControlLoop,DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_FREQUENCY,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVERS

  !Start the creation of the problem solvers
  CALL CMISSSolverTypeInitialise(DynamicSolverNavierStokes,Err)
  CALL CMISSSolverTypeInitialise(NonlinearSolverNavierStokes,Err)
  CALL CMISSSolverTypeInitialise(LinearSolverNavierStokes,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  !Get the dynamic dymamic solver
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(DynamicSolverNavierStokes,DYNAMIC_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set theta
  CALL CMISSSolverDynamicThetaSet(DynamicSolverNavierStokes,DYNAMIC_SOLVER_NAVIER_STOKES_THETA,Err)
!   CALL CMISSSolverDynamicDynamicSet(DynamicSolverNavierStokes,.TRUE.,Err)
  !Get the dynamic nonlinear solver
  CALL CMISSSolverDynamicNonlinearSolverGet(DynamicSolverNavierStokes,NonlinearSolverNavierStokes,Err)
  !Set the nonlinear Jacobian type
!  CALL CMISSSolverNewtonJacobianCalculationTypeSet(NonlinearSolverNavierStokes,CMISSSolverNewtonJacobianAnalyticCalculated,Err)
 !Set the nonlinear Jacobian type
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(NonlinearSolverNavierStokes,CMISSSolverNewtonJacobianFDCalculated,Err)
   !Set the output type
  !CALL CMISSSolverOutputTypeSet(NonlinearSolverNavierStokes,NONLINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(NonlinearSolverNavierStokes,CMISSSolverSolverMatrixOutput,Err)
  !Set the solver settings
  CALL CMISSSolverNewtonAbsoluteToleranceSet(NonlinearSolverNavierStokes,ABSOLUTE_TOLERANCE,Err)
  CALL CMISSSolverNewtonRelativeToleranceSet(NonlinearSolverNavierStokes,RELATIVE_TOLERANCE,Err)
  !Get the dynamic nonlinear linear solver
  CALL CMISSSolverNewtonLinearSolverGet(NonlinearSolverNavierStokes,LinearSolverNavierStokes,Err)
  !Set the output type
  CALL CMISSSolverOutputTypeSet(LinearSolverNavierStokes,LINEAR_SOLVER_NAVIER_STOKES_OUTPUT_TYPE,Err)
  !Set the solver settings
  IF(LINEAR_SOLVER_NAVIER_STOKES_DIRECT_FLAG) THEN
    CALL CMISSSolverLinearTypeSet(LinearSolverNavierStokes,CMISSSolverLinearDirectSolveType,Err)
    CALL CMISSSolverLibraryTypeSet(LinearSolverNavierStokes,CMISSSolverMUMPSLibrary,Err)
  ELSE
    CALL CMISSSolverLinearTypeSet(LinearSolverNavierStokes,CMISSSolverLinearIterativeSolveType,Err)
    CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolverNavierStokes,MAXIMUM_ITERATIONS,Err)
    CALL CMISSSolverLinearIterativeDivergenceToleranceSet(LinearSolverNavierStokes,DIVERGENCE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeRelativeToleranceSet(LinearSolverNavierStokes,RELATIVE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeAbsoluteToleranceSet(LinearSolverNavierStokes,ABSOLUTE_TOLERANCE,Err)
    CALL CMISSSolverLinearIterativeGMRESRestartSet(LinearSolverNavierStokes,RESTART_VALUE,Err)
  ENDIF
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !SOLVER EQUATIONS

  !Start the creation of the problem solver equations
  CALL CMISSSolverTypeInitialise(DynamicSolverNavierStokes,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquationsNavierStokes,Err)

  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the dynamic solver equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverNavierStokesUserNumber,DynamicSolverNavierStokes,Err)
  CALL CMISSSolverSolverEquationsGet(DynamicSolverNavierStokes,SolverEquationsNavierStokes,Err)
  !Set the solver equations sparsity
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsNavierStokes,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquationsNavierStokes,CMISSSolverEquationsFullMatrices,Err)
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquationsNavierStokes,EquationsSetNavierStokes,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !
  !================================================================================================================================
  !

  !BOUNDARY CONDITIONS

  !Start the creation of the equations set boundary conditions for Navier-Stokes
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditionsNavierStokes,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquationsNavierStokes,BoundaryConditionsNavierStokes,Err)

  !Set velocity boundary conditions
  IF(INLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES
      NODE_NUMBER=INLET_WALL_NODES_NAVIER_STOKES(NODE_COUNTER)
      CONDITION=CMISSBoundaryConditionInletWall
      CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
        DO COMPONENT_NUMBER=1,NUMBER_OF_DIMENSIONS
          VALUE=BOUNDARY_CONDITIONS_NAVIER_STOKES(COMPONENT_NUMBER)
          CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes,CMISSFieldUVariableType,1, &
            & CMISSNoGlobalDerivative, &
          & NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  !Set area boundary conditions
  IF(OUTLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
    DO NODE_COUNTER=1,NUMBER_OF_OUTLET_WALL_NODES_NAVIER_STOKES
      NODE_NUMBER=OUTLET_WALL_NODES_NAVIER_STOKES(NODE_COUNTER)
      CONDITION=CMISSBoundaryConditionFixed
      DO COMPONENT_NUMBER=2,NUMBER_OF_DIMENSIONS+1
        VALUE=BOUNDARY_CONDITIONS_NAVIER_STOKES(COMPONENT_NUMBER)
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes,CMISSFieldUVariableType,1, &
          & CMISSNoGlobalDerivative, &
          & NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
      ENDDO
    ENDDO
  ENDIF

  !Set pressure boundary conditions
!  IF(INLET_WALL_NODES_NAVIER_STOKES_FLAG) THEN
!    DO NODE_COUNTER=1,NUMBER_OF_INLET_WALL_NODES_NAVIER_STOKES
!      NODE_NUMBER=INLET_WALL_NODES_NAVIER_STOKES(NODE_COUNTER)
!      CONDITION=CMISSBoundaryConditionInletWall
!      CALL CMISSDecompositionNodeDomainGet(Decomposition,NODE_NUMBER,1,BoundaryNodeDomain,Err)
!      IF(BoundaryNodeDomain==ComputationalNodeNumber) THEN
!        DO COMPONENT_NUMBER=3,NUMBER_OF_DIMENSIONS+2
!          VALUE=BOUNDARY_CONDITIONS_NAVIER_STOKES(COMPONENT_NUMBER)
!          CALL CMISSBoundaryConditionsSetNode(BoundaryConditionsNavierStokes,DependentFieldNavierStokes,CMISSFieldUVariableType,1,CMISSNoGlobalDerivative, &
!          & NODE_NUMBER,COMPONENT_NUMBER,CONDITION,VALUE,Err)
!        ENDDO
!      ENDIF
!    ENDDO
!  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquationsNavierStokes,Err)

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
    CALL CMISSFieldIONodesExport(Fields,"1DTransientNavierStokes","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"1DTransientNavierStokes","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)
    WRITE(*,'(A)') "Field exported!"
  ENDIF
  
  !Finialise CMISS
!  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

END PROGRAM NAVIERSTOKES1DTRANSIENTEXAMPLE
