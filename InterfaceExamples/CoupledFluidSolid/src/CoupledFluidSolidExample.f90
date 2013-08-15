!> \file
!> \author Andreas Hessenthaler
!> \brief This is an example program which solves a weakly coupled FiniteElasticity-ALENavierStokes equation using OpenCMISS calls.
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

!> \example InterfaceExamples/Coupled-FiniteElasticity-ALENavierStokes/src/Coupled-FiniteElasticity-ALENavierStokes.f90
!! This is an example program which solves a weakly coupled FiniteElasticity-ALENavierStokes equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/InterfaceExamples/Coupled-FiniteElasticity-ALENavierStokes/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/InterfaceExamples/Coupled-FiniteElasticity-ALENavierStokes/build-gnu'>Linux GNU Build</a>
!<


    !   2D case
    !
    !                    =>>=  =  =  =  =  =  =
    !                    =>>=  =  =  =  =  =  =
    !                    =>>=  =  =  =  =  =  =
    !   HeightFluid=3    =>>=  =  =  =  =  =  =          FLUID     REGION(2)
    !                    => =  =  =  =  =  =  =
    !                    => =  =  =  =  =  =  =
    !                    x == == == == == == ==
    !                    x=====================
    !   HeightSolid=1    x=====================          SOLID     REGION(1)
    !
    !                           Width=4 


!> Main program
PROGRAM CoupledFluidSolidExample

  USE OPENCMISS
  
#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HeightSolid=2.0_CMISSDP !cm
  REAL(CMISSDP), PARAMETER :: HeightFluid=3.0_CMISSDP !cm
  REAL(CMISSDP), PARAMETER :: WidthSolid=1.0_CMISSDP !cm
  REAL(CMISSDP), PARAMETER :: WidthFluid=5.0_CMISSDP !cm
  REAL(CMISSDP), PARAMETER :: Length=0.0_CMISSDP !cm
  
  REAL(CMISSDP) :: SolidDensity ! kg / cm^3
  REAL(CMISSDP), ALLOCATABLE :: Gravity(:)
  
  INTEGER(CMISSIntg) :: NumberOfFluidNodes=73
  INTEGER(CMISSIntg) :: NumberOfFluidElements=13
  INTEGER(CMISSIntg) :: NumberOfSolidNodes=15
  INTEGER(CMISSIntg) :: NumberOfSolidElements=2
  INTEGER(CMISSIntg) :: NumberOfInterfaceNodes=11
  INTEGER(CMISSIntg) :: NumberOfInterfaceElements=5
  
  INTEGER(CMISSIntg) :: TimeDependenceTypeCheck(2)=(/0,0/)
  
  INTEGER(CMISSIntg), ALLOCATABLE :: ElementIndices(:),NodeIndices(:),NodeIndexFMappings(:),NodeIndexSMappings(:), &
    & SolidElements(:), &
    & SolidElementSDNodes(:),SolidElementHPNodes(:),SolidElementSDNodeMappings(:),SolidElementHPNodeMappings(:),FluidElements(:), &
    & FluidElementSVNodes(:),FluidElementPNodes(:),FluidElementSVNodeMappings(:),FluidElementPNodeMappings(:), &
    & MovingMeshElements(:), &
    & MovingMeshElementNodes(:),MovingMeshBoundaryElementNodes(:),MovingMeshFixedNodes(:),MovingMeshMovingNodes(:), &
    & MovingMeshBoundaryElementNodeMappings(:),MovingMeshFixedNodeMappings(:),MovingMeshMovingNodeMappings(:), &
    & InterfaceElements(:),InterfaceInterfaceElementNodes(:),SolidInterfaceElements(:),SolidInterfaceElementNodes(:), &
    & SolidInterfaceElementNodeMappings(:),FluidInterfaceElements(:),FluidInterfaceElementNodes(:), &
    & FluidInterfaceElementNodeMappings(:), &
    & InterfaceNodeNumbers(:),SolidInterfaceNodeNumbers(:),FluidInterfaceNodeNumbers(:),SolidInterfaceNodeNumberMappings(:), &
    & FluidInterfaceNodeNumberMappings(:),DisplacementBoundaryNodes(:),DisplacementBoundaryNodeMappings(:), &
    & VelocityInletBoundaryNodes(:),VelocityInletBoundaryNodeMappings(:), &
    & VelocityBoundaryNodes(:),VelocityBoundaryNodeMappings(:),PressureBoundaryNodes(:), &
    & PressureBoundaryNodeMappings(:)
  
  REAL(CMISSDP), ALLOCATABLE :: SolidGeometry(:),FluidGeometry(:),InterfaceGeometry(:),InterfaceXiPosition(:), &
    & SolidInterfaceXiPosition(:),FluidInterfaceXiPosition(:),DisplacementBoundaryCondition(:),VelocityBoundaryCondition(:), &
    & PressureBoundaryCondition(:)
    
    
    
    
  ! 2D plate variables TODO Share variables for different geometries
  INTEGER(CMISSIntg), ALLOCATABLE :: SolidNodeNumbers(:),FluidNodeNumbers(:),SolidElementNodes(:,:),FluidElementNodes(:,:), &
    & InterfaceElementNodes(:,:),InterfaceSolidElementNodes(:,:),InterfaceFluidElementNodes(:,:),ConnectedInterfaceNodes(:,:), &
    & InterfaceNodeNumbersForGeometry(:), &
    & NoDisplacementNodes(:),FixedNodes(:),FixedZNodes(:),MovedNodes(:),MovedYNodes(:),InletNodes(:),NoSlipNodes(:), &
    & OutletNodes(:),SlipNodes(:),InterfaceInterfaceNodeInformationNE(:,:),SolidInterfaceNodeInformationNE(:,:), &
    & FluidInterfaceNodeInformationNE(:,:),LagrangeNodes(:)
  
  REAL(CMISSDP), ALLOCATABLE :: SolidGeometryX(:),SolidGeometryY(:),SolidGeometryZ(:), &
    & FluidGeometryX(:),FluidGeometryY(:),FluidGeometryZ(:), &
    & InterfaceGeometryX(:),InterfaceGeometryY(:),InterfaceGeometryZ(:), &
    & SolidXi1(:,:),SolidXi2(:,:),SolidXi3(:,:),FluidXi1(:,:),FluidXi2(:,:),FluidXi3(:,:), &
    & InterfaceInterfaceNodeInformationXi(:,:), SolidInterfaceNodeInformationXi(:,:),FluidInterfaceNodeInformationXi(:,:)
  
  
                                                          
  INTEGER(CMISSIntg), ALLOCATABLE :: IndexArray(:)
  
  
  REAL :: t(2)=(/0.0,0.0/)
  REAL :: e
  

  INTEGER(CMISSIntg), PARAMETER :: SolidCoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FluidCoordinateSystemUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: InterfaceCoordinateSystemUserNumber=4
  
  INTEGER(CMISSIntg), PARAMETER :: SolidRegionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: FluidRegionUserNumber=7
  
!  INTEGER(CMISSIntg), PARAMETER :: Basis1UserNumber=9
!  INTEGER(CMISSIntg), PARAMETER :: Basis2UserNumber=10
!  INTEGER(CMISSIntg), PARAMETER :: Basis3UserNumber=11
!  INTEGER(CMISSIntg), PARAMETER :: PressureBasis1UserNumber=52
!  INTEGER(CMISSIntg), PARAMETER :: PressureBasis2UserNumber=53
!  INTEGER(CMISSIntg), PARAMETER :: PressureBasis3UserNumber=54
  INTEGER(CMISSIntg), PARAMETER :: InterfaceBasisUserNumber=12
  
  INTEGER(CMISSIntg), PARAMETER :: SolidMeshUserNumber=19
  INTEGER(CMISSIntg), PARAMETER :: FluidMeshUserNumber=20
  INTEGER(CMISSIntg), PARAMETER :: InterfaceMeshUserNumber=22
  
  INTEGER(CMISSIntg), PARAMETER :: SolidDecompositionUserNumber=24
  INTEGER(CMISSIntg), PARAMETER :: FluidDecompositionUserNumber=25
  INTEGER(CMISSIntg), PARAMETER :: InterfaceDecompositionUserNumber=27
  
  INTEGER(CMISSIntg), PARAMETER :: SolidGeometricFieldUserNumber=29
  INTEGER(CMISSIntg), PARAMETER :: FluidGeometricFieldUserNumber=30
  INTEGER(CMISSIntg), PARAMETER :: InterfaceGeometricFieldUserNumber=32
  
  INTEGER(CMISSIntg), PARAMETER :: SolidEquationsSetUserNumber=34
  INTEGER(CMISSIntg), PARAMETER :: FluidEquationsSetUserNumber=35
  
  
  
  
  
  
  
  
  
  
  
  
  
  INTEGER(CMISSIntg), PARAMETER :: DependentField1UserNumber=37
  INTEGER(CMISSIntg), PARAMETER :: DependentField2UserNumber=38
  
  INTEGER(CMISSIntg), PARAMETER :: Interface1UserNumber=40
  
  INTEGER(CMISSIntg), PARAMETER :: InterfaceCondition1UserNumber=42
  
  INTEGER(CMISSIntg), PARAMETER :: LagrangeField1UserNumber=44
  
  INTEGER(CMISSIntg), PARAMETER :: CoupledProblemUserNumber=46
  
  INTEGER(CMISSIntg), PARAMETER :: InterfaceMappingBasis1UserNumber=47
  
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetField1UserNumber=49
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetField2UserNumber=50
  
  !======
  INTEGER(CMISSIntg), PARAMETER :: BasisSpaceSolidUserNumber=100
  INTEGER(CMISSIntg), PARAMETER :: BasisDisplacementUserNumber=101
  INTEGER(CMISSIntg), PARAMETER :: BasisHydrostaticPressureUserNumber=102
  INTEGER(CMISSIntg), PARAMETER :: BasisSpaceFluidUserNumber=103
  INTEGER(CMISSIntg), PARAMETER :: BasisVelocityUserNumber=104
  INTEGER(CMISSIntg), PARAMETER :: BasisPressureUserNumber=105
  INTEGER(CMISSIntg), PARAMETER :: FibreFieldUserNumber=106
  INTEGER(CMISSIntg), PARAMETER :: MovingMeshUserNumber=107
  INTEGER(CMISSIntg), PARAMETER :: MovingMeshFieldUserNumber=108
  INTEGER(CMISSIntg), PARAMETER :: MovingMeshEquationsSetUserNumber=109
  INTEGER(CMISSIntg), PARAMETER :: MaterialField1UserNumber=110
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldMovingMeshUserNumber=111
  INTEGER(CMISSIntg), PARAMETER :: SourceField1UserNumber=112
  
  INTEGER(CMISSIntg), PARAMETER :: LinearSolverMovingMeshUserNumber=2!
  INTEGER(CMISSIntg), PARAMETER :: DynamicSolverUserNumber=1!1
  
  INTEGER(CMISSIntg), PARAMETER :: MaterialField2UserNumberMu=1 !???
  INTEGER(CMISSIntg), PARAMETER :: MaterialField2UserNumberRho=2 !???
  INTEGER(CMISSIntg), PARAMETER :: MaterialFieldMovingMeshUserNumberK=1
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldMovingMeshUserNumberK=1
  INTEGER(CMISSIntg), PARAMETER :: MaterialField2UserNumber=117
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldMovingMeshUserNumber=118
  INTEGER(CMISSIntg), PARAMETER :: MaterialFieldMovingMeshUserNumber=119
  INTEGER(CMISSIntg), PARAMETER :: IndependentField2UserNumber=120
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldMovingMeshUserNumber=121
  INTEGER(CMISSIntg), PARAMETER :: LinearSolverMovingMeshEquationsUserNumber=122
  
  
  
  
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS,ARGUMENT_Length,STATUS
  INTEGER(CMISSIntg) :: NumberGlobalElementsX,NumberGlobalElementsY_Fluid,NumberGlobalElementsY_Solid, &
    & NumberGlobalElementsZ,InterpolationTypeInterface,NumberOfGaussXi,NUMBER_OF_NODE_XI,NumberOfDimensions,component_idx, &
    & NumberOfGaussXiSpace,NumberOfGaussXiVelocity,NumberOfGaussXiPressure,arraySize
  CHARACTER(LEN=255) :: COMMAND_ARGUMENT

  INTEGER(CMISSIntg) :: EquationsSet1Index=1
  INTEGER(CMISSIntg) :: EquationsSet2Index=2
  INTEGER(CMISSIntg) :: EquationsSet3Index
  INTEGER(CMISSIntg) :: LinearSolverMovingMeshEquationsIndex=1
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain,LastNodeDomain
  INTEGER(CMISSIntg) :: InterfaceCondition1Index=1
  INTEGER(CMISSIntg) :: InterfaceCondition2Index
  INTEGER(CMISSIntg) :: Mesh1Index=1
  INTEGER(CMISSIntg) :: Mesh2Index=2
  INTEGER(CMISSIntg) :: Mesh3Index
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: x_element_idx,z_element_idx,mesh_local_x_node,mesh_local_z_node,PressureMeshComponent
  INTEGER(CMISSIntg) :: Mesh2ComponentNumberVelocity,ComponentNumber,FirstMesh1ElementNumber
  
  INTEGER(CMISSIntg) :: OutputFrequency
  INTEGER(CMISSIntg) :: DynamicSolver_OutputType
  INTEGER(CMISSIntg) :: NonlinearSolver_OutputType
  INTEGER(CMISSIntg) :: LinearSolver_OutputType
  INTEGER(CMISSIntg) :: LinearSolverMovingMesh_OutputType
  INTEGER(CMISSIntg) :: MaximumIterations,MaxFunctionEvaluations
  INTEGER(CMISSIntg) :: RestartValue
  INTEGER(CMISSIntg) :: MMCondition=0
  
  INTEGER(CMISSIntg) :: EquationsNavierStokesOutput,InterfaceMeshComponentNumber,InterpolationTypeDisplacement, &
    & InterpolationTypeHydrostaticPressure,InterpolationTypePressure,InterpolationTypeSpace,InterpolationTypeVelocity, &
    & MovingMeshEquationsSetIndex,Mesh1ComponentNumberDisplacement,Mesh1ComponentNumberHydrostaticPressure, &
    & Mesh1ComponentNumberSpace, &
    & NodeDomain,Mesh2ComponentNumberPressure,Mesh2ComponentNumberSpace,MeshNumberOfComponents,NumberGlobalElementsY
  
  
  
  
  
  
  ! LOOP INTEGERS
  INTEGER(CMISSIntg) :: NodeNumber,S,F,P,ElementIndex,NodeIndex,ComponentIndex,MM,MaterialSpecification,LocalNodeIndex
  INTEGER(CMISSIntg), PARAMETER :: Blood=1
  INTEGER(CMISSIntg), PARAMETER :: Plate2D=2
  INTEGER(CMISSIntg), PARAMETER :: Plate3D=3
  
  
  
  !check below for units
  REAL(CMISSDP) :: XI2(2),XI3(3)
  REAL(CMISSDP) :: MovingMeshParameterK
  REAL(CMISSDP) :: FluidDynamicViscosity
  REAL(CMISSDP) :: FluidDensity
  REAL(CMISSDP) :: YoungsModulus
  REAL(CMISSDP) :: PoissonsRatio
  REAL(CMISSDP) :: ShearModulus
  REAL(CMISSDP) :: MooneyRivlin1
  REAL(CMISSDP) :: MooneyRivlin2

  REAL(CMISSDP) :: InitialFieldNavierStokes(3)
  REAL(CMISSDP) :: InitialFieldMovingMesh(3)
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_NAVIER_STOKES(3)
  REAL(CMISSDP) :: BOUNDARY_CONDITIONS_MOVING_MESH(3)
  REAL(CMISSDP) :: DivergenceTolerance
  REAL(CMISSDP) :: RelativeTolerance
  REAL(CMISSDP) :: AbsoluteTolerance
  REAL(CMISSDP) :: LinesearchAlpha
  REAL(CMISSDP) :: VALUE

  REAL(CMISSDP) :: StartTime
  REAL(CMISSDP) :: StopTime
  REAL(CMISSDP) :: DynamicSolver_Theta
  REAL(CMISSDP) :: TimeStepSize
  
  REAL(CMISSDP) :: Frac

  LOGICAL :: EXPORT_FIELD_IO
  LOGICAL :: LinearSolverMovingMesh_DirectFlag
  LOGICAL :: FIXED_WALL_NODES_NAVIER_STOKES_FLAG
  LOGICAL :: MOVED_WALL_NODES_NAVIER_STOKES_FLAG
  LOGICAL :: INLET_WALL_NODES_NAVIER_STOKES_FLAG
  LOGICAL :: FIXED_WALL_NODES_MOVING_MESH_FLAG
  LOGICAL :: MOVED_WALL_NODES_MOVING_MESH_FLAG
  LOGICAL :: FileReadDiagnostics=.FALSE.
  LOGICAL :: ExampleFileProgressDiagnostics=.TRUE.
  LOGICAL :: GeometryCheck=.FALSE.
  

  !CMISS variables
  
  TYPE(CMISSBasisType) :: BasisSpaceSolid,BasisDisplacement,BasisHydrostaticPressure
  TYPE(CMISSBasisType) :: BasisSpaceFluid,BasisVelocity,BasisPressure
  TYPE(CMISSBasisType) :: InterfaceBasis1,InterfaceMappingBasis1
  
  TYPE(CMISSNodesType) :: SolidNodes,HydrostaticPressureNodes
  TYPE(CMISSNodesType) :: FluidNodes,PressureNodes,InterfaceNodes,MovingMeshNodes
  
  TYPE(CMISSMeshElementsType) :: Mesh1ElementsSpace,Mesh1ElementsDisplacement,Mesh1ElementsHydrostaticPressure
  TYPE(CMISSMeshElementsType) :: Mesh2ElementsSpace,Mesh2ElementsVelocity,Mesh2ElementsPressure
  TYPE(CMISSMeshElementsType) :: InterfaceMeshElements
  
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions,BoundaryConditionsMovingMesh
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem1,CoordinateSystem2,CoordinateSystem3,CoordinateSystemInterface1, &
    & CoordinateSystemInterface2,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition1,Decomposition2,Decomposition3,InterfaceDecomposition1,InterfaceDecomposition2
  TYPE(CMISSEquationsType) :: Equations1,Equations2,Equations3,EquationsMovingMesh
  TYPE(CMISSEquationsSetType) :: EquationsSet1,EquationsSet2,EquationsSet3,MovingMeshEquationsSet
  TYPE(CMISSFieldType) :: GeometricField1,GeometricField2,GeometricField3,InterfaceGeometricField1,InterfaceGeometricField2, &
    & DependentField1,DependentField2,DependentField3,LagrangeField1,LagrangeField2,EquationsSetField1,EquationsSetField2, &
    & EquationsSetField3,DependentFieldMovingMesh,MaterialFieldMovingMesh,IndependentField2,IndependentFieldMovingMesh
  
  TYPE(CMISSFieldType) :: FibreField,MovingMeshField,MaterialField1,MaterialField2,EquationsSetFieldMovingMesh,SourceField1
  
  TYPE(CMISSFieldsType) :: Fields1,Fields2,FieldsI,InterfaceFields1,InterfaceFields2
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh1,GeneratedMesh2,GeneratedMesh3,InterfaceGeneratedMesh1,InterfaceGeneratedMesh2
  TYPE(CMISSInterfaceType) :: Interface1,Interface2
  TYPE(CMISSInterfaceConditionType) :: InterfaceCondition1,InterfaceCondition2
  TYPE(CMISSInterfaceEquationsType) :: InterfaceEquations1,InterfaceEquations2
  TYPE(CMISSInterfaceMeshConnectivityType) :: InterfaceMeshConnectivity1,InterfaceMeshConnectivity2
  TYPE(CMISSMeshType) :: Mesh1,Mesh2,Mesh3,InterfaceMesh1,InterfaceMesh2
  TYPE(CMISSNodesType) :: Nodes
  TYPE(CMISSProblemType) :: CoupledProblem
  TYPE(CMISSControlLoopType) :: ControlLoop,MovingMeshControlLoop,FENSControlLoop
  TYPE(CMISSRegionType) :: Region1,Region2,Region3,WorldRegion
  TYPE(CMISSSolverType) :: DynamicSolver
  TYPE(CMISSSolverType) :: NonlinearSolver
  TYPE(CMISSSolverType) :: LinearSolver
  TYPE(CMISSSolverType) :: LinearSolverMovingMesh
  TYPE(CMISSSolverType) :: CoupledSolver
  TYPE(CMISSSolverEquationsType) :: CoupledSolverEquations,LinearSolverMovingMeshEquations
  
  
  
  !
  !=================================================================================================================================
  !
  !
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
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

  e=etime(t)
  
  !
  !=================================================================================================================================
  !
  ! C O M M A N D L I N E _ A R G U M E N T S
  
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  IF(NUMBER_OF_ARGUMENTS == 4) THEN
    CALL GET_COMMAND_ARGUMENT(1,COMMAND_ARGUMENT,ARGUMENT_Length,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_Length),*) NumberGlobalElementsX
    IF(NumberGlobalElementsX<=0) CALL HANDLE_ERROR("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(2,COMMAND_ARGUMENT,ARGUMENT_Length,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_Length),*) NumberGlobalElementsY
    IF(NumberGlobalElementsY<=0) CALL HANDLE_ERROR("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,COMMAND_ARGUMENT,ARGUMENT_Length,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 3.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_Length),*) NumberGlobalElementsZ
    IF(NumberGlobalElementsY<0) CALL HANDLE_ERROR("Invalid number of Z elements.")
    CALL GET_COMMAND_ARGUMENT(4,COMMAND_ARGUMENT,ARGUMENT_Length,STATUS)
    IF(STATUS>0) CALL HANDLE_ERROR("Error for command argument 4.")
    READ(COMMAND_ARGUMENT(1:ARGUMENT_Length),*) InterpolationTypeInterface
    IF(InterpolationTypeInterface<=0.OR.InterpolationTypeInterface>4) CALL HANDLE_ERROR("Invalid Interpolation specification.")
  ELSE IF(NUMBER_OF_ARGUMENTS == 0) THEN
    NumberGlobalElementsZ=0 !TODO get rid of this variable
    InterpolationTypeDisplacement=2
    InterpolationTypeSpace=2
    InterpolationTypeVelocity=2
    InterpolationTypePressure=1
    InterpolationTypeInterface=2
  ELSE
    CALL HANDLE_ERROR("Invalid number of arguments.")
  ENDIF
  !
  !=================================================================================================================================
  !
  ! I N I T I A L I S E _ O P E N C M I S S / V A R I A B L E S

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Set error handling mode
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
 
  !Set diganostics for testing
  !CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["SOLVER_MAPPING_CALCULATE         ", &
  !  & "SOLVER_MATRIX_STRUCTURE_CALCULATE"],Err)
  
  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)
  !Set initial values
  InitialFieldNavierStokes(1)=0.1_CMISSDP
  InitialFieldNavierStokes(2)=0.0_CMISSDP
  InitialFieldNavierStokes(3)=0.0_CMISSDP
  InitialFieldMovingMesh(1)=0.0_CMISSDP
  InitialFieldMovingMesh(2)=0.0_CMISSDP
  InitialFieldMovingMesh(3)=0.0_CMISSDP
  !Set output parameter
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LinearSolverMovingMesh_OutputType=CMISS_SOLVER_SOLVER_OUTPUT
  DynamicSolver_OutputType=CMISS_SOLVER_SOLVER_OUTPUT
  LinearSolver_OutputType=CMISS_SOLVER_SOLVER_OUTPUT
  NonlinearSolver_OutputType=CMISS_SOLVER_SOLVER_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EquationsNavierStokesOutput=CMISS_EQUATIONS_NO_OUTPUT
  !Set result output parameter
  OutputFrequency=1
  
  
  
  
  MaterialSpecification=Plate3D
  
  SELECT CASE(MaterialSpecification)
  CASE(Blood)
    !Set solver parameters
    RelativeTolerance=1.0E-12_CMISSDP !default: 1.0E-05_CMISSDP
    AbsoluteTolerance=1.0E-12_CMISSDP !default: 1.0E-10_CMISSDP
    DivergenceTolerance=1.0E30 !default: 1.0E5
    MaximumIterations=100000000 !default: 100000
    MaxFunctionEvaluations=100000
    RestartValue=30 !default: 30
    LinesearchAlpha=1.0_CMISSDP
    MovingMeshParameterK=1.0
    !Set time parameter
    StartTime=0.0_CMISSDP
    StopTime=2.0_CMISSDP
    TimeStepSize=0.125_CMISSDP
    DynamicSolver_Theta=1.0_CMISSDP
    !BLOOD
    FluidDynamicViscosity=4.0E-5_CMISSDP !kg/(cm s)
    FluidDensity=1.050E-3_CMISSDP ! kg / cm^-3
    !ARTERIAL WALL
    SolidDensity=1.160E-3_CMISSDP ! kg / cm^3
    !Young's modulus E
    YoungsModulus=3.0E4_CMISSDP ! 3.0E4_CMISSDP kg / (cm s^2)
    !Poisson's ratio
    PoissonsRatio=0.45
    !Homogenous, isotropic material G=E/(2*(1+poissonsRatio))
    ShearModulus=YoungsModulus/(2.0_CMISSDP*(1+PoissonsRatio)) ! kg / (cm s^2)
    !Neo-Hookean material: c1=G/2 c2=0
    MooneyRivlin1=0.5_CMISSDP*ShearModulus ! kg / (cm s^2)
    MooneyRivlin2=0.0_CMISSDP !
    !==========================================================================
    ElementIndices=(/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
    NodeIndices=(/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28, &
      & 29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53, &
      & 54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77/)
                                                      
    NodeIndexFMappings=(/1,25,2,26,3,4,27,5,28,6,29,30,31,32,33,34,35,36,37,38,7,39,8,40,9, &
      & 10,41,11,42,12,43,44,45,46,47,48,49,50,51,52,13,53,14,54,15,55,16,56,17,57,18,58,59, &
      & 60,61,62,63,64,65,66,67,68,19,69,20,70,21,71,22,72,23,73,24/)
    NodeIndexSMappings=(/1,7,2,8,9,10,3,11,4,12,13,14,5,15,6/)
    SolidElements=(/1,2/)
    SolidElementSDNodes=(/1,2,3,4,5,6,7,8,9,7,8,9,10,11,12,13,14,15/)
    SolidElementHPNodes=(/1,2,3,4,3,4,5,6/)
    SolidElementSDNodeMappings=(/1,7,2,8,9,10,3,11,4,3,11,4,12,13,14,5,15,6/)
    SolidElementHPNodeMappings=(/1,2,3,4,3,4,5,6/)
    
    FluidElements=(/1,2,3,4,5,6,7,8,9,10,11,12,13/)
    FluidElementSVNodes=(/1,2,3,11,12,13,21,22,23, &
                                                              & 3,4,5,13,14,15,23,24,25, &
                                                              & 6,7,8,16,17,18,26,27,28, &
                                                              & 8,9,10,18,19,20,28,29,30, &
                                                              & 21,22,23,31,32,33,41,42,43, &
                                                              & 23,24,25,33,34,35,43,44,45, &
                                                              & 26,27,28,36,37,38,47,48,49, &
                                                              & 28,29,30,38,39,40,49,50,51, &
                                                              & 41,42,43,52,53,54,63,64,65, &
                                                              & 43,44,45,54,55,56,65,66,67, &
                                                              & 45,46,47,56,57,58,67,68,69, &
                                                              & 47,48,49,58,59,60,69,70,71, &
                                                              & 49,50,51,60,61,62,71,72,73/)
    FluidElementPNodes=(/1,2,7,8,2,3,8,9,4,5,10,11,5,6,11,12,7,8,13,14,8,9,14,15,10,11,16,17, &
      & 11,12,17,18,13,14,19,20,14,15,20,21,15,16,21,22,16,17,22,23,17,18,23,24/)
    FluidElementSVNodeMappings=(/1,25,2,29,30,31,7,39,8,2,26,3,31,32,33,8,40,9, &
      & 4,27,5,34,35,36,10,41,11,5,28,6,36,37,38,11,41,12,7,39,8,43,44,45,13,53,14,8,40,9,45,46,47,14,54,15, &
      & 10,41,11,48,49,50,16,56,17,11,42,12,50,51,52,17,57,18,13,53,14,58,59,60,19,69,20,14,54,15,60,61,62,20,70,21, &
      & 15,55,16,62,63,64,21,71,22,16,56,17,64,65,66,22,72,23,17,57,18,66,67,68,23,73,24/)
    FluidElementPNodeMappings=(/1,2,7,8,2,3,8,9,4,5,10,11,5,6,11,12,7,8,13,14,8,9,14,15, &
      & 10,11,16,17,11,12,17,18,13,14,19,20,14,15,20,21,15,16,21,22,16,17,22,23,17,18,23,24/)
    MovingMeshElements=FluidElements(:)
    MovingMeshElementNodes=FluidElementSVNodes(:)
    MovingMeshBoundaryElementNodes=(/1,2,3,4,5,15,25,35,45,46,47,36,26,16,6,7,8,9,10,20,30, &
      & 40,51,62,73,72,71,70,69,68,67,66,65,64,63,52,41,31,21,11/)
    MovingMeshFixedNodes=(/1,2,3,4,5,6,7,8,9,10,20,30,40,51,62,73,72,71,70,69,68,67,66,65,64, &
      & 63,52,41,31,21,11/)
    MovingMeshMovingNodes=(/15,25,35,45,46,47,36,26,16/)
    MovingMeshBoundaryElementNodeMappings=(/1,25,2,26,3,33,9,47,15,55,16,48,10,34, &
      & 4,27,5,28,6,38,12,52,18,68,24,73,23,72,22,71,21,70,20,69,19,58,13,43,7,29/)
    MovingMeshFixedNodeMappings=(/1,25,2,26,3,4,27,5,28,6,38,12,52,18,68,24,73,23,72,22,71,21, &
      & 70,20,69,19,58,13,43,7,29/)
    MovingMeshMovingNodeMappings=(/33,9,47,15,55,16,48,10,34/)
    
    InterfaceElements=(/1,2,3,4,5/)
    InterfaceInterfaceElementNodes=(/1,2,3,3,4,5,5,6,7,7,8,9,9,10,11/)
    SolidInterfaceElements=(/1,2,2,2,1/)
    SolidInterfaceElementNodes=(/1,4,7,7,10,13,13,14,15,15,12,9,9,6,3/)
    SolidInterfaceElementNodeMappings=(/1,8,3,3,12,5,5,15,6,6,14,4,4,10,2/)
    FluidInterfaceElements=(/2,6,11,7,3/)
    FluidInterfaceElementNodes=(/5,15,25,25,35,45,45,46,47,47,36,26,26,16,6/)
    FluidInterfaceElementNodeMappings=(/3,33,9,9,47,15,15,55,16,16,48,10,10,34,4/)
    
    InterfaceNodeNumbers=(/1,2,3,4,5,6,7,8,9,10,11/)
    SolidInterfaceNodeNumbers=(/1,4,7,10,13,14,15,12,9,6,3/)
    FluidInterfaceNodeNumbers=(/5,15,25,35,45,46,47,36,26,16,6/)
    SolidInterfaceNodeNumberMappings=(/1,8,3,12,5,15,6,14,4,10,2/)
    FluidInterfaceNodeNumberMappings=(/3,33,9,47,15,55,16,48,10,34,4/)
    
    ! G E O M E T R I C _ P O S I T I O N S
    SolidGeometry=(/2.0_CMISSDP,2.5_CMISSDP,3.0_CMISSDP,2.0_CMISSDP,2.5_CMISSDP,3.0_CMISSDP, &
      & 2.0_CMISSDP,2.5_CMISSDP,3.0_CMISSDP,2.0_CMISSDP,2.5_CMISSDP,3.0_CMISSDP,2.0_CMISSDP,2.5_CMISSDP,3.0_CMISSDP, &
      !2nd component
      & 0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.5_CMISSDP,0.5_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.5_CMISSDP, &
      & 1.5_CMISSDP,1.5_CMISSDP,2.0_CMISSDP,2.0_CMISSDP,2.0_CMISSDP/)
    FluidGeometry=(/0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,1.5_CMISSDP,2.0_CMISSDP,3.0_CMISSDP, &
      & 3.5_CMISSDP,4.0_CMISSDP,4.5_CMISSDP,5.0_CMISSDP,0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,1.5_CMISSDP,2.0_CMISSDP,3.0_CMISSDP, &
      & 3.5_CMISSDP,4.0_CMISSDP,4.5_CMISSDP,5.0_CMISSDP,0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,1.5_CMISSDP,2.0_CMISSDP,3.0_CMISSDP, &
      & 3.5_CMISSDP,4.0_CMISSDP,4.5_CMISSDP,5.0_CMISSDP,0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,1.5_CMISSDP,2.0_CMISSDP,3.0_CMISSDP, &
      & 3.5_CMISSDP,4.0_CMISSDP,4.5_CMISSDP,5.0_CMISSDP,0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,1.5_CMISSDP,2.0_CMISSDP,2.5_CMISSDP, &
      & 3.0_CMISSDP,3.5_CMISSDP,4.0_CMISSDP,4.5_CMISSDP,5.0_CMISSDP,0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,1.5_CMISSDP,2.0_CMISSDP, &
      & 2.5_CMISSDP,3.0_CMISSDP,3.5_CMISSDP,4.0_CMISSDP,4.5_CMISSDP,5.0_CMISSDP,0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,1.5_CMISSDP, &
      & 2.0_CMISSDP,2.5_CMISSDP,3.0_CMISSDP,3.5_CMISSDP,4.0_CMISSDP,4.5_CMISSDP,5.0_CMISSDP, &
      !2nd component
      & 0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP, &
      & 0.5_CMISSDP,0.5_CMISSDP,0.5_CMISSDP,0.5_CMISSDP,0.5_CMISSDP,0.5_CMISSDP,0.5_CMISSDP,0.5_CMISSDP,0.5_CMISSDP,0.5_CMISSDP, &
      & 1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP, &
      & 1.5_CMISSDP,1.5_CMISSDP,1.5_CMISSDP,1.5_CMISSDP,1.5_CMISSDP,1.5_CMISSDP,1.5_CMISSDP,1.5_CMISSDP,1.5_CMISSDP,1.5_CMISSDP, &
      & 2.0_CMISSDP,2.0_CMISSDP,2.0_CMISSDP,2.0_CMISSDP,2.0_CMISSDP,2.0_CMISSDP,2.0_CMISSDP,2.0_CMISSDP,2.0_CMISSDP,2.0_CMISSDP, &
      & 2.0_CMISSDP,2.5_CMISSDP,2.5_CMISSDP,2.5_CMISSDP,2.5_CMISSDP,2.5_CMISSDP,2.5_CMISSDP,2.5_CMISSDP,2.5_CMISSDP,2.5_CMISSDP, &
      & 2.5_CMISSDP,2.5_CMISSDP,3.0_CMISSDP,3.0_CMISSDP,3.0_CMISSDP,3.0_CMISSDP,3.0_CMISSDP,3.0_CMISSDP,3.0_CMISSDP,3.0_CMISSDP, &
      & 3.0_CMISSDP,3.0_CMISSDP,3.0_CMISSDP/)
    InterfaceGeometry=(/2.0_CMISSDP,2.0_CMISSDP,2.0_CMISSDP,2.0_CMISSDP,2.0_CMISSDP,2.5_CMISSDP, &
      & 3.0_CMISSDP,3.0_CMISSDP,3.0_CMISSDP,3.0_CMISSDP,3.0_CMISSDP, &
      !2nd component
      & 0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,1.5_CMISSDP,2.0_CMISSDP,2.0_CMISSDP,2.0_CMISSDP,1.5_CMISSDP,1.0_CMISSDP,0.5_CMISSDP, &
      & 0.0_CMISSDP/)
    
    ! X I _ P O S I T I O N S for MeshConnectivity
    InterfaceXiPosition=(/0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP, &
      & 0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP, &
      & 0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP, &
      & 0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP, &
      & 0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP/)
    SolidInterfaceXiPosition=(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP, &
      & 0.0_CMISSDP,0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP, &
      !xi2
      & 0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP, &
      & 0.5_CMISSDP,0.0_CMISSDP,1.0_CMISSDP,0.5_CMISSDP,0.0_CMISSDP/)
    FluidInterfaceXiPosition=(/1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP,1.0_CMISSDP, &
      & 1.0_CMISSDP,0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP, &
      !xi2
      & 0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,0.0_CMISSDP,0.5_CMISSDP,1.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,1.0_CMISSDP, &
      & 0.5_CMISSDP,0.0_CMISSDP,1.0_CMISSDP,0.5_CMISSDP,0.0_CMISSDP/)
    
    ! B O U N D A R Y _ C O N D I T I O N S
    DisplacementBoundaryNodes=(/1,2,3/)
    DisplacementBoundaryNodeMappings=(/1,7,2/)
    DisplacementBoundaryCondition=(/2.0_CMISSDP,2.5_CMISSDP,3.0_CMISSDP, &
                                                                 & 0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
    VelocityInletBoundaryNodes=(/11,21,31,41,52/)
    VelocityInletBoundaryNodeMappings=(/29,7,43,13,58/)
    VelocityBoundaryNodes=(/1,2,3,4,5,6,7,8,9,10,11,21,31,41,52, &
                                                              & 63,64,65,66,67,68,69,70,71,72,73/)
    VelocityBoundaryNodeMappings=(/1,25,2,26,3,4,27,5,28,6,29,7,43,13,58,19,69,20,70,21,71,22, &
      & 72,23,73,24/)
    VelocityBoundaryCondition=(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP, &
                                                            & 0.0_CMISSDP,0.0_CMISSDP, &
                                                            & 0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,1.0_CMISSDP/30,2.0_CMISSDP/30, &
                                                            & 3.0_CMISSDP/30, &
                                                            & 4.0_CMISSDP/30,5.0_CMISSDP/30,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP, &
                                                            & 0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP, &
                                                            & 0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP, &
                                                            !2nd component
                                                            & 0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP, &
                                                            & 0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP, &
                                                            & 0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP, &
                                                            & 0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP, &
                                                            & 0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP, &
                                                            & 0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
    PressureBoundaryNodes=(/12,18/)!VelocityBoundaryNodes(:)
    PressureBoundaryNodeMappings=(/12,18/)!VelocityBoundaryNodes(:)
    PressureBoundaryCondition=(/0.0_CMISSDP,0.0_CMISSDP/)
    NumberOfDimensions=2
    
    ALLOCATE(Gravity(NumberOfDimensions))
    Gravity(:)=[0.0_CMISSDP,-981.0_CMISSDP] !in cm s^-2
    
    !===============================================================================================================================
  CASE(Plate2D,Plate3D)!NOTE: USE OF SI UNITS
    !Set solver parameters
    RelativeTolerance=1.0E-05_CMISSDP !default: 1.0E-05_CMISSDP
    AbsoluteTolerance=1.0E-10_CMISSDP !default: 1.0E-10_CMISSDP
    DivergenceTolerance=1.0E5 !default: 1.0E5
    MaximumIterations=100000000 !default: 100000
    MaxFunctionEvaluations=100000
    RestartValue=30 !default: 30
    LinesearchAlpha=1.0_CMISSDP
    MovingMeshParameterK=1.0
    !Set time parameter
    StartTime=0.0_CMISSDP
    StopTime=2.0_CMISSDP
    TimeStepSize=0.05_CMISSDP
    DynamicSolver_Theta=1.0_CMISSDP
    
    !Fluid
    FluidDynamicViscosity=0.1_CMISSDP !kg/(m s)
    FluidDensity=1.0E3_CMISSDP ! kg / m^-3
    
    !Solid
    SolidDensity=1.0E6_CMISSDP ! kg / m^3
    !Young's modulus E
    YoungsModulus=7.0E10_CMISSDP ! N / m^2
    !Poisson's ratio
    PoissonsRatio=0.3
    !Homogenous, isotropic material G=E/(2*(1+poissonsRatio))
    ShearModulus=YoungsModulus/(2.0_CMISSDP*(1+PoissonsRatio)) ! N / m^2
    !Neo-Hookean material: c1=G/2 c2=0
    MooneyRivlin1=0.5_CMISSDP*ShearModulus ! N / m^2
    MooneyRivlin2=0.0_CMISSDP !
    
    !Set geometric dimension n gravity
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      NumberOfDimensions=2
      ALLOCATE(Gravity(NumberOfDimensions))
      Gravity(:)=[0.0_CMISSDP,-9.81_CMISSDP] !in m s^-2
    CASE(Plate3D)
      NumberOfDimensions=3
      ALLOCATE(Gravity(NumberOfDimensions))
      Gravity(:)=[0.0_CMISSDP,-9.81_CMISSDP,0.0_CMISSDP] !in m s^-2
    END SELECT
    
    
    arraySize=0
    PRINT *, 'Reading geometry from files..'
    !Read node numbers
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DSolidNodes.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DSolidNodes.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(SolidNodeNumbers(arraySize))
    DO S=1,arraySize
      READ(1,*) SolidNodeNumbers(S)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Solid nodes: ',SolidNodeNumbers
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DFluidNodes.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DFluidNodes.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(FluidNodeNumbers(arraySize))
    DO S=1,arraySize
      READ(1,*) FluidNodeNumbers(S)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Fluid nodes: ',FluidNodeNumbers
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DInterfaceNodes.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DInterfaceNodes.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(InterfaceNodeNumbersForGeometry(arraySize))                                   !TODO TODO TODO InterfaceNodeNumbers is 1:22 not mixed!!
    DO S=1,arraySize
      READ(1,*) InterfaceNodeNumbersForGeometry(S)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Interface nodes: ',InterfaceNodeNumbersForGeometry
    
    !Read geometry
    IF(NumberOfDimensions==3) THEN
      OPEN (unit=1, file='Plate3DSolidX.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(SolidGeometryX(arraySize))
      DO S=1,arraySize
        READ(1,*) SolidGeometryX(S)
      ENDDO
      arraySize=0
      CLOSE(1)
      IF(FileReadDiagnostics) PRINT *, 'Solid x coordinates: ',SolidGeometryX
    ENDIF
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DSolidY.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DSolidY.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(SolidGeometryY(arraySize))
    DO S=1,arraySize
      READ(1,*) SolidGeometryY(S)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Solid y coordinates: ',SolidGeometryY
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DSolidZ.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DSolidZ.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(SolidGeometryZ(arraySize))
    DO S=1,arraySize
      READ(1,*) SolidGeometryZ(S)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Solid z coordinates: ',SolidGeometryZ
    
    IF(NumberOfDimensions==3) THEN
      OPEN (unit=1, file='Plate3DFluidX.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(FluidGeometryX(arraySize))
      DO S=1,arraySize
        READ(1,*) FluidGeometryX(S)
      ENDDO
      arraySize=0
      CLOSE(1)
      IF(FileReadDiagnostics) PRINT *, 'Fluid x coordinates: ',FluidGeometryX
    ENDIF
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DFluidY.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DFluidY.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(FluidGeometryY(arraySize))
    DO S=1,arraySize
      READ(1,*) FluidGeometryY(S)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Fluid y coordinates: ',FluidGeometryY
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DFluidZ.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DFluidZ.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(FluidGeometryZ(arraySize))
    DO S=1,arraySize
      READ(1,*) FluidGeometryZ(S)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Fluid z coordinates: ',FluidGeometryZ
    
    IF(NumberOfDimensions==3) THEN
      OPEN (unit=1, file='Plate3DInterfaceX.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(InterfaceGeometryX(arraySize))
      DO S=1,arraySize
        READ(1,*) InterfaceGeometryX(S)
      ENDDO
      arraySize=0
      CLOSE(1)
      IF(FileReadDiagnostics) PRINT *, 'Interface x coordinates: ',InterfaceGeometryX
    ENDIF
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DInterfaceY.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DInterfaceY.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(InterfaceGeometryY(arraySize))
    DO S=1,arraySize
      READ(1,*) InterfaceGeometryY(S)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Interface y coordinates: ',InterfaceGeometryY
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DInterfaceZ.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DInterfaceZ.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(InterfaceGeometryZ(arraySize))
    DO S=1,arraySize
      READ(1,*) InterfaceGeometryZ(S)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Interface z coordinates: ',InterfaceGeometryZ
    
    !Read element information
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DSolidElementNodeNumbers.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DSolidElementNodeNumbers.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(SolidElementNodes(arraySize,3**NumberOfDimensions))
    DO S=1,arraySize
      READ(1,*) SolidElementNodes(S,:)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Solid element node numbers: ',SolidElementNodes
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DFluidElementNodeNumbers.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DFluidElementNodeNumbers.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(FluidElementNodes(arraySize,3**NumberOfDimensions))
    DO S=1,arraySize
      READ(1,*) FluidElementNodes(S,:)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Fluid element node numbers: ',FluidElementNodes
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DInterfaceElementNodeNumbers.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DInterfaceElementNodeNumbers.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(InterfaceElementNodes(arraySize,3**(NumberOfDimensions-1)))
    DO S=1,arraySize
      READ(1,*) InterfaceElementNodes(S,:)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Interface element node numbers: ',InterfaceElementNodes
    
    !Read interface nodes for solid/fluid/interface
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      !Read xi position's
      OPEN (unit=1, file='Plate2DSolidXi1.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(SolidXi2(arraySize,3))
      DO S=1,arraySize
        READ(1,*) SolidXi2(S,:)
      ENDDO
      arraySize=0
      CLOSE(1)
      IF(FileReadDiagnostics) PRINT *, 'Solid xi2: ',SolidXi2
      
      OPEN (unit=1, file='Plate2DSolidXi2.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(SolidXi3(arraySize,3))
      DO S=1,arraySize
        READ(1,*) SolidXi3(S,:)
      ENDDO
      arraySize=0
      CLOSE(1)
      IF(FileReadDiagnostics) PRINT *, 'Solid xi3: ',SolidXi3
      
      OPEN (unit=1, file='Plate2DFluidXi1.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(FluidXi2(arraySize,3))
      DO S=1,arraySize
        READ(1,*) FluidXi2(S,:)
      ENDDO
      arraySize=0
      CLOSE(1)
      IF(FileReadDiagnostics) PRINT *, 'Fluid xi2: ',FluidXi2
      
      OPEN (unit=1, file='Plate2DFluidXi2.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(FluidXi3(arraySize,3))
      DO S=1,arraySize
        READ(1,*) FluidXi3(S,:)
      ENDDO
      arraySize=0
      CLOSE(1)
      IF(FileReadDiagnostics) PRINT *, 'Fluid xi3: ',FluidXi3
      
      OPEN (unit=1, file='Plate2DInterfaceSolidElements.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(SolidInterfaceElements(arraySize))
      READ(1,*) SolidInterfaceElements(:)
      arraySize=0
      CLOSE(1)
      IF(FileReadDiagnostics) PRINT *, 'Solid interface element numbers: ',SolidInterfaceElements
      
      OPEN (unit=1, file='Plate2DInterfaceFluidElements.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(FluidInterfaceElements(arraySize))
      READ(1,*) FluidInterfaceElements(:)
      arraySize=0
      CLOSE(1)
      IF(FileReadDiagnostics) PRINT *, 'Fluid interface element numbers: ',FluidInterfaceElements
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DSolidInterfaceNodeInformationNodeElement.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(SolidInterfaceNodeInformationNE(arraySize,2))
      DO S=1,arraySize
        READ(1,*) SolidInterfaceNodeInformationNE(S,:)
      ENDDO
      arraySize=0
      CLOSE(1)
      OPEN (unit=1, file='Plate3DFluidInterfaceNodeInformationNodeElement.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(FluidInterfaceNodeInformationNE(arraySize,2))
      DO S=1,arraySize
        READ(1,*) FluidInterfaceNodeInformationNE(S,:)
      ENDDO
      arraySize=0
      CLOSE(1)
      OPEN (unit=1, file='Plate3DInterfaceInterfaceNodeInformationNodeElement.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(InterfaceInterfaceNodeInformationNE(arraySize,2))
      DO S=1,arraySize
        READ(1,*) InterfaceInterfaceNodeInformationNE(S,:)
      ENDDO
      arraySize=0
      CLOSE(1)
      OPEN (unit=1, file='Plate3DSolidInterfaceNodeInformationXi.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(SolidInterfaceNodeInformationXi(arraySize,3))
      DO S=1,arraySize
        READ(1,*) SolidInterfaceNodeInformationXi(S,:)
      ENDDO
      arraySize=0
      CLOSE(1)
      OPEN (unit=1, file='Plate3DFluidInterfaceNodeInformationXi.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(FluidInterfaceNodeInformationXi(arraySize,3))
      DO S=1,arraySize
        READ(1,*) FluidInterfaceNodeInformationXi(S,:)
      ENDDO
      arraySize=0
      CLOSE(1)
      OPEN (unit=1, file='Plate3DlagrangeNodes.txt', status='old', action='read')
      READ(1, *), arraySize
      IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
      ALLOCATE(LagrangeNodes(arraySize))
      READ(1,*) LagrangeNodes(:)
      arraySize=0
      CLOSE(1)
    ENDSELECT
    
    !Read connected solid/fluid nodes for sorted interface nodes
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DSortedInterfaceNodes.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DSortedInterfaceNodes.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(ConnectedInterfaceNodes(2,arraySize))
    DO S=1,arraySize
      READ(1,*) ConnectedInterfaceNodes(:,S)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Connected nodes: ',ConnectedInterfaceNodes
    
    !Read boundary conditions
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DdisplacementBC.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DdisplacementBC.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(NoDisplacementNodes(arraySize))
    READ(1,*) NoDisplacementNodes(:)
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Solid node numbers with no displacement: ',NoDisplacementNodes
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DfixedNodesBC.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DfixedNodesBC.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(FixedNodes(arraySize))
    READ(1,*) FixedNodes(:)
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'In both coordinates fixed mesh nodes: ',FixedNodes
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DfixedZNodesBC.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DfixedZNodesBC.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(FixedZNodes(arraySize))
    READ(1,*) FixedZNodes(:)
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'In 2nd coordinate fixed mesh nodes: ',FixedZNodes
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DmovedNodesBC.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DmovedNodesBC.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(MovedNodes(arraySize))
    READ(1,*) MovedNodes(:)
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'In both coordinates moving mesh nodes: ',MovedNodes
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DmovedYNodesBC.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DmovedYNodesBC.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(MovedYNodes(arraySize))
    READ(1,*) MovedYNodes(:)
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'In 1st coordinate moving mesh nodes: ',MovedYNodes
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DinletBC.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DinletBC.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(InletNodes(arraySize))
    READ(1,*) InletNodes(:)
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Inlet nodes: ',InletNodes
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DnoSlipBC.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DnoSlipBC.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(NoSlipNodes(arraySize))
    READ(1,*) NoSlipNodes(:)
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'No-slip nodes: ',NoSlipNodes
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DpressureBC.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DpressureBC.txt', status='old', action='read')
    ENDSELECT
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(OutletNodes(arraySize))
    READ(1,*) OutletNodes(:)
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Outlet nodes: ',OutletNodes
    
    SELECT CASE(MaterialSpecification)
    CASE(Plate2D)
      OPEN (unit=1, file='Plate2DslipBC.txt', status='old', action='read')
    CASE(Plate3D)
      OPEN (unit=1, file='Plate3DslipBC.txt', status='old', action='read')
    ENDSELECT
    READ(1,*), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(SlipNodes(arraySize))
    READ(1,*) SlipNodes(:)
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Slip BC nodes: ',SlipNodes
    
    !TODO LagrangeNodes (to remove from dof's)
    
    PRINT *, 'Finished reading geometry.'
    
    NumberOfSolidNodes=SIZE(SolidNodeNumbers)
    NumberOfSolidElements=SIZE(SolidElementNodes(:,1))
    NumberOfFluidNodes=SIZE(FluidNodeNumbers)
    NumberOfFluidElements=SIZE(FluidElementNodes(:,1))
    NumberOfInterfaceNodes=SIZE(InterfaceNodeNumbersForGeometry)
    NumberOfInterfaceElements=SIZE(InterfaceElementNodes(:,1))
    
  END SELECT
  
  !
  !=================================================================================================================================
  !
  ! C O O R D I N A T E _ S Y S T E M S
  
  !Create a new RC coordinate system for the solid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID COORDINATE SYSTEM << == '
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem1,Err)
  CALL CMISSCoordinateSystem_CreateStart(SolidCoordinateSystemUserNumber,CoordinateSystem1,Err)
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem1,NumberOfDimensions,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem1,Err)

  !Create a new RC coordinate system for the fluid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID COORDINATE SYSTEM << == '
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem2,Err)
  CALL CMISSCoordinateSystem_CreateStart(FluidCoordinateSystemUserNumber,CoordinateSystem2,Err)
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystem2,NumberOfDimensions,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem2,Err)
  
  !Create a new RC coordinate system for the interface region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE COORDINATE SYSTEM << == '
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystemInterface1,Err)
  CALL CMISSCoordinateSystem_CreateStart(InterfaceCoordinateSystemUserNumber,CoordinateSystemInterface1,Err)
  CALL CMISSCoordinateSystem_DimensionSet(CoordinateSystemInterface1,NumberOfDimensions,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystemInterface1,Err)
  
  !
  !=================================================================================================================================
  !
  ! R E G I O N S
  
  !Create the solid region and set the regions coordinate system to the RC coordinate system that we have created
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID REGION << == '
  CALL CMISSRegion_Initialise(Region1,Err)
  CALL CMISSRegion_CreateStart(SolidRegionUserNumber,WorldRegion,Region1,Err)
  CALL CMISSRegion_LabelSet(Region1,"SolidRegion",Err)
  CALL CMISSRegion_CoordinateSystemSet(Region1,CoordinateSystem1,Err)
  CALL CMISSRegion_CreateFinish(Region1,Err)

  !Create the solid region and set the regions coordinate system to the RC coordinate system that we have created
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID REGION << == '
  CALL CMISSRegion_Initialise(Region2,Err)
  CALL CMISSRegion_CreateStart(FluidRegionUserNumber,WorldRegion,Region2,Err)
  CALL CMISSRegion_LabelSet(Region2,"FluidRegion",Err)
  CALL CMISSRegion_CoordinateSystemSet(Region2,CoordinateSystem2,Err)
  CALL CMISSRegion_CreateFinish(Region2,Err)
  
  !
  !=================================================================================================================================
  !
  ! B A S I S _ S O L I D
  
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> BASIS FOR SOLID: DISPLACEMENT, HYDROSTATIC PRESSURE << == '
  !Create a basis for the dependent field variable displacement
  CALL CMISSBasis_Initialise(BasisDisplacement,Err)
  CALL CMISSBasis_CreateStart(BasisDisplacementUserNumber,BasisDisplacement,Err)
  SELECT CASE(InterpolationTypeDisplacement)
  CASE(1,2,3,4)
    CALL CMISSBasis_TypeSet(BasisDisplacement,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL CMISSBasis_TypeSet(BasisDisplacement,CMISS_BASIS_SIMPLEX_TYPE,Err)
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  SELECT CASE(InterpolationTypeDisplacement)
  CASE(CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION)
    NumberOfGaussXi=2
    PressureMeshComponent=1
  CASE(CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
    NumberOfGaussXi=3
    PressureMeshComponent=2
    InterpolationTypeHydrostaticPressure=1
  CASE(CMISS_BASIS_CUBIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION)
    NumberOfGaussXi=4
    PressureMeshComponent=2
    InterpolationTypeHydrostaticPressure=2
  CASE DEFAULT
    NumberOfGaussXi=0
    PressureMeshComponent=1
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  CALL CMISSBasis_NumberOfXiSet(BasisDisplacement,NumberOfDimensions,Err)
  IF(NumberOfDimensions==2) THEN
    CALL CMISSBasis_InterpolationXiSet(BasisDisplacement,[InterpolationTypeDisplacement,InterpolationTypeDisplacement],Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisDisplacement,[NumberOfGaussXi,NumberOfGaussXi],Err)
  ELSE
    CALL CMISSBasis_InterpolationXiSet(BasisDisplacement,[InterpolationTypeDisplacement,InterpolationTypeDisplacement, &
      & InterpolationTypeDisplacement],Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisDisplacement,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  ENDIF
  CALL CMISSBasis_CreateFinish(BasisDisplacement,Err)
  
  !Use the displacement basis as a space basis
  BasisSpaceSolid=BasisDisplacement
  
  !Create a basis for the dependent field variable hydrostatic pressure
  CALL CMISSBasis_Initialise(BasisHydrostaticPressure,Err)
  CALL CMISSBasis_CreateStart(BasisHydrostaticPressureUserNumber,BasisHydrostaticPressure,Err)
  SELECT CASE(InterpolationTypeHydrostaticPressure)
  CASE(1,2,3,4)
    CALL CMISSBasis_TypeSet(BasisHydrostaticPressure,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL CMISSBasis_TypeSet(BasisHydrostaticPressure,CMISS_BASIS_SIMPLEX_TYPE,Err)
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  CALL CMISSBasis_NumberOfXiSet(BasisHydrostaticPressure,NumberOfDimensions,Err)
  IF(NumberOfDimensions==2) THEN
    CALL CMISSBasis_InterpolationXiSet(BasisHydrostaticPressure,[InterpolationTypeHydrostaticPressure, &
      & InterpolationTypeHydrostaticPressure],Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisHydrostaticPressure,[NumberOfGaussXi,NumberOfGaussXi],Err)
  ELSE
    CALL CMISSBasis_InterpolationXiSet(BasisHydrostaticPressure,[InterpolationTypeHydrostaticPressure, &
      & InterpolationTypeHydrostaticPressure,InterpolationTypeHydrostaticPressure],Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisHydrostaticPressure,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  ENDIF
  CALL CMISSBasis_CreateFinish(BasisHydrostaticPressure,Err)
  
  !
  !=================================================================================================================================
  !
  ! B A S I S _ F L U I D
  
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> CREATING BASIS FOR FLUID: SPACE, VELOCITY, PRESSURE << == '
  !Create a basis for the fluid domain
  MeshNumberOfComponents=1
  CALL CMISSBasis_Initialise(BasisSpaceFluid,Err)
  CALL CMISSBasis_CreateStart(BasisSpaceFluidUserNumber,BasisSpaceFluid,Err)
  SELECT CASE(InterpolationTypeSpace)
  CASE(1,2,3,4)
    CALL CMISSBasis_TypeSet(BasisSpaceFluid,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL CMISSBasis_TypeSet(BasisSpaceFluid,CMISS_BASIS_SIMPLEX_TYPE,Err)
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  IF(InterpolationTypeSpace==2) THEN
    NumberOfGaussXiSpace=3
  ELSE
    CALL HANDLE_ERROR("Interpolation types other than 2 have not been tested in FiniteElasticity-ALENavierStokes example yet.")
  ENDIF
  CALL CMISSBasis_NumberOfXiSet(BasisSpaceFluid,NumberOfDimensions,Err)
  IF(NumberOfDimensions==2) THEN
    CALL CMISSBasis_InterpolationXiSet(BasisSpaceFluid,(/InterpolationTypeSpace,InterpolationTypeSpace/),Err)
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisSpaceFluid,(/NumberOfGaussXiSpace,NumberOfGaussXiSpace/),Err)
  ELSE
    CALL CMISSBasis_InterpolationXiSet(BasisSpaceFluid, &
      & (/InterpolationTypeSpace,InterpolationTypeSpace,InterpolationTypeSpace/),Err)                         
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisSpaceFluid, &
      & (/NumberOfGaussXiSpace,NumberOfGaussXiSpace,NumberOfGaussXiSpace/),Err)
  ENDIF
  CALL CMISSBasis_CreateFinish(BasisSpaceFluid,Err)
  
  !Create a basis for the dependent field variable velocity
  IF(InterpolationTypeVelocity==InterpolationTypeSpace) THEN
    BasisVelocity=BasisSpaceFluid
  ELSE
    MeshNumberOfComponents=MeshNumberOfComponents+1
    CALL CMISSBasis_Initialise(BasisVelocity,Err)
    CALL CMISSBasis_CreateStart(BasisVelocityUserNumber,BasisVelocity,Err)
    CALL CMISSBasis_TypeSet(BasisVelocity,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    NumberOfGaussXiVelocity=3
    CALL CMISSBasis_NumberOfXiSet(BasisVelocity,NumberOfDimensions,Err)
    IF(NumberOfDimensions==2) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisVelocity,(/InterpolationTypeVelocity,InterpolationTypeVelocity/),Err)
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisVelocity,(/NumberOfGaussXiVelocity,NumberOfGaussXiVelocity/),Err)
    ELSE
      CALL CMISSBasis_InterpolationXiSet(BasisVelocity, &
        & (/InterpolationTypeVelocity,InterpolationTypeVelocity,InterpolationTypeVelocity/),Err)                         
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisVelocity, &
        & (/NumberOfGaussXiVelocity,NumberOfGaussXiVelocity,NumberOfGaussXiVelocity/),Err)
    ENDIF
    CALL CMISSBasis_CreateFinish(BasisVelocity,Err)
  ENDIF
  
  !Create a basis for the dependent field variable pressure
  IF(InterpolationTypePressure==InterpolationTypeSpace) THEN
    BasisPressure=BasisSpaceFluid
  ELSE IF(InterpolationTypePressure==InterpolationTypeVelocity) THEN
    BasisPressure=BasisVelocity
  ELSE
    MeshNumberOfComponents=MeshNumberOfComponents+1
    NumberOfGaussXiPressure=3
    CALL CMISSBasis_Initialise(BasisPressure,Err)
    CALL CMISSBasis_CreateStart(BasisPressureUserNumber,BasisPressure,Err)
    CALL CMISSBasis_TypeSet(BasisPressure,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CALL CMISSBasis_NumberOfXiSet(BasisPressure,NumberOfDimensions,Err)
    IF(NumberOfDimensions==2) THEN
      CALL CMISSBasis_InterpolationXiSet(BasisPressure,(/InterpolationTypePressure,InterpolationTypePressure/),Err)
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisPressure,(/NumberOfGaussXiPressure,NumberOfGaussXiPressure/),Err)
    ELSE
      CALL CMISSBasis_InterpolationXiSet(BasisPressure, &
        & (/InterpolationTypePressure,InterpolationTypePressure,InterpolationTypePressure/),Err)                         
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(BasisPressure, &
        & (/NumberOfGaussXiPressure,NumberOfGaussXiPressure,NumberOfGaussXiPressure/),Err)
    ENDIF
    CALL CMISSBasis_CreateFinish(BasisPressure,Err)
  ENDIF

  !
  !================================================================================================================================
  !
  ! M E S H E S
  
  !Create the solid mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID MESH << == '
  CALL CMISSNodes_Initialise(SolidNodes,Err)
  CALL CMISSMesh_Initialise(Mesh1,Err)
  CALL CMISSNodes_CreateStart(Region1,NumberOfSolidNodes,SolidNodes,Err)
  CALL CMISSNodes_CreateFinish(SolidNodes,Err)
  
  CALL CMISSMesh_CreateStart(SolidMeshUserNumber,Region1,NumberOfDimensions,Mesh1,Err)
  CALL CMISSMesh_NumberOfElementsSet(Mesh1,NumberOfSolidElements,Err)
  CALL CMISSMesh_NumberOfComponentsSet(Mesh1,2,Err)
  
  CALL CMISSMeshElements_Initialise(Mesh1ElementsSpace,Err)
  CALL CMISSMeshElements_Initialise(Mesh1ElementsDisplacement,Err)
  CALL CMISSMeshElements_Initialise(Mesh1ElementsHydrostaticPressure,Err)
  
  Mesh1ComponentNumberSpace=1
  Mesh1ComponentNumberDisplacement=1
  Mesh1ComponentNumberHydrostaticPressure=1
  
  CALL CMISSMeshElements_CreateStart(Mesh1,Mesh1ComponentNumberSpace,BasisSpaceSolid,Mesh1ElementsSpace,Err)
  SELECT CASE(MaterialSpecification)
  CASE(Blood)
    ALLOCATE(IndexArray(9))
    IF(.NOT.ALLOCATED(IndexArray)) CALL HANDLE_ERROR('IndexArray not allocated.')
    IndexArray(:)=(/1,2,3,4,5,6,7,8,9/)
    DO ElementIndex=1,NumberOfSolidElements
      CALL CMISSMeshElements_NodesSet(Mesh1ElementsSpace,SolidElements(ElementIndex), &
        & SolidElementSDNodeMappings(IndexArray),Err)
      IndexArray=IndexArray+9
    ENDDO
    IF(ALLOCATED(IndexArray)) DEALLOCATE(IndexArray)
  CASE(Plate2D,Plate3D)
    DO ElementIndex=1,NumberOfSolidElements
      CALL CMISSMeshElements_NodesSet(Mesh1ElementsSpace,ElementIndex, &
        & SolidElementNodes(ElementIndex,:),Err)
    ENDDO
  END SELECT
  CALL CMISSMeshElements_CreateFinish(Mesh1ElementsSpace,Err)
  
  Mesh1ElementsDisplacement=Mesh1ElementsSpace
  
  Mesh1ComponentNumberHydrostaticPressure=Mesh1ComponentNumberDisplacement+1
  CALL CMISSMeshElements_CreateStart(Mesh1,Mesh1ComponentNumberHydrostaticPressure, &
    & BasisHydrostaticPressure,Mesh1ElementsHydrostaticPressure,Err)
  SELECT CASE(MaterialSpecification)
  CASE(Blood)
    ALLOCATE(IndexArray(4))
    IF(.NOT.ALLOCATED(IndexArray)) CALL HANDLE_ERROR('IndexArray not allocated.')
    IndexArray(:)=(/1,2,3,4/)
    DO ElementIndex=1,NumberOfSolidElements
      CALL CMISSMeshElements_NodesSet(Mesh1ElementsHydrostaticPressure,SolidElements(ElementIndex), &
        & SolidElementHPNodeMappings(IndexArray),Err)
      IndexArray=IndexArray+4
    ENDDO
    IF(ALLOCATED(IndexArray)) DEALLOCATE(IndexArray)
  CASE(Plate2D)
    DO ElementIndex=1,NumberOfSolidElements
      CALL CMISSMeshElements_NodesSet(Mesh1ElementsHydrostaticPressure,ElementIndex, &
        & SolidElementNodes(ElementIndex,(/1,3,7,9/)),Err)
    ENDDO
  CASE(Plate3D)
    DO ElementIndex=1,NumberOfSolidElements
      CALL CMISSMeshElements_NodesSet(Mesh1ElementsHydrostaticPressure,ElementIndex, &
        & SolidElementNodes(ElementIndex,(/1,3,7,9,19,21,25,27/)),Err)
    ENDDO
  END SELECT
  CALL CMISSMeshElements_CreateFinish(Mesh1ElementsHydrostaticPressure,Err)
  CALL CMISSMesh_CreateFinish(Mesh1,Err)
  
  !Create fluid mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID MESH << == '
  CALL CMISSNodes_Initialise(FluidNodes,Err)
  CALL CMISSMesh_Initialise(Mesh2,Err)
  CALL CMISSNodes_CreateStart(Region2,NumberOfFluidNodes,FluidNodes,Err)
  CALL CMISSNodes_CreateFinish(FluidNodes,Err)
  
  CALL CMISSMesh_CreateStart(FluidMeshUserNumber,Region2,NumberOfDimensions,Mesh2,Err)
  CALL CMISSMesh_NumberOfElementsSet(Mesh2,NumberOfFluidElements,Err)
  CALL CMISSMesh_NumberOfComponentsSet(Mesh2,2,Err)
  
  CALL CMISSMeshElements_Initialise(Mesh2ElementsSpace,Err)
  CALL CMISSMeshElements_Initialise(Mesh2ElementsVelocity,Err)
  CALL CMISSMeshElements_Initialise(Mesh2ElementsPressure,Err)
  
  Mesh2ComponentNumberSpace=1
  Mesh2ComponentNumberVelocity=1
  Mesh2ComponentNumberPressure=1
  
  CALL CMISSMeshElements_CreateStart(Mesh2,Mesh2ComponentNumberSpace,BasisSpaceFluid,Mesh2ElementsSpace,Err)
  SELECT CASE(MaterialSpecification)
  CASE(Blood)
    ALLOCATE(IndexArray(9))
    IF(.NOT.ALLOCATED(IndexArray)) CALL HANDLE_ERROR('IndexArray not allocated.')
    IndexArray(:)=(/1,2,3,4,5,6,7,8,9/)
    DO ElementIndex=1,NumberOfFluidElements
      CALL CMISSMeshElements_NodesSet(Mesh2ElementsSpace,FluidElements(ElementIndex), &
        & FluidElementSVNodeMappings(IndexArray),Err)
      IndexArray=IndexArray+9
    ENDDO
    IF(ALLOCATED(IndexArray)) DEALLOCATE(IndexArray)
  CASE(Plate2D,Plate3D)
    DO ElementIndex=1,NumberOfFluidElements
      CALL CMISSMeshElements_NodesSet(Mesh2ElementsSpace,ElementIndex, &
        & FluidElementNodes(ElementIndex,:),Err)
    ENDDO
  END SELECT
  CALL CMISSMeshElements_CreateFinish(Mesh2ElementsSpace,Err)
  
  Mesh2ElementsVelocity=Mesh2ElementsSpace
  
  Mesh2ComponentNumberPressure=Mesh2ComponentNumberVelocity+1
  CALL CMISSMeshElements_CreateStart(Mesh2,Mesh2ComponentNumberPressure, &
    & BasisPressure,Mesh2ElementsPressure,Err)
  SELECT CASE(MaterialSpecification)
  CASE(Blood)
    ALLOCATE(IndexArray(4))
    IF(.NOT.ALLOCATED(IndexArray)) CALL HANDLE_ERROR('IndexArray not allocated.')
    IndexArray(:)=(/1,2,3,4/)
    DO ElementIndex=1,NumberOfFluidElements
      CALL CMISSMeshElements_NodesSet(Mesh2ElementsPressure,FluidElements(ElementIndex), &
        & FluidElementPNodeMappings(IndexArray),Err)
      IndexArray=IndexArray+4
    ENDDO
    IF(ALLOCATED(IndexArray)) DEALLOCATE(IndexArray)
  CASE(Plate2D)
    DO ElementIndex=1,NumberOfFluidElements
      CALL CMISSMeshElements_NodesSet(Mesh2ElementsPressure,ElementIndex, &
        & FluidElementNodes(ElementIndex,(/1,3,7,9/)),Err)
    ENDDO
  CASE(Plate3D)
    DO ElementIndex=1,NumberOfFluidElements
      CALL CMISSMeshElements_NodesSet(Mesh2ElementsPressure,ElementIndex, &
        & FluidElementNodes(ElementIndex,(/1,3,7,9,19,21,25,27/)),Err)
    ENDDO
  END SELECT
  CALL CMISSMeshElements_CreateFinish(Mesh2ElementsPressure,Err)
  CALL CMISSMesh_CreateFinish(Mesh2,Err)
  
  !
  !================================================================================================================================
  !
  ! I N T E R F A C E

  !Create an interface between the two meshes
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE << == '
  CALL CMISSInterface_Initialise(Interface1,Err)
  CALL CMISSInterface_CreateStart(Interface1UserNumber,WorldRegion,Interface1,Err)
  CALL CMISSInterface_LabelSet(Interface1,"Interface1",Err)
  !Add in the two meshes
  CALL CMISSInterface_MeshAdd(Interface1,Mesh1,Mesh1Index,Err)
  CALL CMISSInterface_MeshAdd(Interface1,Mesh2,Mesh2Index,Err)
  CALL CMISSInterface_CoordinateSystemSet(Interface1,CoordinateSystemInterface1,Err)
  CALL CMISSInterface_CreateFinish(Interface1,Err)

  !Create a (bi)-quadratic-Lagrange basis (3D: faces // 2D: lines)
  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> INTERFACE BASIS << == '
  CALL CMISSBasis_Initialise(InterfaceBasis1,Err)
  CALL CMISSBasis_CreateStart(InterfaceBasisUserNumber,InterfaceBasis1,Err)
  CALL CMISSBasis_NumberOfXiSet(InterfaceBasis1,NumberOfDimensions-1,Err)
  IF(NumberOfDimensions==2) THEN
    CALL CMISSBasis_InterpolationXiSet(InterfaceBasis1,[CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  ELSE
    CALL CMISSBasis_InterpolationXiSet(InterfaceBasis1, &
      & [CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  ENDIF
  CALL CMISSBasis_CreateFinish(InterfaceBasis1,Err)

  !Create a (bi)-quadratic-Lagrange basis for the interface mapping (3D: faces // 2D: lines)
  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> INTERFACE MAPPING BASIS << == '
  CALL CMISSBasis_Initialise(InterfaceMappingBasis1,Err)
  CALL CMISSBasis_CreateStart(InterfaceMappingBasis1UserNumber,InterfaceMappingBasis1,Err)
  CALL CMISSBasis_NumberOfXiSet(InterfaceMappingBasis1,NumberOfDimensions-1,Err)
  IF(NumberOfDimensions==2) THEN
    CALL CMISSBasis_InterpolationXiSet(InterfaceMappingBasis1,[CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  ELSE
    CALL CMISSBasis_InterpolationXiSet(InterfaceMappingBasis1, &
      & [CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  ENDIF
  CALL CMISSBasis_CreateFinish(InterfaceMappingBasis1,Err)
  
  !
  !================================================================================================================================
  !
  ! I N T E R F A C E _ M E S H
  
  !Create an interface mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE MESH << == '
  CALL CMISSNodes_Initialise(InterfaceNodes,Err)
  CALL CMISSMesh_Initialise(InterfaceMesh1,Err)
  CALL CMISSNodes_CreateStart(Interface1,NumberOfInterfaceNodes,InterfaceNodes,Err)
  CALL CMISSNodes_CreateFinish(InterfaceNodes,Err)
  
  CALL CMISSMesh_CreateStart(InterfaceMeshUserNumber,Interface1,NumberOfDimensions-1,InterfaceMesh1,Err)
  CALL CMISSMesh_NumberOfElementsSet(InterfaceMesh1,NumberOfInterfaceElements,Err)
  CALL CMISSMesh_NumberOfComponentsSet(InterfaceMesh1,1,Err)
  
  CALL CMISSMeshElements_Initialise(InterfaceMeshElements,Err)
  InterfaceMeshComponentNumber=1
  CALL CMISSMeshElements_CreateStart(InterfaceMesh1,InterfaceMeshComponentNumber,InterfaceBasis1,InterfaceMeshElements,Err)
  SELECT CASE(MaterialSpecification)
  CASE(Blood)
    ALLOCATE(IndexArray(3))
    IF(.NOT.ALLOCATED(IndexArray)) CALL HANDLE_ERROR('IndexArray not allocated.')
    IndexArray(:)=(/1,2,3/)
    DO ElementIndex=1,NumberOfInterfaceElements
      CALL CMISSMeshElements_NodesSet(InterfaceMeshElements,InterfaceElements(ElementIndex), &
        & InterfaceInterfaceElementNodes(IndexArray),Err)
      IndexArray=IndexArray+3
    ENDDO
    IF(ALLOCATED(IndexArray)) DEALLOCATE(IndexArray)
  CASE(Plate2D,Plate3D)
    DO ElementIndex=1,NumberOfInterfaceElements
      CALL CMISSMeshElements_NodesSet(InterfaceMeshElements,ElementIndex, &
        & InterfaceElementNodes(ElementIndex,:),Err)
    ENDDO
  END SELECT
  CALL CMISSMeshElements_CreateFinish(InterfaceMeshElements,Err)
  CALL CMISSMesh_CreateFinish(InterfaceMesh1,Err)

  !
  !================================================================================================================================
  !
  ! M E S H _ C O N N E C T I V I T Y

  !Couple the interface meshes
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE MESH CONNECTIVITY << == '
  IF (NumberOfComputationalNodes/=1) CALL HANDLE_ERROR("MeshConnectivity does not work in parallel yet!")
  CALL CMISSInterfaceMeshConnectivity_Initialise(InterfaceMeshConnectivity1,Err)
  CALL CMISSInterfaceMeshConnectivity_CreateStart(Interface1,InterfaceMesh1,InterfaceMeshConnectivity1,Err)
  CALL CMISSInterfaceMeshConnectivity_BasisSet(InterfaceMeshConnectivity1,InterfaceMappingBasis1,Err)
  !TODO Modify to match interface (mapping) basis interpolation type
  SELECT CASE(InterpolationTypeInterface)
  CASE(1,4)
    NUMBER_OF_NODE_XI=2
  CASE(2)
    NUMBER_OF_NODE_XI=3
  CASE(3)
    NUMBER_OF_NODE_XI=4
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  IF(NumberOfDimensions==2) THEN
    SELECT CASE(MaterialSpecification)
    CASE(Blood)
      DO ElementIndex=1,NumberOfInterfaceElements*NUMBER_OF_NODE_XI,NUMBER_OF_NODE_XI
        !SUBROUTINE CMISSInterfaceMeshConnectivity_ElementXiSetObj(interfaceMeshConnectivity,interfaceElementNumber, &
        !  & coupledMeshIndexNumber,coupledMeshElementNumber,interfaceMeshLocalNodeNumber,interfaceMeshComponentNodeNumber,xi,err)
        
        !Map the interface element to the elements in mesh 1
        CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1, &
          & InterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),Mesh1Index, &
          & SolidInterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),Err)
        
        XI2=[SolidInterfaceXiPosition(ElementIndex), &
          & SolidInterfaceXiPosition(ElementIndex+NumberOfInterfaceElements*NUMBER_OF_NODE_XI)]
        PRINT *, XI2
        CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1, &
          & InterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),Mesh1Index, &
          & SolidInterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),1,1,XI2,Err)
        
        XI2=[SolidInterfaceXiPosition(ElementIndex+1), &
          & SolidInterfaceXiPosition(ElementIndex+1+NumberOfInterfaceElements*NUMBER_OF_NODE_XI)]
        PRINT *, XI2
        CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1, &
          & InterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),Mesh1Index, &
          & SolidInterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),2,1,XI2,Err)
        
        XI2=[SolidInterfaceXiPosition(ElementIndex+2), &
          & SolidInterfaceXiPosition(ElementIndex+2+NumberOfInterfaceElements*NUMBER_OF_NODE_XI)]
        PRINT *, XI2
        CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1, &
          & InterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),Mesh1Index, &
          & SolidInterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),3,1,XI2,Err)
        
        !Map the interface element to the elements in mesh 2
        CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1, &
          & InterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),Mesh2Index, &
          & FluidInterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),Err)
        
        XI2=[FluidInterfaceXiPosition(ElementIndex), &
          & FluidInterfaceXiPosition(ElementIndex+NumberOfInterfaceElements*NUMBER_OF_NODE_XI)]
        CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1, &
          & InterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),Mesh2Index, &
          & FluidInterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),1,1,XI2,Err)
        
        XI2=[FluidInterfaceXiPosition(ElementIndex+1), &
          & FluidInterfaceXiPosition(ElementIndex+1+NumberOfInterfaceElements*NUMBER_OF_NODE_XI)]
        CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1, &
          & InterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),Mesh2Index, &
          & FluidInterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),2,1,XI2,Err)
        
        XI2=[FluidInterfaceXiPosition(ElementIndex+2), &
          & FluidInterfaceXiPosition(ElementIndex+2+NumberOfInterfaceElements*NUMBER_OF_NODE_XI)]
        CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1, &
          & InterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),Mesh2Index, &
          & FluidInterfaceElements((ElementIndex+NUMBER_OF_NODE_XI-1)/NUMBER_OF_NODE_XI),3,1,XI2,Err)
      ENDDO !ElementIndex
      CALL CMISSInterfaceMeshConnectivity_NodeNumberSet(InterfaceMeshConnectivity1,InterfaceNodeNumbers, &
        & Mesh1Index,SolidInterfaceNodeNumberMappings,Mesh2Index,FluidInterfaceNodeNumberMappings,Err)
    CASE(Plate2D)
      DO ElementIndex=1,NumberOfInterfaceElements
        !Map the interface element to the elements in mesh 1
        CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1,ElementIndex,Mesh1Index, &
          & SolidInterfaceElements(ElementIndex),Err)
        DO LocalNodeIndex=1,3**(NumberOfDimensions-1)
          XI2 = [SolidXi2(ElementIndex,LocalNodeIndex),SolidXi3(ElementIndex,LocalNodeIndex)]
          CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1,ElementIndex,Mesh1Index, &
            & SolidInterfaceElements(ElementIndex),LocalNodeIndex,1,XI2,Err)
        ENDDO
        !Map the interface element to the elements in mesh 2
        CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1,ElementIndex,Mesh2Index, &
          & FluidInterfaceElements(ElementIndex),Err)
        DO LocalNodeIndex=1,3**(NumberOfDimensions-1)
          XI2 = [FluidXi2(ElementIndex,LocalNodeIndex),FluidXi3(ElementIndex,LocalNodeIndex)]
          CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1,ElementIndex,Mesh2Index, &
            & FluidInterfaceElements(ElementIndex),LocalNodeIndex,1,XI2,Err)
        ENDDO
      ENDDO !ElementIndex
      CALL CMISSInterfaceMeshConnectivity_NodeNumberSet(InterfaceMeshConnectivity1,InterfaceNodeNumbersForGeometry, &
        & Mesh1Index,ConnectedInterfaceNodes(1,:),Mesh2Index,ConnectedInterfaceNodes(2,:),Err)
    CASE DEFAULT
      CALL HANDLE_ERROR("Invalid material specification for mesh connectivity.")
    END SELECT
  ELSE
  !  CALL HANDLE_ERROR("Only 2D Mesh Connectivity implemented in this example file yet.")
    SELECT CASE(MaterialSpecification)
    CASE(Plate3D)
      LocalNodeIndex=0
      DO ElementIndex=1,NumberOfInterfaceElements*9
        LocalNodeIndex=LocalNodeIndex+1
        !Map the interface element to the elements in the fluid mesh
        CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1, &
          & InterfaceInterfaceNodeInformationNE(ElementIndex,2),Mesh1Index,SolidInterfaceNodeInformationNE(ElementIndex,2),Err)
        XI3 = [SolidInterfaceNodeInformationXi(ElementIndex,1),SolidInterfaceNodeInformationXi(ElementIndex,2), &
          & SolidInterfaceNodeInformationXi(ElementIndex,3)]
        CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1, &
          & InterfaceInterfaceNodeInformationNE(ElementIndex,2),Mesh1Index, &
          & SolidInterfaceNodeInformationNE(ElementIndex,2),LocalNodeIndex,1,XI3,Err)
        !Map the interface element to the elements in the solid mesh
        CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1, &
          & InterfaceInterfaceNodeInformationNE(ElementIndex,2),Mesh2Index, &
          & FluidInterfaceNodeInformationNE(ElementIndex,2),Err)
        XI3 = [FluidInterfaceNodeInformationXi(ElementIndex,1),FluidInterfaceNodeInformationXi(ElementIndex,2), &
          & FluidInterfaceNodeInformationXi(ElementIndex,3)]
        CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1, &
          & InterfaceInterfaceNodeInformationNE(ElementIndex,2),Mesh2Index, &
          & FluidInterfaceNodeInformationNE(ElementIndex,2),LocalNodeIndex,1,XI3,Err)
        IF(LocalNodeIndex==9) LocalNodeIndex=0
      ENDDO !ElementIndex
      CALL CMISSInterfaceMeshConnectivity_NodeNumberSet(InterfaceMeshConnectivity1,InterfaceNodeNumbersForGeometry, &
        & Mesh1Index,ConnectedInterfaceNodes(1,:),Mesh2Index,ConnectedInterfaceNodes(2,:),Err)
    CASE DEFAULT
      CALL HANDLE_ERROR("Invalid material specification for mesh connectivity.")
    END SELECT
  ENDIF
  CALL CMISSInterfaceMeshConnectivity_CreateFinish(InterfaceMeshConnectivity1,Err)

  !
  !================================================================================================================================
  !
  ! D E C O M P O S I T I O N

  !Create a decomposition for the solid mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID MESH DECOMPOSITION << == '
  CALL CMISSDecomposition_Initialise(Decomposition1,Err)
  CALL CMISSDecomposition_CreateStart(SolidDecompositionUserNumber,Mesh1,Decomposition1,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition1,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition1,NumberOfComputationalNodes,Err)
  CALL CMISSDecomposition_CalculateFacesSet(Decomposition1,.TRUE.,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition1,Err)

  !Create a decomposition for the fluid mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID MESH DECOMPOSITION << == '
  CALL CMISSDecomposition_Initialise(Decomposition2,Err)
  CALL CMISSDecomposition_CreateStart(FluidDecompositionUserNumber,Mesh2,Decomposition2,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(Decomposition2,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition2,NumberOfComputationalNodes,Err)
  CALL CMISSDecomposition_CalculateFacesSet(Decomposition2,.TRUE.,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(Decomposition2,Err)
  
  !Create a decomposition for the interface mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE DECOMPOSITION << == '
  CALL CMISSDecomposition_Initialise(InterfaceDecomposition1,Err)
  CALL CMISSDecomposition_CreateStart(InterfaceDecompositionUserNumber,InterfaceMesh1,InterfaceDecomposition1,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(InterfaceDecomposition1,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(InterfaceDecomposition1,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(InterfaceDecomposition1,Err)

  !
  !================================================================================================================================
  !
  ! G E O M E T R I C _ F I E L D S

  !Start to create a default (geometric) field on the solid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID MESH GEOMETRIC FIELD << == '
  CALL CMISSField_Initialise(GeometricField1,Err)
  CALL CMISSField_CreateStart(SolidGeometricFieldUserNumber,Region1,GeometricField1,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField1,Decomposition1,Err)
  !Set the scaling to use
  CALL CMISSField_ScalingTypeSet(GeometricField1,CMISS_FIELD_NO_SCALING,Err)
  CALL CMISSField_VariableLabelSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,"SolidGF",Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NumberOfDimensions==3) CALL CMISSField_ComponentMeshComponentSet(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the first field
  CALL CMISSField_CreateFinish(GeometricField1,Err)

  !Start to create a default (geometric) field on the fluid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID MESH GEOMETRIC FIELD << == '
  CALL CMISSField_Initialise(GeometricField2,Err)
  CALL CMISSField_CreateStart(FluidGeometricFieldUserNumber,Region2,GeometricField2,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField2,Decomposition2,Err)
  !Set the scaling to use
  CALL CMISSField_ScalingTypeSet(GeometricField2,CMISS_FIELD_NO_SCALING,Err)
  CALL CMISSField_VariableLabelSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,"FluidGF",Err)
  !Set the domain to be used by the field components.
  CALL CMISSField_ComponentMeshComponentSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NumberOfDimensions==3) CALL CMISSField_ComponentMeshComponentSet(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the second field
  CALL CMISSField_CreateFinish(GeometricField2,Err)

  !Start to create a default (geometric) field on the Interface
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE GEOMETRIC FIELD << == '
  CALL CMISSField_Initialise(InterfaceGeometricField1,Err)
  CALL CMISSField_CreateStart(InterfaceGeometricFieldUserNumber,Interface1,InterfaceGeometricField1,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(InterfaceGeometricField1,InterfaceDecomposition1,Err)
  CALL CMISSField_VariableLabelSet(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,"InterfaceGF",Err)
  !Set the domain to be used by the field components
  CALL CMISSField_ComponentMeshComponentSet(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NumberOfDimensions==3) CALL CMISSField_ComponentMeshComponentSet(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the first field
  CALL CMISSField_CreateFinish(InterfaceGeometricField1,Err)

  SELECT CASE(MaterialSpecification)
  CASE(Blood)
    !Update the geometric field parameters (solid)
    DO NodeIndex=1,NumberOfSolidNodes
      CALL CMISSDecomposition_NodeDomainGet(Decomposition1,NodeIndexSMappings(NodeIndex),1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,NodeIndexSMappings(NodeIndex),1,SolidGeometry(NodeIndex),Err)
        CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,NodeIndexSMappings(NodeIndex),2,SolidGeometry(NodeIndex+NumberOfSolidNodes),Err)
        IF(NumberGlobalElementsZ/=0) THEN
          CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,NodeIndexSMappings(NodeIndex),3,SolidGeometry(NodeIndex+NumberOfSolidNodes*2),Err)
        ENDIF
      ENDIF
    ENDDO
    !Update the geometric field parameters (fluid)
    DO NodeIndex=1,NumberOfFluidNodes
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NodeIndexFMappings(NodeIndex),1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,NodeIndexFMappings(NodeIndex),1,FluidGeometry(NodeIndex),Err)
        CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,NodeIndexFMappings(NodeIndex),2,FluidGeometry(NodeIndex+NumberOfFluidNodes),Err)
        IF(NumberGlobalElementsZ/=0) THEN
          CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,NodeIndexFMappings(NodeIndex),3,FluidGeometry(NodeIndex+NumberOfFluidNodes*2),Err)
        ENDIF
      ENDIF
    ENDDO
    !Update the geometric field parameters (interface)
    DO NodeIndex=1,NumberOfInterfaceNodes
      CALL CMISSDecomposition_NodeDomainGet(InterfaceDecomposition1,NodeIndex,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSField_ParameterSetUpdateNode(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,NodeIndex,1,InterfaceGeometry(NodeIndex),Err)
        CALL CMISSField_ParameterSetUpdateNode(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,NodeIndex,2,InterfaceGeometry(NodeIndex+NumberOfInterfaceNodes),Err)
        IF(NumberGlobalElementsZ/=0) THEN
          CALL CMISSField_ParameterSetUpdateNode(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,NodeIndex,3,InterfaceGeometry(NodeIndex+2*NumberOfInterfaceNodes),Err)
        ENDIF
      ENDIF
    ENDDO
  CASE(Plate2D,Plate3D)
    !Update the geometric field parameters (solid)
    DO NodeIndex=1,NumberOfSolidNodes
      CALL CMISSDecomposition_NodeDomainGet(Decomposition1,SolidNodeNumbers(NodeIndex),1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,SolidNodeNumbers(NodeIndex),1,SolidGeometryY(NodeIndex),Err)
        CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,SolidNodeNumbers(NodeIndex),2,SolidGeometryZ(NodeIndex),Err)
        IF(NumberOfDimensions==3) THEN
          CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,SolidNodeNumbers(NodeIndex),3,SolidGeometryX(NodeIndex),Err)
        ENDIF
      ENDIF
    ENDDO
    !Update the geometric field parameters (fluid)
    DO NodeIndex=1,NumberOfFluidNodes
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,FluidNodeNumbers(NodeIndex),1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,FluidNodeNumbers(NodeIndex),1,FluidGeometryY(NodeIndex),Err)
        CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,FluidNodeNumbers(NodeIndex),2,FluidGeometryZ(NodeIndex),Err)
        IF(NumberOfDimensions==3) THEN
          CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,FluidNodeNumbers(NodeIndex),3,FluidGeometryX(NodeIndex),Err)
        ENDIF
      ENDIF
    ENDDO
    !Update the geometric field parameters (interface)
    DO NodeIndex=1,NumberOfInterfaceNodes
      CALL CMISSDecomposition_NodeDomainGet(InterfaceDecomposition1,InterfaceNodeNumbersForGeometry(NodeIndex),1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSField_ParameterSetUpdateNode(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,InterfaceNodeNumbersForGeometry(NodeIndex),1,InterfaceGeometryY(NodeIndex),Err)
        CALL CMISSField_ParameterSetUpdateNode(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,InterfaceNodeNumbersForGeometry(NodeIndex),2,InterfaceGeometryZ(NodeIndex),Err)
        IF(NumberOfDimensions==3) THEN
          CALL CMISSField_ParameterSetUpdateNode(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,InterfaceNodeNumbersForGeometry(NodeIndex),3,InterfaceGeometryX(NodeIndex),Err)
        ENDIF
      ENDIF
    ENDDO
  END SELECT
  
  CALL CMISSField_ParameterSetUpdateStart(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateStart(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateStart(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  
  !
  !================================================================================================================================
  !
  ! F I B R E _ F I E L D _ S O L I D
  
  !Create a fibre field and attach it to the geometric field 1
!  CALL CMISSField_Initialise(FibreField,Err)
!  CALL CMISSField_CreateStart(FibreFieldUserNumber,Region1,FibreField,Err)
!  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
!  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition1,Err)
!  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField1,Err)
!  CALL CMISSField_VariableLabelSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,"SolidFF",Err)
!  CALL CMISSField_ScalingTypeSet(FibreField,CMISS_FIELD_NO_SCALING,Err)
!  CALL CMISSField_CreateFinish(FibreField,Err)

  !
  !================================================================================================================================
  !
  ! E Q U A T I O N _ S E T S

  !Create the equations set for the solid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID EQUATION SET << == '
  !Create the equations_set for FiniteElasticity MooneyRivlin
  CALL CMISSField_Initialise(EquationsSetField1,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSet1,Err)
  CALL CMISSEquationsSet_CreateStart(SolidEquationsSetUserNumber,Region1,GeometricField1,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EquationsSetField1UserNumber, &
    & EquationsSetField1,EquationsSet1,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSet1,Err)

  !Create the equations set for the fluid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID EQUATION SET << == '
  !Create the equations set for ALE Navier-Stokes
  CALL CMISSField_Initialise(EquationsSetField2,Err)
  CALL CMISSEquationsSet_Initialise(EquationsSet2,Err)
  CALL CMISSEquationsSet_CreateStart(FluidEquationsSetUserNumber,Region2,GeometricField2,CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS,&
    & CMISS_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE,CMISS_EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE,&
    & EquationsSetField2UserNumber,EquationsSetField2,EquationsSet2,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(EquationsSet2,Err)

  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> MOVING MESH EQUATION SET << == '
  !Create the equations set for the moving mesh
  CALL CMISSEquationsSet_Initialise(MovingMeshEquationsSet,Err)
  CALL CMISSField_Initialise(EquationsSetFieldMovingMesh,Err)
  CALL CMISSEquationsSet_CreateStart(MovingMeshEquationsSetUserNumber,Region2,GeometricField2, &
    & CMISS_EQUATIONS_SET_CLASSICAL_FIELD_CLASS,CMISS_EQUATIONS_SET_LAPLACE_EQUATION_TYPE, &
    & CMISS_EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE,EquationsSetFieldMovingMeshUserNumber,EquationsSetFieldMovingMesh, &
    & MovingMeshEquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(MovingMeshEquationsSet,Err)

  !
  !================================================================================================================================
  !
  ! D E P E N D E N T _ F I E L D S

  !Create the equations set dependent field variables for the first equations set
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID DEPENDENT FIELD << == '
  !Create the dependent field
  CALL CMISSField_Initialise(DependentField1,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet1,DependentField1UserNumber,DependentField1,Err)
  CALL CMISSField_VariableLabelSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,"SolidDF",Err)
  DO component_idx=1,NumberOfDimensions
    CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,component_idx,1,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,component_idx,1,Err)
  ENDDO
  CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,NumberOfDimensions+1, &
    & PressureMeshComponent,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,NumberOfDimensions+1, &
    & PressureMeshComponent,Err)
  IF(PressureMeshComponent==1) THEN
    CALL CMISSField_ComponentInterpolationSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,NumberOfDimensions+1, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,NumberOfDimensions+1, &
      & CMISS_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ELSE
    CALL CMISSField_ComponentInterpolationSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,NumberOfDimensions+1, &
      & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,NumberOfDimensions+1, &
      & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
  ENDIF
  CALL CMISSField_ScalingTypeSet(DependentField1,CMISS_FIELD_NO_SCALING,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet1,Err)
  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,1,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE, &
    & CMISS_FIELD_VALUES_SET_TYPE,2,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  IF (NumberOfDimensions==3) THEN
    CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE, &
      & CMISS_FIELD_VALUES_SET_TYPE,3,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  ENDIF
  CALL CMISSField_ComponentValuesInitialise(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & NumberOfDimensions+1,-14.0_CMISSDP,Err)
  
  CALL CMISSField_ParameterSetUpdateStart(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> SOLID MATERIAL FIELD << == '
  !Create the material field
  CALL CMISSField_Initialise(MaterialField1,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet1,MaterialField1UserNumber,MaterialField1,Err)
  CALL CMISSField_VariableLabelSet(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,"Material1",Err)
  CALL CMISSField_VariableLabelSet(MaterialField1,CMISS_FIELD_V_VARIABLE_TYPE,"SolidDensity",Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet1,Err)

  !Set Mooney-Rivlin constants c10 and c01 (default?2.0 and 6.0) respectively
  CALL CMISSField_ComponentValuesInitialise(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & MooneyRivlin1,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2, &
    & MooneyRivlin2,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField1,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & SolidDensity,Err)

  !IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOURCE FIELD - GRAVITY << == '
  !!!Create the source field with the gravity vector
  !CALL CMISSField_Initialise(SourceField1,Err)
  !PRINT *, 'HERE'
  !CALL CMISSEquationsSet_SourceCreateStart(EquationsSet1,SourceField1UserNumber,SourceField1,Err)
  !PRINT *, 'OR HERE'
  !CALL CMISSField_ScalingTypeSet(SourceField1,CMISS_FIELD_NO_SCALING,Err)
  !PRINT *, 'OR HERE'
  !CALL CMISSEquationsSet_SourceCreateFinish(EquationsSet1,Err)
  !PRINT *, 'OR HERE'
  !DO component_idx=1,NumberOfDimensions
  !  PRINT *, component_idx
  !  CALL CMISSField_ComponentValuesInitialise(SourceField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
  !      & component_idx,Gravity(component_idx),Err)
  !ENDDO
  
  !=========================

  !Create the equations set dependent field variables for dynamic Navier-Stokes
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID DEPENDENT FIELD << == '
  CALL CMISSField_Initialise(DependentField2,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet2,DependentField2UserNumber,DependentField2,Err)
  CALL CMISSField_VariableLabelSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,"FluidDF",Err)
!  Mesh2ComponentNumberSpace=1
!  Mesh2ComponentNumberVelocity=1
!  Mesh2ComponentNumberPressure=2 ! TODO CHECK if this is correct
  !Set the mesh component to be used by the field components.
  DO ComponentNumber=1,NumberOfDimensions
    CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,ComponentNumber,1,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,ComponentNumber,1,Err)
  ENDDO
  CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,NumberOfDimensions+1, & 
    & Mesh2ComponentNumberPressure,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,NumberOfDimensions+1, & 
    & Mesh2ComponentNumberPressure,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet2,Err)
  !Initialise dependent field
  DO ComponentNumber=1,NumberOfDimensions
    CALL CMISSField_ComponentValuesInitialise(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & ComponentNumber,InitialFieldNavierStokes(ComponentNumber),Err)
  ENDDO
  !Initialise pressure component
  CALL CMISSField_ComponentValuesInitialise(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & NumberOfDimensions+1,0.1_CMISSDP,Err)
  
  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> MOVING MESH DEPENDENT FIELD << == '
  !Create the equations set dependent field variables for moving mesh
  CALL CMISSField_Initialise(DependentFieldMovingMesh,Err)
  CALL CMISSEquationsSet_DependentCreateStart(MovingMeshEquationsSet,DependentFieldMovingMeshUserNumber, & 
    & DependentFieldMovingMesh,Err)
  CALL CMISSField_VariableLabelSet(DependentFieldMovingMesh,CMISS_FIELD_U_VARIABLE_TYPE,"MovingMeshDF",Err)
  !Set the mesh component to be used by the field components.
  DO ComponentNumber=1,NumberOfDimensions
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldMovingMesh,CMISS_FIELD_U_VARIABLE_TYPE,ComponentNumber, & 
      & Mesh2ComponentNumberSpace,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentFieldMovingMesh,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,ComponentNumber, & 
      & Mesh2ComponentNumberSpace,Err)
  ENDDO
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSet_DependentCreateFinish(MovingMeshEquationsSet,Err)
  !Initialise dependent field moving mesh
  DO ComponentNumber=1,NumberOfDimensions
    CALL CMISSField_ComponentValuesInitialise(DependentFieldMovingMesh,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & ComponentNumber,InitialFieldMovingMesh(ComponentNumber),Err)
  ENDDO
  
  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> FLUID MATERIAL FIELD << == '
  !Create the equations set materials field variables for dynamic Navier-Stokes
  CALL CMISSField_Initialise(MaterialField2,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet2,MaterialField2UserNumber,MaterialField2,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet2,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & MaterialField2UserNumberMu,FluidDynamicViscosity,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & MaterialField2UserNumberRho,FluidDensity,Err)
    
  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> MATERIAL FIELD MOVING MESH << == '
  !Create the equations set materials field variables for moving mesh
  CALL CMISSField_Initialise(MaterialFieldMovingMesh,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(MovingMeshEquationsSet,MaterialFieldMovingMeshUserNumber, &
  & MaterialFieldMovingMesh,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(MovingMeshEquationsSet,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialFieldMovingMesh,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & MaterialFieldMovingMeshUserNumberK,MovingMeshParameterK,Err)
    
  !
  !================================================================================================================================
  !
  ! I N D E P E N D E N T _ F I E L D // mesh velocity for fluid equations set && mesh stiffness for moving mesh equations set
  
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID INDEPENDENT FIELD << == '
  !Create the equations set independent field variables for ALE Navier-Stokes
  CALL CMISSField_Initialise(IndependentField2,Err)
  CALL CMISSEquationsSet_IndependentCreateStart(EquationsSet2,IndependentField2UserNumber,IndependentField2,Err)
  CALL CMISSField_VariableLabelSet(IndependentField2,CMISS_FIELD_U_VARIABLE_TYPE,"FluidInDF",Err)
  !Set the mesh component to be used by the field components.
  DO ComponentNumber=1,NumberOfDimensions
    CALL CMISSField_ComponentMeshComponentSet(IndependentField2,CMISS_FIELD_U_VARIABLE_TYPE,ComponentNumber, & 
      & Mesh2ComponentNumberSpace,Err)
  ENDDO
  !Finish the equations set independent field variables
  CALL CMISSEquationsSet_IndependentCreateFinish(EquationsSet2,Err)
  
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INDEPENDENT FIELD MOVING MESH << == '
  !Create the equations set independent field variables for moving mesh
  CALL CMISSField_Initialise(IndependentFieldMovingMesh,Err)
  CALL CMISSEquationsSet_IndependentCreateStart(MovingMeshEquationsSet,IndependentFieldMovingMeshUserNumber, & 
    & IndependentFieldMovingMesh,Err)
  CALL CMISSField_VariableLabelSet(IndependentFieldMovingMesh,CMISS_FIELD_U_VARIABLE_TYPE,"MovingMeshInDF",Err)
  !Set the mesh component to be used by the field components.
  DO ComponentNumber=1,NumberOfDimensions
    CALL CMISSField_ComponentMeshComponentSet(IndependentFieldMovingMesh,CMISS_FIELD_U_VARIABLE_TYPE,ComponentNumber, & 
      & Mesh2ComponentNumberSpace,Err)
  ENDDO
  !Finish the equations set independent field variables
  CALL CMISSEquationsSet_IndependentCreateFinish(MovingMeshEquationsSet,Err)
  
  CALL CMISSField_ComponentValuesInitialise(IndependentFieldMovingMesh,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & IndependentFieldMovingMeshUserNumberK,MovingMeshParameterK,Err)
  !
  !================================================================================================================================
  !
  IF(GeometryCheck) THEN
    IF(ExampleFileProgressDiagnostics) WRITE(*,'(A)') "Exporting fields..."
    !Export initial fields
    CALL CMISSFields_Initialise(Fields1,Err)
    CALL CMISSFields_Create(Region1,Fields1,Err)
    CALL CMISSFields_NodesExport(Fields1,"GeometryCheckSolid","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields1,"GeometryCheckSolid","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields1,Err)
    CALL CMISSFields_Initialise(Fields2,Err)
    CALL CMISSFields_Create(Region2,Fields2,Err)
    CALL CMISSFields_NodesExport(Fields2,"GeometryCheckFluid","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(Fields2,"GeometryCheckFluid","FORTRAN",Err)
    CALL CMISSFields_Finalise(Fields2,Err)
    CALL CMISSFields_Initialise(FieldsI,Err)
    CALL CMISSFields_Create(Interface1,FieldsI,Err)
    CALL CMISSFields_NodesExport(FieldsI,"GeometryCheckInterface","FORTRAN",Err)
    CALL CMISSFields_ElementsExport(FieldsI,"GeometryCheckInterface","FORTRAN",Err)
    CALL CMISSFields_Finalise(FieldsI,Err)
    IF(ExampleFileProgressDiagnostics) WRITE(*,'(A)') "Field exported!"
  ENDIF
  
  !
  !================================================================================================================================
  !
  ! E Q U A T I O N S

  !Create the equations set equations for the first equations set
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID EQUATIONS << == '
  CALL CMISSEquations_Initialise(Equations1,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet1,Equations1,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations1,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet1,Err)

  !Create the equations set equations for the second equations set
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID EQUATIONS << == '
  CALL CMISSEquations_Initialise(Equations2,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet2,Equations2,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations2,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations lumping type
  CALL CMISSEquations_LumpingTypeSet(Equations2,CMISS_EQUATIONS_UNLUMPED_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
!  CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EquationsNavierStokesOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet2,Err)
  
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> MOVING MESH EQUATIONS << == '
  !Create the equations set equations
  CALL CMISSEquations_Initialise(EquationsMovingMesh,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(MovingMeshEquationsSet,EquationsMovingMesh,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(EquationsMovingMesh,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(EquationsMovingMesh,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(MovingMeshEquationsSet,Err)

  !
  !================================================================================================================================
  !
  ! I N T E R F A C E _ C O N D I T I O N

  !Create an interface condition between the two meshes
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE CONDITIONS << == '
  CALL CMISSInterfaceCondition_Initialise(InterfaceCondition1,Err)
  CALL CMISSInterfaceCondition_CreateStart(InterfaceCondition1UserNumber,Interface1,InterfaceGeometricField1, &
    & InterfaceCondition1,Err)
  !Specify the method for the interface condition
  CALL CMISSInterfaceCondition_MethodSet(InterfaceCondition1,CMISS_INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,Err)
  !Specify the type of interface condition operator
  CALL CMISSInterfaceCondition_OperatorSet(InterfaceCondition1,CMISS_INTERFACE_CONDITION_SOLID_FLUID_OPERATOR,Err)
  !Add in the dependent variables from the equations sets
  CALL CMISSInterfaceCondition_DependentVariableAdd(InterfaceCondition1,Mesh1Index,EquationsSet1,CMISS_FIELD_U_VARIABLE_TYPE,Err)
  CALL CMISSInterfaceCondition_DependentVariableAdd(InterfaceCondition1,Mesh2Index,EquationsSet2,CMISS_FIELD_U_VARIABLE_TYPE,Err)
  !Finish creating the interface condition
  CALL CMISSInterfaceCondition_CreateFinish(InterfaceCondition1,Err)

  !Create the Lagrange multipliers field
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE LAGRANGE FIELD << == '
  CALL CMISSField_Initialise(LagrangeField1,Err)
  CALL CMISSInterfaceCondition_LagrangeFieldCreateStart(InterfaceCondition1,LagrangeField1UserNumber,LagrangeField1,Err)
  CALL CMISSField_VariableLabelSet(LagrangeField1,CMISS_FIELD_U_VARIABLE_TYPE,"InterfaceLF",Err)
  !Finish the Lagrange multipliers field
  CALL CMISSInterfaceCondition_LagrangeFieldCreateFinish(InterfaceCondition1,Err)  

  !Create the interface condition equations
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE EQUATIONS << == '
  CALL CMISSInterfaceEquations_Initialise(InterfaceEquations1,Err)
  CALL CMISSInterfaceCondition_EquationsCreateStart(InterfaceCondition1,InterfaceEquations1,Err)
  !Set the interface equations sparsity
  CALL CMISSInterfaceEquations_SparsitySet(InterfaceEquations1,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the interface equations output
  CALL CMISSInterfaceEquations_OutputTypeSet(InterfaceEquations1,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !Finish creating the interface equations
  CALL CMISSInterfaceCondition_EquationsCreateFinish(InterfaceCondition1,Err)

  !
  !================================================================================================================================
  !
  ! P R O B L E M
  
  !Start the creation of a coupled problem.
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> COUPLED PROBLEM << == '
  CALL CMISSProblem_Initialise(CoupledProblem,Err)
  CALL CMISSProblem_CreateStart(CoupledProblemUserNumber,CoupledProblem,Err)
  !Set the problem to be a MULTI PHYSICS class - FINITE ELASTICITY NAVIER STOKES type - FINITE ELASTICITY NAVIER STOKES ALE subtype  
  CALL CMISSProblem_SpecificationSet(CoupledProblem,CMISS_PROBLEM_MULTI_PHYSICS_CLASS, &
    & CMISS_PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE,CMISS_PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE,Err)
  
  !Finish the creation of a problem.
  CALL CMISSProblem_CreateFinish(CoupledProblem,Err)  
  
  !
  !================================================================================================================================
  !
  ! C O N T R O L _ L O O P

  !Start the creation of the problem control loop for the coupled problem
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> PROBLEM CONTROL LOOP << == '
  CALL CMISSProblem_ControlLoopCreateStart(CoupledProblem,Err)
  
  !Time loop (main loop)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)  
  CALL CMISSProblem_ControlLoopGet(CoupledProblem,CMISS_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL CMISSControlLoop_LabelSet(ControlLoop,'TimeLoop',Err)
  CALL CMISSControlLoop_TimesSet(ControlLoop,StartTime,StopTime,TimeStepSize,Err)
  CALL CMISSControlLoop_TimeInputSet(ControlLoop,MaterialSpecification,Err)
  CALL CMISSControlLoop_TimeOutputSet(ControlLoop,OutputFrequency,Err)
  
  !Finish creating the problem control loop
  CALL CMISSProblem_ControlLoopCreateFinish(CoupledProblem,Err)
  
  !
  !================================================================================================================================
  !
  ! S O L V E R S

  !Start the creation of the problem solver for the coupled problem
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> PROBLEM SOLVERS << == '
  CALL CMISSProblem_SolversCreateStart(CoupledProblem,Err)
  
  !Linear solver for moving mesh
  CALL CMISSSolver_Initialise(LinearSolverMovingMesh,Err)
  CALL CMISSProblem_SolverGet(CoupledProblem,CMISS_CONTROL_LOOP_NODE, &
    & LinearSolverMovingMeshUserNumber,LinearSolverMovingMesh,Err)
  CALL CMISSSolver_OutputTypeSet(LinearSolverMovingMesh,LinearSolverMovingMesh_OutputType,Err)
  !TODO Additional settings
  
  !Solvers for coupled FiniteElasticity NavierStokes problem
  CALL CMISSSolver_Initialise(DynamicSolver,Err)
  CALL CMISSSolver_Initialise(NonlinearSolver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  !Get the dynamic ALE solver
  CALL CMISSProblem_SolverGet(CoupledProblem,CMISS_CONTROL_LOOP_NODE,DynamicSolverUserNumber, &
    & DynamicSolver,Err)
  CALL CMISSSolver_OutputTypeSet(DynamicSolver,DynamicSolver_OutputType,Err)
  CALL CMISSSolver_DynamicThetaSet(DynamicSolver,DynamicSolver_Theta,Err)
  !Get the dynamic nonlinear solver
  CALL CMISSSolver_DynamicNonlinearSolverGet(DynamicSolver,NonlinearSolver,Err)
  CALL CMISSSolver_NewtonLineSearchTypeSet(NonlinearSolver,CMISS_SOLVER_NEWTON_LINESEARCH_LINEAR,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolver,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED, &
    & Err)
  CALL CMISSSolver_NewtonMaximumFunctionEvaluationsSet(NonlinearSolver,MaxFunctionEvaluations,Err)
  CALL CMISSSolver_OutputTypeSet(NonlinearSolver,NonlinearSolver_OutputType,Err)
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(NonlinearSolver,AbsoluteTolerance,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(NonlinearSolver,RelativeTolerance,Err)
  CALL CMISSSolver_NewtonMaximumIterationsSet(NonlinearSolver,MaximumIterations,Err)
  CALL CMISSSolver_NewtonLineSearchAlphaSet(NonlinearSolver,LinesearchAlpha,Err)
  !Get the dynamic nonlinear linear solver
  CALL CMISSSolver_NewtonLinearSolverGet(NonlinearSolver,LinearSolver,Err)
  IF(.FALSE.) THEN
    CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolver,MaximumIterations,Err)
    CALL CMISSSolver_LinearIterativeAbsoluteToleranceSet(LinearSolver,AbsoluteTolerance,Err)
    CALL CMISSSolver_LinearIterativeRelativeToleranceSet(LinearSolver,RelativeTolerance,Err)
    CALL CMISSSolver_LinearIterativeDivergenceToleranceSet(LinearSolver,DivergenceTolerance,Err)
  ELSE
    CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  ENDIF
  CALL CMISSSolver_OutputTypeSet(LinearSolver,LinearSolver_OutputType,Err)
  
  !Finish the creation of the problem solver
  CALL CMISSProblem_SolversCreateFinish(CoupledProblem,Err)
  
  !
  !================================================================================================================================
  !
  ! S O L V E R _ E Q U A T I O N S

  CALL CMISSProblem_SolverEquationsCreateStart(CoupledProblem,Err)

  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> MOVING MESH SOLVER EQUATIONS << == '
  !Start the creation of the problem solver equations
  CALL CMISSSolver_Initialise(LinearSolverMovingMesh,Err)
  CALL CMISSSolverEquations_Initialise(LinearSolverMovingMeshEquations,Err)
  !Get the linear solver equations
  CALL CMISSProblem_SolverGet(CoupledProblem,CMISS_CONTROL_LOOP_NODE, &
    & LinearSolverMovingMeshUserNumber,LinearSolverMovingMesh,Err)
  CALL CMISSSolver_SolverEquationsGet(LinearSolverMovingMesh,LinearSolverMovingMeshEquations,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(LinearSolverMovingMeshEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(LinearSolverMovingMeshEquations,MovingMeshEquationsSet, &
    & MovingMeshEquationsSetIndex,Err)
    
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLVER EQUATIONS << == '
  CALL CMISSSolver_Initialise(DynamicSolver,Err)
  CALL CMISSSolverEquations_Initialise(CoupledSolverEquations,Err)
  !Get the dynamic solver equations
  CALL CMISSProblem_SolverGet(CoupledProblem,CMISS_CONTROL_LOOP_NODE,DynamicSolverUserNumber, &
    & DynamicSolver,Err)
  CALL CMISSSolver_SolverEquationsGet(DynamicSolver,CoupledSolverEquations,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(CoupledSolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(CoupledSolverEquations,EquationsSet1,EquationsSet1Index,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(CoupledSolverEquations,EquationsSet2,EquationsSet2Index,Err)
  CALL CMISSSolverEquations_InterfaceConditionAdd(CoupledSolverEquations,InterfaceCondition1,InterfaceCondition1Index,Err)
  !Set the time dependence of the interface matrix to determine the interface matrix coefficient in the solver matrix
  ! (basically position in big coupled matrix system)
  CALL CMISSInterfaceMatrices_TimeDependenceTypeSet(InterfaceCondition1, &
    & EquationsSet1Index,.TRUE., &
    & (/CMISS_INTERFACE_MATRIX_STATIC,CMISS_INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC/),Err)
  CALL CMISSInterfaceMatrices_TimeDependenceTypeSet(InterfaceCondition1, &
    & EquationsSet2Index,.TRUE., &
    & (/CMISS_INTERFACE_MATRIX_STATIC,CMISS_INTERFACE_MATRIX_STATIC/),Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(CoupledProblem,Err)

  !
  !================================================================================================================================
  !
  ! B O U N D A R Y _ C O N D I T I O N S

  !Start the creation of the equations set boundary conditions
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> BOUNDARY CONDITIONS << == '
  
  !SUBROUTINE CMISSBoundaryConditions_SetNodeObj(boundaryConditions,field,variableType,versionNumber,derivativeNumber, &
  !  & nodeUserNumber,componentNumber,condition,value,err)
  SELECT CASE(MaterialSpecification)
  CASE(Blood)
    CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
    CALL CMISSSolverEquations_BoundaryConditionsCreateStart(CoupledSolverEquations,BoundaryConditions,Err)
    !No displacement boundary for solid
    DO S=1,SIZE(DisplacementBoundaryNodeMappings)
      NodeNumber=DisplacementBoundaryNodeMappings(S)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition1,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,DisplacementBoundaryCondition(S),Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,2,CMISS_BOUNDARY_CONDITION_FIXED,DisplacementBoundaryCondition(S+SIZE(DisplacementBoundaryNodeMappings)),Err)
      ENDIF
    ENDDO
    
    !Linear inflow profile on left boundary  0 --> 1  (fluid)
    DO F=1,SIZE(VelocityBoundaryNodeMappings)
      NodeNumber=VelocityBoundaryNodeMappings(F)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        MMCondition=CMISS_BOUNDARY_CONDITION_FIXED
        IF(ANY(VelocityInletBoundaryNodeMappings.EQ.NodeNumber)) MMCondition=CMISS_BOUNDARY_CONDITION_FIXED_INLET
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,1,MMCondition,VelocityBoundaryCondition(F),Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,2,MMCondition,VelocityBoundaryCondition(F+SIZE(VelocityBoundaryNodeMappings)),Err)
      ENDIF
    ENDDO
    
    !Zero pressure on outflow boundaries
    DO P=1,SIZE(PressureBoundaryNodeMappings)
      NodeNumber=PressureBoundaryNodeMappings(P)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,NumberOfDimensions+1,CMISS_BOUNDARY_CONDITION_FIXED,PressureBoundaryCondition(P),Err)
      ENDIF
    ENDDO
    !Finish equations set boundary conditions
    CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(CoupledSolverEquations,Err)
    
    !Start the creation of the moving mesh boundary conditions
    CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsMovingMesh,Err)
    CALL CMISSSolverEquations_BoundaryConditionsCreateStart(LinearSolverMovingMeshEquations,BoundaryConditionsMovingMesh,Err)
    !Set which wall nodes are moved and which have a fixed position in space
    DO MM=1,SIZE(MovingMeshBoundaryElementNodeMappings)
      NodeNumber=MovingMeshBoundaryElementNodeMappings(MM)
      IF(COUNT(MovingMeshMovingNodeMappings==NodeNumber)/=0) THEN
        MMCondition=CMISS_BOUNDARY_CONDITION_MOVED_WALL
      ELSE
        MMCondition=CMISS_BOUNDARY_CONDITION_FIXED_WALL
      ENDIF
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh,CMISS_FIELD_U_VARIABLE_TYPE, &
          & 1, CMISS_NO_GLOBAL_DERIV,NodeNumber,1,MMCondition,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh,CMISS_FIELD_U_VARIABLE_TYPE, &
          & 1, CMISS_NO_GLOBAL_DERIV,NodeNumber,2,MMCondition,0.0_CMISSDP,Err)
      ENDIF
    ENDDO
    !Finish moving mesh boundary conditions
    CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(LinearSolverMovingMeshEquations,Err)
  !=================================================================================================================================
  CASE(Plate2D,Plate3D)
    IF(MaterialSpecification==Plate2D) THEN
      PRINT *, 'Boundary conditions for 2D plate reference case.'
    ELSE
      PRINT *, 'Boundary conditions for 3D plate reference case.'
    ENDIF
    CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
    CALL CMISSSolverEquations_BoundaryConditionsCreateStart(CoupledSolverEquations,BoundaryConditions,Err)
    !No displacement boundary for solid
    DO S=1,SIZE(NoDisplacementNodes)
      NodeNumber=NoDisplacementNodes(S)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition1,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,2,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        IF(NumberOfDimensions==3) THEN
          CALL CMISSBoundaryConditions_AddNode(BoundaryConditions,DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,3,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        ENDIF
      ENDIF
    ENDDO
    !Set outlet (zero) pressure nodes
    DO S=1,SIZE(OutletNodes)
      NodeNumber=OutletNodes(S)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,NumberOfDimensions+1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      ENDIF
    ENDDO
    !Set no-slip BC
    DO S=1,SIZE(NoSlipNodes)
      NodeNumber=NoSlipNodes(S)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,2,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        IF(NumberOfDimensions==3) THEN
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
            & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,3,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        ENDIF
      ENDIF
    ENDDO
    !Set slip BC
    DO S=1,SIZE(SlipNodes)
      NodeNumber=SlipNodes(S)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, & ! TODO TODO reference case has slip BC!
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,2,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        IF(NumberOfDimensions==3) THEN
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, & ! TODO TODO reference case has slip BC!
            & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,3,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        ENDIF
      ENDIF
    ENDDO
    !Inlet velocity nodes, small starting velocity in 1st coordinate direction
    DO S=1,SIZE(InletNodes)
      NodeNumber=InletNodes(S)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED_INLET,0.1_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,2,CMISS_BOUNDARY_CONDITION_FIXED_INLET,0.0_CMISSDP,Err)
        IF(NumberOfDimensions==3) THEN
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
            & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,3,CMISS_BOUNDARY_CONDITION_FIXED_INLET,0.0_CMISSDP,Err)
        ENDIF
      ENDIF
    ENDDO
    IF(NumberOfDimensions==2) THEN
      !Remove dof's at nodes where solid displacement and zero velocity is set (first n last interface node)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
        & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
        & 1,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
        & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
        & 1,2,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
        & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
        & NumberOfInterfaceNodes,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
        & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
        & NumberOfInterfaceNodes,2,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ELSE
      DO S=1,SIZE(LagrangeNodes)
        NodeNumber=LagrangeNodes(S)
        CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
            & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
            & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,2,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
            & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,3,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        ENDIF
      ENDDO
    ENDIF
    !Finish equations set boundary conditions
    CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(CoupledSolverEquations,Err)
    
    !Start the creation of the moving mesh boundary conditions
    CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsMovingMesh,Err)
    CALL CMISSSolverEquations_BoundaryConditionsCreateStart(LinearSolverMovingMeshEquations,BoundaryConditionsMovingMesh,Err)
    !Mesh nodes that are moving in 1st coordinate direction; fixed in 2nd coordinate direction
    DO S=1,SIZE(MovedYNodes)
      NodeNumber=MovedYNodes(S)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED_WALL,0.0_CMISSDP,Err)!TODO Maybe move nodes to have higher quality elements
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,2,CMISS_BOUNDARY_CONDITION_FIXED_WALL,0.0_CMISSDP,Err)
        IF(NumberOfDimensions==3) THEN
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
            & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,3,CMISS_BOUNDARY_CONDITION_FIXED_WALL,0.0_CMISSDP,Err)
        ENDIF
      ENDIF
    ENDDO
    !Mesh nodes that are moving wall nodes
    DO S=1,SIZE(MovedNodes)
      NodeNumber=MovedNodes(S)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_MOVED_WALL,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,2,CMISS_BOUNDARY_CONDITION_MOVED_WALL,0.0_CMISSDP,Err)
        IF(NumberOfDimensions==3) THEN
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
            & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,3,CMISS_BOUNDARY_CONDITION_MOVED_WALL,0.0_CMISSDP,Err)
        ENDIF
      ENDIF
    ENDDO
    !Mesh nodes that are fixed in space
    DO S=1,SIZE(FixedNodes)
      NodeNumber=FixedNodes(S)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition2,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED_WALL,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,2,CMISS_BOUNDARY_CONDITION_FIXED_WALL,0.0_CMISSDP,Err)
        IF(NumberOfDimensions==3) THEN
          CALL CMISSBoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
            & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,3,CMISS_BOUNDARY_CONDITION_FIXED_WALL,0.0_CMISSDP,Err)
        ENDIF
      ENDIF
    ENDDO
    !Finish moving mesh boundary conditions
    CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(LinearSolverMovingMeshEquations,Err)
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid material specification for boundary conditions.")
  END SELECT
  

  !
  !================================================================================================================================
  !
  ! S O L V E
  PRINT *, ' == >> SOLVING PROBLEM << == '
  CALL CMISSProblem_Solve(CoupledProblem,Err)
  
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  e=etime(t)

  PRINT *, "Program successfully completed in ",e," seconds."

  STOP
 
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR
     
END PROGRAM CoupledFluidSolidExample
