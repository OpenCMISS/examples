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

!> Main program
PROGRAM FortranExample

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWINCMISS
#endif

  IMPLICIT NONE

  !Test program parameters==========================================================================================================
  !Material specification
  INTEGER(CMISSIntg), PARAMETER :: Plate2D=2
  INTEGER(CMISSIntg), PARAMETER :: Plate3D=3
  
  REAL(CMISSDP) :: SolidDensity
  REAL(CMISSDP), ALLOCATABLE :: Gravity(:)
  
  INTEGER(CMISSIntg) :: NumberOfFluidNodes
  INTEGER(CMISSIntg) :: NumberOfFluidElements
  INTEGER(CMISSIntg) :: NumberOfSolidNodes
  INTEGER(CMISSIntg) :: NumberOfSolidElements
  INTEGER(CMISSIntg) :: NumberOfInterfaceNodes
  INTEGER(CMISSIntg) :: NumberOfInterfaceElements
  
  INTEGER(CMISSIntg), ALLOCATABLE :: SolidInterfaceElements(:),FluidInterfaceElements(:)

  ! 2D plate variables TODO Share variables for different geometries
  INTEGER(CMISSIntg), ALLOCATABLE :: SolidNodeNumbers(:),FluidNodeNumbers(:),SolidElementNodes(:,:),FluidElementNodes(:,:), &
    & InterfaceElementNodes(:,:),ConnectedInterfaceNodes(:,:), &
    & InterfaceNodeNumbersForGeometry(:), &
    & NoDisplacementNodes(:),FixedNodes(:),FixedZNodes(:),MovedNodes(:),MovedYNodes(:),InletNodes(:),NoSlipNodes(:), &
    & OutletNodes(:),SlipNodes(:),InterfaceInterfaceNodeInformationNE(:,:),SolidInterfaceNodeInformationNE(:,:), &
    & FluidInterfaceNodeInformationNE(:,:),LagrangeNodes(:),SlipNodesTop(:),SlipNodesRightLeft(:)

  REAL(CMISSDP), ALLOCATABLE :: SolidGeometryX(:),SolidGeometryY(:),SolidGeometryZ(:), &
    & FluidGeometryX(:),FluidGeometryY(:),FluidGeometryZ(:), &
    & InterfaceGeometryX(:),InterfaceGeometryY(:),InterfaceGeometryZ(:), &
    & SolidXi2(:,:),SolidXi3(:,:),FluidXi2(:,:),FluidXi3(:,:), &
    & SolidInterfaceNodeInformationXi(:,:),FluidInterfaceNodeInformationXi(:,:)

  !variables for timing
  REAL :: e,t(2)
  
  INTEGER(CMISSIntg), PARAMETER :: SolidCoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FluidCoordinateSystemUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: InterfaceCoordinateSystemUserNumber=4
  
  INTEGER(CMISSIntg), PARAMETER :: SolidRegionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: FluidRegionUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: InterfaceUserNumber=40
  
  INTEGER(CMISSIntg), PARAMETER :: SolidMeshUserNumber=19
  INTEGER(CMISSIntg), PARAMETER :: FluidMeshUserNumber=20
  INTEGER(CMISSIntg), PARAMETER :: InterfaceMeshUserNumber=22
  INTEGER(CMISSIntg), PARAMETER :: MovingMeshUserNumber=107
  
  INTEGER(CMISSIntg), PARAMETER :: SolidDecompositionUserNumber=24
  INTEGER(CMISSIntg), PARAMETER :: FluidDecompositionUserNumber=25
  INTEGER(CMISSIntg), PARAMETER :: InterfaceDecompositionUserNumber=27
  
  INTEGER(CMISSIntg), PARAMETER :: SolidGeometricFieldUserNumber=29
  INTEGER(CMISSIntg), PARAMETER :: FluidGeometricFieldUserNumber=30
  INTEGER(CMISSIntg), PARAMETER :: InterfaceGeometricFieldUserNumber=32
  
  INTEGER(CMISSIntg), PARAMETER :: SolidEquationsSetUserNumber=34
  INTEGER(CMISSIntg), PARAMETER :: FluidEquationsSetUserNumber=35
  INTEGER(CMISSIntg), PARAMETER :: InterfaceConditionUserNumber=42
  
  INTEGER(CMISSIntg), PARAMETER :: SolidDependentFieldUserNumber=37
  INTEGER(CMISSIntg), PARAMETER :: FluidDependentFieldUserNumber=38
  INTEGER(CMISSIntg), PARAMETER :: LagrangeFieldUserNumber=44
  
  INTEGER(CMISSIntg), PARAMETER :: CoupledProblemUserNumber=46
  
  INTEGER(CMISSIntg), PARAMETER :: SolidEquationsSetFieldUserNumber=49
  INTEGER(CMISSIntg), PARAMETER :: FluidEquationsSetFieldUserNumber=50
  
  !======
  INTEGER(CMISSIntg), PARAMETER :: BasisSpaceSolidUserNumber=100
  INTEGER(CMISSIntg), PARAMETER :: BasisDisplacementUserNumber=101
  INTEGER(CMISSIntg), PARAMETER :: BasisHydrostaticPressureUserNumber=102
  INTEGER(CMISSIntg), PARAMETER :: BasisSpaceFluidUserNumber=103
  INTEGER(CMISSIntg), PARAMETER :: BasisVelocityUserNumber=104
  INTEGER(CMISSIntg), PARAMETER :: BasisPressureUserNumber=105
  INTEGER(CMISSIntg), PARAMETER :: InterfaceBasisUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: InterfaceMappingBasisUserNumber=47
  
  INTEGER(CMISSIntg), PARAMETER :: FibreFieldUserNumber=106
  INTEGER(CMISSIntg), PARAMETER :: MovingMeshFieldUserNumber=108
  INTEGER(CMISSIntg), PARAMETER :: MovingMeshEquationsSetUserNumber=109
  INTEGER(CMISSIntg), PARAMETER :: FluidMaterialFieldUserNumber=117
  INTEGER(CMISSIntg), PARAMETER :: SolidMaterialFieldUserNumber=110
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldMovingMeshUserNumber=111
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=112
  
  INTEGER(CMISSIntg), PARAMETER :: DynamicSolverIndex=1
  INTEGER(CMISSIntg), PARAMETER :: LinearSolverMovingMeshIndex=2
  
  INTEGER(CMISSIntg), PARAMETER :: FluidMaterialFieldComponentMu=1
  INTEGER(CMISSIntg), PARAMETER :: FluidMaterialFieldComponentRho=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialFieldMovingMeshUserNumberK=1
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldMovingMeshUserNumberK=1
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldMovingMeshUserNumber=118
  INTEGER(CMISSIntg), PARAMETER :: MaterialFieldMovingMeshUserNumber=119
  INTEGER(CMISSIntg), PARAMETER :: IndependentField2UserNumber=120
  INTEGER(CMISSIntg), PARAMETER :: IndependentFieldMovingMeshUserNumber=121
  INTEGER(CMISSIntg), PARAMETER :: LinearSolverMovingMeshEquationsUserNumber=122
  
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_OF_ARGUMENTS
  INTEGER(CMISSIntg) :: InterpolationTypeInterface,NumberOfGaussXi,NUMBER_OF_NODE_XI,NumberOfDimensions,component_idx, &
    & NumberOfGaussXiSpace,NumberOfGaussXiVelocity,NumberOfGaussXiPressure,arraySize

  INTEGER(CMISSIntg) :: SolidEquationsSetIndex=1
  INTEGER(CMISSIntg) :: FluidEquationsSetIndex=2
  INTEGER(CMISSIntg) :: InterfaceConditionIndex=1
  INTEGER(CMISSIntg) :: SolidMeshIndex=1
  INTEGER(CMISSIntg) :: FluidMeshIndex=2
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: PressureMeshComponent
  INTEGER(CMISSIntg) :: FluidMeshComponentNumberVelocity,ComponentNumber
  
  INTEGER(CMISSIntg) :: OutputFrequency
  INTEGER(CMISSIntg) :: DynamicSolver_OutputType
  INTEGER(CMISSIntg) :: NonlinearSolver_OutputType
  INTEGER(CMISSIntg) :: LinearSolver_OutputType
  INTEGER(CMISSIntg) :: LinearSolverMovingMesh_OutputType
  INTEGER(CMISSIntg) :: MaximumIterations,MaxFunctionEvaluations
  INTEGER(CMISSIntg) :: RestartValue
  
  INTEGER(CMISSIntg) :: EquationsNavierStokesOutput,InterfaceMeshComponentNumber,InterpolationTypeDisplacement, &
    & InterpolationTypeHydrostaticPressure,InterpolationTypePressure,InterpolationTypeSpace,InterpolationTypeVelocity, &
    & MovingMeshEquationsSetIndex,Mesh1ComponentNumberDisplacement,Mesh1ComponentNumberHydrostaticPressure, &
    & Mesh1ComponentNumberSpace,Components(3),NodeDomain,Mesh2ComponentNumberPressure,Mesh2ComponentNumberSpace, &
    & MeshNumberOfComponents
  
  ! LOOP INTEGERS
  INTEGER(CMISSIntg) :: NodeNumber,S,ElementIndex,NodeIndex,MaterialSpecification,LocalNodeIndex
  
  !check below for units
  REAL(CMISSDP) :: XI2(2),XI3(3)
  REAL(CMISSDP) :: MovingMeshParameterK
  REAL(CMISSDP) :: FluidDynamicViscosity
  REAL(CMISSDP) :: FluidDensity
  REAL(CMISSDP) :: YoungsModulus
  REAL(CMISSDP) :: PoissonsRatio
  REAL(CMISSDP) :: ShearModulus
  REAL(CMISSDP) :: BulkModulus
  REAL(CMISSDP) :: MooneyRivlin1
  REAL(CMISSDP) :: MooneyRivlin2

  REAL(CMISSDP) :: InitialFieldNavierStokes(3)
  REAL(CMISSDP) :: InitialFieldMovingMesh(3)
  REAL(CMISSDP) :: DivergenceTolerance
  REAL(CMISSDP) :: RelativeTolerance
  REAL(CMISSDP) :: AbsoluteTolerance
  REAL(CMISSDP) :: LinesearchAlpha

  REAL(CMISSDP) :: StartTime
  REAL(CMISSDP) :: StopTime
  REAL(CMISSDP) :: DynamicSolver_Theta
  REAL(CMISSDP) :: TimeStepSize

  LOGICAL :: FileReadDiagnostics=.FALSE.
  LOGICAL :: ExampleFileProgressDiagnostics=.FALSE.
  LOGICAL :: GeometryCheck=.FALSE.
  LOGICAL :: GravityFlag=.FALSE.
  LOGICAL :: CheckWithoutInterfaceCondition=.FALSE.
  LOGICAL :: SetupOutput=.FALSE.
  
  !CMISS variables
  TYPE(CMISSBasisType) :: BasisSpaceSolid,BasisDisplacement,BasisHydrostaticPressure, &
    & BasisSpaceFluid,BasisVelocity,BasisPressure,InterfaceBasis1,InterfaceMappingBasis1
  TYPE(CMISSNodesType) :: SolidNodes,FluidNodes,InterfaceNodes
  TYPE(CMISSMeshElementsType) :: SolidMeshElementsSpace,SolidMeshElementsDisplacement,SolidMeshElementsHydrostaticPressure, &
    & FluidMeshElementsSpace,FluidMeshElementsVelocity,FluidMeshElementsPressure,InterfaceMeshElements
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions,BoundaryConditionsMovingMesh
  TYPE(CMISSCoordinateSystemType) :: SolidCoordinateSystem,FluidCoordinateSystem,InterfaceCoordinateSystem, &
    & WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: SolidDecomposition,FluidDecomposition,InterfaceDecomposition
  TYPE(CMISSEquationsType) :: Equations1,Equations2,EquationsMovingMesh
  TYPE(CMISSEquationsSetType) :: SolidEquationsSet,FluidEquationsSet,MovingMeshEquationsSet
  TYPE(CMISSFieldType) :: GeometricField1,GeometricField2,InterfaceGeometricField1, &
    & DependentField1,DependentField2,LagrangeField1,EquationsSetField1,EquationsSetField2, &
    & DependentFieldMovingMesh,MaterialFieldMovingMesh,IndependentField2,IndependentFieldMovingMesh, &
    & MaterialField1,MaterialField2,EquationsSetFieldMovingMesh,SourceField1
  TYPE(CMISSFieldsType) :: Fields1,Fields2,FieldsI
  TYPE(CMISSInterfaceType) :: Interface1
  TYPE(CMISSInterfaceConditionType) :: InterfaceCondition
  TYPE(CMISSInterfaceEquationsType) :: InterfaceEquations
  TYPE(CMISSInterfaceMeshConnectivityType) :: InterfaceMeshConnectivity1
  TYPE(CMISSMeshType) :: Mesh1,Mesh2,InterfaceMesh1
  TYPE(CMISSProblemType) :: CoupledProblem
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSRegionType) :: Region1,Region2,WorldRegion
  TYPE(CMISSSolverType) :: DynamicSolver,NonlinearSolver,LinearSolver,LinearSolverMovingMesh
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

  !
  !=================================================================================================================================
  !
  ! C O M M A N D L I N E _ A R G U M E N T S
  
  NUMBER_OF_ARGUMENTS = COMMAND_ARGUMENT_COUNT()
  !Not used, read in specifics from files
  IF(NUMBER_OF_ARGUMENTS == 0) THEN
    !quadratic interpolation for fluid velocities, solid displacements
    InterpolationTypeDisplacement=2
    InterpolationTypeSpace=2
    InterpolationTypeVelocity=2
    !linear interpolation for pressures
    InterpolationTypePressure=1
    !Interpolation type on interface matches interpolation type for fluid velocities, solid displacements
    InterpolationTypeInterface=InterpolationTypeDisplacement
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
  !CALL CMISSOutputSetOn("Testing",Err)
  !Set diganostics for testing
  !CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["SOLVER_MAPPING_CALCULATE", &
  !  & "SOLVER_MATRIX_STRUCTURE_CALCULATE"],Err)
  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)
  !
  !=================================================================================================================================
  !
  ! I N I T I A L _ V A L U E S _ & _ P R O B L E M _ C O N S T A N T S
  !Set initial values
  InitialFieldNavierStokes(1)=0.0_CMISSDP
  InitialFieldNavierStokes(2)=0.0_CMISSDP
  InitialFieldNavierStokes(3)=0.0_CMISSDP
  InitialFieldMovingMesh(1)=0.0_CMISSDP
  InitialFieldMovingMesh(2)=0.0_CMISSDP
  InitialFieldMovingMesh(3)=0.0_CMISSDP
  !Set output parameters
  SetupOutput=.TRUE.
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LinearSolverMovingMesh_OutputType=CMISS_SOLVER_NO_OUTPUT
  DynamicSolver_OutputType=CMISS_SOLVER_NO_OUTPUT
  LinearSolver_OutputType=CMISS_SOLVER_NO_OUTPUT
  NonlinearSolver_OutputType=CMISS_SOLVER_NO_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EquationsNavierStokesOutput=CMISS_EQUATIONS_NO_OUTPUT
  !Set result output parameter
  OutputFrequency=1
  
  !Choose 2D or 3D case
  MaterialSpecification=Plate3D
  !Set solver parameters
  RelativeTolerance=1.0E-4_CMISSDP !default: 1.0E-05_CMISSDP
  AbsoluteTolerance=1.0E-4_CMISSDP !default: 1.0E-10_CMISSDP
  DivergenceTolerance=1.0E5 !default: 1.0E5
  MaximumIterations=100000000 !default: 100000
  MaxFunctionEvaluations=100000
  RestartValue=30 !default: 30
  LinesearchAlpha=1.0_CMISSDP
  MovingMeshParameterK=1.0 !default
  DynamicSolver_Theta=1.0_CMISSDP
  StartTime=0.0_CMISSDP
  StopTime=2000.0_CMISSDP
  TimeStepSize=1.0_CMISSDP
  IF(MaterialSpecification==Plate2D) THEN
    StopTime=0.1_CMISSDP
    TimeStepSize=0.05_CMISSDP
  ENDIF
  
  !Material properties
  !NOTE: USE OF SI UNITS unless comment
  !Low density fluid, rubber-like solid
  FluidDynamicViscosity=0.05! kg/(m.s)
  FluidDensity=100! kg/m3
  SolidDensity=300! kg/m3
  YoungsModulus=2.3E4! Pa
  PoissonsRatio=0.49
  !Neo-Hookean material law
  ShearModulus=YoungsModulus/(2.0_CMISSDP*(1.0_CMISSDP+PoissonsRatio)) ! N/m2
  BulkModulus=YoungsModulus/(3.0_CMISSDP*(1.0_CMISSDP-2.0_CMISSDP*PoissonsRatio))
  !MooneyRivlin1=0.5_CMISSDP*ShearModulus ! N/m2
  MooneyRivlin1=0.0595_CMISSDP
  MooneyRivlin2=0.0_CMISSDP
  
  
  !Set geometric dimension n gravity
  SELECT CASE(MaterialSpecification)
  CASE(Plate2D)
    NumberOfDimensions=2
    ALLOCATE(Gravity(NumberOfDimensions))
    Gravity(:)=[0.0_CMISSDP,9.81_CMISSDP] ! m/s2
  CASE(Plate3D)
    NumberOfDimensions=3
    CheckWithoutInterfaceCondition=.FALSE.!if set to true we remove all Lagrange field dofs by setting them as zero dirichlet BC
    IF(CheckWithoutInterfaceCondition) THEN
      GravityFlag=.TRUE.
      ALLOCATE(Gravity(NumberOfDimensions))
      Gravity(:)=[0.0_CMISSDP,9.81_CMISSDP,0.0_CMISSDP] ! m/s2
    ENDIF
  END SELECT
    
  !Reading geometry from files
  arraySize=0
  PRINT *, 'Reading geometry from files..'
  !Read node numbers
  SELECT CASE(MaterialSpecification)
  CASE(Plate2D)
    OPEN (unit=1, file='./Input/2D/Solid/Plate2DSolidNodes.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Solid/Plate3DSolidNodes.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/Fluid/Plate2DFluidNodes.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Fluid/Plate3DFluidNodes.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/Interface/Plate2DInterfaceNodes.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Interface/Plate3DInterfaceNodes.txt', status='old', action='read')
  ENDSELECT
  READ(1, *), arraySize
  IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
  ALLOCATE(InterfaceNodeNumbersForGeometry(arraySize))
  DO S=1,arraySize
    READ(1,*) InterfaceNodeNumbersForGeometry(S)
  ENDDO
  arraySize=0
  CLOSE(1)
  IF(FileReadDiagnostics) PRINT *, 'Interface nodes: ',InterfaceNodeNumbersForGeometry
  
  !Read geometry
  IF(NumberOfDimensions==3) THEN
    OPEN (unit=1, file='./Input/3D/Solid/Plate3DSolidX.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/Solid/Plate2DSolidY.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Solid/Plate3DSolidY.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/Solid/Plate2DSolidZ.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Solid/Plate3DSolidZ.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/3D/Fluid/Plate3DFluidX.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/Fluid/Plate2DFluidY.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Fluid/Plate3DFluidY.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/Fluid/Plate2DFluidZ.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Fluid/Plate3DFluidZ.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/3D/Interface/Plate3DInterfaceX.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/Interface/Plate2DInterfaceY.txt', &
      & status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Interface/Plate3DInterfaceY.txt', &
      & status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/Interface/Plate2DInterfaceZ.txt', &
      & status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Interface/Plate3DInterfaceZ.txt', &
      & status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/Solid/Plate2DSolidElementNodeNumbers.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Solid/Plate3DSolidElementNodeNumbers.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/Fluid/Plate2DFluidElementNodeNumbers.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Fluid/Plate3DFluidElementNodeNumbers.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/Interface/Plate2DInterfaceElementNodeNumbers.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Interface/Plate3DInterfaceElementNodeNumbers.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/Solid/Plate2DSolidXi1.txt', status='old', action='read')
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(SolidXi2(arraySize,3))
    DO S=1,arraySize
      READ(1,*) SolidXi2(S,:)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Solid xi2: ',SolidXi2
    
    OPEN (unit=1, file='./Input/2D/Solid/Plate2DSolidXi2.txt', status='old', action='read')
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(SolidXi3(arraySize,3))
    DO S=1,arraySize
      READ(1,*) SolidXi3(S,:)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Solid xi3: ',SolidXi3
    
    OPEN (unit=1, file='./Input/2D/Fluid/Plate2DFluidXi1.txt', status='old', action='read')
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(FluidXi2(arraySize,3))
    DO S=1,arraySize
      READ(1,*) FluidXi2(S,:)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Fluid xi2: ',FluidXi2
    
    OPEN (unit=1, file='./Input/2D/Fluid/Plate2DFluidXi2.txt', status='old', action='read')
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(FluidXi3(arraySize,3))
    DO S=1,arraySize
      READ(1,*) FluidXi3(S,:)
    ENDDO
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Fluid xi3: ',FluidXi3
    
    OPEN (unit=1, file='./Input/2D/Interface/Plate2DInterfaceSolidElements.txt', status='old', action='read')
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(SolidInterfaceElements(arraySize))
    READ(1,*) SolidInterfaceElements(:)
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Solid interface element numbers: ',SolidInterfaceElements
    
    OPEN (unit=1, file='./Input/2D/Interface/Plate2DInterfaceFluidElements.txt', status='old', action='read')
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(FluidInterfaceElements(arraySize))
    READ(1,*) FluidInterfaceElements(:)
    arraySize=0
    CLOSE(1)
    IF(FileReadDiagnostics) PRINT *, 'Fluid interface element numbers: ',FluidInterfaceElements
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Solid/Plate3DSolidInterfaceNodeInformationNodeElement.txt',status='old', action='read')
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(SolidInterfaceNodeInformationNE(arraySize,2))
    DO S=1,arraySize
      READ(1,*) SolidInterfaceNodeInformationNE(S,:)
    ENDDO
    arraySize=0
    CLOSE(1)
    
    OPEN (unit=1, file='./Input/3D/Fluid/Plate3DFluidInterfaceNodeInformationNodeElement.txt', status='old', action='read')
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(FluidInterfaceNodeInformationNE(arraySize,2))
    DO S=1,arraySize
      READ(1,*) FluidInterfaceNodeInformationNE(S,:)
    ENDDO
    arraySize=0
    CLOSE(1)
    
    OPEN (unit=1, file='./Input/3D/Interface/Plate3DIInterfaceNodeInformationNodeElement.txt',status='old', action='read')
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(InterfaceInterfaceNodeInformationNE(arraySize,2))
    DO S=1,arraySize
      READ(1,*) InterfaceInterfaceNodeInformationNE(S,:)
    ENDDO
    arraySize=0
    CLOSE(1)
    
    OPEN (unit=1, file= './Input/3D/Solid/Plate3DSolidInterfaceNodeInformationXi.txt', status='old', action='read')
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(SolidInterfaceNodeInformationXi(arraySize,3))
    DO S=1,arraySize
      READ(1,*) SolidInterfaceNodeInformationXi(S,:)
    ENDDO
    arraySize=0
    CLOSE(1)
    
    OPEN (unit=1, file= './Input/3D/Fluid/Plate3DFluidInterfaceNodeInformationXi.txt', status='old', action='read')
    READ(1, *), arraySize
    IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
    ALLOCATE(FluidInterfaceNodeInformationXi(arraySize,3))
    DO S=1,arraySize
      READ(1,*) FluidInterfaceNodeInformationXi(S,:)
    ENDDO
    arraySize=0
    CLOSE(1)
    
    OPEN (unit=1, file='./Input/3D/BC/Plate3DlagrangeNodes.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/Interface/Plate2DSortedInterfaceNodes.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/Interface/Plate3DSortedInterfaceNodes.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/BC/Plate2DdisplacementBC.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/BC/Plate3DdisplacementBC.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/BC/Plate2DfixedNodesBC.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/BC/Plate3DfixedNodesBC.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/BC/Plate2DfixedZNodesBC.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/BC/Plate3DfixedZNodesBC.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/BC/Plate2DmovedNodesBC.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/BC/Plate3DmovedNodesBC.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/BC/Plate2DmovedYNodesBC.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/BC/Plate3DmovedYNodesBC.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/BC/Plate2DinletBC.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/BC/Plate3DinletBC.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/BC/Plate2DnoSlipBC.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/BC/Plate3DnoSlipBC.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/BC/Plate2DslipBC.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/BC/Plate3DslipTop.txt', status='old', action='read')
  ENDSELECT
  READ(1, *), arraySize
  IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
  ALLOCATE(SlipNodesTop(arraySize))
  READ(1,*) SlipNodesTop(:)
  arraySize=0
  CLOSE(1)
  IF(FileReadDiagnostics) PRINT *, 'Slip nodes top: ',SlipNodesTop
  
  SELECT CASE(MaterialSpecification)
  CASE(Plate2D)
    OPEN (unit=1, file='./Input/2D/BC/Plate2DslipBC.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/BC/Plate3DslipRightLeft.txt', status='old', action='read')
  ENDSELECT
  READ(1, *), arraySize
  IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
  ALLOCATE(SlipNodesRightLeft(arraySize))
  READ(1,*) SlipNodesRightLeft(:)
  arraySize=0
  CLOSE(1)
  IF(FileReadDiagnostics) PRINT *, 'Slip nodes r/l: ',SlipNodesRightLeft
  
  SELECT CASE(MaterialSpecification)
  CASE(Plate2D)
    OPEN (unit=1, file='./Input/2D/BC/Plate2DpressureBC.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/BC/Plate3DpressureBC.txt', status='old', action='read')
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
    OPEN (unit=1, file='./Input/2D/BC/Plate2DslipBC.txt', status='old', action='read')
  CASE(Plate3D)
    OPEN (unit=1, file='./Input/3D/BC/Plate3DslipBC.txt', status='old', action='read')
  ENDSELECT
  READ(1,*), arraySize
  IF(arraySize==0) CALL HANDLE_ERROR("Incorrect number.")
  ALLOCATE(SlipNodes(arraySize))
  READ(1,*) SlipNodes(:)
  arraySize=0
  CLOSE(1)
  IF(FileReadDiagnostics) PRINT *, 'Slip BC nodes: ',SlipNodes
  
  NumberOfSolidNodes=SIZE(SolidNodeNumbers)
  NumberOfSolidElements=SIZE(SolidElementNodes(:,1))
  NumberOfFluidNodes=SIZE(FluidNodeNumbers)
  NumberOfFluidElements=SIZE(FluidElementNodes(:,1))
  NumberOfInterfaceNodes=SIZE(InterfaceNodeNumbersForGeometry)
  NumberOfInterfaceElements=SIZE(InterfaceElementNodes(:,1))
  
  PRINT *, "Finished reading geometry."
  
  IF(SetupOutput) THEN
1   FORMAT(1X,A,F10.2)
2   FORMAT(1X,A,F10.2,3X,A)
3   FORMAT(1X,A,F10.2,1X,A,F5.2,1X,A,F5.2,1X,A)
4   FORMAT(1X,A,F10.2,1X,A,F5.2,1X,A)
5   FORMAT(1X,A,I7)
    PRINT *," "
    PRINT *,"==========SUMMARY=========="
    PRINT *," "
    PRINT 1,"Start time:",StartTime
    PRINT 1,"Stop time:",StopTime
    PRINT 1,"Time increment:",TimeStepSize
    PRINT *," "
    PRINT *,"FLUID:"
    PRINT 2,"Dynamic viscosity:",FluidDynamicViscosity,"kg/(m.s)"
    PRINT 2,"Density:",FluidDensity,"kg/m3"
    IF(NumberOfDimensions==3) THEN
      PRINT 3,"Domain:",FluidGeometryX(NumberOfFluidNodes),"x",FluidGeometryY(NumberOfFluidNodes),"x", &
        & FluidGeometryZ(NumberOfFluidNodes),"m3"
    ELSE
      PRINT 4,"Domain:",FluidGeometryY(NumberOfFluidNodes),"x",FluidGeometryZ(NumberOfFluidNodes),"m2"
    ENDIF
    PRINT 5, "Number of nodes:",NumberOfFluidNodes
    PRINT *," "
    PRINT *,"SOLID:"
    PRINT 2,"Density:",SolidDensity,"kg/m3"
    PRINT 2,"Young's modulus:",YoungsModulus,"Pa"
    PRINT 1,"Poisson's ratio:",PoissonsRatio
    PRINT 1,"Neo-Hookean constant:",MooneyRivlin1
    IF(NumberOfDimensions==3) THEN
      PRINT 3,"Domain:",SolidGeometryX(NumberOfSolidNodes)-SolidGeometryX(1),"x", &
        & SolidGeometryY(NumberOfSolidNodes)-SolidGeometryY(1),"x", &
        & SolidGeometryZ(NumberOfSolidNodes)-SolidGeometryZ(1),"m3"
    ELSE
      PRINT 4,"Domain:",SolidGeometryY(NumberOfSolidNodes)-SolidGeometryY(1),"x", &
        & SolidGeometryZ(NumberOfSolidNodes)-SolidGeometryZ(1),"m2"
    ENDIF
    PRINT 5, "Number of nodes:",NumberOfSolidNodes
    PRINT *," "
    PRINT *,"INTERFACE:"
    PRINT 5, "Number of interface nodes:",NumberOfInterfaceNodes
    PRINT *," "
    PRINT *,"==========================="
    PRINT *," "
  ENDIF
  !
  !=================================================================================================================================
  !
  ! C O O R D I N A T E _ S Y S T E M S
  !Create a new RC coordinate system for the solid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID COORDINATE SYSTEM << == '
  CALL CMISSCoordinateSystem_Initialise(SolidCoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(SolidCoordinateSystemUserNumber,SolidCoordinateSystem,Err)
  CALL CMISSCoordinateSystem_DimensionSet(SolidCoordinateSystem,NumberOfDimensions,Err)
  CALL CMISSCoordinateSystem_CreateFinish(SolidCoordinateSystem,Err)
  !Create a new RC coordinate system for the fluid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID COORDINATE SYSTEM << == '
  CALL CMISSCoordinateSystem_Initialise(FluidCoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(FluidCoordinateSystemUserNumber,FluidCoordinateSystem,Err)
  CALL CMISSCoordinateSystem_DimensionSet(FluidCoordinateSystem,NumberOfDimensions,Err)
  CALL CMISSCoordinateSystem_CreateFinish(FluidCoordinateSystem,Err)
  !Create a new RC coordinate system for the interface region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE COORDINATE SYSTEM << == '
  CALL CMISSCoordinateSystem_Initialise(InterfaceCoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(InterfaceCoordinateSystemUserNumber,InterfaceCoordinateSystem,Err)
  CALL CMISSCoordinateSystem_DimensionSet(InterfaceCoordinateSystem,NumberOfDimensions,Err)
  CALL CMISSCoordinateSystem_CreateFinish(InterfaceCoordinateSystem,Err)
  !
  !=================================================================================================================================
  !
  ! R E G I O N S
  !Create the solid region and set the regions coordinate system to the RC coordinate system that we have created
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID REGION << == '
  CALL CMISSRegion_Initialise(Region1,Err)
  CALL CMISSRegion_CreateStart(SolidRegionUserNumber,WorldRegion,Region1,Err)
  CALL CMISSRegion_LabelSet(Region1,"SolidRegion",Err)
  CALL CMISSRegion_CoordinateSystemSet(Region1,SolidCoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region1,Err)
  !Create the fluid region and set the regions coordinate system to the RC coordinate system that we have created
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID REGION << == '
  CALL CMISSRegion_Initialise(Region2,Err)
  CALL CMISSRegion_CreateStart(FluidRegionUserNumber,WorldRegion,Region2,Err)
  CALL CMISSRegion_LabelSet(Region2,"FluidRegion",Err)
  CALL CMISSRegion_CoordinateSystemSet(Region2,FluidCoordinateSystem,Err)
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
  CALL CMISSMeshElements_Initialise(SolidMeshElementsSpace,Err)
  CALL CMISSMeshElements_Initialise(SolidMeshElementsDisplacement,Err)
  CALL CMISSMeshElements_Initialise(SolidMeshElementsHydrostaticPressure,Err)
  
  Mesh1ComponentNumberSpace=1
  Mesh1ComponentNumberDisplacement=1
  Mesh1ComponentNumberHydrostaticPressure=1
  
  CALL CMISSMeshElements_CreateStart(Mesh1,Mesh1ComponentNumberSpace,BasisSpaceSolid,SolidMeshElementsSpace,Err)
  DO ElementIndex=1,NumberOfSolidElements
    CALL CMISSMeshElements_NodesSet(SolidMeshElementsSpace,ElementIndex, &
      & SolidElementNodes(ElementIndex,:),Err)
  ENDDO
  CALL CMISSMeshElements_CreateFinish(SolidMeshElementsSpace,Err)
  
  SolidMeshElementsDisplacement=SolidMeshElementsSpace
  Mesh1ComponentNumberHydrostaticPressure=Mesh1ComponentNumberDisplacement+1
  CALL CMISSMeshElements_CreateStart(Mesh1,Mesh1ComponentNumberHydrostaticPressure, &
    & BasisHydrostaticPressure,SolidMeshElementsHydrostaticPressure,Err)
  SELECT CASE(MaterialSpecification)
  CASE(Plate2D)
    DO ElementIndex=1,NumberOfSolidElements
      CALL CMISSMeshElements_NodesSet(SolidMeshElementsHydrostaticPressure,ElementIndex, &
        & SolidElementNodes(ElementIndex,(/1,3,7,9/)),Err)
    ENDDO
  CASE(Plate3D)
    DO ElementIndex=1,NumberOfSolidElements
      CALL CMISSMeshElements_NodesSet(SolidMeshElementsHydrostaticPressure,ElementIndex, &
        & SolidElementNodes(ElementIndex,(/1,3,7,9,19,21,25,27/)),Err)
    ENDDO
  END SELECT
  CALL CMISSMeshElements_CreateFinish(SolidMeshElementsHydrostaticPressure,Err)
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
  CALL CMISSMeshElements_Initialise(FluidMeshElementsSpace,Err)
  CALL CMISSMeshElements_Initialise(FluidMeshElementsVelocity,Err)
  CALL CMISSMeshElements_Initialise(FluidMeshElementsPressure,Err)
  
  Mesh2ComponentNumberSpace=1
  FluidMeshComponentNumberVelocity=1
  Mesh2ComponentNumberPressure=1
  
  CALL CMISSMeshElements_CreateStart(Mesh2,Mesh2ComponentNumberSpace,BasisSpaceFluid,FluidMeshElementsSpace,Err)
  DO ElementIndex=1,NumberOfFluidElements
    CALL CMISSMeshElements_NodesSet(FluidMeshElementsSpace,ElementIndex, &
      & FluidElementNodes(ElementIndex,:),Err)
  ENDDO
  CALL CMISSMeshElements_CreateFinish(FluidMeshElementsSpace,Err)
  
  FluidMeshElementsVelocity=FluidMeshElementsSpace
  Mesh2ComponentNumberPressure=FluidMeshComponentNumberVelocity+1
  CALL CMISSMeshElements_CreateStart(Mesh2,Mesh2ComponentNumberPressure, &
    & BasisPressure,FluidMeshElementsPressure,Err)
  SELECT CASE(MaterialSpecification)
  CASE(Plate2D)
    DO ElementIndex=1,NumberOfFluidElements
      CALL CMISSMeshElements_NodesSet(FluidMeshElementsPressure,ElementIndex, &
        & FluidElementNodes(ElementIndex,(/1,3,7,9/)),Err)
    ENDDO
  CASE(Plate3D)
    DO ElementIndex=1,NumberOfFluidElements
      CALL CMISSMeshElements_NodesSet(FluidMeshElementsPressure,ElementIndex, &
        & FluidElementNodes(ElementIndex,(/1,3,7,9,19,21,25,27/)),Err)
    ENDDO
  END SELECT
  CALL CMISSMeshElements_CreateFinish(FluidMeshElementsPressure,Err)
  CALL CMISSMesh_CreateFinish(Mesh2,Err)
  !
  !================================================================================================================================
  !
  ! I N T E R F A C E
  !Create an interface between the two meshes
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE << == '
  CALL CMISSInterface_Initialise(Interface1,Err)
  CALL CMISSInterface_CreateStart(InterfaceUserNumber,WorldRegion,Interface1,Err)
  CALL CMISSInterface_LabelSet(Interface1,"Interface1",Err)
  !Add in the two meshes
  CALL CMISSInterface_MeshAdd(Interface1,Mesh1,SolidMeshIndex,Err)
  CALL CMISSInterface_MeshAdd(Interface1,Mesh2,FluidMeshIndex,Err)
  CALL CMISSInterface_CoordinateSystemSet(Interface1,InterfaceCoordinateSystem,Err)
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
  CALL CMISSBasis_CreateStart(InterfaceMappingBasisUserNumber,InterfaceMappingBasis1,Err)
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
  DO ElementIndex=1,NumberOfInterfaceElements
    CALL CMISSMeshElements_NodesSet(InterfaceMeshElements,ElementIndex, &
      & InterfaceElementNodes(ElementIndex,:),Err)
  ENDDO
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
    DO ElementIndex=1,NumberOfInterfaceElements
      !Map the interface element to the elements in mesh 1
      CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1,ElementIndex,SolidMeshIndex, &
        & SolidInterfaceElements(ElementIndex),Err)
      DO LocalNodeIndex=1,3**(NumberOfDimensions-1)
        XI2 = [SolidXi2(ElementIndex,LocalNodeIndex),SolidXi3(ElementIndex,LocalNodeIndex)]
        CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1,ElementIndex,SolidMeshIndex, &
          & SolidInterfaceElements(ElementIndex),LocalNodeIndex,1,XI2,Err)
      ENDDO
      !Map the interface element to the elements in mesh 2
      CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1,ElementIndex,FluidMeshIndex, &
        & FluidInterfaceElements(ElementIndex),Err)
      DO LocalNodeIndex=1,3**(NumberOfDimensions-1)
        XI2 = [FluidXi2(ElementIndex,LocalNodeIndex),FluidXi3(ElementIndex,LocalNodeIndex)]
        CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1,ElementIndex,FluidMeshIndex, &
          & FluidInterfaceElements(ElementIndex),LocalNodeIndex,1,XI2,Err)
      ENDDO
    ENDDO !ElementIndex
    CALL CMISSInterfaceMeshConnectivity_NodeNumberSet(InterfaceMeshConnectivity1,InterfaceNodeNumbersForGeometry, &
      & SolidMeshIndex,ConnectedInterfaceNodes(1,:),FluidMeshIndex,ConnectedInterfaceNodes(2,:),Err)
  ELSE
    LocalNodeIndex=0
    DO ElementIndex=1,NumberOfInterfaceElements*9
      LocalNodeIndex=LocalNodeIndex+1
      !Map the interface element to the elements in the solid mesh
      CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1, &
        & InterfaceInterfaceNodeInformationNE(ElementIndex,2),SolidMeshIndex,SolidInterfaceNodeInformationNE(ElementIndex,2),Err)
      XI3 = [SolidInterfaceNodeInformationXi(ElementIndex,1),SolidInterfaceNodeInformationXi(ElementIndex,2), &
        & SolidInterfaceNodeInformationXi(ElementIndex,3)]
   !   SUBROUTINE CMISSInterfaceMeshConnectivity_ElementXiSetObj(interfaceMeshConnectivity,interfaceElementNumber, &
   !     &  coupledMeshIndexNumber,coupledMeshElementNumber,interfaceMeshLocalNodeNumber,interfaceMeshComponentNodeNumber,xi,err)
      CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1, &
        & InterfaceInterfaceNodeInformationNE(ElementIndex,2),SolidMeshIndex, &
        & SolidInterfaceNodeInformationNE(ElementIndex,2),LocalNodeIndex,1,XI3,Err)
      !Map the interface element to the elements in the fluid mesh
      CALL CMISSInterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1, &
        & InterfaceInterfaceNodeInformationNE(ElementIndex,2),FluidMeshIndex, &
        & FluidInterfaceNodeInformationNE(ElementIndex,2),Err)
      XI3 = [FluidInterfaceNodeInformationXi(ElementIndex,1),FluidInterfaceNodeInformationXi(ElementIndex,2), &
        & FluidInterfaceNodeInformationXi(ElementIndex,3)]
      CALL CMISSInterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1, &
        & InterfaceInterfaceNodeInformationNE(ElementIndex,2),FluidMeshIndex, &
        & FluidInterfaceNodeInformationNE(ElementIndex,2),LocalNodeIndex,1,XI3,Err)
      IF(LocalNodeIndex==9) LocalNodeIndex=0
    ENDDO !ElementIndex
    CALL CMISSInterfaceMeshConnectivity_NodeNumberSet(InterfaceMeshConnectivity1,InterfaceNodeNumbersForGeometry, &
      & SolidMeshIndex,ConnectedInterfaceNodes(1,:),FluidMeshIndex,ConnectedInterfaceNodes(2,:),Err)
  ENDIF
  CALL CMISSInterfaceMeshConnectivity_CreateFinish(InterfaceMeshConnectivity1,Err)
  !
  !================================================================================================================================
  !
  ! D E C O M P O S I T I O N
  !Create a decomposition for the solid mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID MESH DECOMPOSITION << == '
  CALL CMISSDecomposition_Initialise(SolidDecomposition,Err)
  CALL CMISSDecomposition_CreateStart(SolidDecompositionUserNumber,Mesh1,SolidDecomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(SolidDecomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(SolidDecomposition,NumberOfComputationalNodes,Err)
  CALL CMISSDecomposition_CalculateFacesSet(SolidDecomposition,.TRUE.,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(SolidDecomposition,Err)
  !Create a decomposition for the fluid mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID MESH DECOMPOSITION << == '
  CALL CMISSDecomposition_Initialise(FluidDecomposition,Err)
  CALL CMISSDecomposition_CreateStart(FluidDecompositionUserNumber,Mesh2,FluidDecomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(FluidDecomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(FluidDecomposition,NumberOfComputationalNodes,Err)
  CALL CMISSDecomposition_CalculateFacesSet(FluidDecomposition,.TRUE.,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(FluidDecomposition,Err)
  !Create a decomposition for the interface mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE DECOMPOSITION << == '
  CALL CMISSDecomposition_Initialise(InterfaceDecomposition,Err)
  CALL CMISSDecomposition_CreateStart(InterfaceDecompositionUserNumber,InterfaceMesh1,InterfaceDecomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecomposition_TypeSet(InterfaceDecomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(InterfaceDecomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMISSDecomposition_CreateFinish(InterfaceDecomposition,Err)
  !
  !================================================================================================================================
  !
  ! G E O M E T R I C _ F I E L D S
  !Start to create a default (geometric) field on the solid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID MESH GEOMETRIC FIELD << == '
  CALL CMISSField_Initialise(GeometricField1,Err)
  CALL CMISSField_CreateStart(SolidGeometricFieldUserNumber,Region1,GeometricField1,Err)
  !Set the decomposition to use
  CALL CMISSField_MeshDecompositionSet(GeometricField1,SolidDecomposition,Err)
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
  CALL CMISSField_MeshDecompositionSet(GeometricField2,FluidDecomposition,Err)
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
  CALL CMISSField_MeshDecompositionSet(InterfaceGeometricField1,InterfaceDecomposition,Err)
  CALL CMISSField_VariableLabelSet(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,"InterfaceGF",Err)
  !Set the domain to be used by the field components
  CALL CMISSField_ComponentMeshComponentSet(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMISSField_ComponentMeshComponentSet(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NumberOfDimensions==3) CALL CMISSField_ComponentMeshComponentSet(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the first field
  CALL CMISSField_CreateFinish(InterfaceGeometricField1,Err)
  !Update the geometric field parameters (solid)
  DO NodeIndex=1,NumberOfSolidNodes
    CALL CMISSDecomposition_NodeDomainGet(SolidDecomposition,SolidNodeNumbers(NodeIndex),1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      IF(NumberOfDimensions==3) THEN
        Components=(/2,3,1/)
      ELSE
        Components=(/1,2,0/)
      ENDIF
      CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
        & 1,CMISS_NO_GLOBAL_DERIV,SolidNodeNumbers(NodeIndex),Components(1),SolidGeometryY(NodeIndex),Err)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
        & 1,CMISS_NO_GLOBAL_DERIV,SolidNodeNumbers(NodeIndex),Components(2),SolidGeometryZ(NodeIndex),Err)
      IF(NumberOfDimensions==3) THEN
        CALL CMISSField_ParameterSetUpdateNode(GeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,SolidNodeNumbers(NodeIndex),Components(3),SolidGeometryX(NodeIndex),Err)
      ENDIF
    ENDIF
  ENDDO
  !Update the geometric field parameters (fluid)
  DO NodeIndex=1,NumberOfFluidNodes
    CALL CMISSDecomposition_NodeDomainGet(FluidDecomposition,FluidNodeNumbers(NodeIndex),1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
        & 1,CMISS_NO_GLOBAL_DERIV,FluidNodeNumbers(NodeIndex),Components(1),FluidGeometryY(NodeIndex),Err)
      CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
        & 1,CMISS_NO_GLOBAL_DERIV,FluidNodeNumbers(NodeIndex),Components(2),FluidGeometryZ(NodeIndex),Err)
      IF(NumberOfDimensions==3) THEN
        CALL CMISSField_ParameterSetUpdateNode(GeometricField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,FluidNodeNumbers(NodeIndex),Components(3),FluidGeometryX(NodeIndex),Err)
      ENDIF
    ENDIF
  ENDDO
  !Update the geometric field parameters (interface)
  DO NodeIndex=1,NumberOfInterfaceNodes
    CALL CMISSDecomposition_NodeDomainGet(InterfaceDecomposition,InterfaceNodeNumbersForGeometry(NodeIndex),1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSField_ParameterSetUpdateNode(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
        & 1,CMISS_NO_GLOBAL_DERIV,InterfaceNodeNumbersForGeometry(NodeIndex),Components(1),InterfaceGeometryY(NodeIndex),Err)
      CALL CMISSField_ParameterSetUpdateNode(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
        & 1,CMISS_NO_GLOBAL_DERIV,InterfaceNodeNumbersForGeometry(NodeIndex),Components(2),InterfaceGeometryZ(NodeIndex),Err)
      IF(NumberOfDimensions==3) THEN
        CALL CMISSField_ParameterSetUpdateNode(InterfaceGeometricField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
          & 1,CMISS_NO_GLOBAL_DERIV,InterfaceNodeNumbersForGeometry(NodeIndex),Components(3),InterfaceGeometryX(NodeIndex),Err)
      ENDIF
    ENDIF
  ENDDO
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
!  CALL CMISSField_MeshDecompositionSet(FibreField,SolidDecomposition,Err)
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
  CALL CMISSEquationsSet_Initialise(SolidEquationsSet,Err)
  CALL CMISSEquationsSet_CreateStart(SolidEquationsSetUserNumber,Region1,GeometricField1,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,SolidEquationsSetFieldUserNumber, &
    & EquationsSetField1,SolidEquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(SolidEquationsSet,Err)
  !Create the equations set for the fluid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID EQUATION SET << == '
  !Create the equations set for ALE Navier-Stokes
  CALL CMISSField_Initialise(EquationsSetField2,Err)
  CALL CMISSEquationsSet_Initialise(FluidEquationsSet,Err)
  CALL CMISSEquationsSet_CreateStart(FluidEquationsSetUserNumber,Region2,GeometricField2,CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS,&
    & CMISS_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE,CMISS_EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE,&
    & FluidEquationsSetFieldUserNumber,EquationsSetField2,FluidEquationsSet,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSet_CreateFinish(FluidEquationsSet,Err)
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
  CALL CMISSEquationsSet_DependentCreateStart(SolidEquationsSet,SolidDependentFieldUserNumber,DependentField1,Err)
  CALL CMISSField_VariableLabelSet(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,"SolidDF",Err)
  CALL CMISSField_VariableLabelSet(DependentField1,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,"SolidDerivatives",Err)
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
  CALL CMISSEquationsSet_DependentCreateFinish(SolidEquationsSet,Err)
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
    & NumberOfDimensions+1,-MooneyRivlin1,Err)!14.0_CMISSDP,Err)
  
  CALL CMISSField_ParameterSetUpdateStart(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(DependentField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)

  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> SOLID MATERIAL FIELD << == '
  !Create the material field
  CALL CMISSField_Initialise(MaterialField1,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(SolidEquationsSet,SolidMaterialFieldUserNumber,MaterialField1,Err)
  CALL CMISSField_VariableLabelSet(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,"Material1",Err)
  CALL CMISSField_VariableLabelSet(MaterialField1,CMISS_FIELD_V_VARIABLE_TYPE,"SolidDensity",Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(SolidEquationsSet,Err)
  !Set Mooney-Rivlin constants c10 and c01 (default?2.0 and 6.0) respectively
  CALL CMISSField_ComponentValuesInitialise(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & MooneyRivlin1,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2, &
    & MooneyRivlin2,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField1,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & SolidDensity,Err)
  IF(GravityFlag) THEN
    IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOURCE FIELD - GRAVITY << == '
    !Create the source field with the gravity vector
    CALL CMISSField_Initialise(SourceField1,Err)
    CALL CMISSEquationsSet_SourceCreateStart(SolidEquationsSet,SourceFieldUserNumber,SourceField1,Err)
    CALL CMISSField_ScalingTypeSet(SourceField1,CMISS_FIELD_NO_SCALING,Err)
    CALL CMISSEquationsSet_SourceCreateFinish(SolidEquationsSet,Err)
    DO component_idx=1,NumberOfDimensions
      CALL CMISSField_ComponentValuesInitialise(SourceField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
          & component_idx,Gravity(component_idx),Err)
    ENDDO
  ENDIF
  !=========================
  !Create the equations set dependent field variables for dynamic Navier-Stokes
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID DEPENDENT FIELD << == '
  CALL CMISSField_Initialise(DependentField2,Err)
  CALL CMISSEquationsSet_DependentCreateStart(FluidEquationsSet,FluidDependentFieldUserNumber,DependentField2,Err)
  CALL CMISSField_VariableLabelSet(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,"FluidDF",Err)
  CALL CMISSField_VariableLabelSet(DependentField2,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,"FluidDerivatives",Err)
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
  CALL CMISSEquationsSet_DependentCreateFinish(FluidEquationsSet,Err)
  !Initialise dependent field
  DO ComponentNumber=1,NumberOfDimensions
    CALL CMISSField_ComponentValuesInitialise(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
      & ComponentNumber,InitialFieldNavierStokes(ComponentNumber),Err)
  ENDDO
  !Initialise pressure component
  CALL CMISSField_ComponentValuesInitialise(DependentField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & NumberOfDimensions+1,0.1_CMISSDP,Err)
  !=========================
  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> MOVING MESH DEPENDENT FIELD << == '
  !Create the equations set dependent field variables for moving mesh
  CALL CMISSField_Initialise(DependentFieldMovingMesh,Err)
  CALL CMISSEquationsSet_DependentCreateStart(MovingMeshEquationsSet,DependentFieldMovingMeshUserNumber, & 
    & DependentFieldMovingMesh,Err)
  CALL CMISSField_VariableLabelSet(DependentFieldMovingMesh,CMISS_FIELD_U_VARIABLE_TYPE,"MovingMeshDF",Err)
  CALL CMISSField_VariableLabelSet(DependentFieldMovingMesh,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,"MovingMeshDerivatives",Err)
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
  CALL CMISSEquationsSet_MaterialsCreateStart(FluidEquationsSet,FluidMaterialFieldUserNumber,MaterialField2,Err)
  !Finish the equations set materials field variables
  CALL CMISSEquationsSet_MaterialsCreateFinish(FluidEquationsSet,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & FluidMaterialFieldComponentMu,FluidDynamicViscosity,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField2,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, & 
    & FluidMaterialFieldComponentRho,FluidDensity,Err)
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
  CALL CMISSEquationsSet_IndependentCreateStart(FluidEquationsSet,IndependentField2UserNumber,IndependentField2,Err)
  CALL CMISSField_VariableLabelSet(IndependentField2,CMISS_FIELD_U_VARIABLE_TYPE,"FluidInDF",Err)
  !Set the mesh component to be used by the field components.
  DO ComponentNumber=1,NumberOfDimensions
    CALL CMISSField_ComponentMeshComponentSet(IndependentField2,CMISS_FIELD_U_VARIABLE_TYPE,ComponentNumber, & 
      & Mesh2ComponentNumberSpace,Err)
  ENDDO
  !Finish the equations set independent field variables
  CALL CMISSEquationsSet_IndependentCreateFinish(FluidEquationsSet,Err)
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INDEPENDENT FIELD MOVING MESH << == '
  !Create the equations set independent field variables for moving mesh
  CALL CMISSField_Initialise(IndependentFieldMovingMesh,Err)
  CALL CMISSEquationsSet_IndependentCreateStart(MovingMeshEquationsSet,IndependentFieldMovingMeshUserNumber, & 
    & IndependentFieldMovingMesh,Err)
  !Set the scaling to use
  CALL CMISSField_ScalingTypeSet(IndependentFieldMovingMesh,CMISS_FIELD_NO_SCALING,Err)
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
  CALL CMISSEquationsSet_EquationsCreateStart(SolidEquationsSet,Equations1,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations1,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL CMISSEquations_OutputTypeSet(Equations1,CMISS_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(SolidEquationsSet,Err)
  !Create the equations set equations for the second equations set
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID EQUATIONS << == '
  CALL CMISSEquations_Initialise(Equations2,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(FluidEquationsSet,Equations2,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquations_SparsityTypeSet(Equations2,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output (NO_OUTPUT/TIMING_OUTPUT/MATRIX_OUTPUT/ELEMENT_MATRIX_OUTPUT)
  CALL CMISSEquations_OutputTypeSet(Equations2,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSet_EquationsCreateFinish(FluidEquationsSet,Err)
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
  CALL CMISSInterfaceCondition_Initialise(InterfaceCondition,Err)
  CALL CMISSInterfaceCondition_CreateStart(InterfaceConditionUserNumber,Interface1,InterfaceGeometricField1, &
    & InterfaceCondition,Err)
  !Specify the method for the interface condition
  CALL CMISSInterfaceCondition_MethodSet(InterfaceCondition,CMISS_INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,Err)
  !Specify the type of interface condition operator
  CALL CMISSInterfaceCondition_OperatorSet(InterfaceCondition,CMISS_INTERFACE_CONDITION_SOLID_FLUID_OPERATOR,Err)
  !Add in the dependent variables from the equations sets
  CALL CMISSInterfaceCondition_DependentVariableAdd(InterfaceCondition,SolidMeshIndex,SolidEquationsSet, &
    & CMISS_FIELD_U_VARIABLE_TYPE,Err)
  CALL CMISSInterfaceCondition_DependentVariableAdd(InterfaceCondition,FluidMeshIndex,FluidEquationsSet, &
    & CMISS_FIELD_U_VARIABLE_TYPE,Err)
  !Finish creating the interface condition
  CALL CMISSInterfaceCondition_CreateFinish(InterfaceCondition,Err)
  !Create the Lagrange multipliers field
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE LAGRANGE FIELD << == '
  CALL CMISSField_Initialise(LagrangeField1,Err)
  CALL CMISSInterfaceCondition_LagrangeFieldCreateStart(InterfaceCondition,LagrangeFieldUserNumber,LagrangeField1,Err)
  CALL CMISSField_VariableLabelSet(LagrangeField1,CMISS_FIELD_U_VARIABLE_TYPE,"InterfaceLF",Err)
  !Finish the Lagrange multipliers field
  CALL CMISSInterfaceCondition_LagrangeFieldCreateFinish(InterfaceCondition,Err)
  DO ComponentNumber=1,NumberOfDimensions
    CALL CMISSField_ComponentValuesInitialise(LagrangeField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
      & ComponentNumber,0.0_CMISSDP,Err)
  ENDDO
  CALL CMISSField_ParameterSetUpdateStart(LagrangeField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  CALL CMISSField_ParameterSetUpdateFinish(LagrangeField1,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,Err)
  !Create the interface condition equations
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE EQUATIONS << == '
  CALL CMISSInterfaceEquations_Initialise(InterfaceEquations,Err)
  CALL CMISSInterfaceCondition_EquationsCreateStart(InterfaceCondition,InterfaceEquations,Err)
  !Set the interface equations sparsity
  CALL CMISSInterfaceEquations_SparsitySet(InterfaceEquations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the interface equations output
  CALL CMISSInterfaceEquations_OutputTypeSet(InterfaceEquations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  !Finish creating the interface equations
  CALL CMISSInterfaceCondition_EquationsCreateFinish(InterfaceCondition,Err)
  !
  !================================================================================================================================
  !
  ! P R O B L E M
  !Start the creation of a coupled problem
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
    & LinearSolverMovingMeshIndex,LinearSolverMovingMesh,Err)
  CALL CMISSSolver_OutputTypeSet(LinearSolverMovingMesh,LinearSolverMovingMesh_OutputType,Err)
  !Solvers for coupled FiniteElasticity NavierStokes problem
  CALL CMISSSolver_Initialise(DynamicSolver,Err)
  CALL CMISSSolver_Initialise(NonlinearSolver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  !Get the dynamic ALE solver
  CALL CMISSProblem_SolverGet(CoupledProblem,CMISS_CONTROL_LOOP_NODE,DynamicSolverIndex,DynamicSolver,Err)
  CALL CMISSSolver_OutputTypeSet(DynamicSolver,DynamicSolver_OutputType,Err)
  CALL CMISSSolver_DynamicThetaSet(DynamicSolver,DynamicSolver_Theta,Err)
  !Get the dynamic nonlinear solver
  CALL CMISSSolver_DynamicNonlinearSolverGet(DynamicSolver,NonlinearSolver,Err)
  CALL CMISSSolver_NewtonLineSearchTypeSet(NonlinearSolver,CMISS_SOLVER_NEWTON_LINESEARCH_LINEAR,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(NonlinearSolver,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL CMISSSolver_NewtonMaximumFunctionEvaluationsSet(NonlinearSolver,MaxFunctionEvaluations,Err)
  CALL CMISSSolver_OutputTypeSet(NonlinearSolver,NonlinearSolver_OutputType,Err)
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(NonlinearSolver,AbsoluteTolerance,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(NonlinearSolver,RelativeTolerance,Err)
  CALL CMISSSolver_NewtonMaximumIterationsSet(NonlinearSolver,MaximumIterations,Err)
  CALL CMISSSolver_NewtonLineSearchAlphaSet(NonlinearSolver,LinesearchAlpha,Err)
  !Get the dynamic nonlinear linear solver
  CALL CMISSSolver_NewtonLinearSolverGet(NonlinearSolver,LinearSolver,Err)
  !Choose type of linear solver
  IF(.FALSE.) THEN
    CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    !NO,JACOBI,BLOCK_JACOBI,SOR,INCOMPLETE_CHOLESKY,INCOMPLETE_LU,ADDITIVE_SCHWARZ
    CALL CMISSSolver_LinearIterativePreconditionerTypeSet(LinearSolver,CMISS_SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER,Err)
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
  CALL CMISSProblem_SolverGet(CoupledProblem,CMISS_CONTROL_LOOP_NODE,LinearSolverMovingMeshIndex,LinearSolverMovingMesh,Err)
  CALL CMISSSolver_SolverEquationsGet(LinearSolverMovingMesh,LinearSolverMovingMeshEquations,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(LinearSolverMovingMeshEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(LinearSolverMovingMeshEquations,MovingMeshEquationsSet,MovingMeshEquationsSetIndex,Err)
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLVER EQUATIONS << == '
  CALL CMISSSolver_Initialise(DynamicSolver,Err)
  CALL CMISSSolverEquations_Initialise(CoupledSolverEquations,Err)
  !Get the dynamic solver equations
  CALL CMISSProblem_SolverGet(CoupledProblem,CMISS_CONTROL_LOOP_NODE,DynamicSolverIndex,DynamicSolver,Err)
  CALL CMISSSolver_SolverEquationsGet(DynamicSolver,CoupledSolverEquations,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(CoupledSolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(CoupledSolverEquations,SolidEquationsSet,SolidEquationsSetIndex,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(CoupledSolverEquations,FluidEquationsSet,FluidEquationsSetIndex,Err)
  CALL CMISSSolverEquations_InterfaceConditionAdd(CoupledSolverEquations,InterfaceCondition,InterfaceConditionIndex,Err)
  !Set the time dependence of the interface matrix to determine the interface matrix coefficient in the solver matrix
  ! (basically position in big coupled matrix system)
  CALL CMISSInterfaceMatrices_TimeDependenceTypeSet(InterfaceCondition,SolidEquationsSetIndex,.TRUE., &
    & (/CMISS_INTERFACE_MATRIX_STATIC,CMISS_INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC/),Err)
  CALL CMISSInterfaceMatrices_TimeDependenceTypeSet(InterfaceCondition,FluidEquationsSetIndex,.TRUE., &
    & (/CMISS_INTERFACE_MATRIX_STATIC,CMISS_INTERFACE_MATRIX_STATIC/),Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblem_SolverEquationsCreateFinish(CoupledProblem,Err)
  !
  !================================================================================================================================
  !
  ! B O U N D A R Y _ C O N D I T I O N S
  !Start the creation of the equations set boundary conditions
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> BOUNDARY CONDITIONS << == '
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
    CALL CMISSDecomposition_NodeDomainGet(SolidDecomposition,NodeNumber,1,NodeDomain,Err)
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
    CALL CMISSDecomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
        & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,NumberOfDimensions+1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO
  !Inlet velocity nodes, small starting velocity in 1st coordinate direction
  DO S=1,SIZE(InletNodes)
    NodeNumber=InletNodes(S)
    CALL CMISSDecomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
        & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED_INLET,0.0_CMISSDP,Err)
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
  !Set no-slip BC
  DO S=1,SIZE(NoSlipNodes)
    NodeNumber=NoSlipNodes(S)
    CALL CMISSDecomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
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
  DO S=1,SIZE(SlipNodesTop)
    NodeNumber=SlipNodesTop(S)
    CALL CMISSDecomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
        & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,NumberOfDimensions,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO
  !Set slip BC
  IF(NumberOfDimensions==3) THEN
    DO S=1,SIZE(SlipNodesRightLeft)
      NodeNumber=SlipNodesRightLeft(S)
      CALL CMISSDecomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
          & CMISS_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      ENDIF
    ENDDO
  ENDIF
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
    IF(CheckWithoutInterfaceCondition) THEN
      DO S=1,SIZE(InterfaceNodeNumbersForGeometry)
        NodeNumber=InterfaceNodeNumbersForGeometry(S)
        CALL CMISSDecomposition_NodeDomainGet(InterfaceDecomposition,NodeNumber,1,NodeDomain,Err)
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
    ELSE
      DO S=1,SIZE(LagrangeNodes)
        NodeNumber=LagrangeNodes(S)
        CALL CMISSDecomposition_NodeDomainGet(InterfaceDecomposition,NodeNumber,1,NodeDomain,Err)
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
  ENDIF
  !Finish equations set boundary conditions
  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(CoupledSolverEquations,Err)
  
  !Start the creation of the moving mesh boundary conditions
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditionsMovingMesh,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(LinearSolverMovingMeshEquations,BoundaryConditionsMovingMesh,Err)
  !Fixed boundary nodes. May be used to move nodes..
  DO S=1,SIZE(MovedYNodes)
    NodeNumber=MovedYNodes(S)
    CALL CMISSDecomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
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
  !Mesh nodes that are moving wall nodes
  DO S=1,SIZE(MovedNodes)
    NodeNumber=MovedNodes(S)
    CALL CMISSDecomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
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
    CALL CMISSDecomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
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
  
  !
  !================================================================================================================================
  !
  ! S O L V E
  
  PRINT *, ' == >> SOLVING PROBLEM << == '
  CALL CMISSProblem_Solve(CoupledProblem,Err)
  
  !Finialise CMISS
  CALL CMISSFinalise(Err)

  e=etime(t)/60.0_CMISSSP
  PRINT *, "Program successfully completed in ",e," minutes."
  
  STOP
  
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)
    
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING
    
    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP
    
  END SUBROUTINE HANDLE_ERROR
  
END PROGRAM FortranExample
