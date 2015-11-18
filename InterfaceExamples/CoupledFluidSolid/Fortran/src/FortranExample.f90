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

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWINCMISS
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters==========================================================================================================
  !Material specification
  INTEGER(CMISSIntg), PARAMETER :: Plate2D=2
  INTEGER(CMISSIntg), PARAMETER :: Plate3D=3
  
  REAL(CMISSRP) :: SolidDensity
  REAL(CMISSRP), ALLOCATABLE :: Gravity(:)
  
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

  REAL(CMISSRP), ALLOCATABLE :: SolidGeometryX(:),SolidGeometryY(:),SolidGeometryZ(:), &
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
  REAL(CMISSRP) :: XI2(2),XI3(3)
  REAL(CMISSRP) :: MovingMeshParameterK
  REAL(CMISSRP) :: FluidDynamicViscosity
  REAL(CMISSRP) :: FluidDensity
  REAL(CMISSRP) :: YoungsModulus
  REAL(CMISSRP) :: PoissonsRatio
  REAL(CMISSRP) :: ShearModulus
  REAL(CMISSRP) :: BulkModulus
  REAL(CMISSRP) :: MooneyRivlin1
  REAL(CMISSRP) :: MooneyRivlin2

  REAL(CMISSRP) :: InitialFieldNavierStokes(3)
  REAL(CMISSRP) :: InitialFieldMovingMesh(3)
  REAL(CMISSRP) :: DivergenceTolerance
  REAL(CMISSRP) :: RelativeTolerance
  REAL(CMISSRP) :: AbsoluteTolerance
  REAL(CMISSRP) :: LinesearchAlpha

  REAL(CMISSRP) :: StartTime
  REAL(CMISSRP) :: StopTime
  REAL(CMISSRP) :: DynamicSolver_Theta
  REAL(CMISSRP) :: TimeStepSize

  LOGICAL :: FileReadDiagnostics=.FALSE.
  LOGICAL :: ExampleFileProgressDiagnostics=.FALSE.
  LOGICAL :: GeometryCheck=.FALSE.
  LOGICAL :: GravityFlag=.FALSE.
  LOGICAL :: CheckWithoutInterfaceCondition=.FALSE.
  LOGICAL :: SetupOutput=.FALSE.
  
  !CMISS variables
  TYPE(cmfe_BasisType) :: BasisSpaceSolid,BasisDisplacement,BasisHydrostaticPressure, &
    & BasisSpaceFluid,BasisVelocity,BasisPressure,InterfaceBasis1,InterfaceMappingBasis1
  TYPE(cmfe_NodesType) :: SolidNodes,FluidNodes,InterfaceNodes
  TYPE(cmfe_MeshElementsType) :: SolidMeshElementsSpace,SolidMeshElementsDisplacement,SolidMeshElementsHydrostaticPressure, &
    & FluidMeshElementsSpace,FluidMeshElementsVelocity,FluidMeshElementsPressure,InterfaceMeshElements
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions,BoundaryConditionsMovingMesh
  TYPE(cmfe_CoordinateSystemType) :: SolidCoordinateSystem,FluidCoordinateSystem,InterfaceCoordinateSystem, &
    & WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: SolidDecomposition,FluidDecomposition,InterfaceDecomposition
  TYPE(cmfe_EquationsType) :: Equations1,Equations2,EquationsMovingMesh
  TYPE(cmfe_EquationsSetType) :: SolidEquationsSet,FluidEquationsSet,MovingMeshEquationsSet
  TYPE(cmfe_FieldType) :: GeometricField1,GeometricField2,InterfaceGeometricField1, &
    & DependentField1,DependentField2,LagrangeField1,EquationsSetField1,EquationsSetField2, &
    & DependentFieldMovingMesh,MaterialFieldMovingMesh,IndependentField2,IndependentFieldMovingMesh, &
    & MaterialField1,MaterialField2,EquationsSetFieldMovingMesh,SourceField1
  TYPE(cmfe_FieldsType) :: Fields1,Fields2,FieldsI
  TYPE(cmfe_InterfaceType) :: Interface1
  TYPE(cmfe_InterfaceConditionType) :: InterfaceCondition
  TYPE(cmfe_InterfaceEquationsType) :: InterfaceEquations
  TYPE(cmfe_InterfaceMeshConnectivityType) :: InterfaceMeshConnectivity1
  TYPE(cmfe_MeshType) :: Mesh1,Mesh2,InterfaceMesh1
  TYPE(cmfe_ProblemType) :: CoupledProblem
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_RegionType) :: Region1,Region2,WorldRegion
  TYPE(cmfe_SolverType) :: DynamicSolver,NonlinearSolver,LinearSolver,LinearSolverMovingMesh
  TYPE(cmfe_SolverEquationsType) :: CoupledSolverEquations,LinearSolverMovingMeshEquations
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
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  !Set error handling mode
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  !CALL cmfe_OutputSetOn("Testing",Err)
  !Set diganostics for testing
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["SOLVER_MAPPING_CALCULATE", &
  !  & "SOLVER_MATRIX_STRUCTURE_CALCULATE"],Err)
  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)
  !
  !=================================================================================================================================
  !
  ! I N I T I A L _ V A L U E S _ & _ P R O B L E M _ C O N S T A N T S
  !Set initial values
  InitialFieldNavierStokes(1)=0.0_CMISSRP
  InitialFieldNavierStokes(2)=0.0_CMISSRP
  InitialFieldNavierStokes(3)=0.0_CMISSRP
  InitialFieldMovingMesh(1)=0.0_CMISSRP
  InitialFieldMovingMesh(2)=0.0_CMISSRP
  InitialFieldMovingMesh(3)=0.0_CMISSRP
  !Set output parameters
  SetupOutput=.TRUE.
  !(NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
  LinearSolverMovingMesh_OutputType=CMFE_SOLVER_NO_OUTPUT
  DynamicSolver_OutputType=CMFE_SOLVER_NO_OUTPUT
  LinearSolver_OutputType=CMFE_SOLVER_NO_OUTPUT
  NonlinearSolver_OutputType=CMFE_SOLVER_NO_OUTPUT
  !(NoOutput/TimingOutput/MatrixOutput/ElementOutput)
  EquationsNavierStokesOutput=CMFE_EQUATIONS_NO_OUTPUT
  !Set result output parameter
  OutputFrequency=1
  
  !Choose 2D or 3D case
  MaterialSpecification=Plate3D
  !Set solver parameters
  RelativeTolerance=1.0E-4_CMISSRP !default: 1.0E-05_CMISSRP
  AbsoluteTolerance=1.0E-4_CMISSRP !default: 1.0E-10_CMISSRP
  DivergenceTolerance=1.0E5 !default: 1.0E5
  MaximumIterations=100000000 !default: 100000
  MaxFunctionEvaluations=100000
  RestartValue=30 !default: 30
  LinesearchAlpha=1.0_CMISSRP
  MovingMeshParameterK=1.0 !default
  DynamicSolver_Theta=1.0_CMISSRP
  StartTime=0.0_CMISSRP
  StopTime=2000.0_CMISSRP
  TimeStepSize=1.0_CMISSRP
  IF(MaterialSpecification==Plate2D) THEN
    StopTime=0.1_CMISSRP
    TimeStepSize=0.05_CMISSRP
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
  ShearModulus=YoungsModulus/(2.0_CMISSRP*(1.0_CMISSRP+PoissonsRatio)) ! N/m2
  BulkModulus=YoungsModulus/(3.0_CMISSRP*(1.0_CMISSRP-2.0_CMISSRP*PoissonsRatio))
  !MooneyRivlin1=0.5_CMISSRP*ShearModulus ! N/m2
  MooneyRivlin1=0.0595_CMISSRP
  MooneyRivlin2=0.0_CMISSRP
   
  !Set geometric dimension n gravity
  SELECT CASE(MaterialSpecification)
  CASE(Plate2D)
    NumberOfDimensions=2
    ALLOCATE(Gravity(NumberOfDimensions))
    Gravity(:)=[0.0_CMISSRP,9.81_CMISSRP] ! m/s2
  CASE(Plate3D)
    NumberOfDimensions=3
    CheckWithoutInterfaceCondition=.FALSE.!if set to true we remove all Lagrange field dofs by setting them as zero dirichlet BC
    IF(CheckWithoutInterfaceCondition) THEN
      GravityFlag=.TRUE.
      ALLOCATE(Gravity(NumberOfDimensions))
      Gravity(:)=[0.0_CMISSRP,9.81_CMISSRP,0.0_CMISSRP] ! m/s2
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
  CALL cmfe_CoordinateSystem_Initialise(SolidCoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(SolidCoordinateSystemUserNumber,SolidCoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(SolidCoordinateSystem,NumberOfDimensions,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(SolidCoordinateSystem,Err)
  !Create a new RC coordinate system for the fluid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID COORDINATE SYSTEM << == '
  CALL cmfe_CoordinateSystem_Initialise(FluidCoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(FluidCoordinateSystemUserNumber,FluidCoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(FluidCoordinateSystem,NumberOfDimensions,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(FluidCoordinateSystem,Err)
  !Create a new RC coordinate system for the interface region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE COORDINATE SYSTEM << == '
  CALL cmfe_CoordinateSystem_Initialise(InterfaceCoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(InterfaceCoordinateSystemUserNumber,InterfaceCoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(InterfaceCoordinateSystem,NumberOfDimensions,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(InterfaceCoordinateSystem,Err)
  !
  !=================================================================================================================================
  !
  ! R E G I O N S
  !Create the solid region and set the regions coordinate system to the RC coordinate system that we have created
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID REGION << == '
  CALL cmfe_Region_Initialise(Region1,Err)
  CALL cmfe_Region_CreateStart(SolidRegionUserNumber,WorldRegion,Region1,Err)
  CALL cmfe_Region_LabelSet(Region1,"SolidRegion",Err)
  CALL cmfe_Region_CoordinateSystemSet(Region1,SolidCoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region1,Err)
  !Create the fluid region and set the regions coordinate system to the RC coordinate system that we have created
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID REGION << == '
  CALL cmfe_Region_Initialise(Region2,Err)
  CALL cmfe_Region_CreateStart(FluidRegionUserNumber,WorldRegion,Region2,Err)
  CALL cmfe_Region_LabelSet(Region2,"FluidRegion",Err)
  CALL cmfe_Region_CoordinateSystemSet(Region2,FluidCoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region2,Err)
  !
  !=================================================================================================================================
  !
  ! B A S I S _ S O L I D
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> BASIS FOR SOLID: DISPLACEMENT, HYDROSTATIC PRESSURE << == '
  !Create a basis for the dependent field variable displacement
  CALL cmfe_Basis_Initialise(BasisDisplacement,Err)
  CALL cmfe_Basis_CreateStart(BasisDisplacementUserNumber,BasisDisplacement,Err)
  SELECT CASE(InterpolationTypeDisplacement)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(BasisDisplacement,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(BasisDisplacement,CMFE_BASIS_SIMPLEX_TYPE,Err)
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  SELECT CASE(InterpolationTypeDisplacement)
  CASE(CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION)
    NumberOfGaussXi=2
    PressureMeshComponent=1
  CASE(CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
    NumberOfGaussXi=3
    PressureMeshComponent=2
    InterpolationTypeHydrostaticPressure=1
  CASE(CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION)
    NumberOfGaussXi=4
    PressureMeshComponent=2
    InterpolationTypeHydrostaticPressure=2
  CASE DEFAULT
    NumberOfGaussXi=0
    PressureMeshComponent=1
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  CALL cmfe_Basis_NumberOfXiSet(BasisDisplacement,NumberOfDimensions,Err)
  IF(NumberOfDimensions==2) THEN
    CALL cmfe_Basis_InterpolationXiSet(BasisDisplacement,[InterpolationTypeDisplacement,InterpolationTypeDisplacement],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisDisplacement,[NumberOfGaussXi,NumberOfGaussXi],Err)
  ELSE
    CALL cmfe_Basis_InterpolationXiSet(BasisDisplacement,[InterpolationTypeDisplacement,InterpolationTypeDisplacement, &
      & InterpolationTypeDisplacement],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisDisplacement,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  ENDIF
  CALL cmfe_Basis_CreateFinish(BasisDisplacement,Err)
  
  !Use the displacement basis as a space basis
  BasisSpaceSolid=BasisDisplacement
  
  !Create a basis for the dependent field variable hydrostatic pressure
  CALL cmfe_Basis_Initialise(BasisHydrostaticPressure,Err)
  CALL cmfe_Basis_CreateStart(BasisHydrostaticPressureUserNumber,BasisHydrostaticPressure,Err)
  SELECT CASE(InterpolationTypeHydrostaticPressure)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(BasisHydrostaticPressure,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(BasisHydrostaticPressure,CMFE_BASIS_SIMPLEX_TYPE,Err)
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  CALL cmfe_Basis_NumberOfXiSet(BasisHydrostaticPressure,NumberOfDimensions,Err)
  IF(NumberOfDimensions==2) THEN
    CALL cmfe_Basis_InterpolationXiSet(BasisHydrostaticPressure,[InterpolationTypeHydrostaticPressure, &
      & InterpolationTypeHydrostaticPressure],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisHydrostaticPressure,[NumberOfGaussXi,NumberOfGaussXi],Err)
  ELSE
    CALL cmfe_Basis_InterpolationXiSet(BasisHydrostaticPressure,[InterpolationTypeHydrostaticPressure, &
      & InterpolationTypeHydrostaticPressure,InterpolationTypeHydrostaticPressure],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisHydrostaticPressure,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  ENDIF
  CALL cmfe_Basis_CreateFinish(BasisHydrostaticPressure,Err)
  !
  !=================================================================================================================================
  !
  ! B A S I S _ F L U I D
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> CREATING BASIS FOR FLUID: SPACE, VELOCITY, PRESSURE << == '
  !Create a basis for the fluid domain
  MeshNumberOfComponents=1
  CALL cmfe_Basis_Initialise(BasisSpaceFluid,Err)
  CALL cmfe_Basis_CreateStart(BasisSpaceFluidUserNumber,BasisSpaceFluid,Err)
  SELECT CASE(InterpolationTypeSpace)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(BasisSpaceFluid,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(BasisSpaceFluid,CMFE_BASIS_SIMPLEX_TYPE,Err)
  CASE DEFAULT
    CALL HANDLE_ERROR("Invalid interpolation type.")
  END SELECT
  IF(InterpolationTypeSpace==2) THEN
    NumberOfGaussXiSpace=3
  ELSE
    CALL HANDLE_ERROR("Interpolation types other than 2 have not been tested in FiniteElasticity-ALENavierStokes example yet.")
  ENDIF
  CALL cmfe_Basis_NumberOfXiSet(BasisSpaceFluid,NumberOfDimensions,Err)
  IF(NumberOfDimensions==2) THEN
    CALL cmfe_Basis_InterpolationXiSet(BasisSpaceFluid,[InterpolationTypeSpace,InterpolationTypeSpace],Err)
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisSpaceFluid,[NumberOfGaussXiSpace,NumberOfGaussXiSpace],Err)
  ELSE
    CALL cmfe_Basis_InterpolationXiSet(BasisSpaceFluid, &
      & [InterpolationTypeSpace,InterpolationTypeSpace,InterpolationTypeSpace],Err)                         
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisSpaceFluid, &
      & [NumberOfGaussXiSpace,NumberOfGaussXiSpace,NumberOfGaussXiSpace],Err)
  ENDIF
  CALL cmfe_Basis_CreateFinish(BasisSpaceFluid,Err)
  !Create a basis for the dependent field variable velocity
  IF(InterpolationTypeVelocity==InterpolationTypeSpace) THEN
    BasisVelocity=BasisSpaceFluid
  ELSE
    MeshNumberOfComponents=MeshNumberOfComponents+1
    CALL cmfe_Basis_Initialise(BasisVelocity,Err)
    CALL cmfe_Basis_CreateStart(BasisVelocityUserNumber,BasisVelocity,Err)
    CALL cmfe_Basis_TypeSet(BasisVelocity,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    NumberOfGaussXiVelocity=3
    CALL cmfe_Basis_NumberOfXiSet(BasisVelocity,NumberOfDimensions,Err)
    IF(NumberOfDimensions==2) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisVelocity,[InterpolationTypeVelocity,InterpolationTypeVelocity],Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisVelocity,[NumberOfGaussXiVelocity,NumberOfGaussXiVelocity],Err)
    ELSE
      CALL cmfe_Basis_InterpolationXiSet(BasisVelocity, &
        & [InterpolationTypeVelocity,InterpolationTypeVelocity,InterpolationTypeVelocity],Err)                         
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisVelocity, &
        & [NumberOfGaussXiVelocity,NumberOfGaussXiVelocity,NumberOfGaussXiVelocity],Err)
    ENDIF
    CALL cmfe_Basis_CreateFinish(BasisVelocity,Err)
  ENDIF
  !Create a basis for the dependent field variable pressure
  IF(InterpolationTypePressure==InterpolationTypeSpace) THEN
    BasisPressure=BasisSpaceFluid
  ELSE IF(InterpolationTypePressure==InterpolationTypeVelocity) THEN
    BasisPressure=BasisVelocity
  ELSE
    MeshNumberOfComponents=MeshNumberOfComponents+1
    NumberOfGaussXiPressure=3
    CALL cmfe_Basis_Initialise(BasisPressure,Err)
    CALL cmfe_Basis_CreateStart(BasisPressureUserNumber,BasisPressure,Err)
    CALL cmfe_Basis_TypeSet(BasisPressure,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CALL cmfe_Basis_NumberOfXiSet(BasisPressure,NumberOfDimensions,Err)
    IF(NumberOfDimensions==2) THEN
      CALL cmfe_Basis_InterpolationXiSet(BasisPressure,[InterpolationTypePressure,InterpolationTypePressure],Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisPressure,[NumberOfGaussXiPressure,NumberOfGaussXiPressure],Err)
    ELSE
      CALL cmfe_Basis_InterpolationXiSet(BasisPressure, &
        & [InterpolationTypePressure,InterpolationTypePressure,InterpolationTypePressure],Err)                         
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(BasisPressure, &
        & [NumberOfGaussXiPressure,NumberOfGaussXiPressure,NumberOfGaussXiPressure],Err)
    ENDIF
    CALL cmfe_Basis_CreateFinish(BasisPressure,Err)
  ENDIF
  !
  !================================================================================================================================
  !
  ! M E S H E S
  !Create the solid mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID MESH << == '
  CALL cmfe_Nodes_Initialise(SolidNodes,Err)
  CALL cmfe_Mesh_Initialise(Mesh1,Err)
  CALL cmfe_Nodes_CreateStart(Region1,NumberOfSolidNodes,SolidNodes,Err)
  CALL cmfe_Nodes_CreateFinish(SolidNodes,Err)
  CALL cmfe_Mesh_CreateStart(SolidMeshUserNumber,Region1,NumberOfDimensions,Mesh1,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh1,NumberOfSolidElements,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh1,2,Err)
  CALL cmfe_MeshElements_Initialise(SolidMeshElementsSpace,Err)
  CALL cmfe_MeshElements_Initialise(SolidMeshElementsDisplacement,Err)
  CALL cmfe_MeshElements_Initialise(SolidMeshElementsHydrostaticPressure,Err)
  
  Mesh1ComponentNumberSpace=1
  Mesh1ComponentNumberDisplacement=1
  Mesh1ComponentNumberHydrostaticPressure=1
  
  CALL cmfe_MeshElements_CreateStart(Mesh1,Mesh1ComponentNumberSpace,BasisSpaceSolid,SolidMeshElementsSpace,Err)
  DO ElementIndex=1,NumberOfSolidElements
    CALL cmfe_MeshElements_NodesSet(SolidMeshElementsSpace,ElementIndex, &
      & SolidElementNodes(ElementIndex,:),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(SolidMeshElementsSpace,Err)
  
  SolidMeshElementsDisplacement=SolidMeshElementsSpace
  Mesh1ComponentNumberHydrostaticPressure=Mesh1ComponentNumberDisplacement+1
  CALL cmfe_MeshElements_CreateStart(Mesh1,Mesh1ComponentNumberHydrostaticPressure, &
    & BasisHydrostaticPressure,SolidMeshElementsHydrostaticPressure,Err)
  SELECT CASE(MaterialSpecification)
  CASE(Plate2D)
    DO ElementIndex=1,NumberOfSolidElements
      CALL cmfe_MeshElements_NodesSet(SolidMeshElementsHydrostaticPressure,ElementIndex, &
        & SolidElementNodes(ElementIndex,[1,3,7,9]),Err)
    ENDDO
  CASE(Plate3D)
    DO ElementIndex=1,NumberOfSolidElements
      CALL cmfe_MeshElements_NodesSet(SolidMeshElementsHydrostaticPressure,ElementIndex, &
        & SolidElementNodes(ElementIndex,[1,3,7,9,19,21,25,27]),Err)
    ENDDO
  END SELECT
  CALL cmfe_MeshElements_CreateFinish(SolidMeshElementsHydrostaticPressure,Err)
  CALL cmfe_Mesh_CreateFinish(Mesh1,Err)
  !Create fluid mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID MESH << == '
  CALL cmfe_Nodes_Initialise(FluidNodes,Err)
  CALL cmfe_Mesh_Initialise(Mesh2,Err)
  CALL cmfe_Nodes_CreateStart(Region2,NumberOfFluidNodes,FluidNodes,Err)
  CALL cmfe_Nodes_CreateFinish(FluidNodes,Err)
  CALL cmfe_Mesh_CreateStart(FluidMeshUserNumber,Region2,NumberOfDimensions,Mesh2,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh2,NumberOfFluidElements,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh2,2,Err)
  CALL cmfe_MeshElements_Initialise(FluidMeshElementsSpace,Err)
  CALL cmfe_MeshElements_Initialise(FluidMeshElementsVelocity,Err)
  CALL cmfe_MeshElements_Initialise(FluidMeshElementsPressure,Err)
  
  Mesh2ComponentNumberSpace=1
  FluidMeshComponentNumberVelocity=1
  Mesh2ComponentNumberPressure=1
  
  CALL cmfe_MeshElements_CreateStart(Mesh2,Mesh2ComponentNumberSpace,BasisSpaceFluid,FluidMeshElementsSpace,Err)
  DO ElementIndex=1,NumberOfFluidElements
    CALL cmfe_MeshElements_NodesSet(FluidMeshElementsSpace,ElementIndex, &
      & FluidElementNodes(ElementIndex,:),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(FluidMeshElementsSpace,Err)
  
  FluidMeshElementsVelocity=FluidMeshElementsSpace
  Mesh2ComponentNumberPressure=FluidMeshComponentNumberVelocity+1
  CALL cmfe_MeshElements_CreateStart(Mesh2,Mesh2ComponentNumberPressure, &
    & BasisPressure,FluidMeshElementsPressure,Err)
  SELECT CASE(MaterialSpecification)
  CASE(Plate2D)
    DO ElementIndex=1,NumberOfFluidElements
      CALL cmfe_MeshElements_NodesSet(FluidMeshElementsPressure,ElementIndex, &
        & FluidElementNodes(ElementIndex,[1,3,7,9]),Err)
    ENDDO
  CASE(Plate3D)
    DO ElementIndex=1,NumberOfFluidElements
      CALL cmfe_MeshElements_NodesSet(FluidMeshElementsPressure,ElementIndex, &
        & FluidElementNodes(ElementIndex,[1,3,7,9,19,21,25,27]),Err)
    ENDDO
  END SELECT
  CALL cmfe_MeshElements_CreateFinish(FluidMeshElementsPressure,Err)
  CALL cmfe_Mesh_CreateFinish(Mesh2,Err)
  !
  !================================================================================================================================
  !
  ! I N T E R F A C E
  !Create an interface between the two meshes
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE << == '
  CALL cmfe_Interface_Initialise(Interface1,Err)
  CALL cmfe_Interface_CreateStart(InterfaceUserNumber,WorldRegion,Interface1,Err)
  CALL cmfe_Interface_LabelSet(Interface1,"Interface1",Err)
  !Add in the two meshes
  CALL cmfe_Interface_MeshAdd(Interface1,Mesh1,SolidMeshIndex,Err)
  CALL cmfe_Interface_MeshAdd(Interface1,Mesh2,FluidMeshIndex,Err)
  CALL cmfe_Interface_CoordinateSystemSet(Interface1,InterfaceCoordinateSystem,Err)
  CALL cmfe_Interface_CreateFinish(Interface1,Err)
  !Create a (bi)-quadratic-Lagrange basis (3D: faces // 2D: lines)
  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> INTERFACE BASIS << == '
  CALL cmfe_Basis_Initialise(InterfaceBasis1,Err)
  CALL cmfe_Basis_CreateStart(InterfaceBasisUserNumber,InterfaceBasis1,Err)
  CALL cmfe_Basis_NumberOfXiSet(InterfaceBasis1,NumberOfDimensions-1,Err)
  IF(NumberOfDimensions==2) THEN
    CALL cmfe_Basis_InterpolationXiSet(InterfaceBasis1,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  ELSE
    CALL cmfe_Basis_InterpolationXiSet(InterfaceBasis1, &
      & [CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  ENDIF
  CALL cmfe_Basis_CreateFinish(InterfaceBasis1,Err)
  !Create a (bi)-quadratic-Lagrange basis for the interface mapping (3D: faces // 2D: lines)
  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> INTERFACE MAPPING BASIS << == '
  CALL cmfe_Basis_Initialise(InterfaceMappingBasis1,Err)
  CALL cmfe_Basis_CreateStart(InterfaceMappingBasisUserNumber,InterfaceMappingBasis1,Err)
  CALL cmfe_Basis_NumberOfXiSet(InterfaceMappingBasis1,NumberOfDimensions-1,Err)
  IF(NumberOfDimensions==2) THEN
    CALL cmfe_Basis_InterpolationXiSet(InterfaceMappingBasis1,[CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  ELSE
    CALL cmfe_Basis_InterpolationXiSet(InterfaceMappingBasis1, &
      & [CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION],Err)
  ENDIF
  CALL cmfe_Basis_CreateFinish(InterfaceMappingBasis1,Err)
  !
  !================================================================================================================================
  !
  ! I N T E R F A C E _ M E S H
  !Create an interface mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE MESH << == '
  CALL cmfe_Nodes_Initialise(InterfaceNodes,Err)
  CALL cmfe_Mesh_Initialise(InterfaceMesh1,Err)
  CALL cmfe_Nodes_CreateStart(Interface1,NumberOfInterfaceNodes,InterfaceNodes,Err)
  CALL cmfe_Nodes_CreateFinish(InterfaceNodes,Err)
  CALL cmfe_Mesh_CreateStart(InterfaceMeshUserNumber,Interface1,NumberOfDimensions-1,InterfaceMesh1,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(InterfaceMesh1,NumberOfInterfaceElements,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(InterfaceMesh1,1,Err)
  CALL cmfe_MeshElements_Initialise(InterfaceMeshElements,Err)
  InterfaceMeshComponentNumber=1
  CALL cmfe_MeshElements_CreateStart(InterfaceMesh1,InterfaceMeshComponentNumber,InterfaceBasis1,InterfaceMeshElements,Err)
  DO ElementIndex=1,NumberOfInterfaceElements
    CALL cmfe_MeshElements_NodesSet(InterfaceMeshElements,ElementIndex, &
      & InterfaceElementNodes(ElementIndex,:),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(InterfaceMeshElements,Err)
  CALL cmfe_Mesh_CreateFinish(InterfaceMesh1,Err)
  !
  !================================================================================================================================
  !
  ! M E S H _ C O N N E C T I V I T Y
  !Couple the interface meshes
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE MESH CONNECTIVITY << == '
  IF (NumberOfComputationalNodes/=1) CALL HANDLE_ERROR("MeshConnectivity does not work in parallel yet!")
  CALL cmfe_InterfaceMeshConnectivity_Initialise(InterfaceMeshConnectivity1,Err)
  CALL cmfe_InterfaceMeshConnectivity_CreateStart(Interface1,InterfaceMesh1,InterfaceMeshConnectivity1,Err)
  CALL cmfe_InterfaceMeshConnectivity_BasisSet(InterfaceMeshConnectivity1,InterfaceMappingBasis1,Err)
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
      CALL cmfe_InterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1,ElementIndex,SolidMeshIndex, &
        & SolidInterfaceElements(ElementIndex),Err)
      DO LocalNodeIndex=1,3**(NumberOfDimensions-1)
        XI2 = [SolidXi2(ElementIndex,LocalNodeIndex),SolidXi3(ElementIndex,LocalNodeIndex)]
        CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1,ElementIndex,SolidMeshIndex, &
          & SolidInterfaceElements(ElementIndex),LocalNodeIndex,1,XI2,Err)
      ENDDO
      !Map the interface element to the elements in mesh 2
      CALL cmfe_InterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1,ElementIndex,FluidMeshIndex, &
        & FluidInterfaceElements(ElementIndex),Err)
      DO LocalNodeIndex=1,3**(NumberOfDimensions-1)
        XI2 = [FluidXi2(ElementIndex,LocalNodeIndex),FluidXi3(ElementIndex,LocalNodeIndex)]
        CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1,ElementIndex,FluidMeshIndex, &
          & FluidInterfaceElements(ElementIndex),LocalNodeIndex,1,XI2,Err)
      ENDDO
    ENDDO !ElementIndex
    CALL cmfe_InterfaceMeshConnectivity_NodeNumberSet(InterfaceMeshConnectivity1,InterfaceNodeNumbersForGeometry, &
      & SolidMeshIndex,ConnectedInterfaceNodes(1,:),FluidMeshIndex,ConnectedInterfaceNodes(2,:),Err)
  ELSE
    LocalNodeIndex=0
    DO ElementIndex=1,NumberOfInterfaceElements*9
      LocalNodeIndex=LocalNodeIndex+1
      !Map the interface element to the elements in the solid mesh
      CALL cmfe_InterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1, &
        & InterfaceInterfaceNodeInformationNE(ElementIndex,2),SolidMeshIndex,SolidInterfaceNodeInformationNE(ElementIndex,2),Err)
      XI3 = [SolidInterfaceNodeInformationXi(ElementIndex,1),SolidInterfaceNodeInformationXi(ElementIndex,2), &
        & SolidInterfaceNodeInformationXi(ElementIndex,3)]
   !   SUBROUTINE cmfe_InterfaceMeshConnectivity_ElementXiSetObj(interfaceMeshConnectivity,interfaceElementNumber, &
   !     &  coupledMeshIndexNumber,coupledMeshElementNumber,interfaceMeshLocalNodeNumber,interfaceMeshComponentNodeNumber,xi,err)
      CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1, &
        & InterfaceInterfaceNodeInformationNE(ElementIndex,2),SolidMeshIndex, &
        & SolidInterfaceNodeInformationNE(ElementIndex,2),LocalNodeIndex,1,XI3,Err)
      !Map the interface element to the elements in the fluid mesh
      CALL cmfe_InterfaceMeshConnectivity_ElementNumberSet(InterfaceMeshConnectivity1, &
        & InterfaceInterfaceNodeInformationNE(ElementIndex,2),FluidMeshIndex, &
        & FluidInterfaceNodeInformationNE(ElementIndex,2),Err)
      XI3 = [FluidInterfaceNodeInformationXi(ElementIndex,1),FluidInterfaceNodeInformationXi(ElementIndex,2), &
        & FluidInterfaceNodeInformationXi(ElementIndex,3)]
      CALL cmfe_InterfaceMeshConnectivity_ElementXiSet(InterfaceMeshConnectivity1, &
        & InterfaceInterfaceNodeInformationNE(ElementIndex,2),FluidMeshIndex, &
        & FluidInterfaceNodeInformationNE(ElementIndex,2),LocalNodeIndex,1,XI3,Err)
      IF(LocalNodeIndex==9) LocalNodeIndex=0
    ENDDO !ElementIndex
    CALL cmfe_InterfaceMeshConnectivity_NodeNumberSet(InterfaceMeshConnectivity1,InterfaceNodeNumbersForGeometry, &
      & SolidMeshIndex,ConnectedInterfaceNodes(1,:),FluidMeshIndex,ConnectedInterfaceNodes(2,:),Err)
  ENDIF
  CALL cmfe_InterfaceMeshConnectivity_CreateFinish(InterfaceMeshConnectivity1,Err)
  !
  !================================================================================================================================
  !
  ! D E C O M P O S I T I O N
  !Create a decomposition for the solid mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID MESH DECOMPOSITION << == '
  CALL cmfe_Decomposition_Initialise(SolidDecomposition,Err)
  CALL cmfe_Decomposition_CreateStart(SolidDecompositionUserNumber,Mesh1,SolidDecomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(SolidDecomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(SolidDecomposition,NumberOfComputationalNodes,Err)
  CALL cmfe_Decomposition_CalculateFacesSet(SolidDecomposition,.TRUE.,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(SolidDecomposition,Err)
  !Create a decomposition for the fluid mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID MESH DECOMPOSITION << == '
  CALL cmfe_Decomposition_Initialise(FluidDecomposition,Err)
  CALL cmfe_Decomposition_CreateStart(FluidDecompositionUserNumber,Mesh2,FluidDecomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(FluidDecomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(FluidDecomposition,NumberOfComputationalNodes,Err)
  CALL cmfe_Decomposition_CalculateFacesSet(FluidDecomposition,.TRUE.,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(FluidDecomposition,Err)
  !Create a decomposition for the interface mesh
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE DECOMPOSITION << == '
  CALL cmfe_Decomposition_Initialise(InterfaceDecomposition,Err)
  CALL cmfe_Decomposition_CreateStart(InterfaceDecompositionUserNumber,InterfaceMesh1,InterfaceDecomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(InterfaceDecomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(InterfaceDecomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(InterfaceDecomposition,Err)
  !
  !================================================================================================================================
  !
  ! G E O M E T R I C _ F I E L D S
  !Start to create a default (geometric) field on the solid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID MESH GEOMETRIC FIELD << == '
  CALL cmfe_Field_Initialise(GeometricField1,Err)
  CALL cmfe_Field_CreateStart(SolidGeometricFieldUserNumber,Region1,GeometricField1,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField1,SolidDecomposition,Err)
  !Set the scaling to use
  CALL cmfe_Field_ScalingTypeSet(GeometricField1,CMFE_FIELD_NO_SCALING,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,"SolidGF",Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NumberOfDimensions==3) CALL cmfe_Field_ComponentMeshComponentSet(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the first field
  CALL cmfe_Field_CreateFinish(GeometricField1,Err)
  !Start to create a default (geometric) field on the fluid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID MESH GEOMETRIC FIELD << == '
  CALL cmfe_Field_Initialise(GeometricField2,Err)
  CALL cmfe_Field_CreateStart(FluidGeometricFieldUserNumber,Region2,GeometricField2,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField2,FluidDecomposition,Err)
  !Set the scaling to use
  CALL cmfe_Field_ScalingTypeSet(GeometricField2,CMFE_FIELD_NO_SCALING,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,"FluidGF",Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NumberOfDimensions==3) CALL cmfe_Field_ComponentMeshComponentSet(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the second field
  CALL cmfe_Field_CreateFinish(GeometricField2,Err)
  !Start to create a default (geometric) field on the Interface
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE GEOMETRIC FIELD << == '
  CALL cmfe_Field_Initialise(InterfaceGeometricField1,Err)
  CALL cmfe_Field_CreateStart(InterfaceGeometricFieldUserNumber,Interface1,InterfaceGeometricField1,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(InterfaceGeometricField1,InterfaceDecomposition,Err)
  CALL cmfe_Field_VariableLabelSet(InterfaceGeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,"InterfaceGF",Err)
  !Set the domain to be used by the field components
  CALL cmfe_Field_ComponentMeshComponentSet(InterfaceGeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(InterfaceGeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  IF(NumberOfDimensions==3) CALL cmfe_Field_ComponentMeshComponentSet(InterfaceGeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the first field
  CALL cmfe_Field_CreateFinish(InterfaceGeometricField1,Err)
  !Update the geometric field parameters (solid)
  DO NodeIndex=1,NumberOfSolidNodes
    CALL cmfe_Decomposition_NodeDomainGet(SolidDecomposition,SolidNodeNumbers(NodeIndex),1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      IF(NumberOfDimensions==3) THEN
        Components=[2,3,1]
      ELSE
        Components=[1,2,0]
      ENDIF
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,SolidNodeNumbers(NodeIndex),Components(1),SolidGeometryY(NodeIndex),Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,SolidNodeNumbers(NodeIndex),Components(2),SolidGeometryZ(NodeIndex),Err)
      IF(NumberOfDimensions==3) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
          & 1,CMFE_NO_GLOBAL_DERIV,SolidNodeNumbers(NodeIndex),Components(3),SolidGeometryX(NodeIndex),Err)
      ENDIF
    ENDIF
  ENDDO
  !Update the geometric field parameters (fluid)
  DO NodeIndex=1,NumberOfFluidNodes
    CALL cmfe_Decomposition_NodeDomainGet(FluidDecomposition,FluidNodeNumbers(NodeIndex),1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,FluidNodeNumbers(NodeIndex),Components(1),FluidGeometryY(NodeIndex),Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,FluidNodeNumbers(NodeIndex),Components(2),FluidGeometryZ(NodeIndex),Err)
      IF(NumberOfDimensions==3) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
          & 1,CMFE_NO_GLOBAL_DERIV,FluidNodeNumbers(NodeIndex),Components(3),FluidGeometryX(NodeIndex),Err)
      ENDIF
    ENDIF
  ENDDO
  !Update the geometric field parameters (interface)
  DO NodeIndex=1,NumberOfInterfaceNodes
    CALL cmfe_Decomposition_NodeDomainGet(InterfaceDecomposition,InterfaceNodeNumbersForGeometry(NodeIndex),1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_Field_ParameterSetUpdateNode(InterfaceGeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,InterfaceNodeNumbersForGeometry(NodeIndex),Components(1),InterfaceGeometryY(NodeIndex),Err)
      CALL cmfe_Field_ParameterSetUpdateNode(InterfaceGeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
        & 1,CMFE_NO_GLOBAL_DERIV,InterfaceNodeNumbersForGeometry(NodeIndex),Components(2),InterfaceGeometryZ(NodeIndex),Err)
      IF(NumberOfDimensions==3) THEN
        CALL cmfe_Field_ParameterSetUpdateNode(InterfaceGeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
          & 1,CMFE_NO_GLOBAL_DERIV,InterfaceNodeNumbersForGeometry(NodeIndex),Components(3),InterfaceGeometryX(NodeIndex),Err)
      ENDIF
    ENDIF
  ENDDO
  CALL cmfe_Field_ParameterSetUpdateStart(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateStart(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricField2,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateStart(InterfaceGeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(InterfaceGeometricField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  !
  !================================================================================================================================
  !
  ! F I B R E _ F I E L D _ S O L I D
  !Create a fibre field and attach it to the geometric field 1
!  CALL cmfe_Field_Initialise(FibreField,Err)
!  CALL cmfe_Field_CreateStart(FibreFieldUserNumber,Region1,FibreField,Err)
!  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
!  CALL cmfe_Field_MeshDecompositionSet(FibreField,SolidDecomposition,Err)
!  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricField1,Err)
!  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"SolidFF",Err)
!  CALL cmfe_Field_ScalingTypeSet(FibreField,CMFE_FIELD_NO_SCALING,Err)
!  CALL cmfe_Field_CreateFinish(FibreField,Err)
  !
  !================================================================================================================================
  !
  ! E Q U A T I O N _ S E T S
  !Create the equations set for the solid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID EQUATION SET << == '
  !Create the equations_set for FiniteElasticity MooneyRivlin
  CALL cmfe_Field_Initialise(EquationsSetField1,Err)
  CALL cmfe_EquationsSet_Initialise(SolidEquationsSet,Err)
  CALL cmfe_EquationsSet_CreateStart(SolidEquationsSetUserNumber,Region1,GeometricField1,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],SolidEquationsSetFieldUserNumber, &
    & EquationsSetField1,SolidEquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(SolidEquationsSet,Err)
  !Create the equations set for the fluid region
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID EQUATION SET << == '
  !Create the equations set for ALE Navier-Stokes
  CALL cmfe_Field_Initialise(EquationsSetField2,Err)
  CALL cmfe_EquationsSet_Initialise(FluidEquationsSet,Err)
  CALL cmfe_EquationsSet_CreateStart(FluidEquationsSetUserNumber,Region2,GeometricField2, &
    & [CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS,CMFE_EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE],FluidEquationsSetFieldUserNumber,EquationsSetField2,FluidEquationsSet,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(FluidEquationsSet,Err)
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> MOVING MESH EQUATION SET << == '
  !Create the equations set for the moving mesh
  CALL cmfe_EquationsSet_Initialise(MovingMeshEquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetFieldMovingMesh,Err)
  CALL cmfe_EquationsSet_CreateStart(MovingMeshEquationsSetUserNumber,Region2,GeometricField2, &
    & [CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS,CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE],EquationsSetFieldMovingMeshUserNumber,EquationsSetFieldMovingMesh, &
    & MovingMeshEquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(MovingMeshEquationsSet,Err)
  !
  !================================================================================================================================
  !
  ! D E P E N D E N T _ F I E L D S
  !Create the equations set dependent field variables for the first equations set
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID DEPENDENT FIELD << == '
  !Create the dependent field
  CALL cmfe_Field_Initialise(DependentField1,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(SolidEquationsSet,SolidDependentFieldUserNumber,DependentField1,Err)
  CALL cmfe_Field_VariableLabelSet(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,"SolidDF",Err)
  CALL cmfe_Field_VariableLabelSet(DependentField1,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"SolidDerivatives",Err)
  DO component_idx=1,NumberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,component_idx,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField1,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,component_idx,1,Err)
  ENDDO
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,NumberOfDimensions+1, &
    & PressureMeshComponent,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField1,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,NumberOfDimensions+1, &
    & PressureMeshComponent,Err)
  IF(PressureMeshComponent==1) THEN
    CALL cmfe_Field_ComponentInterpolationSet(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,NumberOfDimensions+1, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(DependentField1,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,NumberOfDimensions+1, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ELSE
    CALL cmfe_Field_ComponentInterpolationSet(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,NumberOfDimensions+1, &
      & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(DependentField1,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,NumberOfDimensions+1, &
      & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  ENDIF
  CALL cmfe_Field_ScalingTypeSet(DependentField1,CMFE_FIELD_NO_SCALING,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(SolidEquationsSet,Err)
  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,2,DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  IF (NumberOfDimensions==3) THEN
    CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField1,CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE,3,DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
  ENDIF
  CALL cmfe_Field_ComponentValuesInitialise(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & NumberOfDimensions+1,-MooneyRivlin1,Err)!14.0_CMISSRP,Err)
  
  CALL cmfe_Field_ParameterSetUpdateStart(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> SOLID MATERIAL FIELD << == '
  !Create the material field
  CALL cmfe_Field_Initialise(MaterialField1,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(SolidEquationsSet,SolidMaterialFieldUserNumber,MaterialField1,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField1,CMFE_FIELD_U_VARIABLE_TYPE,"Material1",Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField1,CMFE_FIELD_V_VARIABLE_TYPE,"SolidDensity",Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(SolidEquationsSet,Err)
  !Set Mooney-Rivlin constants c10 and c01 (default?2.0 and 6.0) respectively
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & MooneyRivlin1,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
    & MooneyRivlin2,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField1,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & SolidDensity,Err)
  IF(GravityFlag) THEN
    IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOURCE FIELD - GRAVITY << == '
    !Create the source field with the gravity vector
    CALL cmfe_Field_Initialise(SourceField1,Err)
    CALL cmfe_EquationsSet_SourceCreateStart(SolidEquationsSet,SourceFieldUserNumber,SourceField1,Err)
    CALL cmfe_Field_ScalingTypeSet(SourceField1,CMFE_FIELD_NO_SCALING,Err)
    CALL cmfe_EquationsSet_SourceCreateFinish(SolidEquationsSet,Err)
    DO component_idx=1,NumberOfDimensions
      CALL cmfe_Field_ComponentValuesInitialise(SourceField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
          & component_idx,Gravity(component_idx),Err)
    ENDDO
  ENDIF
  !=========================
  !Create the equations set dependent field variables for dynamic Navier-Stokes
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID DEPENDENT FIELD << == '
  CALL cmfe_Field_Initialise(DependentField2,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(FluidEquationsSet,FluidDependentFieldUserNumber,DependentField2,Err)
  CALL cmfe_Field_VariableLabelSet(DependentField2,CMFE_FIELD_U_VARIABLE_TYPE,"FluidDF",Err)
  CALL cmfe_Field_VariableLabelSet(DependentField2,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"FluidDerivatives",Err)
  !Set the mesh component to be used by the field components.
  DO ComponentNumber=1,NumberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField2,CMFE_FIELD_U_VARIABLE_TYPE,ComponentNumber,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField2,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,ComponentNumber,1,Err)
  ENDDO
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField2,CMFE_FIELD_U_VARIABLE_TYPE,NumberOfDimensions+1, & 
    & Mesh2ComponentNumberPressure,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField2,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,NumberOfDimensions+1, & 
    & Mesh2ComponentNumberPressure,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(FluidEquationsSet,Err)
  !Initialise dependent field
  DO ComponentNumber=1,NumberOfDimensions
    CALL cmfe_Field_ComponentValuesInitialise(DependentField2,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
      & ComponentNumber,InitialFieldNavierStokes(ComponentNumber),Err)
  ENDDO
  !Initialise pressure component
  CALL cmfe_Field_ComponentValuesInitialise(DependentField2,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & NumberOfDimensions+1,0.1_CMISSRP,Err)
  !=========================
  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> MOVING MESH DEPENDENT FIELD << == '
  !Create the equations set dependent field variables for moving mesh
  CALL cmfe_Field_Initialise(DependentFieldMovingMesh,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(MovingMeshEquationsSet,DependentFieldMovingMeshUserNumber, & 
    & DependentFieldMovingMesh,Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE,"MovingMeshDF",Err)
  CALL cmfe_Field_VariableLabelSet(DependentFieldMovingMesh,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,"MovingMeshDerivatives",Err)
  !Set the mesh component to be used by the field components.
  DO ComponentNumber=1,NumberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE,ComponentNumber, & 
      & Mesh2ComponentNumberSpace,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentFieldMovingMesh,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,ComponentNumber, & 
      & Mesh2ComponentNumberSpace,Err)
  ENDDO
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(MovingMeshEquationsSet,Err)
  !Initialise dependent field moving mesh
  DO ComponentNumber=1,NumberOfDimensions
    CALL cmfe_Field_ComponentValuesInitialise(DependentFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
      & ComponentNumber,InitialFieldMovingMesh(ComponentNumber),Err)
  ENDDO
  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> FLUID MATERIAL FIELD << == '
  !Create the equations set materials field variables for dynamic Navier-Stokes
  CALL cmfe_Field_Initialise(MaterialField2,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(FluidEquationsSet,FluidMaterialFieldUserNumber,MaterialField2,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(FluidEquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField2,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & FluidMaterialFieldComponentMu,FluidDynamicViscosity,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField2,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & FluidMaterialFieldComponentRho,FluidDensity,Err)
  IF(ExampleFileProgressDiagnostics) PRINT *, '    >> MATERIAL FIELD MOVING MESH << == '
  !Create the equations set materials field variables for moving mesh
  CALL cmfe_Field_Initialise(MaterialFieldMovingMesh,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(MovingMeshEquationsSet,MaterialFieldMovingMeshUserNumber, &
  & MaterialFieldMovingMesh,Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(MovingMeshEquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & MaterialFieldMovingMeshUserNumberK,MovingMeshParameterK,Err)
  !
  !================================================================================================================================
  !
  ! I N D E P E N D E N T _ F I E L D // mesh velocity for fluid equations set && mesh stiffness for moving mesh equations set
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID INDEPENDENT FIELD << == '
  !Create the equations set independent field variables for ALE Navier-Stokes
  CALL cmfe_Field_Initialise(IndependentField2,Err)
  CALL cmfe_EquationsSet_IndependentCreateStart(FluidEquationsSet,IndependentField2UserNumber,IndependentField2,Err)
  CALL cmfe_Field_VariableLabelSet(IndependentField2,CMFE_FIELD_U_VARIABLE_TYPE,"FluidInDF",Err)
  !Set the mesh component to be used by the field components.
  DO ComponentNumber=1,NumberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(IndependentField2,CMFE_FIELD_U_VARIABLE_TYPE,ComponentNumber, & 
      & Mesh2ComponentNumberSpace,Err)
  ENDDO
  !Finish the equations set independent field variables
  CALL cmfe_EquationsSet_IndependentCreateFinish(FluidEquationsSet,Err)
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INDEPENDENT FIELD MOVING MESH << == '
  !Create the equations set independent field variables for moving mesh
  CALL cmfe_Field_Initialise(IndependentFieldMovingMesh,Err)
  CALL cmfe_EquationsSet_IndependentCreateStart(MovingMeshEquationsSet,IndependentFieldMovingMeshUserNumber, & 
    & IndependentFieldMovingMesh,Err)
  !Set the scaling to use
  CALL cmfe_Field_ScalingTypeSet(IndependentFieldMovingMesh,CMFE_FIELD_NO_SCALING,Err)
  CALL cmfe_Field_VariableLabelSet(IndependentFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE,"MovingMeshInDF",Err)
  !Set the mesh component to be used by the field components.
  DO ComponentNumber=1,NumberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(IndependentFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE,ComponentNumber, & 
      & Mesh2ComponentNumberSpace,Err)
  ENDDO
  !Finish the equations set independent field variables
  CALL cmfe_EquationsSet_IndependentCreateFinish(MovingMeshEquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(IndependentFieldMovingMesh,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
    & IndependentFieldMovingMeshUserNumberK,MovingMeshParameterK,Err)
  !
  !================================================================================================================================
  !
  IF(GeometryCheck) THEN
    IF(ExampleFileProgressDiagnostics) WRITE(*,'(A)') "Exporting fields..."
    !Export initial fields
    CALL cmfe_Fields_Initialise(Fields1,Err)
    CALL cmfe_Fields_Create(Region1,Fields1,Err)
    CALL cmfe_Fields_NodesExport(Fields1,"GeometryCheckSolid","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields1,"GeometryCheckSolid","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields1,Err)
    CALL cmfe_Fields_Initialise(Fields2,Err)
    CALL cmfe_Fields_Create(Region2,Fields2,Err)
    CALL cmfe_Fields_NodesExport(Fields2,"GeometryCheckFluid","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields2,"GeometryCheckFluid","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields2,Err)
    CALL cmfe_Fields_Initialise(FieldsI,Err)
    CALL cmfe_Fields_Create(Interface1,FieldsI,Err)
    CALL cmfe_Fields_NodesExport(FieldsI,"GeometryCheckInterface","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(FieldsI,"GeometryCheckInterface","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(FieldsI,Err)
    IF(ExampleFileProgressDiagnostics) WRITE(*,'(A)') "Field exported!"
  ENDIF
  !
  !================================================================================================================================
  !
  ! E Q U A T I O N S
  !Create the equations set equations for the first equations set
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLID EQUATIONS << == '
  CALL cmfe_Equations_Initialise(Equations1,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(SolidEquationsSet,Equations1,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations1,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(Equations1,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations1,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations1,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
  !CALL cmfe_Equations_OutputTypeSet(Equations1,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(SolidEquationsSet,Err)
  !Create the equations set equations for the second equations set
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> FLUID EQUATIONS << == '
  CALL cmfe_Equations_Initialise(Equations2,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(FluidEquationsSet,Equations2,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations2,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output (NO_OUTPUT/TIMING_OUTPUT/MATRIX_OUTPUT/ELEMENT_MATRIX_OUTPUT)
  CALL cmfe_Equations_OutputTypeSet(Equations2,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(FluidEquationsSet,Err)
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> MOVING MESH EQUATIONS << == '
  !Create the equations set equations
  CALL cmfe_Equations_Initialise(EquationsMovingMesh,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(MovingMeshEquationsSet,EquationsMovingMesh,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(EquationsMovingMesh,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(EquationsMovingMesh,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(MovingMeshEquationsSet,Err)
  !
  !================================================================================================================================
  !
  ! I N T E R F A C E _ C O N D I T I O N
  !Create an interface condition between the two meshes
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE CONDITIONS << == '
  CALL cmfe_InterfaceCondition_Initialise(InterfaceCondition,Err)
  CALL cmfe_InterfaceCondition_CreateStart(InterfaceConditionUserNumber,Interface1,InterfaceGeometricField1, &
    & InterfaceCondition,Err)
  !Specify the method for the interface condition
  CALL cmfe_InterfaceCondition_MethodSet(InterfaceCondition,CMFE_INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,Err)
  !Specify the type of interface condition operator
  CALL cmfe_InterfaceCondition_OperatorSet(InterfaceCondition,CMFE_INTERFACE_CONDITION_SOLID_FLUID_OPERATOR,Err)
  !Add in the dependent variables from the equations sets
  CALL cmfe_InterfaceCondition_DependentVariableAdd(InterfaceCondition,SolidMeshIndex,SolidEquationsSet, &
    & CMFE_FIELD_U_VARIABLE_TYPE,Err)
  CALL cmfe_InterfaceCondition_DependentVariableAdd(InterfaceCondition,FluidMeshIndex,FluidEquationsSet, &
    & CMFE_FIELD_U_VARIABLE_TYPE,Err)
  !Finish creating the interface condition
  CALL cmfe_InterfaceCondition_CreateFinish(InterfaceCondition,Err)
  !Create the Lagrange multipliers field
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE LAGRANGE FIELD << == '
  CALL cmfe_Field_Initialise(LagrangeField1,Err)
  CALL cmfe_InterfaceCondition_LagrangeFieldCreateStart(InterfaceCondition,LagrangeFieldUserNumber,LagrangeField1,Err)
  CALL cmfe_Field_VariableLabelSet(LagrangeField1,CMFE_FIELD_U_VARIABLE_TYPE,"InterfaceLF",Err)
  !Finish the Lagrange multipliers field
  CALL cmfe_InterfaceCondition_LagrangeFieldCreateFinish(InterfaceCondition,Err)
  DO ComponentNumber=1,NumberOfDimensions
    CALL cmfe_Field_ComponentValuesInitialise(LagrangeField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
      & ComponentNumber,0.0_CMISSRP,Err)
  ENDDO
  CALL cmfe_Field_ParameterSetUpdateStart(LagrangeField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(LagrangeField1,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  !Create the interface condition equations
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> INTERFACE EQUATIONS << == '
  CALL cmfe_InterfaceEquations_Initialise(InterfaceEquations,Err)
  CALL cmfe_InterfaceCondition_EquationsCreateStart(InterfaceCondition,InterfaceEquations,Err)
  !Set the interface equations sparsity
  CALL cmfe_InterfaceEquations_SparsitySet(InterfaceEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the interface equations output
  CALL cmfe_InterfaceEquations_OutputTypeSet(InterfaceEquations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !Finish creating the interface equations
  CALL cmfe_InterfaceCondition_EquationsCreateFinish(InterfaceCondition,Err)
  !
  !================================================================================================================================
  !
  ! P R O B L E M
  !Start the creation of a coupled problem
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> COUPLED PROBLEM << == '
  CALL cmfe_Problem_Initialise(CoupledProblem,Err)
  CALL cmfe_Problem_CreateStart(CoupledProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
    & CMFE_PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE,CMFE_PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE], &
    & CoupledProblem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(CoupledProblem,Err)  
  !
  !================================================================================================================================
  !
  ! C O N T R O L _ L O O P
  !Start the creation of the problem control loop for the coupled problem
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> PROBLEM CONTROL LOOP << == '
  CALL cmfe_Problem_ControlLoopCreateStart(CoupledProblem,Err)
  !Time loop (main loop)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)  
  CALL cmfe_Problem_ControlLoopGet(CoupledProblem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL cmfe_ControlLoop_LabelSet(ControlLoop,'TimeLoop',Err)
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,StartTime,StopTime,TimeStepSize,Err)
  CALL cmfe_ControlLoop_TimeInputSet(ControlLoop,MaterialSpecification,Err)
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,OutputFrequency,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(CoupledProblem,Err)
  !
  !================================================================================================================================
  !
  ! S O L V E R S
  !Start the creation of the problem solver for the coupled problem
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> PROBLEM SOLVERS << == '
  CALL cmfe_Problem_SolversCreateStart(CoupledProblem,Err)
  !Linear solver for moving mesh
  CALL cmfe_Solver_Initialise(LinearSolverMovingMesh,Err)
  CALL cmfe_Problem_SolverGet(CoupledProblem,CMFE_CONTROL_LOOP_NODE, &
    & LinearSolverMovingMeshIndex,LinearSolverMovingMesh,Err)
  CALL cmfe_Solver_OutputTypeSet(LinearSolverMovingMesh,LinearSolverMovingMesh_OutputType,Err)
  !Solvers for coupled FiniteElasticity NavierStokes problem
  CALL cmfe_Solver_Initialise(DynamicSolver,Err)
  CALL cmfe_Solver_Initialise(NonlinearSolver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  !Get the dynamic ALE solver
  CALL cmfe_Problem_SolverGet(CoupledProblem,CMFE_CONTROL_LOOP_NODE,DynamicSolverIndex,DynamicSolver,Err)
  CALL cmfe_Solver_OutputTypeSet(DynamicSolver,DynamicSolver_OutputType,Err)
  CALL cmfe_Solver_DynamicThetaSet(DynamicSolver,DynamicSolver_Theta,Err)
  !Get the dynamic nonlinear solver
  CALL cmfe_Solver_DynamicNonlinearSolverGet(DynamicSolver,NonlinearSolver,Err)
  CALL cmfe_Solver_NewtonLineSearchTypeSet(NonlinearSolver,CMFE_SOLVER_NEWTON_LINESEARCH_LINEAR,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(NonlinearSolver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL cmfe_Solver_NewtonMaximumFunctionEvaluationsSet(NonlinearSolver,MaxFunctionEvaluations,Err)
  CALL cmfe_Solver_OutputTypeSet(NonlinearSolver,NonlinearSolver_OutputType,Err)
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(NonlinearSolver,AbsoluteTolerance,Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(NonlinearSolver,RelativeTolerance,Err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(NonlinearSolver,MaximumIterations,Err)
  CALL cmfe_Solver_NewtonLineSearchAlphaSet(NonlinearSolver,LinesearchAlpha,Err)
  !Get the dynamic nonlinear linear solver
  CALL cmfe_Solver_NewtonLinearSolverGet(NonlinearSolver,LinearSolver,Err)
  !Choose type of linear solver
  IF(.FALSE.) THEN
    CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
    !NO,JACOBI,BLOCK_JACOBI,SOR,INCOMPLETE_CHOLESKY,INCOMPLETE_LU,ADDITIVE_SCHWARZ
    CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(LinearSolver,CMFE_SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER,Err)
    CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolver,MaximumIterations,Err)
    CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(LinearSolver,AbsoluteTolerance,Err)
    CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(LinearSolver,RelativeTolerance,Err)
    CALL cmfe_Solver_LinearIterativeDivergenceToleranceSet(LinearSolver,DivergenceTolerance,Err)
  ELSE
    CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  ENDIF
  CALL cmfe_Solver_OutputTypeSet(LinearSolver,LinearSolver_OutputType,Err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(CoupledProblem,Err)
  !
  !================================================================================================================================
  !
  ! S O L V E R _ E Q U A T I O N S
  CALL cmfe_Problem_SolverEquationsCreateStart(CoupledProblem,Err)
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> MOVING MESH SOLVER EQUATIONS << == '
  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(LinearSolverMovingMesh,Err)
  CALL cmfe_SolverEquations_Initialise(LinearSolverMovingMeshEquations,Err)
  !Get the linear solver equations
  CALL cmfe_Problem_SolverGet(CoupledProblem,CMFE_CONTROL_LOOP_NODE,LinearSolverMovingMeshIndex,LinearSolverMovingMesh,Err)
  CALL cmfe_Solver_SolverEquationsGet(LinearSolverMovingMesh,LinearSolverMovingMeshEquations,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(LinearSolverMovingMeshEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(LinearSolverMovingMeshEquations,MovingMeshEquationsSet,MovingMeshEquationsSetIndex,Err)
  IF(ExampleFileProgressDiagnostics) PRINT *, ' == >> SOLVER EQUATIONS << == '
  CALL cmfe_Solver_Initialise(DynamicSolver,Err)
  CALL cmfe_SolverEquations_Initialise(CoupledSolverEquations,Err)
  !Get the dynamic solver equations
  CALL cmfe_Problem_SolverGet(CoupledProblem,CMFE_CONTROL_LOOP_NODE,DynamicSolverIndex,DynamicSolver,Err)
  CALL cmfe_Solver_SolverEquationsGet(DynamicSolver,CoupledSolverEquations,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(CoupledSolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(CoupledSolverEquations,SolidEquationsSet,SolidEquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(CoupledSolverEquations,FluidEquationsSet,FluidEquationsSetIndex,Err)
  CALL cmfe_SolverEquations_InterfaceConditionAdd(CoupledSolverEquations,InterfaceCondition,InterfaceConditionIndex,Err)
  !Set the time dependence of the interface matrix to determine the interface matrix coefficient in the solver matrix
  ! (basically position in big coupled matrix system)
  CALL cmfe_InterfaceMatrices_TimeDependenceTypeSet(InterfaceCondition,SolidEquationsSetIndex,.TRUE., &
    & [CMFE_INTERFACE_MATRIX_STATIC,CMFE_INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC],Err)
  CALL cmfe_InterfaceMatrices_TimeDependenceTypeSet(InterfaceCondition,FluidEquationsSetIndex,.TRUE., &
    & [CMFE_INTERFACE_MATRIX_STATIC,CMFE_INTERFACE_MATRIX_STATIC],Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(CoupledProblem,Err)
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
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(CoupledSolverEquations,BoundaryConditions,Err)
  !No displacement boundary for solid
  DO S=1,SIZE(NoDisplacementNodes)
    NodeNumber=NoDisplacementNodes(S)
    CALL cmfe_Decomposition_NodeDomainGet(SolidDecomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
      CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
      IF(NumberOfDimensions==3) THEN
        CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField1,CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
      ENDIF
    ENDIF
  ENDDO
  !Set outlet (zero) pressure nodes
  DO S=1,SIZE(OutletNodes)
    NodeNumber=OutletNodes(S)
    CALL cmfe_Decomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,NumberOfDimensions+1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO
  !Inlet velocity nodes, small starting velocity in 1st coordinate direction
  DO S=1,SIZE(InletNodes)
    NodeNumber=InletNodes(S)
    CALL cmfe_Decomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,1,CMFE_BOUNDARY_CONDITION_FIXED_INLET,0.0_CMISSRP,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,2,CMFE_BOUNDARY_CONDITION_FIXED_INLET,0.0_CMISSRP,Err)
      IF(NumberOfDimensions==3) THEN
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
          & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,3,CMFE_BOUNDARY_CONDITION_FIXED_INLET,0.0_CMISSRP,Err)
      ENDIF
    ENDIF
  ENDDO
  !Set no-slip BC
  DO S=1,SIZE(NoSlipNodes)
    NodeNumber=NoSlipNodes(S)
    CALL cmfe_Decomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
      IF(NumberOfDimensions==3) THEN
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
          & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
      ENDIF
    ENDIF
  ENDDO
  !Set slip BC
  DO S=1,SIZE(SlipNodesTop)
    NodeNumber=SlipNodesTop(S)
    CALL cmfe_Decomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,NumberOfDimensions,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO
  !Set slip BC
  IF(NumberOfDimensions==3) THEN
    DO S=1,SIZE(SlipNodesRightLeft)
      NodeNumber=SlipNodesRightLeft(S)
      CALL cmfe_Decomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField2, &
          & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
      ENDIF
    ENDDO
  ENDIF
  IF(NumberOfDimensions==2) THEN
    !Remove dof's at nodes where solid displacement and zero velocity is set (first n last interface node)
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
      & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
      & 1,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
      & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
      & 1,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
      & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
      & NumberOfInterfaceNodes,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
      & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
      & NumberOfInterfaceNodes,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
  ELSE
    IF(CheckWithoutInterfaceCondition) THEN
      DO S=1,SIZE(InterfaceNodeNumbersForGeometry)
        NodeNumber=InterfaceNodeNumbersForGeometry(S)
        CALL cmfe_Decomposition_NodeDomainGet(InterfaceDecomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
            & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
            & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
            & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
        ENDIF
      ENDDO
    ELSE
      DO S=1,SIZE(LagrangeNodes)
        NodeNumber=LagrangeNodes(S)
        CALL cmfe_Decomposition_NodeDomainGet(InterfaceDecomposition,NodeNumber,1,NodeDomain,Err)
        IF(NodeDomain==ComputationalNodeNumber) THEN
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
            & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
            & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
          CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,LagrangeField1, &
            & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
            & NodeNumber,3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
        ENDIF
      ENDDO
    ENDIF
  ENDIF
  !Finish equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(CoupledSolverEquations,Err)
  
  !Start the creation of the moving mesh boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditionsMovingMesh,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(LinearSolverMovingMeshEquations,BoundaryConditionsMovingMesh,Err)
  !Fixed boundary nodes. May be used to move nodes..
  DO S=1,SIZE(MovedYNodes)
    NodeNumber=MovedYNodes(S)
    CALL cmfe_Decomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,1,CMFE_BOUNDARY_CONDITION_FIXED_WALL,0.0_CMISSRP,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,2,CMFE_BOUNDARY_CONDITION_FIXED_WALL,0.0_CMISSRP,Err)
      IF(NumberOfDimensions==3) THEN
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
          & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,3,CMFE_BOUNDARY_CONDITION_FIXED_WALL,0.0_CMISSRP,Err)
      ENDIF
    ENDIF
  ENDDO
  !Mesh nodes that are moving wall nodes
  DO S=1,SIZE(MovedNodes)
    NodeNumber=MovedNodes(S)
    CALL cmfe_Decomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,1,CMFE_BOUNDARY_CONDITION_MOVED_WALL,0.0_CMISSRP,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,2,CMFE_BOUNDARY_CONDITION_MOVED_WALL,0.0_CMISSRP,Err)
      IF(NumberOfDimensions==3) THEN
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
          & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,3,CMFE_BOUNDARY_CONDITION_MOVED_WALL,0.0_CMISSRP,Err)
      ENDIF
    ENDIF
  ENDDO
  !Mesh nodes that are fixed in space
  DO S=1,SIZE(FixedNodes)
    NodeNumber=FixedNodes(S)
    CALL cmfe_Decomposition_NodeDomainGet(FluidDecomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,1,CMFE_BOUNDARY_CONDITION_FIXED_WALL,0.0_CMISSRP,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
        & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
        & NodeNumber,2,CMFE_BOUNDARY_CONDITION_FIXED_WALL,0.0_CMISSRP,Err)
      IF(NumberOfDimensions==3) THEN
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditionsMovingMesh,DependentFieldMovingMesh, &
          & CMFE_FIELD_U_VARIABLE_TYPE,1,1, &
          & NodeNumber,3,CMFE_BOUNDARY_CONDITION_FIXED_WALL,0.0_CMISSRP,Err)
      ENDIF
    ENDIF
  ENDDO
  !Finish moving mesh boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(LinearSolverMovingMeshEquations,Err)
  
  !
  !================================================================================================================================
  !
  ! S O L V E
  
  PRINT *, ' == >> SOLVING PROBLEM << == '
  CALL cmfe_Problem_Solve(CoupledProblem,Err)
  
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  e=etime(t)/60.0_CMISSRP
  PRINT *, "Program successfully completed in ",e," minutes."
  
  STOP
  
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)
    
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING
    
    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP
    
  END SUBROUTINE HANDLE_ERROR
  
END PROGRAM FortranExample
