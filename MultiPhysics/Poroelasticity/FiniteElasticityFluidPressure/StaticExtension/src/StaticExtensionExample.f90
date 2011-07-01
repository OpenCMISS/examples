!> \file
!> \author Adam Reeve
!> \brief This is an example program to solve a coupled finite elasticity and fluid pressure problem using OpenCMISS calls.
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

!> \example MultiPhysics/Poroelasticity/FiniteElasticityFluidPressure/StaticExtension/src/StaticExtensionExample.f90
!! Example program to solve a coupled finite elasticity and fluid pressure problem using OpenCMISS calls.
!<

!> Main program
PROGRAM POROELASTICITYEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: Height=10.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: Width=10.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: Length=10.0_CMISSDP
  INTEGER(CMISSIntg), PARAMETER :: GeometricInterpolationType=CMISSBasisQuadraticLagrangeInterpolation
  INTEGER(CMISSIntg), PARAMETER :: PressureInterpolationType=CMISSBasisLinearLagrangeInterpolation
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussXi=3
  INTEGER(CMISSIntg) :: NumberElements
  INTEGER(CMISSIntg) :: NumberGlobalXElements=2
  INTEGER(CMISSIntg) :: NumberGlobalYElements=2
  INTEGER(CMISSIntg) :: NumberGlobalZElements=2
  REAL(CMISSDP) :: FixedWidth=12.0_CMISSDP
  REAL(CMISSDP) :: FluidPressureBC=0.5E-3_CMISSDP
  REAL(CMISSDP) :: FluidPressureBC2=0.0_CMISSDP
  REAL(CMISSDP) :: InitialPressure
  INTEGER(CMISSIntg) :: Increments=10
  LOGICAL :: FixRightSide=.TRUE.
  REAL (CMISSDP) :: K1=2.0E-3_CMISSDP
  REAL (CMISSDP) :: K2=33.0E-6_CMISSDP
  REAL (CMISSDP) :: K=2.20E-1_CMISSDP
  REAL (CMISSDP) :: M=2.18E-1_CMISSDP
  REAL (CMISSDP) :: b=1.0_CMISSDP
  REAL (CMISSDP) :: p_0=0.0_CMISSDP
  REAL (CMISSDP) :: permeability=1.0E-3_CMISSDP

  !Object user numbers

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricMeshComponent=1
  INTEGER(CMISSIntg), PARAMETER :: PressureMeshComponent=2
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FibreFieldUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldFluidMaterialUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldSolidMaterialUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: SolidEquationsSetFieldUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: FluidEquationsSetFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: FluidEquationsSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolidEquationsSetUserNumber=2

  !Program variables
  INTEGER(CMISSIntg),ALLOCATABLE :: BottomSurfaceNodes(:),TopSurfaceNodes(:),LeftSurfaceNodes(:), &
    & RightSurfaceNodes(:),FrontSurfaceNodes(:),BackSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: PressureBottomSurfaceNodes(:),PressureTopSurfaceNodes(:),PressureLeftSurfaceNodes(:), &
    & PressureRightSurfaceNodes(:),PressureFrontSurfaceNodes(:),PressureBackSurfaceNodes(:)
  INTEGER(CMISSIntg) :: TopNormalXi,LeftNormalXi,RightNormalXi,FrontNormalXi,BackNormalXi,BottomNormalXi
  INTEGER(CMISSIntg) :: node_idx,NodeNumber,NodeDomain
  INTEGER(CMISSIntg) :: NumberOfDimensions,component_idx
  CHARACTER(LEN=255) :: Filename
  INTEGER(CMISSIntg) :: NumberOfArguments,ArgumentLength,ArgStatus
  CHARACTER(LEN=255) :: CommandArgument

  !CMISS variables

  TYPE(CMISSBasisType) :: GeometricBasis,PressureBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: SolidEquations,FluidEquations
  TYPE(CMISSEquationsSetType) :: SolidEquationsSet,FluidEquationsSet
  TYPE(CMISSFieldType) :: GeometricField,SolidEquationsSetField,FluidEquationsSetField,DependentField, &
    & FibreField,SolidMaterialField,FluidMaterialField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations

  !Generic CMISS variables

  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: Err

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)
  CALL CMISSDiagnosticsSetOn(CMISSInDiagType,[1,2,3,4,5],"Diagnostics",[ &
      & "EQUATIONS_SET_LOAD_INCREMENT_APPLY   " &
      ],Err)

  !Read in arguments and overwrite default values
  NumberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(NumberOfArguments >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(CommandArgument(1:ArgumentLength),*) FluidPressureBC
  ENDIF
  IF(NumberOfArguments >= 2) THEN
    CALL GET_COMMAND_ARGUMENT(2,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(CommandArgument(1:ArgumentLength),*) FixedWidth
  ENDIF
  IF(NumberOfArguments >= 3) THEN
    CALL GET_COMMAND_ARGUMENT(3,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 3.")
    READ(CommandArgument(1:ArgumentLength),*) Increments
    IF(Increments<1) CALL HANDLE_ERROR("Invalid number of loading increments.")
  ENDIF
  IF(NumberOfArguments >= 4) THEN
    CALL GET_COMMAND_ARGUMENT(4,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 4.")
    READ(CommandArgument(1:ArgumentLength),*) NumberElements
    IF(NumberElements<1) CALL HANDLE_ERROR("Invalid number of elements.")
    NumberGlobalXElements=NumberElements
    NumberGlobalYElements=NumberElements
    NumberGlobalZElements=NumberElements
  ENDIF
  IF(FixedWidth<1.0E-5_CMISSDP) FixRightSide=.FALSE.
  WRITE(*,*) "Fluid pressure BC:", FluidPressureBC
  IF (FixRightSide) THEN
    WRITE(*,*) "Fixed width:", FixedWidth
  ELSE
    WRITE(*,*) "Width not fixed"
  ENDIF
  WRITE(*,*) "Increments:", Increments
  WRITE(*,*) "Elements:", NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  InitialPressure=(FluidPressureBC+FluidPressureBC2)*0.5_CMISSDP/REAL(Increments)

  WRITE(Filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "StaticExtension",NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements,GeometricInterpolationType
  CALL CMISSOutputSetOn(Filename,Err)

  !Get the computational nodes information
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !Create a 3D rectangular cartesian coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !Define basis function and second basis function for pressure etc
  CALL CMISSBasisTypeInitialise(GeometricBasis,Err)
  CALL CMISSBasisTypeInitialise(PressureBasis,Err)
  CALL CMISSBasisCreateStart(GeometricBasisUserNumber,GeometricBasis,Err)
  CALL CMISSBasisCreateStart(PressureBasisUserNumber,PressureBasis,Err)
  SELECT CASE(GeometricInterpolationType)
  CASE(1,2,3,4)
    CALL CMISSBasisTypeSet(GeometricBasis,CMISSBasisLagrangeHermiteTPType,Err)
    CALL CMISSBasisTypeSet(PressureBasis,CMISSBasisLagrangeHermiteTPType,Err)
  CASE(7,8,9)
    CALL CMISSBasisTypeSet(GeometricBasis,CMISSBasisSimplexType,Err)
    CALL CMISSBasisTypeSet(PressureBasis,CMISSBasisSimplexType,Err)
  END SELECT
  IF(NumberGlobalZElements==0) THEN
    NumberOfDimensions=2
    CALL CMISSBasisNumberOfXiSet(GeometricBasis,2,Err)
    CALL CMISSBasisNumberOfXiSet(PressureBasis,2,Err)
    CALL CMISSBasisInterpolationXiSet(GeometricBasis,[GeometricInterpolationType,GeometricInterpolationType],Err)
    CALL CMISSBasisInterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      !Finite elasticity doesn't support different numbers of gauss points for different components
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(GeometricBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ELSE
    NumberOfDimensions=3
    CALL CMISSBasisNumberOfXiSet(GeometricBasis,3,Err)
    CALL CMISSBasisNumberOfXiSet(PressureBasis,3,Err)
    CALL CMISSBasisInterpolationXiSet(GeometricBasis,[GeometricInterpolationType,GeometricInterpolationType, &
        & GeometricInterpolationType],Err)
    CALL CMISSBasisInterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType, &
        & PressureInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(GeometricBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ENDIF
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(GeometricBasis,.TRUE.,Err)
  CALL CMISSBasisQuadratureLocalFaceGaussEvaluateSet(PressureBasis,.TRUE.,Err)
  CALL CMISSBasisCreateFinish(GeometricBasis,Err)
  CALL CMISSBasisCreateFinish(PressureBasis,Err)

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the bases
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,[GeometricBasis,PressureBasis],Err)
  !Define the mesh on the region
  IF(NumberGlobalXElements==0) THEN
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[Width,Height],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[Width,Height,Length],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
      & NumberGlobalZElements],Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL CMISSMeshTypeInitialise(Mesh,Err)
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionMeshComponentSet(Decomposition,1,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  CALL CMISSDecompositionCalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

  !Create a fibre field and attach it to the geometric field
  CALL CMISSFieldTypeInitialise(FibreField,Err)
  CALL CMISSFieldCreateStart(FibreFieldUserNumber,Region,FibreField,Err)
  CALL CMISSFieldTypeSet(FibreField,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreField,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSFieldVariableLabelSet(FibreField,CMISSFieldUVariableType,"Fibre",Err)
  CALL CMISSFieldCreateFinish(FibreField,Err)

  !Create the equations sets
  CALL CMISSEquationsSetTypeInitialise(SolidEquationsSet,Err)
  CALL CMISSFieldTypeInitialise(SolidEquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(SolidEquationsSetUserNumber,Region,FibreField,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetElasticityFluidPressureStaticSubtype, &
    & SolidEquationsSetFieldUserNumber,SolidEquationsSetField,SolidEquationsSet,Err)
  CALL CMISSEquationsSetCreateFinish(SolidEquationsSet,Err)

  CALL CMISSEquationsSetTypeInitialise(FluidEquationsSet,Err)
  CALL CMISSFieldTypeInitialise(FluidEquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(FluidEquationsSetUserNumber,Region,FibreField,CMISSEquationsSetFluidMechanicsClass, &
    & CMISSEquationsSetDarcyPressureEquationType,CMISSEquationsSetElasticityFluidPressureStaticSubtype, &
    & FluidEquationsSetFieldUserNumber,FluidEquationsSetField,FluidEquationsSet,Err)
  CALL CMISSEquationsSetCreateFinish(FluidEquationsSet,Err)

  !Create the dependent field
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(FluidEquationsSet,FieldDependentUserNumber,DependentField,Err)
  DO component_idx=1,NumberOfDimensions
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,component_idx,GeometricMeshComponent,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,component_idx,GeometricMeshComponent,Err)
  ENDDO
  DO component_idx=1,1
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldVVariableType,component_idx,PressureMeshComponent,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelVDelNVariableType,component_idx,PressureMeshComponent,Err)
  ENDDO
  CALL CMISSEquationsSetDependentCreateFinish(FluidEquationsSet,Err)
  CALL CMISSEquationsSetDependentCreateStart(SolidEquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSetDependentCreateFinish(SolidEquationsSet,Err)

  !Create the material field
  CALL CMISSFieldTypeInitialise(FluidMaterialField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(FluidEquationsSet,FieldFluidMaterialUserNumber,FluidMaterialField,Err)
  CALL CMISSFieldVariableLabelSet(FluidMaterialField,CMISSFieldUVariableType,"FluidMaterial",Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(FluidEquationsSet,Err)

  CALL CMISSFieldTypeInitialise(SolidMaterialField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(SolidEquationsSet,FieldSolidMaterialUserNumber,SolidMaterialField,Err)
  CALL CMISSFieldVariableLabelSet(SolidMaterialField,CMISSFieldUVariableType,"SolidMaterial",Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(SolidEquationsSet,Err)

  !Set material constants
  !Solid constitutive law
  CALL CMISSFieldComponentValuesInitialise(SolidMaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,K1,Err)
  CALL CMISSFieldComponentValuesInitialise(SolidMaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,K2,Err)
  CALL CMISSFieldComponentValuesInitialise(SolidMaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,K,Err)
  CALL CMISSFieldComponentValuesInitialise(SolidMaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,M,Err)
  CALL CMISSFieldComponentValuesInitialise(SolidMaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,b,Err)
  CALL CMISSFieldComponentValuesInitialise(SolidMaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,6,p_0,Err)

  !Permeability tensor
  CALL CMISSFieldComponentValuesInitialise(FluidMaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,permeability,Err)
  CALL CMISSFieldComponentValuesInitialise(FluidMaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(FluidMaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(FluidMaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,permeability,Err)
  CALL CMISSFieldComponentValuesInitialise(FluidMaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(FluidMaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,6,permeability,Err)
  !Density
  CALL CMISSFieldComponentValuesInitialise(FluidMaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,7,1000.0_CMISSDP,Err)

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(FluidEquations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(FluidEquationsSet,FluidEquations,Err)
  CALL CMISSEquationsSparsityTypeSet(FluidEquations,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(FluidEquations,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(FluidEquationsSet,Err)

  CALL CMISSEquationsTypeInitialise(SolidEquations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(SolidEquationsSet,SolidEquations,Err)
  CALL CMISSEquationsSparsityTypeSet(SolidEquations,CMISSEquationsSparseMatrices,Err)
  CALL CMISSEquationsOutputTypeSet(SolidEquations,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsSetEquationsCreateFinish(SolidEquationsSet,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 1,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 2,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,Err)
  CALL CMISSFieldParametersToFieldParametersComponentCopy(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
    & 3,DependentField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,Err)
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldVVariableType,CMISSFieldValuesSetType,1,InitialPressure,Err)

  !Define the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemMultiPhysicsClass,CMISSProblemFiniteElasticityFluidPressureType, &
    & CMISSProblemStandardElasticityFluidPressureSubtype,Err)
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  CALL CMISSProblemControlLoopGet(Problem,[CMISSControlLoopNode],ControlLoop,Err)
  CALL CMISSControlLoopMaximumIterationsSet(ControlLoop,Increments,Err)
  CALL CMISSControlLoopOutputTypeSet(ControlLoop,CMISSControlLoopProgressOutput,Err)
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  CALL CMISSSolverNewtonJacobianCalculationTypeSet(Solver,CMISSSolverNewtonJacobianFDCalculated,Err)
  CALL CMISSSolverNewtonMaximumIterationsSet(Solver,200,Err)
  CALL CMISSSolverNewtonMaximumFunctionEvaluationsSet(Solver,10000,Err)
  CALL CMISSSolverNewtonAbsoluteToleranceSet(Solver,1.0E-10_CMISSDP,Err)
  CALL CMISSSolverNewtonRelativeToleranceSet(Solver,1.0E-9_CMISSDP,Err)
  CALL CMISSSolverNewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearIterativeSolveType,Err)
  CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolver,1000,Err)
  CALL CMISSProblemSolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,SolidEquationsSet,EquationsSetIndex,Err)
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,FluidEquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Set the boundary conditions for the solver equations
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Get geometric nodes
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMISSGeneratedMeshRegularBottomSurface,BottomSurfaceNodes,BottomNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMISSGeneratedMeshRegularTopSurface,TopSurfaceNodes,TopNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMISSGeneratedMeshRegularLeftSurface,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMISSGeneratedMeshRegularRightSurface,RightSurfaceNodes,RightNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMISSGeneratedMeshRegularFrontSurface,FrontSurfaceNodes,FrontNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMISSGeneratedMeshRegularBackSurface,BackSurfaceNodes,BackNormalXi,Err)
  !Get pressure nodes
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMISSGeneratedMeshRegularBottomSurface,PressureBottomSurfaceNodes,BottomNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMISSGeneratedMeshRegularTopSurface,PressureTopSurfaceNodes,TopNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMISSGeneratedMeshRegularLeftSurface,PressureLeftSurfaceNodes,LeftNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMISSGeneratedMeshRegularRightSurface,PressureRightSurfaceNodes,RightNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMISSGeneratedMeshRegularFrontSurface,PressureFrontSurfaceNodes,FrontNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMISSGeneratedMeshRegularBackSurface,PressureBackSurfaceNodes,BackNormalXi,Err)

  !Set boundary conditions on geometry
  !Set x=0 nodes to no x displacment in x
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,GeometricMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NodeNumber,1, &
        & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
    ENDIF
  ENDDO
  !Set y=0 nodes to no y displacement
  DO node_idx=1,SIZE(FrontSurfaceNodes,1)
    NodeNumber=FrontSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,GeometricMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NodeNumber,2, &
        & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
    ENDIF
  ENDDO
  !Set z=0 nodes to no z displacement
  DO node_idx=1,SIZE(BottomSurfaceNodes,1)
    NodeNumber=BottomSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,GeometricMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NodeNumber,3, &
        & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
    ENDIF
  ENDDO
  !Fix right surface nodes
  IF (FixRightSide) THEN
    DO node_idx=1,SIZE(RightSurfaceNodes,1)
      NodeNumber=RightSurfaceNodes(node_idx)
      CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,GeometricMeshComponent,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NodeNumber,1, &
          & CMISSBoundaryConditionFixedIncremented,FixedWidth,Err)
      ENDIF
    ENDDO
  ENDIF

  !Set boundary conditions on fluid pressure
  DO node_idx=1,SIZE(PressureLeftSurfaceNodes,1)
    NodeNumber=PressureLeftSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,PressureMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldVVariableType,1,1,NodeNumber,1, &
        & CMISSBoundaryConditionFixedIncremented,FluidPressureBC,Err)
    ENDIF
  ENDDO
  DO node_idx=1,SIZE(PressureRightSurfaceNodes,1)
    NodeNumber=PressureRightSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,PressureMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldVVariableType,1,1,NodeNumber,1, &
        & CMISSBoundaryConditionFixed,FluidPressureBC2,Err)
    ENDIF
  ENDDO

  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblemSolve(Problem,Err)

  !Output solution
  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"StaticExtension","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"StaticExtension","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP
  END SUBROUTINE HANDLE_ERROR

END PROGRAM POROELASTICITYEXAMPLE
