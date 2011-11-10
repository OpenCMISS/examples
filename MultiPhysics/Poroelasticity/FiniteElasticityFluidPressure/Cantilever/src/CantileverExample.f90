!> \file
!> \author Adam Reeve
!> \brief This is an example program to solve a coupled finite elasticity and fluid pressure problem with gravity loading using OpenCMISS calls.
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

!> \example MultiPhysics/Poroelasticity/FiniteElasticityFluidPressure/Cantilever/src/CantileverExample.f90
!! Example program to solve a coupled finite elasticity and fluid pressure problem with gravity loading using OpenCMISS calls.
!<

!> Main program
PROGRAM COUPLEDCANTILEVER

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  !Units system
  !length: m
  !force: N
  !mass: kg

  REAL(CMISSDP), PARAMETER :: Height=0.045_CMISSDP
  REAL(CMISSDP), PARAMETER :: Width=0.065_CMISSDP
  REAL(CMISSDP), PARAMETER :: Length=0.045_CMISSDP
  INTEGER(CMISSIntg), PARAMETER :: GeometricInterpolationType=CMISSBasisQuadraticLagrangeInterpolation
  INTEGER(CMISSIntg), PARAMETER :: PressureInterpolationType=CMISSBasisLinearLagrangeInterpolation
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussXi=3
  REAL(CMISSDP), PARAMETER :: Density=1.0E3_CMISSDP
  REAL(CMISSDP), PARAMETER :: FluidDensity=1.0E3_CMISSDP
  REAL(CMISSDP), PARAMETER :: Gravity(3)=[0.0_CMISSDP,0.0_CMISSDP,-9.81_CMISSDP]
  INTEGER(CMISSIntg) :: Increments=5
  REAL(CMISSDP) :: FluidPressureBC=30.0_CMISSDP
  REAL (CMISSDP) :: K1=3.0E3_CMISSDP
  REAL (CMISSDP) :: K2=50.0_CMISSDP
  REAL (CMISSDP) :: K=3.20E4_CMISSDP
  REAL (CMISSDP) :: M=3.18E4_CMISSDP
  REAL (CMISSDP) :: b=1.0_CMISSDP
  REAL (CMISSDP) :: p_0=0.0_CMISSDP
  REAL (CMISSDP) :: permeability=1.0E-3_CMISSDP

  INTEGER(CMISSIntg) :: NumberGlobalXElements=3
  INTEGER(CMISSIntg) :: NumberGlobalYElements=2
  INTEGER(CMISSIntg) :: NumberGlobalZElements=2

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricMeshComponent=1
  INTEGER(CMISSIntg), PARAMETER :: PressureMeshComponent=2
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: FieldSourceUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: SolidEquationsSetFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: FluidEquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: FluidEquationsSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: SolidEquationsSetUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program variables

  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx,component_idx
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:),PressureLeftSurfaceNodes(:)
  INTEGER(CMISSIntg) :: LeftNormalXi
  REAL(CMISSDP) :: Value
  INTEGER(CMISSIntg) :: NumberOfArguments,ArgumentLength,ArgStatus
  CHARACTER(LEN=255) :: CommandArgument

  !CMISS variables
  TYPE(CMISSBasisType) :: GeometricBasis,PressureBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: SolidEquations,FluidEquations
  TYPE(CMISSEquationsSetType) :: SolidEquationsSet,FluidEquationsSet
  TYPE(CMISSFieldType) :: GeometricField,FibreField,MaterialField,DependentField, &
      & SourceField,SolidEquationsSetField,FluidEquationsSetField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations
  TYPE(CMISSControlLoopType) :: ControlLoop

  !Generic CMISS variables
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

  !Intialise cmiss
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  !Read in arguments and overwrite default values
  NumberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(NumberOfArguments >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(CommandArgument(1:ArgumentLength),*) FluidPressureBC
  ENDIF

  CALL CMISSErrorHandlingModeSet(CMISSTrapError,Err)

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains=NumberOfComputationalNodes

  !Create a 3D rectangular cartesian coordinate system
  CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL CMISSRegionTypeInitialise(Region,Err)
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegionCreateFinish(Region,Err)

  !Define basis function for displacement
  CALL CMISSBasisTypeInitialise(GeometricBasis,Err)
  CALL CMISSBasisCreateStart(GeometricBasisUserNumber,GeometricBasis,Err)
  SELECT CASE(GeometricInterpolationType)
  CASE(1,2,3,4)
    CALL CMISSBasisTypeSet(GeometricBasis,CMISSBasisLagrangeHermiteTPType,Err)
  CASE(7,8,9)
    CALL CMISSBasisTypeSet(GeometricBasis,CMISSBasisSimplexType,Err)
  END SELECT
  CALL CMISSBasisNumberOfXiSet(GeometricBasis,3,Err)
  CALL CMISSBasisInterpolationXiSet(GeometricBasis,[GeometricInterpolationType,GeometricInterpolationType, &
      & GeometricInterpolationType],Err)
  IF(NumberOfGaussXi>0) THEN
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(GeometricBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  ENDIF
  CALL CMISSBasisCreateFinish(GeometricBasis,Err)

  !Basis for pressure
  CALL CMISSBasisTypeInitialise(PressureBasis,Err)
  CALL CMISSBasisCreateStart(PressureBasisUserNumber,PressureBasis,Err)
  SELECT CASE(PressureInterpolationType)
  CASE(1,2,3,4)
    CALL CMISSBasisTypeSet(PressureBasis,CMISSBasisLagrangeHermiteTPType,Err)
  CASE(7,8,9)
    CALL CMISSBasisTypeSet(PressureBasis,CMISSBasisSimplexType,Err)
  END SELECT
  CALL CMISSBasisNumberOfXiSet(PressureBasis,3,Err)
  CALL CMISSBasisInterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType, &
      & PressureInterpolationType],Err)
  IF(NumberOfGaussXi>0) THEN
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  ENDIF
  CALL CMISSBasisCreateFinish(PressureBasis,Err)

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,[GeometricBasis,PressureBasis],Err)
  !Define the mesh on the region
  IF(NumberGlobalXElements==0) THEN
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[Width,Height],Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,[Width,Length,Height],Err)
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
  CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSFieldVariableLabelSet(GeometricField,CMISSFieldUVariableType,"Geometry",Err)
  CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,3,Err)
  CALL CMISSFieldScalingTypeSet(GeometricField,CMISSFieldArithmeticMeanScaling,Err)
  CALL CMISSFieldCreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

  !Create a fibre field and attach it to the geometric field
  CALL CMISSFieldTypeInitialise(FibreField,Err)
  CALL CMISSFieldCreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSFieldTypeSet(FibreField,CMISSFieldFibreType,Err)
  CALL CMISSFieldMeshDecompositionSet(FibreField,Decomposition,Err)
  CALL CMISSFieldGeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSFieldVariableLabelSet(FibreField,CMISSFieldUVariableType,"Fibre",Err)
  CALL CMISSFieldCreateFinish(FibreField,Err)

  !Create the equations sets
  CALL CMISSEquationsSetTypeInitialise(SolidEquationsSet,Err)
  CALL CMISSFieldTypeInitialise(SolidEquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(SolidEquationsSetUserNumber,Region,FibreField,CMISSEquationsSetElasticityClass, &
    & CMISSEquationsSetFiniteElasticityType,CMISSEquationsSetElasticityFluidPressureStaticINRIASubtype, &
    & SolidEquationsSetFieldUserNumber,SolidEquationsSetField,SolidEquationsSet,Err)
  CALL CMISSEquationsSetCreateFinish(SolidEquationsSet,Err)

  CALL CMISSEquationsSetTypeInitialise(FluidEquationsSet,Err)
  CALL CMISSFieldTypeInitialise(FluidEquationsSetField,Err)
  CALL CMISSEquationsSetCreateStart(FluidEquationsSetUserNumber,Region,FibreField,CMISSEquationsSetFluidMechanicsClass, &
    & CMISSEquationsSetDarcyPressureEquationType,CMISSEquationsSetElasticityFluidPressureQuadraticSubtype, &
    & FluidEquationsSetFieldUserNumber,FluidEquationsSetField,FluidEquationsSet,Err)
  CALL CMISSEquationsSetCreateFinish(FluidEquationsSet,Err)

  !Create the dependent field
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(SolidEquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSFieldVariableLabelSet(DependentField,CMISSFieldUVariableType,"Dependent",Err)
  DO component_idx=1,3
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,component_idx,GeometricMeshComponent,Err)
    CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,component_idx,GeometricMeshComponent,Err)
  ENDDO
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldVVariableType,1,PressureMeshComponent,Err)
  CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelVDelNVariableType,1,PressureMeshComponent,Err)
  CALL CMISSEquationsSetDependentCreateFinish(SolidEquationsSet,Err)
  CALL CMISSEquationsSetDependentCreateStart(FluidEquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSetDependentCreateFinish(FluidEquationsSet,Err)

  !Create the material field
  CALL CMISSFieldTypeInitialise(MaterialField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(SolidEquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL CMISSFieldVariableLabelSet(MaterialField,CMISSFieldUVariableType,"SolidMaterial",Err)
  CALL CMISSFieldVariableLabelSet(MaterialField,CMISSFieldVVariableType,"Density",Err)
  CALL CMISSFieldVariableLabelSet(MaterialField,CMISSFieldU1VariableType,"FluidMaterial",Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(SolidEquationsSet,Err)

  CALL CMISSEquationsSetMaterialsCreateStart(FluidEquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL CMISSEquationsSetMaterialsCreateFinish(FluidEquationsSet,Err)

  !Set material constants
  !Solid constitutive law
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,K1,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,2,K2,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,3,K,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,4,M,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,5,b,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,6,p_0,Err)
  !Solid density
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldVVariableType,CMISSFieldValuesSetType,1,Density,Err)

  !Permeability tensor
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldU1VariableType,CMISSFieldValuesSetType,1,permeability,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldU1VariableType,CMISSFieldValuesSetType,2,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldU1VariableType,CMISSFieldValuesSetType,3,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldU1VariableType,CMISSFieldValuesSetType,4,permeability,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldU1VariableType,CMISSFieldValuesSetType,5,0.0_CMISSDP,Err)
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldU1VariableType,CMISSFieldValuesSetType,6,permeability,Err)
  !Fluid Density
  CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldU1VariableType,CMISSFieldValuesSetType,7,FluidDensity,Err)

  !Create the source field with the gravity vector
  CALL CMISSFieldTypeInitialise(SourceField,Err)
  CALL CMISSEquationsSetSourceCreateStart(SolidEquationsSet,FieldSourceUserNumber,SourceField,Err)
  CALL CMISSEquationsSetSourceCreateFinish(SolidEquationsSet,Err)
  DO component_idx=1,3
    CALL CMISSFieldComponentValuesInitialise(SourceField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
        & component_idx,Gravity(component_idx),Err)
  ENDDO

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
  CALL CMISSFieldComponentValuesInitialise(DependentField,CMISSFieldVVariableType,CMISSFieldValuesSetType,1,FluidPressureBC,Err)

  !Define the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemMultiPhysicsClass,CMISSProblemFiniteElasticityFluidPressureType, &
    & CMISSProblemStandardElasticityFluidPressureSubtype,Err)
  CALL CMISSProblemCreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  CALL CMISSControlLoopTypeInitialise(ControlLoop,Err)
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
  CALL CMISSSolverNewtonAbsoluteToleranceSet(Solver,1.0E-6_CMISSDP,Err)
  CALL CMISSSolverNewtonRelativeToleranceSet(Solver,1.0E-7_CMISSDP,Err)
  CALL CMISSSolverNewtonMaximumIterationsSet(Solver,200,Err)
  CALL CMISSSolverNewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearDirectSolveType,Err)
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

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSSolverEquationsBoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,GeometricMeshComponent,CMISSGeneratedMeshRegularLeftSurface, &
      & LeftSurfaceNodes,LeftNormalXi,Err)
  CALL CMISSGeneratedMeshSurfaceGet(GeneratedMesh,PressureMeshComponent,CMISSGeneratedMeshRegularLeftSurface, &
      & PressureLeftSurfaceNodes,LeftNormalXi,Err)

  !Fix x=0 nodes
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      DO component_idx=1,3
        CALL CMISSFieldParameterSetGetNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,1,1,NodeNumber, &
          & component_idx,Value,Err)
        CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldUVariableType,1,1,NodeNumber, &
          & component_idx,CMISSBoundaryConditionFixed,Value,Err)
      ENDDO
    ENDIF
  ENDDO
  DO node_idx=1,SIZE(PressureLeftSurfaceNodes,1)
    NodeNumber=PressureLeftSurfaceNodes(node_idx)
    CALL CMISSDecompositionNodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,DependentField,CMISSFieldVVariableType,1,1,NodeNumber,1, &
        & CMISSBoundaryConditionFixedIncremented,FluidPressureBC,Err)
    ENDIF
  ENDDO

  CALL CMISSSolverEquationsBoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblemSolve(Problem,Err)

  !Output solution
  CALL CMISSFieldsTypeInitialise(Fields,Err)
  CALL CMISSFieldsTypeCreate(Region,Fields,Err)
  CALL CMISSFieldIONodesExport(Fields,"CoupledCantilever","FORTRAN",Err)
  CALL CMISSFieldIOElementsExport(Fields,"CoupledCantilever","FORTRAN",Err)
  CALL CMISSFieldsTypeFinalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program completed."

  STOP

CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP
  END SUBROUTINE HANDLE_ERROR

END PROGRAM COUPLEDCANTILEVER

