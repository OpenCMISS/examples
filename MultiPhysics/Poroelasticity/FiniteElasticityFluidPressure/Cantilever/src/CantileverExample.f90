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
  INTEGER(CMISSIntg), PARAMETER :: GeometricInterpolationType=CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: PressureInterpolationType=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussXi=3
  REAL(CMISSDP), PARAMETER :: Density=1.0E3_CMISSDP
  REAL(CMISSDP), PARAMETER :: FluidDensity=1.0E3_CMISSDP
  REAL(CMISSDP), PARAMETER :: Gravity(3)=[0.0_CMISSDP,0.0_CMISSDP,-9.81_CMISSDP]
  INTEGER(CMISSIntg) :: Increments=4
  REAL(CMISSDP) :: FluidPressureBC=500.0_CMISSDP
  REAL (CMISSDP) :: C0=2.5E3_CMISSDP
  REAL (CMISSDP) :: C1=2.0_CMISSDP
  REAL (CMISSDP) :: C2=0.1_CMISSDP
  REAL (CMISSDP) :: permeability=1.0E-3_CMISSDP
  REAL (CMISSDP) :: porosity=0.2_CMISSDP

  INTEGER(CMISSIntg) :: NumberGlobalXElements=3
  INTEGER(CMISSIntg) :: NumberGlobalYElements=2
  INTEGER(CMISSIntg) :: NumberGlobalZElements=2

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricMeshComponent=2
  INTEGER(CMISSIntg), PARAMETER :: PressureMeshComponent=1
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

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains=NumberOfComputationalNodes

  !Create a 3D rectangular cartesian coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Define basis function for displacement
  CALL CMISSBasis_Initialise(GeometricBasis,Err)
  CALL CMISSBasis_CreateStart(GeometricBasisUserNumber,GeometricBasis,Err)
  SELECT CASE(GeometricInterpolationType)
  CASE(1,2,3,4)
    CALL CMISSBasis_TypeSet(GeometricBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL CMISSBasis_TypeSet(GeometricBasis,CMISS_BASIS_SIMPLEX_TYPE,Err)
  END SELECT
  CALL CMISSBasis_NumberOfXiSet(GeometricBasis,3,Err)
  CALL CMISSBasis_InterpolationXiSet(GeometricBasis,[GeometricInterpolationType,GeometricInterpolationType, &
      & GeometricInterpolationType],Err)
  IF(NumberOfGaussXi>0) THEN
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(GeometricBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  ENDIF
  CALL CMISSBasis_CreateFinish(GeometricBasis,Err)

  !Basis for pressure
  CALL CMISSBasis_Initialise(PressureBasis,Err)
  CALL CMISSBasis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
  SELECT CASE(PressureInterpolationType)
  CASE(1,2,3,4)
    CALL CMISSBasis_TypeSet(PressureBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL CMISSBasis_TypeSet(PressureBasis,CMISS_BASIS_SIMPLEX_TYPE,Err)
  END SELECT
  CALL CMISSBasis_NumberOfXiSet(PressureBasis,3,Err)
  CALL CMISSBasis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType, &
      & PressureInterpolationType],Err)
  IF(NumberOfGaussXi>0) THEN
    CALL CMISSBasis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  ENDIF
  CALL CMISSBasis_CreateFinish(PressureBasis,Err)

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,[PressureBasis,GeometricBasis],Err)
  !Define the mesh on the region
  IF(NumberGlobalXElements==0) THEN
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[Width,Height],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[Width,Length,Height],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
      & NumberGlobalZElements],Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL CMISSMesh_Initialise(Mesh,Err)
  CALL CMISSGeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL CMISSDecomposition_Initialise(Decomposition,Err)
  CALL CMISSDecomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL CMISSDecomposition_TypeSet(Decomposition,CMISS_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL CMISSDecomposition_MeshComponentSet(Decomposition,GeometricMeshComponent,Err)
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSField_VariableLabelSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL CMISSField_NumberOfComponentsSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,3,Err)
  CALL CMISSField_ScalingTypeSet(GeometricField,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
  DO component_idx=1,3
    CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
  ENDDO
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Create a fibre field and attach it to the geometric field
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSField_VariableLabelSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  DO component_idx=1,3
    CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
  ENDDO
  CALL CMISSField_CreateFinish(FibreField,Err)

  !Create the equations sets
  CALL CMISSEquationsSet_Initialise(SolidEquationsSet,Err)
  CALL CMISSField_Initialise(SolidEquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(SolidEquationsSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_SUBTYPE, &
    & SolidEquationsSetFieldUserNumber,SolidEquationsSetField,SolidEquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(SolidEquationsSet,Err)

  CALL CMISSEquationsSet_Initialise(FluidEquationsSet,Err)
  CALL CMISSField_Initialise(FluidEquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(FluidEquationsSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & CMISS_EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE,CMISS_EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_SUBTYPE, &
    & FluidEquationsSetFieldUserNumber,FluidEquationsSetField,FluidEquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(FluidEquationsSet,Err)

  !Create the dependent field
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(SolidEquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSField_VariableLabelSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  DO component_idx=1,3
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
      & GeometricMeshComponent,Err)
  ENDDO
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_V_VARIABLE_TYPE,1,PressureMeshComponent,Err)
  CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,1,PressureMeshComponent,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(SolidEquationsSet,Err)
  CALL CMISSEquationsSet_DependentCreateStart(FluidEquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(FluidEquationsSet,Err)

  !Create the material field
  CALL CMISSField_Initialise(MaterialField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(SolidEquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL CMISSField_VariableLabelSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,"SolidMaterial",Err)
  CALL CMISSField_VariableLabelSet(MaterialField,CMISS_FIELD_V_VARIABLE_TYPE,"Density",Err)
  CALL CMISSField_VariableLabelSet(MaterialField,CMISS_FIELD_U1_VARIABLE_TYPE,"FluidMaterial",Err)
  DO component_idx=1,4
    CALL CMISSField_ComponentMeshComponentSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
  ENDDO
  DO component_idx=1,1
    CALL CMISSField_ComponentMeshComponentSet(MaterialField,CMISS_FIELD_V_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
  ENDDO
  DO component_idx=1,8
    CALL CMISSField_ComponentMeshComponentSet(MaterialField,CMISS_FIELD_U1_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
  ENDDO
  CALL CMISSEquationsSet_MaterialsCreateFinish(SolidEquationsSet,Err)

  CALL CMISSEquationsSet_MaterialsCreateStart(FluidEquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(FluidEquationsSet,Err)

  !Set material constants
  !Solid constitutive law
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,C0,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,C1,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,C2,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4, &
    & (1.0_CMISSDP-porosity),Err)
  !Solid density
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Density,Err)

  !Permeability tensor
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,permeability, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,0.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,0.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,permeability, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,5,0.0_CMISSDP, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,6,permeability, &
    & Err)
  !Fluid Density
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,7,FluidDensity, &
    & Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,8,porosity,Err)

  !Create the source field with the gravity vector
  CALL CMISSField_Initialise(SourceField,Err)
  CALL CMISSEquationsSet_SourceCreateStart(SolidEquationsSet,FieldSourceUserNumber,SourceField,Err)
  DO component_idx=1,3
    CALL CMISSField_ComponentMeshComponentSet(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
  ENDDO
  CALL CMISSEquationsSet_SourceCreateFinish(SolidEquationsSet,Err)
  DO component_idx=1,3
    CALL CMISSField_ComponentValuesInitialise(SourceField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
        & component_idx,Gravity(component_idx),Err)
  ENDDO

  !Create the equations set equations
  CALL CMISSEquations_Initialise(FluidEquations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(FluidEquationsSet,FluidEquations,Err)
  CALL CMISSEquations_SparsityTypeSet(FluidEquations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(FluidEquations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(FluidEquationsSet,Err)

  CALL CMISSEquations_Initialise(SolidEquations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(SolidEquationsSet,SolidEquations,Err)
  CALL CMISSEquations_SparsityTypeSet(SolidEquations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(SolidEquations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(SolidEquationsSet,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1, &
    & FluidPressureBC,Err)

  !Define the problem
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_MULTI_PHYSICS_CLASS, &
    & CMISS_PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE, &
    & CMISS_PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE,Err)
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopGet(Problem,[CMISS_CONTROL_LOOP_NODE],ControlLoop,Err)
  CALL CMISSControlLoop_MaximumIterationsSet(ControlLoop,Increments,Err)
  CALL CMISSControlLoop_OutputTypeSet(ControlLoop,CMISS_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(Solver,CMISS_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(Solver,1.0E-6_CMISSDP,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(Solver,1.0E-7_CMISSDP,Err)
  CALL CMISSSolver_NewtonMaximumIterationsSet(Solver,200,Err)
  CALL CMISSSolver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquations_SparsityTypeSet(SolverEquations,CMISS_SOLVER_SPARSE_MATRICES,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,SolidEquationsSet,EquationsSetIndex,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,FluidEquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,GeometricMeshComponent,CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE, &
      & LeftSurfaceNodes,LeftNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,PressureMeshComponent,CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE, &
      & PressureLeftSurfaceNodes,LeftNormalXi,Err)

  !Fix x=0 nodes
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,GeometricMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      DO component_idx=1,3
        CALL CMISSField_ParameterSetGetNode(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,1,NodeNumber, &
          & component_idx,Value,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
          & component_idx,CMISS_BOUNDARY_CONDITION_FIXED,Value,Err)
      ENDDO
    ENDIF
  ENDDO
  DO node_idx=1,SIZE(PressureLeftSurfaceNodes,1)
    NodeNumber=PressureLeftSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,PressureMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_V_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED,FluidPressureBC,Err)
    ENDIF
  ENDDO

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblem_Solve(Problem,Err)

  !Output solution
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"CoupledCantilever","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"CoupledCantilever","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

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

