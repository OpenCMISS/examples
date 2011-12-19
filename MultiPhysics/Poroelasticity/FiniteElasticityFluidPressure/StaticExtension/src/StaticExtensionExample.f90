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
  INTEGER(CMISSIntg), PARAMETER :: GeometricInterpolationType=CMISS_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: PressureInterpolationType=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
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
  INTEGER(CMISSIntg), PARAMETER :: GeometricMeshComponent=2
  INTEGER(CMISSIntg), PARAMETER :: PressureMeshComponent=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FibreFieldUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=4
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
    & FibreField,MaterialField
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
  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)
  CALL CMISSDiagnosticsSetOn(CMISS_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",[ &
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
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Define basis function and second basis function for pressure etc
  CALL CMISSBasis_Initialise(GeometricBasis,Err)
  CALL CMISSBasis_Initialise(PressureBasis,Err)
  CALL CMISSBasis_CreateStart(GeometricBasisUserNumber,GeometricBasis,Err)
  CALL CMISSBasis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
  SELECT CASE(GeometricInterpolationType)
  CASE(1,2,3,4)
    CALL CMISSBasis_TypeSet(GeometricBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CALL CMISSBasis_TypeSet(PressureBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL CMISSBasis_TypeSet(GeometricBasis,CMISS_BASIS_SIMPLEX_TYPE,Err)
    CALL CMISSBasis_TypeSet(PressureBasis,CMISS_BASIS_SIMPLEX_TYPE,Err)
  END SELECT
  IF(NumberGlobalZElements==0) THEN
    NumberOfDimensions=2
    CALL CMISSBasis_NumberOfXiSet(GeometricBasis,2,Err)
    CALL CMISSBasis_NumberOfXiSet(PressureBasis,2,Err)
    CALL CMISSBasis_InterpolationXiSet(GeometricBasis,[GeometricInterpolationType,GeometricInterpolationType],Err)
    CALL CMISSBasis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      !Finite elasticity doesn't support different numbers of gauss points for different components
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(GeometricBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ELSE
    NumberOfDimensions=3
    CALL CMISSBasis_NumberOfXiSet(GeometricBasis,3,Err)
    CALL CMISSBasis_NumberOfXiSet(PressureBasis,3,Err)
    CALL CMISSBasis_InterpolationXiSet(GeometricBasis,[GeometricInterpolationType,GeometricInterpolationType, &
        & GeometricInterpolationType],Err)
    CALL CMISSBasis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType, &
        & PressureInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(GeometricBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ENDIF
  CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(GeometricBasis,.TRUE.,Err)
  CALL CMISSBasis_QuadratureLocalFaceGaussEvaluateSet(PressureBasis,.TRUE.,Err)
  CALL CMISSBasis_CreateFinish(GeometricBasis,Err)
  CALL CMISSBasis_CreateFinish(PressureBasis,Err)

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the bases
  CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,[PressureBasis,GeometricBasis],Err)
  !Define the mesh on the region
  IF(NumberGlobalXElements==0) THEN
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[Width,Height],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[Width,Height,Length],Err)
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
  CALL CMISSDecomposition_MeshComponentSet(Decomposition,GeometricMeshComponent,Err)
  CALL CMISSDecomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  CALL CMISSDecomposition_CalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  DO component_idx=1,NumberOfDimensions
    CALL CMISSField_ComponentMeshComponentSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
  END DO
  CALL CMISSField_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL CMISSGeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Create a fibre field and attach it to the geometric field
  CALL CMISSField_Initialise(FibreField,Err)
  CALL CMISSField_CreateStart(FibreFieldUserNumber,Region,FibreField,Err)
  CALL CMISSField_TypeSet(FibreField,CMISS_FIELD_FIBRE_TYPE,Err)
  CALL CMISSField_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL CMISSField_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL CMISSField_VariableLabelSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  DO component_idx=1,NumberOfDimensions
    CALL CMISSField_ComponentMeshComponentSet(FibreField,CMISS_FIELD_U_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
  END DO
  CALL CMISSField_CreateFinish(FibreField,Err)

  !Create the equations sets
  CALL CMISSEquationsSet_Initialise(SolidEquationsSet,Err)
  CALL CMISSField_Initialise(SolidEquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(SolidEquationsSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_ELASTICITY_FLUID_PRES_STATIC_INRIA_SUBTYPE, &
    & SolidEquationsSetFieldUserNumber,SolidEquationsSetField,SolidEquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(SolidEquationsSet,Err)

  CALL CMISSEquationsSet_Initialise(FluidEquationsSet,Err)
  CALL CMISSField_Initialise(FluidEquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(FluidEquationsSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & CMISS_EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE,CMISS_EQUATIONS_SET_ELASTICITY_FLUID_PRES_STATIC_INRIA_SUBTYPE, &
    & FluidEquationsSetFieldUserNumber,FluidEquationsSetField,FluidEquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(FluidEquationsSet,Err)

  !Create the dependent field
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(FluidEquationsSet,FieldDependentUserNumber,DependentField,Err)
  DO component_idx=1,NumberOfDimensions
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
      & GeometricMeshComponent,Err)
  ENDDO
  DO component_idx=1,1
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_V_VARIABLE_TYPE,component_idx,PressureMeshComponent,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELVDELN_VARIABLE_TYPE,component_idx, &
      & PressureMeshComponent,Err)
  ENDDO
  CALL CMISSEquationsSet_DependentCreateFinish(FluidEquationsSet,Err)
  CALL CMISSEquationsSet_DependentCreateStart(SolidEquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateFinish(SolidEquationsSet,Err)

  !Create the material field
  CALL CMISSField_Initialise(MaterialField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(FluidEquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL CMISSField_VariableLabelSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,"SolidMaterial",Err)
  CALL CMISSField_VariableLabelSet(MaterialField,CMISS_FIELD_V_VARIABLE_TYPE,"SolidDensity",Err)
  CALL CMISSField_VariableLabelSet(MaterialField,CMISS_FIELD_U1_VARIABLE_TYPE,"FluidMaterial",Err)
  DO component_idx=1,NumberOfDimensions
    CALL CMISSField_ComponentMeshComponentSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
  END DO
  CALL CMISSEquationsSet_MaterialsCreateFinish(FluidEquationsSet,Err)

  CALL CMISSEquationsSet_MaterialsCreateStart(SolidEquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(SolidEquationsSet,Err)

  !Set material constants
  !Solid constitutive law
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,K1,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,K2,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,K,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,M,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,5,b,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,6,p_0,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_V_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)

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
  !Density
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U1_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,7, &
    & 1000.0_CMISSDP,Err)

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
    & InitialPressure,Err)

  !Define the problem
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_MULTI_PHYSICS_CLASS, &
    & CMISS_PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE, &
    & CMISS_PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE,Err)
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSControlLoop_Initialise(ControlLoop,Err)
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
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
  CALL CMISSSolver_NewtonMaximumIterationsSet(Solver,200,Err)
  CALL CMISSSolver_NewtonMaximumFunctionEvaluationsSet(Solver,10000,Err)
  CALL CMISSSolver_NewtonAbsoluteToleranceSet(Solver,1.0E-10_CMISSDP,Err)
  CALL CMISSSolver_NewtonRelativeToleranceSet(Solver,1.0E-9_CMISSDP,Err)
  CALL CMISSSolver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
  CALL CMISSSolver_LinearIterativeMaximumIterationsSet(LinearSolver,1000,Err)
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

  !Set the boundary conditions for the solver equations
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Get geometric nodes
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMISS_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMISS_GENERATED_MESH_REGULAR_TOP_SURFACE,TopSurfaceNodes,TopNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMISS_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMISS_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMISS_GENERATED_MESH_REGULAR_BACK_SURFACE,BackSurfaceNodes,BackNormalXi,Err)
  !Get pressure nodes
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMISS_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,PressureBottomSurfaceNodes,BottomNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMISS_GENERATED_MESH_REGULAR_TOP_SURFACE,PressureTopSurfaceNodes,TopNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE,PressureLeftSurfaceNodes,LeftNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMISS_GENERATED_MESH_REGULAR_RIGHT_SURFACE,PressureRightSurfaceNodes,RightNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMISS_GENERATED_MESH_REGULAR_FRONT_SURFACE,PressureFrontSurfaceNodes,FrontNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMISS_GENERATED_MESH_REGULAR_BACK_SURFACE,PressureBackSurfaceNodes,BackNormalXi,Err)

  !Set boundary conditions on geometry
  !Set x=0 nodes to no x displacment in x
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,GeometricMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO
  !Set y=0 nodes to no y displacement
  DO node_idx=1,SIZE(FrontSurfaceNodes,1)
    NodeNumber=FrontSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,GeometricMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO
  !Set z=0 nodes to no z displacement
  DO node_idx=1,SIZE(BottomSurfaceNodes,1)
    NodeNumber=BottomSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,GeometricMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO
  !Fix right surface nodes
  IF (FixRightSide) THEN
    DO node_idx=1,SIZE(RightSurfaceNodes,1)
      NodeNumber=RightSurfaceNodes(node_idx)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,GeometricMeshComponent,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
          & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED,FixedWidth,Err)
      ENDIF
    ENDDO
  ENDIF

  !Set boundary conditions on fluid pressure
  DO node_idx=1,SIZE(PressureLeftSurfaceNodes,1)
    NodeNumber=PressureLeftSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,PressureMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_V_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED_INCREMENTED,FluidPressureBC,Err)
    ENDIF
  ENDDO
  DO node_idx=1,SIZE(PressureRightSurfaceNodes,1)
    NodeNumber=PressureRightSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,PressureMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_V_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED,FluidPressureBC2,Err)
    ENDIF
  ENDDO

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblem_Solve(Problem,Err)

  !Output solution
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"StaticExtension","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"StaticExtension","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

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
