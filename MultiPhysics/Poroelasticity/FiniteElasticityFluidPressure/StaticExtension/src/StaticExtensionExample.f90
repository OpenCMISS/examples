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

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters

  REAL(CMISSRP), PARAMETER :: Height=10.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: Width=10.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: Length=10.0_CMISSRP
  INTEGER(CMISSIntg), PARAMETER :: GeometricInterpolationType=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: PressureInterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussXi=3
  INTEGER(CMISSIntg) :: NumberElements
  INTEGER(CMISSIntg) :: NumberGlobalXElements=2
  INTEGER(CMISSIntg) :: NumberGlobalYElements=2
  INTEGER(CMISSIntg) :: NumberGlobalZElements=2
  REAL(CMISSRP) :: FixedWidth=12.0_CMISSRP
  REAL(CMISSRP) :: FluidPressureBC=0.5E-3_CMISSRP
  REAL(CMISSRP) :: FluidPressureBC2=0.0_CMISSRP
  REAL(CMISSRP) :: InitialPressure
  INTEGER(CMISSIntg) :: Increments=10
  LOGICAL :: FixRightSide=.TRUE.
  REAL (CMISSRP) :: K1=2.0E-3_CMISSRP
  REAL (CMISSRP) :: K2=33.0E-6_CMISSRP
  REAL (CMISSRP) :: K=2.20E-1_CMISSRP
  REAL (CMISSRP) :: M=2.18E-1_CMISSRP
  REAL (CMISSRP) :: b=1.0_CMISSRP
  REAL (CMISSRP) :: p_0=0.0_CMISSRP
  REAL (CMISSRP) :: permeability=1.0E-3_CMISSRP

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

  TYPE(cmfe_BasisType) :: GeometricBasis,PressureBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: SolidEquations,FluidEquations
  TYPE(cmfe_EquationsSetType) :: SolidEquationsSet,FluidEquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,SolidEquationsSetField,FluidEquationsSetField,DependentField, &
    & FibreField,MaterialField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations

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
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",[ &
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
  IF(FixedWidth<1.0E-5_CMISSRP) FixRightSide=.FALSE.
  WRITE(*,*) "Fluid pressure BC:", FluidPressureBC
  IF (FixRightSide) THEN
    WRITE(*,*) "Fixed width:", FixedWidth
  ELSE
    WRITE(*,*) "Width not fixed"
  ENDIF
  WRITE(*,*) "Increments:", Increments
  WRITE(*,*) "Elements:", NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  InitialPressure=(FluidPressureBC+FluidPressureBC2)*0.5_CMISSRP/REAL(Increments)

  WRITE(Filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "StaticExtension",NumberGlobalXElements,NumberGlobalYElements, &
    & NumberGlobalZElements,GeometricInterpolationType
  CALL cmfe_OutputSetOn(Filename,Err)

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !Create a 3D rectangular cartesian coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Define basis function and second basis function for pressure etc
  CALL cmfe_Basis_Initialise(GeometricBasis,Err)
  CALL cmfe_Basis_Initialise(PressureBasis,Err)
  CALL cmfe_Basis_CreateStart(GeometricBasisUserNumber,GeometricBasis,Err)
  CALL cmfe_Basis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
  SELECT CASE(GeometricInterpolationType)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(GeometricBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(GeometricBasis,CMFE_BASIS_SIMPLEX_TYPE,Err)
    CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_SIMPLEX_TYPE,Err)
  END SELECT
  IF(NumberGlobalZElements==0) THEN
    NumberOfDimensions=2
    CALL cmfe_Basis_NumberOfXiSet(GeometricBasis,2,Err)
    CALL cmfe_Basis_NumberOfXiSet(PressureBasis,2,Err)
    CALL cmfe_Basis_InterpolationXiSet(GeometricBasis,[GeometricInterpolationType,GeometricInterpolationType],Err)
    CALL cmfe_Basis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      !Finite elasticity doesn't support different numbers of gauss points for different components
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(GeometricBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ELSE
    NumberOfDimensions=3
    CALL cmfe_Basis_NumberOfXiSet(GeometricBasis,3,Err)
    CALL cmfe_Basis_NumberOfXiSet(PressureBasis,3,Err)
    CALL cmfe_Basis_InterpolationXiSet(GeometricBasis,[GeometricInterpolationType,GeometricInterpolationType, &
        & GeometricInterpolationType],Err)
    CALL cmfe_Basis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType, &
        & PressureInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(GeometricBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ENDIF
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(GeometricBasis,.TRUE.,Err)
  CALL cmfe_Basis_QuadratureLocalFaceGaussEvaluateSet(PressureBasis,.TRUE.,Err)
  CALL cmfe_Basis_CreateFinish(GeometricBasis,Err)
  CALL cmfe_Basis_CreateFinish(PressureBasis,Err)

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the bases
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,[PressureBasis,GeometricBasis],Err)
  !Define the mesh on the region
  IF(NumberGlobalXElements==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[Width,Height],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[Width,Height,Length],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements, &
      & NumberGlobalZElements],Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_MeshComponentSet(Decomposition,GeometricMeshComponent,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  CALL cmfe_Decomposition_CalculateFacesSet(Decomposition,.TRUE.,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  DO component_idx=1,NumberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
  END DO
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FibreFieldUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  DO component_idx=1,NumberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
  END DO
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !Create the equations sets
  CALL cmfe_EquationsSet_Initialise(SolidEquationsSet,Err)
  CALL cmfe_Field_Initialise(SolidEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(SolidEquationsSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_ELASTICITY_FLUID_PRES_STATIC_INRIA_SUBTYPE], &
    & SolidEquationsSetFieldUserNumber,SolidEquationsSetField,SolidEquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(SolidEquationsSet,Err)

  CALL cmfe_EquationsSet_Initialise(FluidEquationsSet,Err)
  CALL cmfe_Field_Initialise(FluidEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(FluidEquationsSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_FLUID_MECHANICS_CLASS, &
    & CMFE_EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE,CMFE_EQUATIONS_SET_ELASTICITY_FLUID_PRES_STATIC_INRIA_SUBTYPE], &
    & FluidEquationsSetFieldUserNumber,FluidEquationsSetField,FluidEquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(FluidEquationsSet,Err)

  !Create the dependent field
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(FluidEquationsSet,FieldDependentUserNumber,DependentField,Err)
  DO component_idx=1,NumberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
      & GeometricMeshComponent,Err)
  ENDDO
  DO component_idx=1,1
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_V_VARIABLE_TYPE,component_idx,PressureMeshComponent,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELVDELN_VARIABLE_TYPE,component_idx, &
      & PressureMeshComponent,Err)
  ENDDO
  CALL cmfe_EquationsSet_DependentCreateFinish(FluidEquationsSet,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(SolidEquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(SolidEquationsSet,Err)

  !Create the material field
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(FluidEquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,"SolidMaterial",Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_V_VARIABLE_TYPE,"SolidDensity",Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_U1_VARIABLE_TYPE,"FluidMaterial",Err)
  DO component_idx=1,NumberOfDimensions
    CALL cmfe_Field_ComponentMeshComponentSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,component_idx,GeometricMeshComponent,Err)
  END DO
  CALL cmfe_EquationsSet_MaterialsCreateFinish(FluidEquationsSet,Err)

  CALL cmfe_EquationsSet_MaterialsCreateStart(SolidEquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(SolidEquationsSet,Err)

  !Set material constants
  !Solid constitutive law
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,K1,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,K2,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,K,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,M,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,b,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6,p_0,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP,Err)

  !Permeability tensor
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,permeability, &
    & Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,0.0_CMISSRP, &
    & Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,0.0_CMISSRP, &
    & Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4,permeability, &
    & Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,5,0.0_CMISSRP, &
    & Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,6,permeability, &
    & Err)
  !Density
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U1_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,7, &
    & 1000.0_CMISSRP,Err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(FluidEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(FluidEquationsSet,FluidEquations,Err)
  CALL cmfe_Equations_SparsityTypeSet(FluidEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(FluidEquations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(FluidEquationsSet,Err)

  CALL cmfe_Equations_Initialise(SolidEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(SolidEquationsSet,SolidEquations,Err)
  CALL cmfe_Equations_SparsityTypeSet(SolidEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(SolidEquations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(SolidEquationsSet,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
    & InitialPressure,Err)

  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
    & CMFE_PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE,CMFE_PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,[CMFE_CONTROL_LOOP_NODE],ControlLoop,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoop,Increments,Err)
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoop,CMFE_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(Solver,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
  CALL cmfe_Solver_NewtonMaximumIterationsSet(Solver,200,Err)
  CALL cmfe_Solver_NewtonMaximumFunctionEvaluationsSet(Solver,10000,Err)
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(Solver,1.0E-10_CMISSRP,Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(Solver,1.0E-9_CMISSRP,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,SolidEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,FluidEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Set the boundary conditions for the solver equations
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  !Get geometric nodes
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMFE_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMFE_GENERATED_MESH_REGULAR_TOP_SURFACE,TopSurfaceNodes,TopNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,GeometricMeshComponent, &
      & CMFE_GENERATED_MESH_REGULAR_BACK_SURFACE,BackSurfaceNodes,BackNormalXi,Err)
  !Get pressure nodes
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMFE_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,PressureBottomSurfaceNodes,BottomNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMFE_GENERATED_MESH_REGULAR_TOP_SURFACE,PressureTopSurfaceNodes,TopNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,PressureLeftSurfaceNodes,LeftNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,PressureRightSurfaceNodes,RightNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE,PressureFrontSurfaceNodes,FrontNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,PressureMeshComponent, &
      & CMFE_GENERATED_MESH_REGULAR_BACK_SURFACE,PressureBackSurfaceNodes,BackNormalXi,Err)

  !Set boundary conditions on geometry
  !Set x=0 nodes to no x displacment in x
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,GeometricMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO
  !Set y=0 nodes to no y displacement
  DO node_idx=1,SIZE(FrontSurfaceNodes,1)
    NodeNumber=FrontSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,GeometricMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO
  !Set z=0 nodes to no z displacement
  DO node_idx=1,SIZE(BottomSurfaceNodes,1)
    NodeNumber=BottomSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,GeometricMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    ENDIF
  ENDDO
  !Fix right surface nodes
  IF (FixRightSide) THEN
    DO node_idx=1,SIZE(RightSurfaceNodes,1)
      NodeNumber=RightSurfaceNodes(node_idx)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,GeometricMeshComponent,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
          & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,FixedWidth,Err)
      ENDIF
    ENDDO
  ENDIF

  !Set boundary conditions on fluid pressure
  DO node_idx=1,SIZE(PressureLeftSurfaceNodes,1)
    NodeNumber=PressureLeftSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,PressureMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_V_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,FluidPressureBC,Err)
    ENDIF
  ENDDO
  DO node_idx=1,SIZE(PressureRightSurfaceNodes,1)
    NodeNumber=PressureRightSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,PressureMeshComponent,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_V_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED,FluidPressureBC2,Err)
    ENDIF
  ENDDO

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !Output solution
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"StaticExtension","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"StaticExtension","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP
  END SUBROUTINE HANDLE_ERROR

END PROGRAM POROELASTICITYEXAMPLE
