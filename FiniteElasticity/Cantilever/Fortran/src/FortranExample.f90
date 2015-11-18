!> \file
!> \author Adam Reeve
!> \brief This is an example program to solve a finite elasticity equation using OpenCMISS calls.
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
!> Contributor(s): Adam Reeve
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

!> \example FiniteElasticity/Cantilever/src/CantileverExample.f90
!! Example program to solve a finite elasticity equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/Cantilever/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/Cantilever/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM CANTILEVEREXAMPLE

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
  REAL(CMISSRP), PARAMETER :: Width=60.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: Length=40.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: Height=40.0_CMISSRP
  INTEGER(CMISSIntg) :: DisplacementInterpolationType
  INTEGER(CMISSIntg) :: PressureInterpolationType
  INTEGER(CMISSIntg) :: PressureMeshComponent
  INTEGER(CMISSIntg) :: NumberOfGaussXi
  INTEGER(CMISSIntg) :: ScalingType
  REAL(CMISSRP), PARAMETER :: Density=9.0E-4_CMISSRP !in g mm^-3
  REAL(CMISSRP), PARAMETER :: Gravity(3)=[0.0_CMISSRP,0.0_CMISSRP,-9.8_CMISSRP] !in m s^-2
  INTEGER(CMISSIntg), PARAMETER :: NumberOfLoadIncrements=2

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DisplacementBasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: FieldSourceUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program variables
  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx,component_idx,deriv_idx
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg) :: LeftNormalXi
  INTEGER(CMISSIntg) :: NumberOfArguments,ArgumentLength,ArgStatus
  CHARACTER(LEN=255) :: CommandArgument

  !CMISS variables
  TYPE(cmfe_BasisType) :: DisplacementBasis,PressureBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,FibreField,MaterialField,DependentField,SourceField,EquationsSetField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_ControlLoopType) :: ControlLoop

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
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  CALL cmfe_OutputSetOn("Cantilever",Err)

  !Read in arguments and overwrite default values
  !Usage: CantileverExample [Displacement Interpolation Type] [X elements] [Y elements] [Z elements] [Scaling Type]
  !Defaults:
  DisplacementInterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  NumberGlobalXElements=3
  NumberGlobalYElements=2
  NumberGlobalZElements=2
  ScalingType=CMFE_FIELD_ARITHMETIC_MEAN_SCALING

  NumberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(NumberOfArguments >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 1.")
    READ(CommandArgument(1:ArgumentLength),*) DisplacementInterpolationType
  ENDIF
  IF(NumberOfArguments >= 2) THEN
    CALL GET_COMMAND_ARGUMENT(2,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 2.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalXElements
    IF(NumberGlobalXElements<1) CALL HANDLE_ERROR("Invalid number of X elements.")
  ENDIF
  IF(NumberOfArguments >= 3) THEN
    CALL GET_COMMAND_ARGUMENT(3,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 3.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalYElements
    IF(NumberGlobalYElements<1) CALL HANDLE_ERROR("Invalid number of Y elements.")
  ENDIF
  IF(NumberOfArguments >= 4) THEN
    CALL GET_COMMAND_ARGUMENT(4,CommandArgument,ArgumentLength,ArgStatus)
    IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 4.")
    READ(CommandArgument(1:ArgumentLength),*) NumberGlobalZElements
    IF(NumberGlobalZElements<1) CALL HANDLE_ERROR("Invalid number of Z elements.")
  ENDIF
  IF(DisplacementInterpolationType==CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
    IF(NumberOfArguments >= 5) THEN
      CALL GET_COMMAND_ARGUMENT(5,CommandArgument,ArgumentLength,ArgStatus)
      IF(ArgStatus>0) CALL HANDLE_ERROR("Error for command argument 5.")
      READ(CommandArgument(1:ArgumentLength),*) ScalingType
      IF(ScalingType<0.OR.ScalingType>5) CALL HANDLE_ERROR("Invalid scaling type.")
    ENDIF
  ELSE
    ScalingType=CMFE_FIELD_NO_SCALING
  ENDIF
  SELECT CASE(DisplacementInterpolationType)
  CASE(CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION)
    NumberOfGaussXi=2
    PressureMeshComponent=1
  CASE(CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION)
    NumberOfGaussXi=3
    PressureMeshComponent=2
    PressureInterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  CASE(CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION)
    NumberOfGaussXi=4
    PressureMeshComponent=2
    !Should generally use quadratic interpolation but use linear to match CMISS example 5e
    PressureInterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  CASE DEFAULT
    NumberOfGaussXi=0
    PressureMeshComponent=1
  END SELECT
  WRITE(*,'("Interpolation: ", i3)') DisplacementInterpolationType
  WRITE(*,'("Elements: ", 3 i3)') NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  WRITE(*,'("Scaling type: ", i3)') ScalingType

  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberOfDomains=NumberOfComputationalNodes

  !Create a 3D rectangular cartesian coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_LabelSet(Region,"Region",Err)
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL cmfe_Region_CreateFinish(Region,Err)

  !Define basis function for displacement
  CALL cmfe_Basis_Initialise(DisplacementBasis,Err)
  CALL cmfe_Basis_CreateStart(DisplacementBasisUserNumber,DisplacementBasis,Err)
  SELECT CASE(DisplacementInterpolationType)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(DisplacementBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(DisplacementBasis,CMFE_BASIS_SIMPLEX_TYPE,Err)
  END SELECT
  CALL cmfe_Basis_NumberOfXiSet(DisplacementBasis,3,Err)
  CALL cmfe_Basis_InterpolationXiSet(DisplacementBasis,[DisplacementInterpolationType,DisplacementInterpolationType, &
      & DisplacementInterpolationType],Err)
  IF(NumberOfGaussXi>0) THEN
    CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(DisplacementBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
  ENDIF
  CALL cmfe_Basis_CreateFinish(DisplacementBasis,Err)

  IF(PressureMeshComponent/=1) THEN
    !Basis for pressure
    CALL cmfe_Basis_Initialise(PressureBasis,Err)
    CALL cmfe_Basis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
    SELECT CASE(PressureInterpolationType)
    CASE(1,2,3,4)
      CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CASE(7,8,9)
      CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_SIMPLEX_TYPE,Err)
    END SELECT
    CALL cmfe_Basis_NumberOfXiSet(PressureBasis,3,Err)
    CALL cmfe_Basis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType, &
        & PressureInterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
    CALL cmfe_Basis_CreateFinish(PressureBasis,Err)
  ENDIF

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  IF(PressureMeshComponent==1) THEN
    CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,[DisplacementBasis],Err)
  ELSE
    CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,[DisplacementBasis,PressureBasis],Err)
  ENDIF
  !Define the mesh on the region
  IF(NumberGlobalXElements==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[Width,Height],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[Width,Length,Height],Err)
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
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfDomains,Err)
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_ScalingTypeSet(GeometricField,ScalingType,Err)
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)

  !Create a fibre field and attach it to the geometric field
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FieldFibreUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_ScalingTypeSet(FibreField,ScalingType,Err)
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !Create the equations_set
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  DO component_idx=1,3
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,component_idx,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,component_idx,1,Err)
  ENDDO
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,PressureMeshComponent,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,PressureMeshComponent,Err)
  IF(PressureMeshComponent==1) THEN
    CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION, &
      & Err)
    CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  ELSE
    CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  ENDIF
  CALL cmfe_Field_ScalingTypeSet(DependentField,ScalingType,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_V_VARIABLE_TYPE,"Density",Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,2.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,6.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_V_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Density,Err)

  !Create the source field with the gravity vector
  CALL cmfe_Field_Initialise(SourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(EquationsSet,FieldSourceUserNumber,SourceField,Err)
  CALL cmfe_Field_ScalingTypeSet(SourceField,ScalingType,Err)
  CALL cmfe_EquationsSet_SourceCreateFinish(EquationsSet,Err)
  DO component_idx=1,3
    CALL cmfe_Field_ComponentValuesInitialise(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & component_idx,Gravity(component_idx),Err)
  ENDDO

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,4, &
    & -14.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateStart(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_NO_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(ControlLoop,NumberOfLoadIncrements,Err)
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(Solver,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
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
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)

  !Fix x=0 nodes in x, y and z
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      DO component_idx=1,3
        CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber, &
          & component_idx,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
        IF(DisplacementInterpolationType==CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
          DO deriv_idx=3,8
            CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,deriv_idx, &
              & NodeNumber, &
              & component_idx,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !Output solution
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"Cantilever","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"Cantilever","FORTRAN",Err)
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

END PROGRAM CANTILEVEREXAMPLE

