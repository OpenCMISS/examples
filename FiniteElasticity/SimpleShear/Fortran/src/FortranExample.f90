!> \file
!> \author Chris Bradley
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
!> Contributor(s): Kumar Mithraratne, Adam Reeve
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

!> \example FiniteElasticity/SimpleShear/src/FortranExample.f90
!! Example program to solve a finite elasticity equation using OpenCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/SimpleShear/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/SimpleShear/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM SIMPLESHEAREXAMPLE

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

  REAL(CMISSRP) :: VALUE

  !Test program parameters

  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP
!  INTEGER(CMISSIntg), PARAMETER :: InterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: InterpolationType=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: PressureInterpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
!  LOGICAL, PARAMETER :: UsePressureBasis=.TRUE.
  LOGICAL, PARAMETER :: UsePressureBasis=.FALSE.
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussXi=3

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldGeometryUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: FieldFibreUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: FieldMaterialUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: FieldDependentUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: EquationSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program types

  !Program variables

  INTEGER(CMISSIntg) :: NumberGlobalXElements,NumberGlobalYElements,NumberGlobalZElements
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,NumberOfDomains,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: NodeNumber,NodeDomain,node_idx
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: FrontSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: BackSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: BottomSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: TopSurfaceNodes(:)
  INTEGER(CMISSIntg) :: LeftNormalXi,RightNormalXi,FrontNormalXi,BackNormalXi,BottomNormalXi,TopNormalXi

  INTEGER(CMISSIntg), PARAMETER :: NUMBER_OF_COMPONENTS = 3 !nearly incompressible

  !CMISS variables
  TYPE(cmfe_BasisType) :: Basis, PressureBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,FibreField,MaterialField,DependentField,EquationsSetField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver,LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations

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

  !Intialise cmiss
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  !Set all diganostic levels on for testing
  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_RESIDUAL_EVALUATE"],Err)

  CALL cmfe_OutputSetOn("SimpleShear",Err)
  
  !Get the number of computational nodes and this computational node number
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberGlobalXElements=2
  NumberGlobalYElements=2
  NumberGlobalZElements=2
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

  !Define geometric basis
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  SELECT CASE(InterpolationType)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_SIMPLEX_TYPE,Err)
  END SELECT
  IF(NumberGlobalZElements==0) THEN
    CALL cmfe_Basis_NumberOfXiSet(Basis,2,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis,[InterpolationType,InterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ELSE
    CALL cmfe_Basis_NumberOfXiSet(Basis,3,Err)
    CALL cmfe_Basis_InterpolationXiSet(Basis,[InterpolationType,InterpolationType,InterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(Basis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ENDIF
  CALL cmfe_Basis_CreateFinish(Basis,Err)

  !Define pressure basis
  IF(UsePressureBasis) THEN
    CALL cmfe_Basis_Initialise(PressureBasis,Err)
    CALL cmfe_Basis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
    SELECT CASE(PressureInterpolationType)
    CASE(1,2,3,4)
      CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CASE(7,8,9)
      CALL cmfe_Basis_TypeSet(PressureBasis,CMFE_BASIS_SIMPLEX_TYPE,Err)
    END SELECT
    IF(NumberGlobalZElements==0) THEN
      CALL cmfe_Basis_NumberOfXiSet(PressureBasis,2,Err)
      CALL cmfe_Basis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType],Err)
      IF(NumberOfGaussXi>0) THEN
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
      ENDIF
    ELSE
      CALL cmfe_Basis_NumberOfXiSet(PressureBasis,3,Err)
      CALL cmfe_Basis_InterpolationXiSet(PressureBasis, &
        & [PressureInterpolationType,PressureInterpolationType,PressureInterpolationType],Err)
      IF(NumberOfGaussXi>0) THEN
        CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
      ENDIF
    ENDIF
    CALL cmfe_Basis_CreateFinish(PressureBasis,Err)
  ENDIF

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  IF(UsePressureBasis) THEN
    CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,[Basis,PressureBasis],Err)
  ELSE
    CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,[Basis],Err)
  ENDIF
  !Define the mesh on the region
  IF(NumberGlobalXElements==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
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

  !Create a field to put the geometry (default is geometry)
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL cmfe_Field_VariableLabelSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL cmfe_Field_ScalingTypeSet(GeometricField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
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
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !Create the dependent field
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_Field_CreateStart(FieldDependentUserNumber,Region,DependentField,Err)
  CALL cmfe_Field_TypeSet(DependentField,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DependentField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DependentField,GeometricField,Err)
  CALL cmfe_Field_DependentTypeSet(DependentField,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DependentField,2,Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  CALL cmfe_Field_NumberOfComponentsSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,Err)
  IF(UsePressureBasis) THEN
    !Set the pressure to be nodally based and use the second mesh component if required
    CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentInterpolationSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,2,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)
  END IF
  CALL cmfe_Field_CreateFinish(DependentField,Err)

  !Create the material field
  CALL cmfe_Field_Initialise(MaterialField,Err)
  CALL cmfe_Field_CreateStart(FieldMaterialUserNumber,Region,MaterialField,Err)
  CALL cmfe_Field_TypeSet(MaterialField,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MaterialField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(MaterialField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MaterialField,1,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL cmfe_Field_NumberOfComponentsSet(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,3,Err)
  CALL cmfe_Field_CreateFinish(MaterialField,Err)

  !Create the equations_set
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE], &
!    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, &
!    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field 
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 0.5 and 0.0 respectively. Third value is kappa (bulk modulus ???)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.5_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,0.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MaterialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,10000.0_CMISSRP,Err)

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs 
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)

  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_NO_SUBTYPE],Problem,Err)
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
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
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_TOP_SURFACE,TopSurfaceNodes,TopNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi, &
    & Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,FrontNormalXi,Err)
  CALL cmfe_GeneratedMesh_SurfaceGet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_BACK_SURFACE,BackSurfaceNodes,BackNormalXi,Err)

  ! set z=WIDTH nodes to 10% x-displacement, no displacement in y- and z-direction
  DO node_idx=1,SIZE(TopSurfaceNodes,1)
    NodeNumber=TopSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      ! x-direction
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,1,&
        & VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE+0.1_CMISSRP*WIDTH,Err)
      ! y-direction
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,2,&
        & VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      ! z-direction
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,3,&
        & VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  ENDDO

  ! set z=0 nodes to no displacement in x-, y- and z-direction
  DO node_idx=1,SIZE(BottomSurfaceNodes,1)
    NodeNumber=BottomSurfaceNodes(node_idx)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      ! x-direction
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,1,&
        & VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      ! y-direction
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,2,&
        & VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
      ! z-direction
      CALL cmfe_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1,NodeNumber,3,&
        & VALUE,Err)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMFE_BOUNDARY_CONDITION_FIXED,VALUE,Err)
    ENDIF
  ENDDO

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !Output solution
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"SimpleShear","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"SimpleShear","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM SIMPLESHEAREXAMPLE

