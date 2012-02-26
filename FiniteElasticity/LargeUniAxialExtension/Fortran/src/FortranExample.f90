!> \file
!> \author Chris Bradley
!> \brief This is an example program to solve a finite elasticity equation using openCMISS calls.
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
!> The Original Code is openCMISS
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

!> \example FiniteElasticity/LargeUniAxialExtension/src/LargeUniAxialExtensionExample.f90
!! Example program to solve a finite elasticity equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/LargeUniAxialExtension/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/FiniteElasticity/LargeUniAxialExtension/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM LARGEUNIAXIALEXTENSIONEXAMPLE

  USE OPENCMISS
  USE MPI

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=1.0_CMISSDP
  INTEGER(CMISSIntg), PARAMETER :: InterpolationType=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  INTEGER(CMISSIntg), PARAMETER :: PressureInterpolationType=CMISS_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  LOGICAL, PARAMETER :: UsePressureBasis=.FALSE.
  INTEGER(CMISSIntg), PARAMETER :: NumberOfGaussXi=2

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
  INTEGER(CMISSIntg),ALLOCATABLE :: BottomSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: LeftSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: RightSurfaceNodes(:)
  INTEGER(CMISSIntg),ALLOCATABLE :: FrontSurfaceNodes(:)
  INTEGER(CMISSIntg) :: BottomNormalXi,LeftNormalXi,RightNormalXi,BackNormalXi

  !CMISS variables
  TYPE(CMISSBasisType) :: Basis, PressureBasis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem, WorldCoordinateSystem
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,FibreField,MaterialField,DependentField,EquationsSetField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver,LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations

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
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  CALL CMISSErrorHandlingModeSet(CMISS_ERRORS_TRAP_ERROR,Err)

  !Set all diganostic levels on for testing
  !CALL CMISSDiagnosticsSetOn(CMISS_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["PROBLEM_RESIDUAL_EVALUATE"],Err)

  !Get the number of computational nodes and this computational node number
  CALL CMISSComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMISSComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  NumberGlobalXElements=5
  NumberGlobalYElements=5
  NumberGlobalZElements=5
  NumberOfDomains=NumberOfComputationalNodes

  !Create a 3D rectangular cartesian coordinate system
  CALL CMISSCoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMISSCoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Create a region and assign the coordinate system to the region
  CALL CMISSRegion_Initialise(Region,Err)
  CALL CMISSRegion_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL CMISSRegion_LabelSet(Region,"Region",Err)
  CALL CMISSRegion_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMISSRegion_CreateFinish(Region,Err)

  !Define geometric basis
  CALL CMISSBasis_Initialise(Basis,Err)
  CALL CMISSBasis_CreateStart(BasisUserNumber,Basis,Err)
  SELECT CASE(InterpolationType)
  CASE(1,2,3,4)
    CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CASE(7,8,9)
    CALL CMISSBasis_TypeSet(Basis,CMISS_BASIS_SIMPLEX_TYPE,Err)
  END SELECT
  IF(NumberGlobalZElements==0) THEN
    CALL CMISSBasis_NumberOfXiSet(Basis,2,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis,[InterpolationType,InterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis,[NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ELSE
    CALL CMISSBasis_NumberOfXiSet(Basis,3,Err)
    CALL CMISSBasis_InterpolationXiSet(Basis,[InterpolationType,InterpolationType,InterpolationType],Err)
    IF(NumberOfGaussXi>0) THEN
      CALL CMISSBasis_QuadratureNumberOfGaussXiSet(Basis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
    ENDIF
  ENDIF
  CALL CMISSBasis_CreateFinish(Basis,Err)

  !Define pressure basis
  IF(UsePressureBasis) THEN
    CALL CMISSBasis_Initialise(PressureBasis,Err)
    CALL CMISSBasis_CreateStart(PressureBasisUserNumber,PressureBasis,Err)
    SELECT CASE(PressureInterpolationType)
    CASE(1,2,3,4)
      CALL CMISSBasis_TypeSet(PressureBasis,CMISS_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
    CASE(7,8,9)
      CALL CMISSBasis_TypeSet(PressureBasis,CMISS_BASIS_SIMPLEX_TYPE,Err)
    END SELECT
    IF(NumberGlobalZElements==0) THEN
      CALL CMISSBasis_NumberOfXiSet(PressureBasis,2,Err)
      CALL CMISSBasis_InterpolationXiSet(PressureBasis,[PressureInterpolationType,PressureInterpolationType],Err)
      IF(NumberOfGaussXi>0) THEN
        CALL CMISSBasis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi],Err)
      ENDIF
    ELSE
      CALL CMISSBasis_NumberOfXiSet(PressureBasis,3,Err)
      CALL CMISSBasis_InterpolationXiSet(PressureBasis, &
        & [PressureInterpolationType,PressureInterpolationType,PressureInterpolationType],Err)
      IF(NumberOfGaussXi>0) THEN
        CALL CMISSBasis_QuadratureNumberOfGaussXiSet(PressureBasis,[NumberOfGaussXi,NumberOfGaussXi,NumberOfGaussXi],Err)
      ENDIF
    ENDIF
    CALL CMISSBasis_CreateFinish(PressureBasis,Err)
  ENDIF

  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL CMISSGeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMesh_TypeSet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  IF(UsePressureBasis) THEN
    CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,[Basis,PressureBasis],Err)
  ELSE
    CALL CMISSGeneratedMesh_BasisSet(GeneratedMesh,[Basis],Err)
  ENDIF
  !Define the mesh on the region
  IF(NumberGlobalXElements==0) THEN
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
    CALL CMISSGeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NumberGlobalXElements,NumberGlobalYElements],Err)
  ELSE
    CALL CMISSGeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
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
  CALL CMISSDecomposition_CreateFinish(Decomposition,Err)

  !Create a field to put the geometry (defualt is geometry)
  CALL CMISSField_Initialise(GeometricField,Err)
  CALL CMISSField_CreateStart(FieldGeometryUserNumber,Region,GeometricField,Err)
  CALL CMISSField_MeshDecompositionSet(GeometricField,Decomposition,Err)
  CALL CMISSField_VariableLabelSet(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
  CALL CMISSField_ScalingTypeSet(GeometricField,CMISS_FIELD_ARITHMETIC_MEAN_SCALING,Err)
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
  CALL CMISSField_CreateFinish(FibreField,Err)

  !Create the equations_set
  CALL CMISSField_Initialise(EquationsSetField,Err)
  CALL CMISSEquationsSet_CreateStart(EquationSetUserNumber,Region,FibreField,CMISS_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMISS_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMISS_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EquationsSetFieldUserNumber, &
      & EquationsSetField, &
    & EquationsSet,Err)
  CALL CMISSEquationsSet_CreateFinish(EquationsSet,Err)

  !Create the dependent field
  CALL CMISSField_Initialise(DependentField,Err)
  CALL CMISSEquationsSet_DependentCreateStart(EquationsSet,FieldDependentUserNumber,DependentField,Err)
  CALL CMISSField_VariableLabelSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
  IF(UsePressureBasis) THEN
    !Set the pressure to be nodally based and use the second mesh component if required
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentInterpolationSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4, &
      & CMISS_FIELD_NODE_BASED_INTERPOLATION,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,4,2,Err)
    CALL CMISSField_ComponentMeshComponentSet(DependentField,CMISS_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)
  END IF
  CALL CMISSEquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the material field
  CALL CMISSField_Initialise(MaterialField,Err)
  CALL CMISSEquationsSet_MaterialsCreateStart(EquationsSet,FieldMaterialUserNumber,MaterialField,Err)
  CALL CMISSField_VariableLabelSet(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,"Material",Err)
  CALL CMISSEquationsSet_MaterialsCreateFinish(EquationsSet,Err)

  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,2.0_CMISSDP,Err)
  CALL CMISSField_ComponentValuesInitialise(MaterialField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,6.0_CMISSDP,Err)

  !Create the equations set equations
  CALL CMISSEquations_Initialise(Equations,Err)
  CALL CMISSEquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  CALL CMISSEquations_SparsityTypeSet(Equations,CMISS_EQUATIONS_SPARSE_MATRICES,Err)
  CALL CMISSEquations_OutputTypeSet(Equations,CMISS_EQUATIONS_NO_OUTPUT,Err)
  CALL CMISSEquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 1,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,1,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 2,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,2,Err)
  CALL CMISSField_ParametersToFieldParametersComponentCopy(GeometricField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE, &
    & 3,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,3,Err)
  CALL CMISSField_ComponentValuesInitialise(DependentField,CMISS_FIELD_U_VARIABLE_TYPE,CMISS_FIELD_VALUES_SET_TYPE,4,-8.0_CMISSDP, &
    & Err)

  !Define the problem
  CALL CMISSProblem_Initialise(Problem,Err)
  CALL CMISSProblem_CreateStart(ProblemUserNumber,Problem,Err)
  CALL CMISSProblem_SpecificationSet(Problem,CMISS_PROBLEM_ELASTICITY_CLASS,CMISS_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMISS_PROBLEM_NO_SUBTYPE,Err)
  CALL CMISSProblem_CreateFinish(Problem,Err)

  !Create the problem control loop
  CALL CMISSProblem_ControlLoopCreateStart(Problem,Err)
  CALL CMISSProblem_ControlLoopCreateFinish(Problem,Err)

  !Create the problem solvers
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolver_Initialise(LinearSolver,Err)
  CALL CMISSProblem_SolversCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_OutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  CALL CMISSSolver_NewtonJacobianCalculationTypeSet(Solver,CMISS_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL CMISSSolver_NewtonLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolver_LinearTypeSet(LinearSolver,CMISS_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL CMISSProblem_SolversCreateFinish(Problem,Err)

  !Create the problem solver equations
  CALL CMISSSolver_Initialise(Solver,Err)
  CALL CMISSSolverEquations_Initialise(SolverEquations,Err)
  CALL CMISSProblem_SolverEquationsCreateStart(Problem,Err)
  CALL CMISSProblem_SolverGet(Problem,CMISS_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMISSSolver_SolverEquationsGet(Solver,SolverEquations,Err)
  CALL CMISSSolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  CALL CMISSProblem_SolverEquationsCreateFinish(Problem,Err)

  !Prescribe boundary conditions (absolute nodal parameters)
  CALL CMISSBoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMISSSolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)

  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_BOTTOM_SURFACE,BottomSurfaceNodes,BottomNormalXi, &
    & Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_LEFT_SURFACE,LeftSurfaceNodes,LeftNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_RIGHT_SURFACE,RightSurfaceNodes,RightNormalXi,Err)
  CALL CMISSGeneratedMesh_SurfaceGet(GeneratedMesh,CMISS_GENERATED_MESH_REGULAR_FRONT_SURFACE,FrontSurfaceNodes,BackNormalXi,Err)

  !Set x=0 nodes to no x displacment in x. Set x=WIDTH nodes to 10% x displacement
  DO node_idx=1,SIZE(LeftSurfaceNodes,1)
    NodeNumber=LeftSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO

  DO node_idx=1,SIZE(RightSurfaceNodes,1)
    NodeNumber=RightSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,1, &
        & CMISS_BOUNDARY_CONDITION_FIXED,1.1_CMISSDP*WIDTH,Err)
    ENDIF
  ENDDO

  !Set y=0 nodes to no y displacement
  DO node_idx=1,SIZE(FrontSurfaceNodes,1)
    NodeNumber=FrontSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,2, &
        & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO

  !Set z=0 nodes to no z displacement
  DO node_idx=1,SIZE(BottomSurfaceNodes,1)
    NodeNumber=BottomSurfaceNodes(node_idx)
    CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1,1,NodeNumber,3, &
        & CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
    ENDIF
  ENDDO

  IF(InterpolationType==CMISS_BASIS_CUBIC_HERMITE_INTERPOLATION) THEN
    !Fix x derivatives at x=0 and x=1 in xi2 and xi3
    DO node_idx=1,SIZE(LeftSurfaceNodes,1)
      NodeNumber=LeftSurfaceNodes(node_idx)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S3, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2_S3, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      ENDIF
    ENDDO
    DO node_idx=1,SIZE(RightSurfaceNodes,1)
      NodeNumber=RightSurfaceNodes(node_idx)
      CALL CMISSDecomposition_NodeDomainGet(Decomposition,NodeNumber,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S3, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
        CALL CMISSBoundaryConditions_SetNode(BoundaryConditions,DependentField,CMISS_FIELD_U_VARIABLE_TYPE,1, &
          & CMISS_GLOBAL_DERIV_S2_S3, &
          & NodeNumber,1,CMISS_BOUNDARY_CONDITION_FIXED,0.0_CMISSDP,Err)
      ENDIF
    ENDDO
  ENDIF

  CALL CMISSSolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve problem
  CALL CMISSProblem_Solve(Problem,Err)

  !Output solution
  CALL CMISSFields_Initialise(Fields,Err)
  CALL CMISSFields_Create(Region,Fields,Err)
  CALL CMISSFields_NodesExport(Fields,"LargeUniaxialExtension","FORTRAN",Err)
  CALL CMISSFields_ElementsExport(Fields,"LargeUniaxialExtension","FORTRAN",Err)
  CALL CMISSFields_Finalise(Fields,Err)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."

  STOP

END PROGRAM LARGEUNIAXIALEXTENSIONEXAMPLE

